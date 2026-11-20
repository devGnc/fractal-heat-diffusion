#include "thermal_diffusion.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_DAYS 60          // mode sans fenêtre : durée maximale simulée
#define TOL_STEADY 1e-5      // régime établi : |P_ext − Q_int| < TOL · Q_int

/* ------------------------------------------------------------------------ */
/* Mode sans fenêtre : simulation complète enregistrée depuis l'itération 0  */
/* ------------------------------------------------------------------------ */

// n_iter = 0 : jusqu'au régime établi (au plus MAX_DAYS jours simulés)
static int run_headless(int mode, int depth, long n_iter, const char *out) {
    Sim s = {0};
    if (!sim_alloc(&s)) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        return 1;
    }
    s.mode = mode;
    s.depth = depth;
    sim_reset(&s);

    FILE *f = fopen(out, "w");
    if (!f) {
        perror(out);
        sim_free(&s);
        return 1;
    }
    csv_header(f, &s);
    csv_row(f, &s);

    char name[16];
    geom_name(&s, name, sizeof name);
    long max_iter = n_iter ? n_iter : (long)(MAX_DAYS * 86400.0 / s.dt);
    const long per_hour = (long)lround(3600.0 / s.dt);
    int steady = 0;
    while (s.iteration < max_iter) {
        step_physics(&s);
        if (s.iteration % per_hour == 0 || s.iteration == max_iter) {
            compute_stats(&s);
            csv_row(f, &s);
            if (s.iteration % (24 * per_hour) == 0) {
                fprintf(stderr, "  jour %3.0f  T_air = %6.3f °C  P_ext = %7.2f W/m\n",
                        s.sim_time / 86400.0, s.stat_T_air, s.stat_p_ext);
            }
            if (fabs(s.stat_p_ext - s.q_int) < TOL_STEADY * s.q_int) {
                steady = 1;
                if (!n_iter) break;
            }
        }
    }
    fclose(f);

    double dT = s.stat_T_air - T_EXT;
    printf("%s : %s après %.1f jours → %s\n", name,
           steady ? "régime établi" : "RÉGIME NON ÉTABLI", s.sim_time / 86400.0, out);
    printf("  T_air moyenne = %.3f °C (max %.3f), T_mur moyenne = %.3f °C\n",
           s.stat_T_air, s.stat_T_air_max, s.stat_T_wall);
    printf("  Bilan : Q_int = %.2f, P_mur = %.2f, P_ext = %.2f W/m (écart %+.4f %%)\n",
           s.q_int, s.stat_p_mur, s.stat_p_ext, 100.0 * (s.stat_p_ext - s.q_int) / s.q_int);
    printf("  Conductance G = P_ext / (T_air moyenne − T_ext) = %.3f W/(m·K)\n",
           s.stat_p_ext / dT);
    sim_free(&s);
    return !steady;
}

// Régime établi résolu directement (gradient conjugué), sans transitoire.
// Écrit un CSV au même format, avec une seule ligne (le régime établi).
static int run_steady(int mode, int depth, const char *out) {
    Sim s = {0};
    if (!sim_alloc(&s)) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        return 1;
    }
    s.mode = mode;
    s.depth = depth;
    sim_reset(&s);

    char name[16];
    geom_name(&s, name, sizeof name);
    clock_t t0 = clock();
    int it = solve_steady(&s, 1e-10, 200000);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    compute_stats(&s);
    if (it < 0) {
        fprintf(stderr, "%s : le gradient conjugué n'a pas convergé\n", name);
        sim_free(&s);
        return 1;
    }

    FILE *f = fopen(out, "w");
    if (!f) {
        perror(out);
        sim_free(&s);
        return 1;
    }
    csv_header(f, &s);
    csv_row(f, &s);
    fclose(f);

    double dT = s.stat_T_air - T_EXT;
    printf("%s : régime établi en %d itérations (%.1f s, grille %d × %d, Δx = %g cm) → %s\n",
           name, it, secs, GRID_SIZE, GRID_SIZE, 100.0 * DX, out);
    printf("  T_air moyenne = %.3f °C (max %.3f), T_mur moyenne = %.3f °C\n",
           s.stat_T_air, s.stat_T_air_max, s.stat_T_wall);
    printf("  Bilan : Q_int = %.2f, P_mur = %.2f, P_ext = %.2f W/m (écart %+.5f %%)\n",
           s.q_int, s.stat_p_mur, s.stat_p_ext, 100.0 * (s.stat_p_ext - s.q_int) / s.q_int);
    printf("  Conductance G = P_ext / (T_air moyenne − T_ext) = %.3f W/(m·K)\n",
           s.stat_p_ext / dT);
    sim_free(&s);
    return 0;
}

static void usage(const char *prog) {
    fprintf(stderr,
            "Usage :\n"
            "  %s                                       fenêtre interactive\n"
            "  %s --headless <geometrie> [iter] [csv]   simulation dans le temps, sans fenêtre\n"
            "      iter : nombre d'itérations (défaut 0 : jusqu'au régime établi)\n"
            "      csv  : fichier de sortie (défaut res_<geometrie>.csv)\n"
            "  %s --regime <geometrie> [csv]            régime établi seul (gradient conjugué)\n"
            "      csv  : fichier de sortie (défaut regime_<geometrie>.csv)\n"
            "  %s --capture <image.bmp> [profondeur] [vue] [heures]\n"
            "      image de la fenêtre sans l'ouvrir (vue : temp ou flux ;\n"
            "      défaut : profondeur 3, temp, 168 h)\n"
            "  geometrie : carre, koch0 ... koch%d (fractal = koch%d)\n",
            prog, prog, prog, prog, KOCH_DEPTH_MAX, KOCH_DEPTH_REF);
}

// Lit une géométrie (carre, koch<n>, fractal) ; renvoie 0 si inconnue
static int parse_geometry(const char *g, int *mode, int *depth) {
    *mode = 1;
    *depth = KOCH_DEPTH_REF;
    if (strcmp(g, "carre") == 0) *mode = 0;
    else if (strcmp(g, "fractal") == 0) return 1;
    else if (strncmp(g, "koch", 4) == 0 && g[4] >= '0' && g[4] <= '0' + KOCH_DEPTH_MAX
             && g[5] == '\0') *depth = g[4] - '0';
    else return 0;
    return 1;
}

int main(int argc, char **argv) {
    if (argc == 1) return run_gui();

    if (strcmp(argv[1], "--capture") == 0 && argc >= 3 && argc <= 6) {
        int depth = argc >= 4 ? atoi(argv[3]) : KOCH_DEPTH_REF + 1;
        int view = argc >= 5 && strcmp(argv[4], "flux") == 0 ? VIEW_FLUX : VIEW_TEMP;
        double hours = argc >= 6 ? atof(argv[5]) : 168.0;
        if (depth < 0 || depth > KOCH_DEPTH_MAX || hours < 0) {
            usage(argv[0]);
            return 1;
        }
        return capture_gui(argv[2], depth, view, hours);
    }

    int mode, depth;
    if (argc < 3 || !parse_geometry(argv[2], &mode, &depth)) {
        usage(argv[0]);
        return 1;
    }
    char name[16];
    if (mode == 0) snprintf(name, sizeof name, "carre");
    else snprintf(name, sizeof name, "koch%d", depth);

    if (strcmp(argv[1], "--regime") == 0 && argc <= 4) {
        char out[64];
        snprintf(out, sizeof out, "regime_%s.csv", name);
        return run_steady(mode, depth, argc == 4 ? argv[3] : out);
    }

    if (strcmp(argv[1], "--headless") != 0 || argc > 5) {
        usage(argv[0]);
        return 1;
    }
    long n_iter = 0;
    if (argc >= 4) {
        char *end;
        n_iter = strtol(argv[3], &end, 10);
        if (*end != '\0' || n_iter < 0) {
            fprintf(stderr, "Nombre d'itérations invalide : %s\n", argv[3]);
            return 1;
        }
    }
    char out[64];
    snprintf(out, sizeof out, "res_%s.csv", name);
    return run_headless(mode, depth, n_iter, argc == 5 ? argv[4] : out);
}
