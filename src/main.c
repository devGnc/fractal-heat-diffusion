#include "thermal_diffusion.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define CSV_HEADER "Iter,Time_h,T_air,T_mur,P_int,P_mur,P_ext\n"
#define MAX_DAYS 60          // mode sans fenêtre : durée maximale simulée
#define TOL_STEADY 1e-5      // régime établi : |P_ext − Q_int| < TOL · Q_int

static int sim_alloc(Sim *s) {
    int N = GRID_SIZE * GRID_SIZE;
    s->T = calloc(N, sizeof(double));
    s->Tnext = calloc(N, sizeof(double));
    s->lambda = calloc(N, sizeof(double));
    s->cap = calloc(N, sizeof(double));
    s->kx = calloc(N, sizeof(double));
    s->ky = calloc(N, sizeof(double));
    s->g_air = calloc(N, sizeof(double));
    s->g_ext = calloc(N, sizeof(double));
    s->zones = calloc(N, 1);
    s->pixels = malloc(N * sizeof(uint32_t));
    s->poly = malloc(POLY_MAX * sizeof(Vec2));
    s->q_int = Q_INT;
    s->r_si = R_SI;
    s->r_se = R_SE;
    return s->T && s->Tnext && s->lambda && s->cap && s->kx && s->ky &&
           s->g_air && s->g_ext && s->zones && s->pixels && s->poly;
}

static void sim_free(Sim *s) {
    free(s->T); free(s->Tnext); free(s->lambda); free(s->cap);
    free(s->kx); free(s->ky); free(s->g_air); free(s->g_ext);
    free(s->zones); free(s->pixels); free(s->poly);
}

// (Ré)initialise la géométrie du mode courant et remet le temps à zéro.
static void sim_reset(Sim *s) {
    s->iteration = 0;
    s->sim_time = 0.0;
    init_geometry(s);
    compute_stats(s);
}

static void geom_name(const Sim *s, char *buf, size_t n) {
    if (s->mode == 0) snprintf(buf, n, "carre");
    else snprintf(buf, n, "koch%d", s->depth);
}

// En-tête : une ligne de commentaire décrivant la géométrie (lue par les
// scripts avec pandas.read_csv(..., comment="#")), puis les noms de colonnes.
static void csv_header(FILE *f, const Sim *s) {
    char name[16];
    geom_name(s, name, sizeof name);
    fprintf(f, "# geometrie=%s profondeur=%d aire_m2=%.4f perimetre_m=%.4f "
               "px_air=%d px_mur=%d faces_air_mur=%d faces_mur_ext=%d "
               "dx=%g dt=%g epaisseur=%g t_ext=%g q_int=%g\n",
            name, s->mode ? s->depth : -1, s->area, s->perimeter,
            s->n_int, s->n_wall, s->n_faces_int_wall, s->n_faces_wall_ext,
            (double)DX, (double)DT, (double)WALL_THICKNESS, T_EXT, s->q_int);
    fprintf(f, CSV_HEADER);
}

static void csv_row(FILE *f, const Sim *s) {
    fprintf(f, "%ld,%.4f,%.5f,%.5f,%.3f,%.4f,%.4f\n",
            s->iteration, s->sim_time / 3600.0, s->T_air, s->stat_T_wall,
            s->q_int, s->stat_p_mur, s->stat_p_ext);
}

static void stop_recording(Sim *s) {
    if (!s->recording) return;
    fclose(s->csv);
    s->csv = NULL;
    s->recording = 0;
    printf("Enregistrement CSV arrêté.\n");
}

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
    check_stability(&s);

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
    long max_iter = n_iter ? n_iter : (long)(MAX_DAYS * 86400.0 / DT);
    int steady = 0;
    while (s.iteration < max_iter) {
        step_physics(&s);
        if (s.iteration % STATS_EVERY == 0 || s.iteration == max_iter) {
            compute_stats(&s);
            csv_row(f, &s);
            if (s.iteration % (24 * STATS_EVERY) == 0) {
                fprintf(stderr, "  jour %3.0f  T_air = %6.3f °C  P_ext = %7.2f W/m\n",
                        s.sim_time / 86400.0, s.T_air, s.stat_p_ext);
            }
            if (fabs(s.stat_p_ext - s.q_int) < TOL_STEADY * s.q_int) {
                steady = 1;
                if (!n_iter) break;
            }
        }
    }
    fclose(f);

    double dT = s.T_air - T_EXT;
    printf("%s : %s après %.1f jours → %s\n", name,
           steady ? "régime établi" : "RÉGIME NON ÉTABLI", s.sim_time / 86400.0, out);
    printf("  T_air = %.3f °C, T_mur moyenne = %.3f °C\n", s.T_air, s.stat_T_wall);
    printf("  Bilan : Q_int = %.2f, P_mur = %.2f, P_ext = %.2f W/m (écart %+.4f %%)\n",
           s.q_int, s.stat_p_mur, s.stat_p_ext, 100.0 * (s.stat_p_ext - s.q_int) / s.q_int);
    printf("  Conductance de l'enveloppe G = P_ext / (T_air − T_ext) = %.3f W/(m·K)\n",
           s.stat_p_ext / dT);
    sim_free(&s);
    return !steady;
}

/* ------------------------------------------------------------------------ */
/* Mode interactif SDL                                                      */
/* ------------------------------------------------------------------------ */

static int run_gui(void) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL Init fail: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow(
        "Diffusion thermique et géométries fractales",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_SIZE, WIN_SIZE, 0
    );
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);

    Sim s = {0};
    if (!sim_alloc(&s)) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        return 1;
    }
    s.tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
                              SDL_TEXTUREACCESS_STREAMING, GRID_SIZE, GRID_SIZE);

    printf("\n============================================\n");
    printf("  Diffusion thermique et géométries fractales\n");
    printf("============================================\n\n");
    printf("MODÈLE PHYSIQUE (par mètre de hauteur) :\n");
    printf("  - Pièce de %.0f m², air brassé (une seule température)\n", ROOM_AREA);
    printf("  - Murs en béton de %.0f cm : λ = %.2f W/(m·K), ρc = %.1e J/(m³·K)\n",
           100.0 * WALL_THICKNESS, LAMBDA_WALL, RHOC_WALL);
    printf("  - Surfaces : R_si = %.2f, R_se = %.2f m²·K/W\n", R_SI, R_SE);
    printf("  - Apports internes : %.0f W/m ; extérieur à %.0f °C\n", Q_INT, T_EXT);
    printf("  - Grille : %d × %d, Δx = %.0f cm, Δt = %.0f s\n\n",
           GRID_SIZE, GRID_SIZE, 100.0 * DX, (double)DT);

    s.mode = 0;
    s.depth = KOCH_DEPTH_REF;
    sim_reset(&s);
    check_stability(&s);

    printf("\nCOMMANDES :\n");
    printf("  ESPACE  - Changer de mode (carré ↔ fractale)\n");
    printf("  0 à %d   - Flocon de Koch de cette profondeur\n", KOCH_DEPTH_MAX);
    printf("  S       - Démarrer/arrêter enregistrement CSV\n");
    printf("  R       - Réinitialiser la simulation\n");
    printf("  Q       - Quitter\n");
    printf("Couleurs : de %.0f °C (bleu) à %.0f °C (jaune) ; l'air est éclairci.\n",
           T_EXT, T_DISPLAY_MAX);
    printf("Pour des données exploitables : --headless (voir README)\n");
    printf("============================================\n\n");

    int running = 1;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = 0;
            }
            if (e.type == SDL_KEYDOWN) {
                SDL_Keycode k = e.key.keysym.sym;
                if (k == SDLK_SPACE || k == SDLK_r ||
                    (k >= SDLK_0 && k <= SDLK_0 + KOCH_DEPTH_MAX)) {
                    // Un CSV ne doit contenir qu'une seule simulation
                    stop_recording(&s);
                    if (k == SDLK_SPACE) s.mode = !s.mode;
                    if (k >= SDLK_0 && k <= SDLK_0 + KOCH_DEPTH_MAX) {
                        s.mode = 1;
                        s.depth = k - SDLK_0;
                    }
                    sim_reset(&s);
                }
                else if (k == SDLK_s) {
                    if (s.recording) {
                        stop_recording(&s);
                    } else {
                        char gname[16], name[64];
                        geom_name(&s, gname, sizeof gname);
                        snprintf(name, sizeof name, "res_%s.csv", gname);
                        s.csv = fopen(name, "w");
                        if (s.csv) {
                            csv_header(s.csv, &s);
                            s.recording = 1;
                            printf("Enregistrement vers %s (à partir de t = %.1f h)\n",
                                   name, s.sim_time / 3600.0);
                        }
                    }
                }
                else if (k == SDLK_q) {
                    running = 0;
                }
            }
        }

        for (int k = 0; k < ITER_PER_FRAME; k++) {
            step_physics(&s);
        }

        if (s.iteration % STATS_EVERY == 0) {
            compute_stats(&s);
            if (s.recording) {
                csv_row(s.csv, &s);
                fflush(s.csv);
            }
        }

        SDL_RenderClear(ren);
        render(ren, &s);
        SDL_RenderPresent(ren);

        char gname[16], title[256];
        if (s.mode) snprintf(gname, sizeof gname, "KOCH %d", s.depth);
        else snprintf(gname, sizeof gname, "CARRÉ");
        snprintf(title, sizeof title,
                 "%s | t = %.1f h | T_air = %.2f °C | P_ext = %.0f W/m | G = %.1f W/(m·K) %s",
                 gname, s.sim_time / 3600.0, s.T_air, s.stat_p_ext,
                 s.T_air > T_EXT ? s.stat_p_ext / (s.T_air - T_EXT) : 0.0,
                 s.recording ? "[REC]" : "");
        SDL_SetWindowTitle(win, title);

        SDL_Delay(10);
    }

    stop_recording(&s);
    sim_free(&s);
    SDL_DestroyTexture(s.tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();

    return 0;
}

static void usage(const char *prog) {
    fprintf(stderr,
            "Usage :\n"
            "  %s                                       fenêtre interactive\n"
            "  %s --headless <geometrie> [iter] [csv]   simulation sans fenêtre\n"
            "      geometrie : carre, koch0 ... koch%d (fractal = koch%d)\n"
            "      iter      : nombre d'itérations (défaut 0 : jusqu'au régime établi)\n"
            "      csv       : fichier de sortie (défaut res_<geometrie>.csv)\n",
            prog, prog, KOCH_DEPTH_MAX, KOCH_DEPTH_REF);
}

int main(int argc, char **argv) {
    if (argc == 1) return run_gui();

    if (strcmp(argv[1], "--headless") != 0 || argc < 3 || argc > 5) {
        usage(argv[0]);
        return 1;
    }

    int mode = 1, depth = KOCH_DEPTH_REF;
    const char *g = argv[2];
    if (strcmp(g, "carre") == 0) {
        mode = 0;
    } else if (strcmp(g, "fractal") == 0) {
        depth = KOCH_DEPTH_REF;
    } else if (strncmp(g, "koch", 4) == 0 && g[4] >= '0' && g[4] <= '0' + KOCH_DEPTH_MAX
               && g[5] == '\0') {
        depth = g[4] - '0';
    } else {
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

    char default_out[64];
    if (mode == 0) snprintf(default_out, sizeof default_out, "res_carre.csv");
    else snprintf(default_out, sizeof default_out, "res_koch%d.csv", depth);
    const char *out = argc == 5 ? argv[4] : default_out;

    return run_headless(mode, depth, n_iter, out);
}
