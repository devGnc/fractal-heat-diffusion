/*
 * Test de validation : pièce circulaire comparée à la solution exacte.
 *
 * Géométrie : air (nœud unique à T_air imposée) dans un cercle de rayon Ra,
 * mur de béton de Ra à Rb = Ra + WALL_THICKNESS, extérieur à T_EXT. En régime
 * établi, la puissance par mètre de hauteur est (trois résistances en série)
 *
 *     P = (T_air − T_ext) / ( R_si/(2π Ra) + ln(Rb/Ra)/(2π λ) + R_se/(2π Rb) )
 *
 * Le test réutilise le code de la simulation (geometry_from_polygon,
 * compute_conductances, step_physics, compute_stats). Il est compilé à
 * plusieurs résolutions (make test) : Δx = 5 cm / s. Un cercle contient
 * toutes les orientations de mur : il mesure le biais de grille des surfaces
 * en escalier.
 *
 * Cas :
 *   - conduction   : R_si = R_se = 0, seule la conduction dans le mur compte ;
 *   - cond_brute   : idem, sans corrections de surface ;
 *   - surf_brutes  : R_si et R_se de la norme, sans corrections de surface ;
 *   - surfaces     : R_si et R_se, avec corrections (comme la simulation).
 * R_surf est la résistance des deux surfaces déduite de la simulation (en
 * retranchant la conduction, exacte), rapportée à sa valeur théorique.
 */
#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define R_ROOM   1.6                // m : rayon de la pièce
#define N_POLY   720                // sommets du polygone qui approche le cercle
#define T_AIR    30.0               // °C, imposée
#define CHECK_EVERY 500
#define MAX_ITER    4000000L
#define TOL_P       1e-9            // variation relative de P entre deux contrôles
#define TOL_BILAN   1e-7            // |P_ext − P_mur| / P_mur

typedef struct {
    const char *nom;
    double r_si, r_se;
    int raw;
} Cas;

static const Cas CAS[] = {
    { "cond_brute",  0.0,  0.0,  1 },
    { "conduction",  0.0,  0.0,  0 },
    { "surf_brutes", R_SI, R_SE, 1 },
    { "surfaces",    R_SI, R_SE, 0 },
};

static int run_case(const Cas *c) {
    const int G = GRID_SIZE, N = G * G;
    Sim s = {0};
    s.T = malloc(N * sizeof(double));
    s.Tnext = malloc(N * sizeof(double));
    s.lambda = calloc(N, sizeof(double));
    s.cap = calloc(N, sizeof(double));
    s.kx = calloc(N, sizeof(double));
    s.ky = calloc(N, sizeof(double));
    s.g_air = calloc(N, sizeof(double));
    s.g_ext = calloc(N, sizeof(double));
    s.zones = calloc(N, 1);
    s.poly = malloc(N_POLY * sizeof(Vec2));
    if (!s.T || !s.Tnext || !s.lambda || !s.cap || !s.kx || !s.ky ||
        !s.g_air || !s.g_ext || !s.zones || !s.poly) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        return 1;
    }

    // Polygone régulier d'aire π R² (même aire que le cercle de rayon R_ROOM)
    const double ctr = G * DX / 2.0;
    const double rp = R_ROOM * sqrt(2.0 * M_PI / (N_POLY * sin(2.0 * M_PI / N_POLY)));
    for (int k = 0; k < N_POLY; k++) {
        double t = 2.0 * M_PI * (k + 0.5) / N_POLY;
        s.poly[k] = (Vec2){ ctr + rp * cos(t), ctr + rp * sin(t) };
    }
    s.n_poly = N_POLY;
    s.r_si = c->r_si;
    s.r_se = c->r_se;
    s.raw_surfaces = c->raw;
    s.air_fixed = 1;
    geometry_from_polygon(&s);  // s'arrête si l'air touche l'extérieur
    s.T_air = T_AIR;

    const double Ra = R_ROOM, Rb = R_ROOM + WALL_THICKNESS, dT = T_AIR - T_EXT;
    const double r_cond = log(Rb / Ra) / (2.0 * M_PI * LAMBDA_WALL);
    const double r_surf = c->r_si / (2.0 * M_PI * Ra) + c->r_se / (2.0 * M_PI * Rb);
    const double p_exact = dT / (r_cond + r_surf);

    // Départ sur le profil exact : le régime établi n'en dépend pas, mais on
    // l'atteint bien plus vite.
    const double T_si = T_AIR - p_exact * c->r_si / (2.0 * M_PI * Ra);
    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            double r = hypot((x + 0.5) * DX - ctr, (y + 0.5) * DX - ctr);
            if (s.zones[i] == ZONE_INT) s.T[i] = T_AIR;
            else if (s.zones[i] == ZONE_WALL)
                s.T[i] = T_si - p_exact * log(fmin(fmax(r, Ra), Rb) / Ra)
                                 / (2.0 * M_PI * LAMBDA_WALL);
        }
    }

    double p_prev = 0.0;
    int converged = 0;
    while (s.iteration < MAX_ITER) {
        step_physics(&s);
        if (s.iteration % CHECK_EVERY) continue;
        compute_stats(&s);
        double p = s.stat_p_mur;
        if (fabs(p - p_prev) < TOL_P * p && fabs(s.stat_p_ext - p) < TOL_BILAN * p) {
            converged = 1;
            break;
        }
        p_prev = p;
    }

    double p = s.stat_p_mur;
    double err = 100.0 * (p - p_exact) / p_exact;
    double bilan = 100.0 * (s.stat_p_ext - p) / p;
    double r_surf_ratio = r_surf > 0.0 ? (dT / p - r_cond) / r_surf : NAN;
    printf("%-12s %d %4d %10.3f %10.3f %+8.3f %+9.6f %8.3f %8ld%s\n",
           c->nom, (int)lround(0.05 / DX), G, p, p_exact, err, bilan, r_surf_ratio,
           s.iteration, converged ? "" : "  NON CONVERGÉ");

    free(s.T); free(s.Tnext); free(s.lambda); free(s.cap); free(s.kx); free(s.ky);
    free(s.g_air); free(s.g_ext); free(s.zones); free(s.poly);
    return !converged;
}

int main(void) {
    int fail = 0;
    for (size_t k = 0; k < sizeof CAS / sizeof CAS[0]; k++) {
        fail |= run_case(&CAS[k]);
        fflush(stdout);
    }
    return fail;
}
