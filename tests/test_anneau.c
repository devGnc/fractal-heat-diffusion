/*
 * Test de validation : pièce circulaire comparée à la solution exacte.
 *
 * Géométrie : source de puissance Q (disque de rayon r_s au centre), air de
 * conductivité effective λ_a jusqu'au rayon Ra, mur de béton de Ra à
 * Rb = Ra + WALL_THICKNESS, extérieur à T_EXT. En régime établi, la conduction
 * radiale donne le profil exact (q = Q / 2π) :
 *
 *   T(Ra) = T_ext + Q · ( R_si/(2π Ra) + ln(Rb/Ra)/(2π λ_m) + R_se/(2π Rb) )
 *   T(r)  = T(Ra) + (q/λ_a) ln(Ra/r)                     pour r_s ≤ r ≤ Ra
 *   T(r)  = T(r_s) + (q/2λ_a) (1 − r²/r_s²)              dans la source
 *
 * On compare la température moyenne de l'air (simulée et exacte, intégrée
 * sur le disque de rayon Ra), via G = Q / (T_air − T_ext).
 *
 * Le test réutilise le code de la simulation (geometry_from_polygon,
 * compute_conductances, step_physics, compute_stats). Il est compilé à
 * plusieurs résolutions (make test) : Δx = 5 cm / s. Un cercle contient
 * toutes les orientations de mur : il mesure le biais de grille des surfaces
 * en escalier.
 *
 * Cas :
 *   - conduction   : R_si = R_se = 0, avec corrections de surface ;
 *   - surf_brutes  : R_si et R_se de la norme, sans corrections de surface ;
 *   - surfaces     : R_si et R_se, avec corrections (comme la simulation).
 */
#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define R_ROOM   1.6                // m : rayon de la pièce
#define N_POLY   720                // sommets du polygone qui approche le cercle
#define CHECK_EVERY 200
#define MAX_ITER    400000L
#define TOL_T       1e-6            // variation de T_air entre deux contrôles (°C)
#define TOL_BILAN   1e-6            // |P_ext − Q| / Q

typedef struct {
    const char *nom;
    double r_si, r_se;
    int raw;
} Cas;

static const Cas CAS[] = {
    { "conduction",  0.0,  0.0,  0 },
    { "surf_brutes", R_SI, R_SE, 1 },
    { "surfaces",    R_SI, R_SE, 0 },
};

// Profil exact T(r)
typedef struct { double Q, rs, Ra, Rb, la, Ta; } Exact;

static double exact_T(const Exact *e, double r) {
    double q = e->Q / (2.0 * M_PI);
    if (r >= e->rs) return e->Ta + q / e->la * log(e->Ra / fmax(r, 1e-12));
    double Ts = e->Ta + q / e->la * log(e->Ra / e->rs);
    return Ts + q / (2.0 * e->la) * (1.0 - r * r / (e->rs * e->rs));
}

// Moyenne sur le disque de rayon Ra (méthode de Simpson)
static double exact_mean(const Exact *e) {
    const int n = 20000;
    double h = e->Ra / n, sum = 0.0;
    for (int k = 0; k <= n; k++) {
        double r = k * h, w = (k == 0 || k == n) ? 1 : (k % 2 ? 4 : 2);
        sum += w * exact_T(e, r) * r;
    }
    return 2.0 / (e->Ra * e->Ra) * sum * h / 3.0;
}

static int run_case(const Cas *c) {
    const int G = GRID_SIZE;
    Sim s = {0};
    if (!sim_alloc(&s)) {
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
    geometry_from_polygon(&s);  // s'arrête si l'air touche l'extérieur

    // Solution exacte, avec le rayon de source effectivement maillé
    Exact e = { s.q_int, sqrt(s.n_src * DX * DX / M_PI), R_ROOM, R_ROOM + WALL_THICKNESS,
                LAMBDA_AIR, 0.0 };
    const double r_wall = c->r_si / (2.0 * M_PI * e.Ra)
                        + log(e.Rb / e.Ra) / (2.0 * M_PI * LAMBDA_WALL)
                        + c->r_se / (2.0 * M_PI * e.Rb);
    e.Ta = T_EXT + e.Q * r_wall;
    const double T_mean = exact_mean(&e);
    const double G_exact = e.Q / (T_mean - T_EXT);

    // Départ sur le profil exact : le régime établi n'en dépend pas, mais on
    // l'atteint bien plus vite.
    const double T_si = T_EXT + e.Q * (r_wall - c->r_si / (2.0 * M_PI * e.Ra));
    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            double r = hypot((x + 0.5) * DX - ctr, (y + 0.5) * DX - ctr);
            if (s.zones[i] == ZONE_INT) s.T[i] = exact_T(&e, r);
            else if (s.zones[i] == ZONE_WALL)
                s.T[i] = T_si - e.Q * log(fmin(fmax(r, e.Ra), e.Rb) / e.Ra)
                                 / (2.0 * M_PI * LAMBDA_WALL);
        }
    }

    double prev = 0.0;
    int converged = 0;
    while (s.iteration < MAX_ITER) {
        step_physics(&s);
        if (s.iteration % CHECK_EVERY) continue;
        compute_stats(&s);
        if (fabs(s.stat_T_air - prev) < TOL_T &&
            fabs(s.stat_p_ext - s.q_int) < TOL_BILAN * s.q_int) {
            converged = 1;
            break;
        }
        prev = s.stat_T_air;
    }

    double G_sim = s.q_int / (s.stat_T_air - T_EXT);
    printf("%-12s %d %4d %9.4f %9.4f %9.3f %9.3f %+8.3f %+10.6f %8ld%s\n",
           c->nom, (int)lround(0.05 / DX), G, s.stat_T_air, T_mean, G_sim, G_exact,
           100.0 * (G_sim - G_exact) / G_exact,
           100.0 * (s.stat_p_ext - s.q_int) / s.q_int,
           s.iteration, converged ? "" : "  NON CONVERGÉ");

    sim_free(&s);
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
