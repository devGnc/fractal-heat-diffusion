/*
 * Test de validation : anneau circulaire comparé à la solution exacte.
 *
 * Géométrie : source (disque de rayon R1 à SOURCE_TEMP), air jusqu'au rayon Ra,
 * mur de Ra à Rb, extérieur (thermostat à AMBIENT_TEMP) au-delà. En régime
 * établi, la conduction radiale donne une puissance par unité d'épaisseur
 *
 *     P = 2π (T_source − T_ext) / ( ln(Ra/R1) / α_air + ln(Rb/Ra) / α_mur )
 *
 * (résistances thermiques des deux couches en série ; capacité 1 par pixel,
 * donc α joue ici le rôle de la conductivité λ).
 *
 * Le test réutilise le code de la simulation (compute_conductances,
 * step_physics, compute_stats, draw_line, flood_fill_safe). Il est compilé à
 * plusieurs résolutions (make test) : SOURCE_RADIUS = 8·s, mur au rayon 32·s.
 * Un cercle contient toutes les orientations de mur, il mesure donc aussi le
 * biais de grille des murs obliques.
 *
 * Conventions de rayon (positions des faces entre cases) :
 *   - mur fin (1 px, tracé avec draw_line comme les flocons) : centres des
 *     cases de mur au rayon R ; Ra = R − 1/2, Rb = R + 1/2 ;
 *   - mur épais : air si d ≤ R, mur si R < d ≤ R + e ; Ra = R + 1/2,
 *     Rb = R + e + 1/2.
 * Ces conventions sont exactes à O(Δx) près : c'est la convergence quand la
 * résolution augmente qui valide le schéma, pas l'erreur à une résolution.
 */
#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if COOLING_FACTOR != 0
#error "Le test se compile avec -DCOOLING_FACTOR=0 (pas de pertes de Newton)"
#endif

#define S       (SOURCE_RADIUS / 8)   // facteur de résolution
#define R1      SOURCE_RADIUS
#define R_MUR   (32 * S)
#define E_EPAIS (4 * S)               // mur épais : épaisseur physique constante

_Static_assert(SOURCE_RADIUS % 8 == 0, "SOURCE_RADIUS doit valoir 8·s");
_Static_assert(GRID_SIZE / 2 > R_MUR + E_EPAIS + 2, "grille trop petite");

#define CHECK_EVERY 2000
#define MAX_ITER    4000000L
#define TOL_P       1e-8    // variation relative de P entre deux contrôles
#define TOL_BILAN   1e-6    // |P_ext − P_source| / P_source

typedef struct {
    const char *nom;
    int epais;          // 0 : mur d'un pixel (draw_line), 1 : mur épais
    double a_air, a_mur;
} Cas;

static const Cas CAS[] = {
    { "homogene",       0, 0.20, 0.20  },  // α_mur = α_air : milieu uniforme
    { "mur_fin",        0, 0.20, 0.02  },  // paramètres de la v2
    { "mur_fin_isol",   0, 0.20, 0.002 },  // le mur domine la résistance
    { "mur_epais",      1, 0.20, 0.02  },  // e = 4·s pixels
};

static void build_ring(Sim *s, const Cas *c, double *Ra, double *Rb) {
    const int G = GRID_SIZE, cx = G / 2, cy = G / 2;
    for (int i = 0; i < G * G; i++) s->zones[i] = ZONE_EXT;

    if (!c->epais) {
        // Polygone régulier à segments d'environ 5 px, comme les côtés d'un
        // flocon de profondeur 3, tracé avec la même fonction draw_line.
        int n = (int)lround(2.0 * M_PI * R_MUR / 5.0);
        for (int k = 0; k < n; k++) {
            double t0 = 2.0 * M_PI * k / n, t1 = 2.0 * M_PI * (k + 1) / n;
            Point p0 = { (int)lround(cx + R_MUR * cos(t0)), (int)lround(cy + R_MUR * sin(t0)) };
            Point p1 = { (int)lround(cx + R_MUR * cos(t1)), (int)lround(cy + R_MUR * sin(t1)) };
            draw_line(s->zones, p0, p1, ZONE_WALL);
        }
        flood_fill_safe(s->zones, cx, cy);
        // Rayon moyen du polygone (entre apothème et rayon circonscrit)
        double r = R_MUR * (1.0 + cos(M_PI / n)) / 2.0;
        *Ra = r - 0.5;
        *Rb = r + 0.5;
    } else {
        for (int y = 0; y < G; y++) {
            for (int x = 0; x < G; x++) {
                double d = hypot(x - cx, y - cy);
                if (d <= R_MUR) s->zones[y * G + x] = ZONE_INT;
                else if (d <= R_MUR + E_EPAIS) s->zones[y * G + x] = ZONE_WALL;
            }
        }
        *Ra = R_MUR + 0.5;
        *Rb = R_MUR + E_EPAIS + 0.5;
    }

    for (int i = 0; i < G * G; i++) {
        s->alpha[i] = s->zones[i] == ZONE_INT ? c->a_air
                    : s->zones[i] == ZONE_WALL ? c->a_mur : 0.0;
    }
}

static int run_case(const Cas *c) {
    const int G = GRID_SIZE, N = G * G;
    Sim s = {0};
    s.T = malloc(N * sizeof(double));
    s.Tnext = malloc(N * sizeof(double));
    s.alpha = calloc(N, sizeof(double));
    s.kx = calloc(N, sizeof(double));
    s.ky = calloc(N, sizeof(double));
    s.zones = calloc(N, 1);
    if (!s.T || !s.Tnext || !s.alpha || !s.kx || !s.ky || !s.zones) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        return 1;
    }

    double Ra, Rb;
    build_ring(&s, c, &Ra, &Rb);
    compute_conductances(&s);  // s'arrête si l'air touche l'extérieur

    const double dT = SOURCE_TEMP - AMBIENT_TEMP;
    const double r_air = log(Ra / R1) / c->a_air, r_mur = log(Rb / Ra) / c->a_mur;
    const double p_exact = 2.0 * M_PI * dT / (r_air + r_mur);

    // Départ sur le profil exact : le régime établi ne dépend pas de la
    // condition initiale, mais on l'atteint bien plus vite.
    const double q = p_exact / (2.0 * M_PI);
    const double T_a = SOURCE_TEMP - q * r_air;
    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            double d = hypot(x - G / 2, y - G / 2), t;
            if (s.zones[i] == ZONE_EXT) t = AMBIENT_TEMP;
            else if (in_source(x, y)) t = SOURCE_TEMP;
            else if (d < Ra) t = SOURCE_TEMP - q * log(fmax(d, R1) / R1) / c->a_air;
            else t = T_a - q * log(fmin(d, Rb) / Ra) / c->a_mur;
            s.T[i] = t;
            s.Tnext[i] = AMBIENT_TEMP;
        }
    }

    double p_prev = 0.0;
    int converged = 0;
    while (s.iteration < MAX_ITER) {
        step_physics(&s);
        if (s.iteration % CHECK_EVERY) continue;
        compute_stats(&s);
        double p = s.stat_power;
        if (fabs(p - p_prev) < TOL_P * p && fabs(s.stat_pext - p) < TOL_BILAN * p) {
            converged = 1;
            break;
        }
        p_prev = p;
    }

    double p = s.stat_power;
    double err = 100.0 * (p - p_exact) / p_exact;
    double bilan = 100.0 * (s.stat_pext - p) / p;
    // Résistance du mur déduite de la simulation, en retranchant celle de l'air
    // (exacte) : rapport à la valeur théorique, 1 si le mur est bien représenté.
    double r_mur_num = 2.0 * M_PI * dT / p - r_air;
    printf("%-13s %d %4d %8.2f %8.2f %11.4f %11.4f %+8.3f %+8.4f %9.3f %8ld%s\n",
           c->nom, S, G, Ra, Rb, p, p_exact, err, bilan,
           c->a_mur == c->a_air ? NAN : r_mur_num / r_mur,
           s.iteration, converged ? "" : "  NON CONVERGÉ");

    free(s.T); free(s.Tnext); free(s.alpha); free(s.kx); free(s.ky); free(s.zones);
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
