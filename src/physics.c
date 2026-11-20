#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Volumes finis, par mètre de hauteur. Chaque case a une capacité
 * C_i = ρc·Δx² (J/(m·K)) ; chaque face, de longueur Δx, une conductance G_f
 * (W/(m·K)) égale à l'inverse de la somme des résistances traversées :
 *
 *   air/air, mur/mur : deux demi-cases      G = Δx / (Δx/2λ_i + Δx/2λ_j)
 *   air/mur : demi-case, surface, demi-case G = Δx / (Δx/2λ_a + R_si + Δx/2λ_m)
 *   mur/ext : demi-case + surface           G = Δx / (Δx/2λ_m + R_se)
 *
 * Corrections de surface (voir surface_geometry dans geometry.c) : la face
 * reçoit en plus la résistance delta/λ de l'épaisseur de mur que les cases
 * entières ne représentent pas, et ses résistances de surface sont multipliées
 * par le facteur d'escalier f = |n_x| + |n_y| :
 *   G = Δx / (Δx/2λ_a + (R_s + delta/λ_m)·f + Δx/2λ_m)
 * La conductance totale d'une surface en marches est alors celle de sa vraie
 * longueur, à sa vraie position.
 *
 * Équation semi-discrète, avec A_x et A_y les échanges par les faces
 * verticales et horizontales, et S la source :
 *   C dT/dt = A_x T + A_y T + S
 */
// f_out : facteur d'escalier de la face (1 hors des surfaces)
static double face_conductance(const Sim *s, int i, int j, double *f_out) {
    uint8_t zi = s->zones[i], zj = s->zones[j];
    *f_out = 1.0;
    if (zi > zj) {  // ordre : EXT < WALL < INT
        int t = i; i = j; j = t;
        uint8_t u = zi; zi = zj; zj = u;
    }
    if (zi == ZONE_EXT && zj != ZONE_WALL) return 0.0;  // ext/ext ; air/ext interdit
    if (zi == zj) {
        double a = s->lambda[i], b = s->lambda[j];
        return DX / (DX / (2.0 * a) + DX / (2.0 * b));
    }
    // Surface : milieu de la face (m) et corrections
    const int G = GRID_SIZE;
    double xm = ((i % G + j % G) / 2.0 + 0.5) * DX;
    double ym = ((i / G + j / G) / 2.0 + 0.5) * DX;
    int outer = zi == ZONE_EXT;
    int w = outer ? j : i;  // case de mur
    double f = 1.0, delta = 0.0, lam = s->lambda[w];
    if (!s->raw_surfaces) surface_geometry(s, xm, ym, outer, &f, &delta);
    *f_out = f;
    double r = (outer ? s->r_se : s->r_si) + delta / lam;
    // Garde-fou : la demi-case ne peut pas être annulée par un excès de mur
    double denom = fmax(r * f + DX / (2.0 * lam), DX / (8.0 * lam));
    if (!outer) denom += DX / (2.0 * s->lambda[j]);  // demi-case d'air
    return DX / denom;
}

void compute_conductances(Sim *s) {
    const int G = GRID_SIZE;
    s->n_int = s->n_wall = s->n_faces_int_wall = s->n_faces_wall_ext = 0;
    s->len_int = s->len_ext = 0.0;
    for (int i = 0; i < G * G; i++) s->g_ext[i] = 0.0;

    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            double f[2] = { 1.0, 1.0 };
            s->kx[i] = (x + 1 < G) ? face_conductance(s, i, i + 1, &f[0]) : 0.0;
            s->ky[i] = (y + 1 < G) ? face_conductance(s, i, i + G, &f[1]) : 0.0;

            if (s->zones[i] == ZONE_INT) s->n_int++;
            if (s->zones[i] == ZONE_WALL) s->n_wall++;

            // Chaque face vue une fois : voisins est et sud
            int nb[2] = { x + 1 < G ? i + 1 : -1, y + 1 < G ? i + G : -1 };
            double k[2] = { s->kx[i], s->ky[i] };
            for (int m = 0; m < 2; m++) {
                if (nb[m] < 0) continue;
                int j = nb[m];
                uint8_t a = s->zones[i], b = s->zones[j];
                if ((a == ZONE_INT && b == ZONE_EXT) || (a == ZONE_EXT && b == ZONE_INT)) {
                    fprintf(stderr, "Erreur : l'air touche l'extérieur en (%d, %d)\n", x, y);
                    exit(1);
                }
                if ((a == ZONE_WALL && b == ZONE_INT) || (a == ZONE_INT && b == ZONE_WALL)) {
                    s->n_faces_int_wall++;
                    s->len_int += DX / f[m];
                }
                if ((a == ZONE_WALL && b == ZONE_EXT) || (a == ZONE_EXT && b == ZONE_WALL)) {
                    s->g_ext[a == ZONE_WALL ? i : j] += k[m];
                    s->n_faces_wall_ext++;
                    s->len_ext += DX / f[m];
                }
            }
        }
    }
}

/*
 * Algorithme de Thomas : résout en O(n) le système tridiagonal
 *   a_k x_{k−1} + b_k x_k + c_k x_{k+1} = d_k   (k = 0 … n−1)
 * par élimination de Gauss sans pivot (la matrice est à diagonale dominante).
 * La solution remplace d ; c est modifié.
 */
static void thomas(int n, const double *a, const double *b, double *c, double *d) {
    c[0] /= b[0];
    d[0] /= b[0];
    for (int k = 1; k < n; k++) {
        double m = b[k] - a[k] * c[k - 1];
        c[k] /= m;
        d[k] = (d[k] - a[k] * d[k - 1]) / m;
    }
    for (int k = n - 2; k >= 0; k--) d[k] -= c[k] * d[k + 1];
}

// Échanges par les faces horizontales (direction y) au temps n : (A_y Tⁿ)_i
static inline double flux_y(const Sim *s, int i) {
    const int G = GRID_SIZE;
    const double *T = s->T;
    return s->ky[i - G] * (T[i - G] - T[i]) + s->ky[i] * (T[i + G] - T[i]);
}

/*
 * Un pas de temps par le schéma de Douglas-Rachford (directions alternées),
 * avec M = C/Δt :
 *   (M − A_x) T* = (M + A_y) Tⁿ + S         lignes : implicite en x
 *   (M − A_y) Tⁿ⁺¹ = M T* − A_y Tⁿ           colonnes : implicite en y
 * Chaque étape est un système tridiagonal par ligne ou par colonne, résolu
 * par l'algorithme de Thomas : un pas coûte O(N) pour N cases. Le schéma est
 * stable pour tout Δt, et son régime établi est exactement celui de
 * l'équation discrète : avec Tⁿ⁺¹ = Tⁿ = T, la seconde étape donne T* = T,
 * puis la première (A_x + A_y) T + S = 0.
 * Les cases extérieures (thermostat) sont des lignes « T = T_EXT ».
 */
void step_physics(Sim *s) {
    const int G = GRID_SIZE;
    double *a = s->work, *b = a + G, *c = b + G, *d = c + G;
    const double dt = s->dt;

    // Étape 1 : chaque ligne, implicite en x
    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            if (s->zones[i] == ZONE_EXT) {
                a[x] = c[x] = 0.0; b[x] = 1.0; d[x] = s->T[i];
                continue;
            }
            double M = s->cap[i] / dt, kl = s->kx[i - 1], kr = s->kx[i];
            a[x] = -kl;
            c[x] = -kr;
            b[x] = M + kl + kr;
            d[x] = M * s->T[i] + flux_y(s, i) + s->src[i];
        }
        thomas(G, a, b, c, d);
        for (int x = 0; x < G; x++) s->Tstar[y * G + x] = d[x];
    }

    // Étape 2 : chaque colonne, implicite en y. Le second membre utilise Tⁿ :
    // on le calcule pour toute la colonne avant de l'écraser par Tⁿ⁺¹.
    for (int x = 0; x < G; x++) {
        for (int y = 0; y < G; y++) {
            int i = y * G + x;
            if (s->zones[i] == ZONE_EXT) {
                a[y] = c[y] = 0.0; b[y] = 1.0; d[y] = s->T[i];
                continue;
            }
            double M = s->cap[i] / dt, kd = s->ky[i - G], ku = s->ky[i];
            a[y] = -kd;
            c[y] = -ku;
            b[y] = M + kd + ku;
            d[y] = M * s->Tstar[i] - flux_y(s, i);
        }
        thomas(G, a, b, c, d);
        for (int y = 0; y < G; y++) s->T[y * G + x] = d[y];
    }

    s->iteration++;
    s->sim_time += dt;
}

/*
 * Régime établi, résolu directement : (A_x + A_y) T + S = 0, sans simuler le
 * transitoire. Les cases extérieures sont fixées à T_EXT ; pour les autres,
 * le système K T = b est symétrique défini positif :
 *   (K v)_i = Σ_faces G_f (v_i − v_j)   (v_j = 0 pour une case extérieure)
 *   b_i     = S_i + Σ_faces extérieures G_f · T_EXT
 * Gradient conjugué préconditionné par la diagonale (Jacobi) : chaque
 * itération coûte O(N). Arrêt quand ‖r‖ ≤ tol · ‖b‖. Renvoie le nombre
 * d'itérations (−1 si max_iter est atteint) ; le résultat est dans s->T.
 */
static void apply_K(const Sim *s, const double *v, double *out) {
    const int G = GRID_SIZE;
    for (int i = 0; i < G * G; i++) {
        if (s->zones[i] == ZONE_EXT) { out[i] = 0.0; continue; }
        out[i] = s->kx[i - 1] * (v[i] - v[i - 1]) + s->kx[i] * (v[i] - v[i + 1])
               + s->ky[i - G] * (v[i] - v[i - G]) + s->ky[i] * (v[i] - v[i + G]);
    }
}

int solve_steady(Sim *s, double tol, int max_iter) {
    const int G = GRID_SIZE, N = G * G;
    double *r = malloc(N * sizeof(double)), *z = malloc(N * sizeof(double));
    double *p = malloc(N * sizeof(double)), *q = malloc(N * sizeof(double));
    double *x = malloc(N * sizeof(double)), *dinv = malloc(N * sizeof(double));
    if (!r || !z || !p || !q || !x || !dinv) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        exit(1);
    }

    // Inconnue x = T − T_EXT : nulle sur les cases extérieures, b_i = S_i
    double bnorm = 0.0;
    for (int i = 0; i < N; i++) {
        x[i] = s->zones[i] == ZONE_EXT ? 0.0 : s->T[i] - T_EXT;
        dinv[i] = s->zones[i] == ZONE_EXT ? 0.0
                : 1.0 / (s->kx[i - 1] + s->kx[i] + s->ky[i - G] + s->ky[i]);
        bnorm += s->src[i] * s->src[i];
    }
    bnorm = sqrt(bnorm);

    apply_K(s, x, q);
    double rz = 0.0;
    for (int i = 0; i < N; i++) {
        r[i] = s->zones[i] == ZONE_EXT ? 0.0 : s->src[i] - q[i];
        z[i] = dinv[i] * r[i];
        p[i] = z[i];
        rz += r[i] * z[i];
    }

    int it = 0, ok = 0;
    for (; it < max_iter; it++) {
        double rr = 0.0;
        for (int i = 0; i < N; i++) rr += r[i] * r[i];
        if (sqrt(rr) <= tol * bnorm) { ok = 1; break; }
        apply_K(s, p, q);
        double pq = 0.0;
        for (int i = 0; i < N; i++) pq += p[i] * q[i];
        double alpha = rz / pq, rz_new = 0.0;
        for (int i = 0; i < N; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * q[i];
            z[i] = dinv[i] * r[i];
            rz_new += r[i] * z[i];
        }
        double beta = rz_new / rz;
        rz = rz_new;
        for (int i = 0; i < N; i++) p[i] = z[i] + beta * p[i];
    }

    for (int i = 0; i < N; i++) s->T[i] = T_EXT + x[i];
    free(r); free(z); free(p); free(q); free(x); free(dinv);
    return ok ? it : -1;
}

void compute_stats(Sim *s) {
    const int G = GRID_SIZE;
    double sum_air = 0.0, max_air = -INFINITY, sum_wall = 0.0, p_mur = 0.0, p_ext = 0.0;
    for (int i = 0; i < G * G; i++) {
        uint8_t z = s->zones[i];
        if (z == ZONE_INT) {
            sum_air += s->T[i];
            if (s->T[i] > max_air) max_air = s->T[i];
        } else if (z == ZONE_WALL) {
            sum_wall += s->T[i];
            p_ext += s->g_ext[i] * (s->T[i] - T_EXT);
            // Faces vers l'air (les bords de grille sont toujours extérieurs)
            int nb[4] = { i - 1, i + 1, i - G, i + G };
            double k[4] = { s->kx[i - 1], s->kx[i], s->ky[i - G], s->ky[i] };
            for (int m = 0; m < 4; m++) {
                if (s->zones[nb[m]] == ZONE_INT) p_mur += k[m] * (s->T[nb[m]] - s->T[i]);
            }
        }
    }
    s->stat_T_air = s->n_int ? sum_air / s->n_int : 0.0;
    s->stat_T_air_max = max_air;
    s->stat_T_wall = s->n_wall ? sum_wall / s->n_wall : 0.0;
    s->stat_p_mur = p_mur;
    s->stat_p_ext = p_ext;
}
