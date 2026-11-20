#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * Volumes finis, par mètre de hauteur. Chaque case de mur a une capacité
 * C_i = ρc·Δx² (J/(m·K)) ; chaque face, de longueur Δx, une conductance G_f
 * (W/(m·K)) égale à l'inverse de la somme des résistances traversées :
 *
 *   mur/mur : deux demi-cases en série   G = Δx / (Δx/2λ_i + Δx/2λ_j)
 *   air/mur : surface + demi-case        G = Δx / (R_si + Δx/2λ)
 *   mur/ext : demi-case + surface        G = Δx / (Δx/2λ + R_se)
 *
 * Corrections de surface (voir surface_geometry dans geometry.c) : la face
 * reçoit en plus la résistance delta/λ de l'épaisseur de mur que les cases
 * entières ne représentent pas, et ses résistances de surface sont multipliées
 * par le facteur d'escalier f = |n_x| + |n_y| :
 *   G = Δx / ((R_s + delta/λ)·f + Δx/2λ)
 * La conductance totale d'une surface en marches est alors celle de sa vraie
 * longueur, à sa vraie position.
 *
 * Schéma explicite (Euler) :
 *   C_i (T_i' − T_i) / Δt = Σ_faces G_f (T_j − T_i)
 *   C_air (T_air' − T_air) / Δt = Q_int − Σ_faces air/mur G_f (T_air − T_j)
 *
 * Le flux d'une case vers sa voisine est l'opposé du flux inverse : l'énergie
 * C_air·T_air + Σ C_i·T_i ne varie que par Q_int et par les pertes vers
 * l'extérieur.
 */
static double face_conductance(const Sim *s, int i, int j) {
    uint8_t zi = s->zones[i], zj = s->zones[j];
    if (zi > zj) {  // ordre : EXT < WALL < INT
        int t = i; i = j; j = t;
        uint8_t u = zi; zi = zj; zj = u;
    }
    if (zi == zj && zi != ZONE_WALL) return 0.0;  // ext/ext, air/air (un seul nœud)
    if (zi == ZONE_EXT && zj == ZONE_INT) return 0.0;  // interdit, détecté plus bas
    if (zi == ZONE_EXT || zj == ZONE_INT) {
        // Surface : milieu de la face (m) et facteur d'escalier
        const int G = GRID_SIZE;
        double xm = ((i % G + j % G) / 2.0 + 0.5) * DX;
        double ym = ((i / G + j / G) / 2.0 + 0.5) * DX;
        int outer = zi == ZONE_EXT;
        int w = outer ? j : i;  // case de mur
        double f = 1.0, delta = 0.0, lam = s->lambda[w];
        if (!s->raw_surfaces) surface_geometry(s, xm, ym, outer, &f, &delta);
        double r = (outer ? s->r_se : s->r_si) + delta / lam;
        // Garde-fou : la demi-case ne peut pas être annulée par un excès de mur
        double denom = fmax(r * f + DX / (2.0 * lam), DX / (8.0 * lam));
        return DX / denom;
    }
    return DX / (DX / (2.0 * s->lambda[i]) + DX / (2.0 * s->lambda[j]));
}

void compute_conductances(Sim *s) {
    const int G = GRID_SIZE;
    s->n_int = s->n_wall = s->n_faces_int_wall = s->n_faces_wall_ext = 0;

    for (int i = 0; i < G * G; i++) s->g_air[i] = s->g_ext[i] = 0.0;

    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            s->kx[i] = (x + 1 < G) ? face_conductance(s, i, i + 1) : 0.0;
            s->ky[i] = (y + 1 < G) ? face_conductance(s, i, i + G) : 0.0;

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
                if (a == ZONE_WALL && b == ZONE_INT) { s->g_air[i] += k[m]; s->n_faces_int_wall++; }
                if (a == ZONE_INT && b == ZONE_WALL) { s->g_air[j] += k[m]; s->n_faces_int_wall++; }
                if (a == ZONE_WALL && b == ZONE_EXT) { s->g_ext[i] += k[m]; s->n_faces_wall_ext++; }
                if (a == ZONE_EXT && b == ZONE_WALL) { s->g_ext[j] += k[m]; s->n_faces_wall_ext++; }
            }
        }
    }
}

void check_stability(const Sim *s) {
    const int G = GRID_SIZE;
    // Case de mur : DT · Σ G_f / C_i ≤ 1. Nœud d'air : DT · Σ G_air / C_air ≤ 1.
    double worst = 0.0, g_air_tot = 0.0;
    for (int y = 1; y < G - 1; y++) {
        for (int x = 1; x < G - 1; x++) {
            int i = y * G + x;
            if (s->zones[i] != ZONE_WALL) continue;
            double sum = s->kx[i] + s->kx[i - 1] + s->ky[i] + s->ky[i - G];
            worst = fmax(worst, DT * sum / s->cap[i]);
            g_air_tot += s->g_air[i];
        }
    }
    double air = s->C_air > 0.0 ? DT * g_air_tot / s->C_air : 0.0;

    printf("Stabilité (doit rester ≤ 1) : murs %.3f, air %.3f (Δt = %g s, Δx = %g m)\n",
           worst, air, (double)DT, (double)DX);
    if (worst > 1.0 || (!s->air_fixed && air > 1.0)) {
        fprintf(stderr, "Erreur : schéma instable, réduire DT\n");
        exit(1);
    }
}

void step_physics(Sim *s) {
    const double *T = s->T, *kx = s->kx, *ky = s->ky;
    const int G = GRID_SIZE;
    double p_mur = 0.0;  // puissance de l'air vers les murs

    for (int y = 1; y < G - 1; y++) {
        for (int x = 1; x < G - 1; x++) {
            int i = y * G + x;
            if (s->zones[i] != ZONE_WALL) continue;
            // Les cases d'air contiennent T_air : le flux venant de l'air est
            // compris dans cette somme.
            double Ti = T[i];
            double flux = kx[i]     * (T[i + 1] - Ti)
                        + kx[i - 1] * (T[i - 1] - Ti)
                        + ky[i]     * (T[i + G] - Ti)
                        + ky[i - G] * (T[i - G] - Ti);
            s->Tnext[i] = Ti + DT * flux / s->cap[i];
            p_mur += s->g_air[i] * (s->T_air - Ti);
        }
    }

    if (!s->air_fixed) s->T_air += DT * (s->q_int - p_mur) / s->C_air;
    for (int i = 0; i < G * G; i++) {
        if (s->zones[i] == ZONE_INT) s->Tnext[i] = s->T_air;
    }

    double *tmp = s->T;
    s->T = s->Tnext;
    s->Tnext = tmp;

    s->iteration++;
    s->sim_time += DT;
}

void compute_stats(Sim *s) {
    double sum = 0.0, p_mur = 0.0, p_ext = 0.0;
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        if (s->zones[i] != ZONE_WALL) continue;
        sum += s->T[i];
        p_mur += s->g_air[i] * (s->T_air - s->T[i]);
        p_ext += s->g_ext[i] * (s->T[i] - T_EXT);
    }
    s->stat_T_wall = s->n_wall ? sum / s->n_wall : 0.0;
    s->stat_p_mur = p_mur;
    s->stat_p_ext = p_ext;
}
