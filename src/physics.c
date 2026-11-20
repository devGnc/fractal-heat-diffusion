#include "thermal_diffusion.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Schéma conservatif (volumes finis) pour ∂T/∂t = ∇·(α∇T) :
 *
 *   T_i ← T_i + DT · Σ_faces k_f (T_j − T_i) − DT · COOLING_FACTOR (T_i − T_ext)
 *
 * k_f est la conductance de la face entre i et j. Elle vaut la moyenne
 * harmonique 2·α_i·α_j / (α_i + α_j) : deux demi-cases en série, comme des
 * résistances thermiques. Le flux de i vers j est l'opposé du flux de j vers i,
 * donc l'énergie est conservée (contrairement à α_i·ΔT, où α_i est pris
 * au centre de la case).
 *
 * L'extérieur est un thermostat parfait à AMBIENT_TEMP : la face entre un
 * mur et l'extérieur n'a que la demi-case de mur en résistance, soit
 * k_f = 2·α_mur.
 */
static float face_conductance(const Sim *s, int i, int j) {
    int ext_i = s->zones[i] == ZONE_EXT, ext_j = s->zones[j] == ZONE_EXT;
    if (ext_i && ext_j) return 0.0f;
    if (ext_i) return 2.0f * s->alpha[j];
    if (ext_j) return 2.0f * s->alpha[i];
    float a = s->alpha[i], b = s->alpha[j];
    return (a + b > 0.0f) ? 2.0f * a * b / (a + b) : 0.0f;
}

void compute_conductances(Sim *s) {
    s->n_int = s->n_wall = s->n_faces_int_wall = s->n_faces_wall_ext = 0;
    
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            int i = y * GRID_SIZE + x;
            s->kx[i] = (x + 1 < GRID_SIZE) ? face_conductance(s, i, i + 1) : 0.0f;
            s->ky[i] = (y + 1 < GRID_SIZE) ? face_conductance(s, i, i + GRID_SIZE) : 0.0f;
            
            if (s->zones[i] == ZONE_INT) s->n_int++;
            if (s->zones[i] == ZONE_WALL) s->n_wall++;
            
            // Comptage des faces (chaque face vue une fois : voisin est et sud)
            int nb[2] = { x + 1 < GRID_SIZE ? i + 1 : -1, y + 1 < GRID_SIZE ? i + GRID_SIZE : -1 };
            for (int k = 0; k < 2; k++) {
                if (nb[k] < 0) continue;
                uint8_t a = s->zones[i], b = s->zones[nb[k]];
                if ((a == ZONE_INT && b == ZONE_WALL) || (a == ZONE_WALL && b == ZONE_INT))
                    s->n_faces_int_wall++;
                if ((a == ZONE_WALL && b == ZONE_EXT) || (a == ZONE_EXT && b == ZONE_WALL))
                    s->n_faces_wall_ext++;
                // L'air ne doit jamais toucher directement l'extérieur
                if ((a == ZONE_INT && b == ZONE_EXT) || (a == ZONE_EXT && b == ZONE_INT)) {
                    fprintf(stderr, "Erreur : fuite air/extérieur en (%d, %d)\n", x, y);
                    exit(1);
                }
            }
        }
    }
}

void check_stability(Sim *s) {
    // Le schéma explicite reste stable (et respecte le principe du maximum)
    // si DT · Σ k_f ≤ 1 pour chaque case.
    float max_sum = 0.0f;
    for (int y = 1; y < GRID_SIZE - 1; y++) {
        for (int x = 1; x < GRID_SIZE - 1; x++) {
            int i = y * GRID_SIZE + x;
            if (s->zones[i] == ZONE_EXT) continue;
            float sum = s->kx[i] + s->kx[i - 1] + s->ky[i] + s->ky[i - GRID_SIZE];
            if (sum > max_sum) max_sum = sum;
        }
    }
    
    float cfl = DT * max_sum;
    printf("\n============================================\n");
    printf("VÉRIFICATION DE LA STABILITÉ NUMÉRIQUE\n");
    printf("============================================\n");
    printf("max Σk_f     = %.5f\n", max_sum);
    printf("ΔT           = %.5f\n", DT);
    printf("Δx           = 1 (pixel)\n");
    printf("CFL = ΔT·Σk_f = %.5f\n", cfl);
    printf("Condition    : CFL ≤ 1\n");
    
    if (cfl > 1.0f) {
        printf("❌ ERREUR : Condition CFL violée (%.5f > 1)\n", cfl);
        printf("   → Réduire ΔT ou α_max\n");
        printf("============================================\n");
        exit(1);
    } else {
        printf("✅ Condition CFL respectée (marge = ×%.1f)\n", 1.0f / cfl);
        printf("============================================\n\n");
    }
}

void step_physics(Sim *s) {
    const float *T = s->T, *kx = s->kx, *ky = s->ky;
    const int G = GRID_SIZE;
    
    for (int y = 1; y < G - 1; y++) {
        for (int x = 1; x < G - 1; x++) {
            int i = y * G + x;
            
            if (s->zones[i] == ZONE_EXT) {
                s->Tnext[i] = AMBIENT_TEMP;
                continue;
            }
            
            float T_C = T[i];
            float flux = kx[i]     * (T[i + 1] - T_C)
                       + kx[i - 1] * (T[i - 1] - T_C)
                       + ky[i]     * (T[i + G] - T_C)
                       + ky[i - G] * (T[i - G] - T_C);
            
            float cooling = COOLING_FACTOR * (T_C - AMBIENT_TEMP);
            
            s->Tnext[i] = T_C + DT * (flux - cooling);
        }
    }
    
    // Remise de la source à SOURCE_TEMP. L'énergie réinjectée pendant ce pas
    // (capacité thermique 1 par pixel) donne la puissance du radiateur.
    int r = SOURCE_RADIUS, c = GRID_SIZE / 2;
    double injected = 0.0;
    for (int y = c - r; y <= c + r; y++) {
        for (int x = c - r; x <= c + r; x++) {
            if (in_source(x, y)) {
                int i = y * GRID_SIZE + x;
                injected += SOURCE_TEMP - s->Tnext[i];
                s->Tnext[i] = SOURCE_TEMP;
            }
        }
    }
    s->stat_power = (float)(injected / DT);
    
    float *tmp = s->T;
    s->T = s->Tnext;
    s->Tnext = tmp;
    
    s->iteration++;
    s->sim_time += DT;
}

// Statistiques sur l'air intérieur uniquement : les murs (plus nombreux pour la
// fractale) et la source (fixée à SOURCE_TEMP) fausseraient la comparaison.
void compute_stats(Sim *s) {
    double sum = 0.0, sum_sq = 0.0;
    float min_temp = 1e6f;
    int count = 0;
    
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            int i = y * GRID_SIZE + x;
            if (s->zones[i] != ZONE_INT || in_source(x, y)) continue;
            double t = s->T[i];
            sum += t;
            sum_sq += t * t;
            if (t < min_temp) min_temp = (float)t;
            count++;
        }
    }
    
    if (count > 0) {
        double avg = sum / count;
        double variance = sum_sq / count - avg * avg;
        s->stat_avg = (float)avg;
        s->stat_min = min_temp;
        s->stat_std = (float)(variance > 0 ? sqrt(variance) : 0.0);
    }
    
    // Bilan de puissance : sortie par les faces mur/extérieur et pertes de Newton
    double p_ext = 0.0, p_cool = 0.0;
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            int i = y * GRID_SIZE + x;
            uint8_t z = s->zones[i];
            if (z == ZONE_EXT) continue;
            p_cool += COOLING_FACTOR * (s->T[i] - AMBIENT_TEMP);
            if (z != ZONE_WALL) continue;
            int nb[4] = { i + 1, i - 1, i + GRID_SIZE, i - GRID_SIZE };
            float k[4] = { s->kx[i], s->kx[i - 1], s->ky[i], s->ky[i - GRID_SIZE] };
            for (int n = 0; n < 4; n++) {
                if (s->zones[nb[n]] == ZONE_EXT)
                    p_ext += k[n] * (s->T[i] - AMBIENT_TEMP);
            }
        }
    }
    s->stat_pext = (float)p_ext;
    s->stat_pcool = (float)p_cool;
}
