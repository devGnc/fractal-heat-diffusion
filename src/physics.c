#include "thermal_diffusion.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void check_stability(Sim *s) {
    float max_alpha = 0.0f;
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        if (s->alpha[i] > max_alpha) max_alpha = s->alpha[i];
    }
    
    float cfl = DT * max_alpha;
    printf("\n============================================\n");
    printf("VÉRIFICATION DE LA STABILITÉ NUMÉRIQUE\n");
    printf("============================================\n");
    printf("α_max        = %.5f\n", max_alpha);
    printf("ΔT           = %.5f\n", DT);
    printf("Δx           = 1 (pixel)\n");
    printf("CFL = α·ΔT/Δx² = %.5f\n", cfl);
    printf("Condition    : CFL ≤ 0.25\n");
    
    if (cfl > 0.25f) {
        printf("❌ ERREUR : Condition CFL violée (%.5f > 0.25)\n", cfl);
        printf("   → Réduire ΔT ou α_max\n");
        printf("============================================\n");
        exit(1);
    } else {
        printf("✅ Condition CFL respectée (marge = ×%.1f)\n", 0.25f / cfl);
        printf("============================================\n\n");
    }
}

void step_physics(Sim *s) {
    for (int y = 1; y < GRID_SIZE - 1; y++) {
        for (int x = 1; x < GRID_SIZE - 1; x++) {
            int i = y * GRID_SIZE + x;
            uint8_t zone = s->zones[i];
            
            if (zone == ZONE_EXT) {
                s->Tnext[i] = AMBIENT_TEMP;
                continue;
            }
            
            float T_C = s->T[i];
            
            float T_E = s->T[i + 1];
            float T_W = s->T[i - 1];
            float T_N = s->T[i - GRID_SIZE];
            float T_S = s->T[i + GRID_SIZE];
            
            float laplacian = T_E + T_W + T_N + T_S - 4.0f * T_C;
            
            float cooling = COOLING_FACTOR * (T_C - AMBIENT_TEMP);
            
            float alpha = s->alpha[i];
            
            s->Tnext[i] = T_C + alpha * laplacian * DT - cooling * DT;
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
}
