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
    int N = GRID_SIZE * GRID_SIZE;
    
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
    
    int cx = GRID_SIZE / 2, cy = GRID_SIZE / 2, rsrc = 15;
    for (int dy = -rsrc; dy <= rsrc; dy++) {
        for (int dx = -rsrc; dx <= rsrc; dx++) {
            int x = cx + dx, y = cy + dy;
            if (x >= 0 && x < GRID_SIZE && y >= 0 && y < GRID_SIZE &&
                dx * dx + dy * dy <= rsrc * rsrc) {
                s->Tnext[y * GRID_SIZE + x] = SOURCE_TEMP;
            }
        }
    }
    
    float *tmp = s->T;
    s->T = s->Tnext;
    s->Tnext = tmp;
    
    s->iteration++;
    s->sim_time += DT;
}

void compute_stats(Sim *s) {
    double sum = 0.0, sum_sq = 0.0;
    float max_temp = -1e6;
    int count = 0;
    
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        if (s->zones[i] == ZONE_INT || s->zones[i] == ZONE_WALL) {
            float t = s->T[i];
            sum += t;
            sum_sq += t * t;
            if (t > max_temp) max_temp = t;
            count++;
        }
    }
    
    if (count > 0) {
        s->stat_avg = sum / count;
        s->stat_max = max_temp;
        double variance = (sum_sq / count) - (s->stat_avg * s->stat_avg);
        s->stat_std = (variance > 0) ? sqrt(variance) : 0.0;
    }
}