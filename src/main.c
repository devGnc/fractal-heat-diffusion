#include "thermal_diffusion.h"
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL Init fail\n");
        return 1;
    }
    
    SDL_Window *win = SDL_CreateWindow(
        "Diffusion Thermique - Pièce d'Air - TIPE MPI",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WIN_SIZE, WIN_SIZE, 0
    );
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED);
    
    Sim s = {0};
    int N = GRID_SIZE * GRID_SIZE;
    s.T = calloc(N, sizeof(float));
    s.Tnext = calloc(N, sizeof(float));
    s.alpha = calloc(N, sizeof(float));
    s.zones = calloc(N, 1);
    s.pixels = malloc(N * sizeof(uint32_t));
    s.tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
                              SDL_TEXTUREACCESS_STREAMING, GRID_SIZE, GRID_SIZE);
    
    printf("\n");
    printf("============================================\n");
    printf("  TIPE - Diffusion thermique\n");
    printf("  Géométries fractales vs. classiques\n");
    printf("============================================\n\n");
    
    printf("MODÈLE PHYSIQUE :\n");
    printf("  - Air intérieur : α = %.2f (diffusion rapide)\n", DIFF_AIR);
    printf("  - Murs isolants : α = %.2f (diffusion lente)\n", DIFF_WALL);
    printf("  - Ratio α_air/α_mur = %.1f\n", DIFF_AIR / DIFF_WALL);
    printf("  - Extérieur : T = %.0f°C (Dirichlet)\n", AMBIENT_TEMP);
    printf("  - Source : T = %.0f°C (radiateur)\n\n", SOURCE_TEMP);
    
    printf("[1] Création de la fractale (référence)...\n");
    s.mode = 1;
    init_geometry(&s);
    
    printf("\n[2] Création du carré iso-aire...\n");
    s.mode = 0;
    s.iteration = 0;
    s.sim_time = 0.0;
    init_geometry(&s);
    
    check_stability(&s);
    
    printf("COMMANDES :\n");
    printf("  ESPACE  - Changer de mode (carré ↔ fractale)\n");
    printf("  S       - Démarrer/arrêter enregistrement CSV\n");
    printf("  R       - Réinitialiser la simulation\n");
    printf("  Q       - Quitter\n");
    printf("============================================\n\n");
    
    int running = 1;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) {
                running = 0;
            }
            if (e.type == SDL_KEYDOWN) {
                if (e.key.keysym.sym == SDLK_SPACE) {
                    s.mode = !s.mode;
                    s.iteration = 0;
                    s.sim_time = 0;
                    init_geometry(&s);
                }
                else if (e.key.keysym.sym == SDLK_r) {
                    s.iteration = 0;
                    s.sim_time = 0;
                    init_geometry(&s);
                }
                else if (e.key.keysym.sym == SDLK_s) {
                    if (s.recording) {
                        fclose(s.csv);
                        s.recording = 0;
                        printf("Enregistrement CSV arrêté.\n");
                    } else {
                        char name[64];
                        sprintf(name, "res_%s.csv", s.mode ? "fractal" : "carre");
                        s.csv = fopen(name, "w");
                        if (s.csv) {
                            fprintf(s.csv, "Iter,Time,Mean,Std,Max\n");
                            s.recording = 1;
                            printf("Enregistrement vers %s\n", name);
                        }
                    }
                }
                else if (e.key.keysym.sym == SDLK_q) {
                    running = 0;
                }
            }
        }
        
        for (int k = 0; k < ITER_PER_FRAME; k++) {
            step_physics(&s);
        }
        
        if (s.iteration % 200 == 0) {
            compute_stats(&s);
        }
        
        if (s.recording && s.iteration % 200 == 0) {
            fprintf(s.csv, "%ld,%.4f,%.4f,%.4f,%.4f\n",
                    s.iteration, s.sim_time, s.stat_avg, s.stat_std, s.stat_max);
            fflush(s.csv);
        }
        
        SDL_RenderClear(ren);
        render(ren, &s);
        SDL_RenderPresent(ren);
        
        char title[256];
        sprintf(title, "Mode:%s | Iter:%ld | T_moy:%.2f°C | T_max:%.2f°C %s",
                s.mode ? "FRACTALE" : "CARRÉ",
                s.iteration, s.stat_avg, s.stat_max,
                s.recording ? "[REC]" : "");
        SDL_SetWindowTitle(win, title);
        
        SDL_Delay(10);
    }
    
    if (s.recording && s.csv) fclose(s.csv);
    free(s.T);
    free(s.Tnext);
    free(s.alpha);
    free(s.zones);
    free(s.pixels);
    SDL_DestroyTexture(s.tex);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    
    return 0;
}