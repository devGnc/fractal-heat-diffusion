#include "thermal_diffusion.h"
#include <stdlib.h>
#include <string.h>

#define CSV_HEADER "Iter,Time,Mean,Std,Min,P_source\n"

static int sim_alloc(Sim *s) {
    int N = GRID_SIZE * GRID_SIZE;
    s->T = calloc(N, sizeof(float));
    s->Tnext = calloc(N, sizeof(float));
    s->alpha = calloc(N, sizeof(float));
    s->zones = calloc(N, 1);
    s->pixels = malloc(N * sizeof(uint32_t));
    if (!s->T || !s->Tnext || !s->alpha || !s->zones || !s->pixels) return 0;
    for (int i = 0; i < N; i++) s->Tnext[i] = AMBIENT_TEMP;
    return 1;
}

static void sim_free(Sim *s) {
    free(s->T);
    free(s->Tnext);
    free(s->alpha);
    free(s->zones);
    free(s->pixels);
}

// (Ré)initialise la géométrie du mode courant et remet le temps à zéro.
static void sim_reset(Sim *s) {
    s->iteration = 0;
    s->sim_time = 0.0;
    s->stat_power = 0.0f;
    init_geometry(s);
    compute_stats(s);
}

// Construit d'abord la fractale (aire de référence du carré iso-aire),
// puis la géométrie demandée.
static void sim_setup(Sim *s, int mode) {
    if (mode == 0) {
        s->mode = 1;
        init_geometry(s);
    }
    s->mode = mode;
    sim_reset(s);
}

static void csv_row(FILE *f, const Sim *s) {
    fprintf(f, "%ld,%.4f,%.4f,%.4f,%.4f,%.4f\n",
            s->iteration, s->sim_time, s->stat_avg, s->stat_std,
            s->stat_min, s->stat_power);
}

static const char *mode_name(int mode) {
    return mode ? "fractal" : "carre";
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

static int run_headless(int mode, long n_iter, const char *out) {
    Sim s = {0};
    if (!sim_alloc(&s)) {
        fprintf(stderr, "Allocation mémoire impossible\n");
        return 1;
    }
    sim_setup(&s, mode);
    check_stability(&s);

    FILE *f = fopen(out, "w");
    if (!f) {
        perror(out);
        sim_free(&s);
        return 1;
    }
    fprintf(f, CSV_HEADER);
    csv_row(f, &s);

    printf("Simulation %s : %ld itérations → %s\n", mode_name(mode), n_iter, out);
    long tick = n_iter / 10 > 0 ? n_iter / 10 : 1;
    while (s.iteration < n_iter) {
        step_physics(&s);
        if (s.iteration % STATS_EVERY == 0 || s.iteration == n_iter) {
            compute_stats(&s);
            csv_row(f, &s);
        }
        if (s.iteration % tick == 0) {
            fprintf(stderr, "  %3ld %%  t = %8.1f  T_moy = %6.2f °C  P = %8.2f\n",
                    100 * s.iteration / n_iter, s.sim_time, s.stat_avg, s.stat_power);
        }
    }
    fclose(f);

    printf("Final : T_moy = %.2f °C, écart-type = %.2f, T_min = %.2f °C, P_source = %.2f\n",
           s.stat_avg, s.stat_std, s.stat_min, s.stat_power);
    sim_free(&s);
    return 0;
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
        "Diffusion Thermique - Pièce d'Air - TIPE MPI",
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

    sim_setup(&s, 0);
    check_stability(&s);

    printf("COMMANDES :\n");
    printf("  ESPACE  - Changer de mode (carré ↔ fractale)\n");
    printf("  S       - Démarrer/arrêter enregistrement CSV\n");
    printf("  R       - Réinitialiser la simulation\n");
    printf("  Q       - Quitter\n");
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
                if (k == SDLK_SPACE || k == SDLK_r) {
                    // Un CSV ne doit contenir qu'une seule simulation
                    stop_recording(&s);
                    if (k == SDLK_SPACE) s.mode = !s.mode;
                    sim_reset(&s);
                }
                else if (k == SDLK_s) {
                    if (s.recording) {
                        stop_recording(&s);
                    } else {
                        char name[64];
                        snprintf(name, sizeof name, "res_%s.csv", mode_name(s.mode));
                        s.csv = fopen(name, "w");
                        if (s.csv) {
                            fprintf(s.csv, CSV_HEADER);
                            s.recording = 1;
                            printf("Enregistrement vers %s (à partir de l'itération %ld)\n",
                                   name, s.iteration);
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

        char title[256];
        snprintf(title, sizeof title,
                 "Mode:%s | Iter:%ld | T_moy:%.2f°C | P:%.1f %s",
                 s.mode ? "FRACTALE" : "CARRÉ",
                 s.iteration, s.stat_avg, s.stat_power,
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
            "  %s                                          fenêtre interactive\n"
            "  %s --headless <carre|fractal> [iter] [csv]  simulation sans fenêtre\n"
            "      iter : nombre d'itérations (défaut 300000)\n"
            "      csv  : fichier de sortie (défaut res_<mode>.csv)\n",
            prog, prog);
}

int main(int argc, char **argv) {
    if (argc == 1) return run_gui();

    if (strcmp(argv[1], "--headless") != 0 || argc < 3 || argc > 5) {
        usage(argv[0]);
        return 1;
    }

    int mode;
    if (strcmp(argv[2], "carre") == 0) mode = 0;
    else if (strcmp(argv[2], "fractal") == 0) mode = 1;
    else {
        usage(argv[0]);
        return 1;
    }

    long n_iter = 300000;
    if (argc >= 4) {
        char *end;
        n_iter = strtol(argv[3], &end, 10);
        if (*end != '\0' || n_iter <= 0) {
            fprintf(stderr, "Nombre d'itérations invalide : %s\n", argv[3]);
            return 1;
        }
    }

    char default_out[64];
    snprintf(default_out, sizeof default_out, "res_%s.csv", mode_name(mode));
    const char *out = argc == 5 ? argv[4] : default_out;

    return run_headless(mode, n_iter, out);
}
