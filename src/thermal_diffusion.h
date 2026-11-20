#ifndef THERMAL_DIFFUSION_H
#define THERMAL_DIFFUSION_H

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdint.h>

#define GRID_SIZE 300
#define WIN_SIZE 900

enum { ZONE_EXT = 0, ZONE_WALL = 1, ZONE_INT = 2 };

#define SOURCE_TEMP 100.0f
#define AMBIENT_TEMP 20.0f
#define WALL_TEMP_INIT 20.0f

#define DIFF_AIR 0.20f
#define DIFF_WALL 0.02f

#define COOLING_FACTOR 0.00002f

#define DT 0.05f
#define ITER_PER_FRAME 80

// Statistiques (et ligne CSV) toutes les STATS_EVERY itérations
#define STATS_EVERY 400
_Static_assert(STATS_EVERY % ITER_PER_FRAME == 0,
               "STATS_EVERY doit être un multiple de ITER_PER_FRAME");

// Source chauffante : disque centré, maintenu à SOURCE_TEMP
#define SOURCE_RADIUS 15

static inline int in_source(int x, int y) {
    int dx = x - GRID_SIZE / 2, dy = y - GRID_SIZE / 2;
    return dx * dx + dy * dy <= SOURCE_RADIUS * SOURCE_RADIUS;
}

typedef struct { int x, y; } Point;

typedef struct {
    float *T, *Tnext;
    float *alpha;
    uint8_t *zones;
    uint32_t *pixels;
    SDL_Texture *tex;
    
    long iteration;
    double sim_time;
    
    int mode;
    
    FILE *csv;
    int recording;
    
    // Statistiques sur l'air intérieur hors source (murs exclus)
    float stat_avg, stat_std, stat_min;
    // Puissance fournie par la source au dernier pas (°C·px / unité de temps)
    float stat_power;
    
    int count_int_square;
    int count_int_fractal;
} Sim;

// Fonctions de géométrie
void init_geometry(Sim *s);
void draw_line(uint8_t *zones, Point p1, Point p2, uint8_t val);
void koch(uint8_t *zones, Point a, Point b, int depth);
void flood_fill_safe(uint8_t *zones, int sx, int sy);
int count_interior_pixels(uint8_t *zones);

// Fonctions de physique
void step_physics(Sim *s);
void check_stability(Sim *s);
void compute_stats(Sim *s);

// Fonctions de rendu
void render(SDL_Renderer *ren, Sim *s);
uint32_t get_color(float T, uint8_t zone);

#endif