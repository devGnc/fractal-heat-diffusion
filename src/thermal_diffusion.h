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

// Stabilité du schéma explicite : DT · (somme des conductances d'une case) ≤ 1,
// soit DT ≤ 1,25 dans l'air (α = 0,2). DT = 0,5 garde une marge de 2,5.
#define DT 0.5f
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

// Profondeur de Koch de référence : c'est sur elle que le carré est iso-aire,
// et toutes les profondeurs ont la même aire théorique qu'elle.
#define KOCH_DEPTH_REF 2
#define KOCH_DEPTH_MAX 4

typedef struct { int x, y; } Point;
typedef struct { double x, y; } Vec2;

typedef struct {
    float *T, *Tnext;
    float *alpha;
    // Conductances des faces : kx[i] entre i et i+1, ky[i] entre i et i+GRID_SIZE
    float *kx, *ky;
    uint8_t *zones;
    uint32_t *pixels;
    SDL_Texture *tex;
    
    long iteration;
    double sim_time;
    
    int mode;   // 0 : carré, 1 : flocon de Koch
    int depth;  // profondeur de Koch (mode 1)
    
    FILE *csv;
    int recording;
    
    // Statistiques sur l'air intérieur hors source (murs exclus)
    float stat_avg, stat_std, stat_min;
    // Puissance fournie par la source au dernier pas (°C·px / unité de temps)
    float stat_power;
    // Bilan : puissance sortant par les faces mur/extérieur, pertes de Newton.
    // En régime établi : stat_power = stat_pext + stat_pcool.
    float stat_pext, stat_pcool;
    
    int count_int_square;
    int count_int_fractal;  // aire du flocon de profondeur KOCH_DEPTH_REF
    
    // Géométrie courante
    int n_int, n_wall;      // pixels d'air et de mur
    int n_faces_int_wall;   // faces air/mur
    int n_faces_wall_ext;   // faces mur/extérieur
} Sim;

// Fonctions de géométrie
void init_geometry(Sim *s);
void draw_line(uint8_t *zones, Point p1, Point p2, uint8_t val);
void koch(uint8_t *zones, Vec2 a, Vec2 b, int depth);
void flood_fill_safe(uint8_t *zones, int sx, int sy);
int count_interior_pixels(uint8_t *zones);

// Fonctions de physique
void compute_conductances(Sim *s);
void step_physics(Sim *s);
void check_stability(Sim *s);
void compute_stats(Sim *s);

// Fonctions de rendu
void render(SDL_Renderer *ren, Sim *s);
uint32_t get_color(float T, uint8_t zone);

#endif