#ifndef THERMAL_DIFFUSION_H
#define THERMAL_DIFFUSION_H

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdint.h>

/*
 * Modèle (v3), en unités physiques. Coupe horizontale d'une pièce : les
 * grandeurs sont données par mètre de hauteur (puissances en W/m).
 *
 *   - Air intérieur brassé par la convection : une seule température T_air,
 *     de capacité ρc_air · aire.
 *   - Murs en béton résolus sur la grille : ρc · ∂T/∂t = ∇·(λ∇T).
 *   - Échanges superficiels (ISO 6946) : résistance R_si entre l'air et le
 *     mur, R_se entre le mur et l'extérieur.
 *   - Extérieur : thermostat à T_EXT.
 *   - Apports internes (occupants, équipements) : puissance Q_INT injectée
 *     dans l'air. En régime établi, elle ressort entièrement par les murs.
 *
 * Grandeur étudiée : la conductance de l'enveloppe G = P / (T_air − T_ext)
 * en W/(m·K) (par mètre de hauteur). Plus G est grande, plus la pièce évacue
 * de chaleur pour un même écart de température.
 */

// Grille. GRID_SIZE, DX, DT et WALL_THICKNESS peuvent être redéfinis à la
// compilation (-D...) : c'est ce que fait le test de l'anneau (make test).
#ifndef GRID_SIZE
#define GRID_SIZE 300
#endif
#ifndef DX
#define DX 0.035                // m : côté d'un pixel (domaine de 10,5 m × 10,5 m)
#endif

enum { ZONE_EXT = 0, ZONE_WALL = 1, ZONE_INT = 2 };

// Matériaux (valeurs usuelles)
#define LAMBDA_WALL 1.75        // W/(m·K) : béton
#define RHOC_WALL   2.3e6       // J/(m³·K) : béton (2 300 kg/m³ × 1 000 J/(kg·K))
#define RHOC_AIR    1206.0      // J/(m³·K) : air (1,2 kg/m³ × 1 005 J/(kg·K))
#ifndef WALL_THICKNESS
#define WALL_THICKNESS 0.20     // m : 4 pixels
#endif

// Résistances superficielles, flux horizontal (ISO 6946), en m²·K/W
#define R_SI 0.13
#define R_SE 0.04

// Conditions
#define T_EXT  20.0             // °C
#define T_INIT 20.0             // °C : air et murs au départ
#define Q_INT  800.0            // W/m : soit 2 kW pour 2,5 m de haut (≈ 55 W/m²)

// Même aire intérieure pour toutes les géométries
#define ROOM_AREA 36.0          // m² (carré de 6 m de côté)

// Pas de temps. Stabilité du schéma explicite : DT · Σ conductances ≤ capacité,
// pour chaque case de mur et pour le nœud d'air (vérifié au démarrage).
#ifndef DT
#define DT 30.0                 // s
#endif

// Statistiques (et ligne CSV) toutes les STATS_EVERY itérations
#define STATS_EVERY 120         // 1 h simulée

// Échelles de couleurs de la fenêtre
#define T_DISPLAY_MAX 30.0      // °C (le bas de l'échelle est T_EXT)
#define FLUX_DISPLAY_MAX 60.0   // W/m²

// Profondeur de Koch par défaut (« fractal ») et profondeur maximale
#define KOCH_DEPTH_REF 2
#define KOCH_DEPTH_MAX 4
#define POLY_MAX (3 * 256)      // 3 · 4^KOCH_DEPTH_MAX sommets

typedef struct { int x, y; } Point;
typedef struct { double x, y; } Vec2;

typedef struct {
    double *T, *Tnext;      // température de chaque case (°C) ; air : T_air
    double *lambda;         // conductivité (W/(m·K)) des cases de mur
    double *cap;            // capacité de chaque case de mur : ρc·Δx² (J/(m·K))
    // Conductances des faces (W/(m·K)) : kx[i] entre i et i+1, ky[i] entre
    // i et i+GRID_SIZE
    double *kx, *ky;
    // Pour chaque case de mur : somme des conductances vers l'air, vers
    // l'extérieur (bilans de puissance)
    double *g_air, *g_ext;
    uint8_t *zones;
    uint32_t *pixels;       // image de la grille (fenêtre)

    // Contour intérieur de la pièce (m)
    Vec2 *poly;
    int n_poly;
    double area, perimeter; // du polygone (m², m)

    // Air brassé
    double T_air, C_air;
    int air_fixed;          // 1 : T_air imposée (test) ; 0 : chauffée par q_int
    double q_int;           // W/m
    double r_si, r_se;      // m²·K/W
    int raw_surfaces;       // 1 : sans corrections de surface (test uniquement)

    long iteration;
    double sim_time;        // s

    int mode;               // 0 : carré, 1 : flocon de Koch
    int depth;              // profondeur de Koch (mode 1)

    FILE *csv;
    int recording;

    // Statistiques
    double stat_T_wall;     // température moyenne des murs
    double stat_p_mur;      // puissance de l'air vers les murs (W/m)
    double stat_p_ext;      // puissance des murs vers l'extérieur (W/m)

    // Géométrie courante
    int n_int, n_wall;      // pixels d'air et de mur
    int n_faces_int_wall;   // faces air/mur
    int n_faces_wall_ext;   // faces mur/extérieur
} Sim;

// Géométrie
void init_geometry(Sim *s);
void geometry_from_polygon(Sim *s);
void surface_geometry(const Sim *s, double x, double y, int outer, double *f, double *delta);

// Physique
void compute_conductances(Sim *s);
void step_physics(Sim *s);
void check_stability(const Sim *s);
void compute_stats(Sim *s);

// Entrées-sorties (io.c)
int sim_alloc(Sim *s);
void sim_free(Sim *s);
void sim_reset(Sim *s);
void geom_name(const Sim *s, char *buf, size_t n);
void csv_header(FILE *f, const Sim *s);
void csv_row(FILE *f, const Sim *s);

// Rendu (render.c)
enum { VIEW_TEMP = 0, VIEW_FLUX = 1 };
double wall_flux(const Sim *s, int i);
SDL_Color temp_color(double T);
SDL_Color flux_color(double q);
void fill_pixels(Sim *s, int view);

// Fenêtre interactive (ui.c)
int run_gui(void);
int capture_gui(const char *path, int depth, int view, double hours);

#endif
