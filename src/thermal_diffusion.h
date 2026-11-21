#ifndef THERMAL_DIFFUSION_H
#define THERMAL_DIFFUSION_H

#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdint.h>

/*
 * Modèle (v4), en unités physiques. Coupe horizontale d'une pièce : les
 * grandeurs sont données par mètre de hauteur (puissances en W/m).
 *
 *   - Source de chaleur (radiateur, appareil) : disque au centre de la pièce,
 *     qui injecte la puissance Q_INT dans l'air.
 *   - Air intérieur résolu sur la grille. La convection le brasse : on la
 *     représente par une conductivité effective LAMBDA_AIR, bien plus grande
 *     que celle de l'air immobile.
 *   - Murs en béton : ρc · ∂T/∂t = ∇·(λ∇T).
 *   - Échanges superficiels (ISO 6946) : résistance R_si entre l'air et le
 *     mur, R_se entre le mur et l'extérieur.
 *   - Extérieur : thermostat à T_EXT.
 *
 * Grandeur étudiée : la conductance de l'enveloppe G = P / (T_air − T_ext)
 * en W/(m·K), où T_air est la température moyenne de l'air. Plus G est grande,
 * plus la pièce évacue de chaleur pour un même écart de température.
 */

// Grille. GRID_SIZE, DX et WALL_THICKNESS peuvent être redéfinis à la
// compilation (-D...) : c'est ce que font make data et make test.
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
// Air brassé : conductivité effective. Avec α = λ/ρc ≈ 0,08 m²/s, la chaleur
// se répartit dans une pièce de 3 m de rayon en R²/α ≈ 2 minutes, l'ordre de
// grandeur du brassage d'une pièce chauffée. Air immobile : 0,026 W/(m·K).
#define LAMBDA_AIR  100.0       // W/(m·K)
#ifndef WALL_THICKNESS
#define WALL_THICKNESS 0.20     // m : valeur par défaut (option --epaisseur)
#endif

// Résistances superficielles, flux horizontal (ISO 6946), en m²·K/W
#define R_SI 0.13
#define R_SE 0.04

// Conditions
#define T_EXT  20.0             // °C
#define T_INIT 20.0             // °C : air et murs au départ
#define Q_INT  800.0            // W/m : soit 2 kW pour 2,5 m de haut (≈ 55 W/m²)
#define SOURCE_R 0.30           // m : rayon de la source

// Même aire intérieure pour toutes les géométries : valeur par défaut
// (option --aire)
#define ROOM_AREA 36.0          // m² (carré de 6 m de côté)

// Pas de temps par défaut. Le schéma est implicite (directions alternées) :
// stable pour tout pas de temps, la fenêtre le fait varier avec la vitesse.
#define DT 30.0                 // s

// Échelles de couleurs de la fenêtre
#define T_DISPLAY_MAX 34.0      // °C (le bas de l'échelle est T_EXT)
#define FLUX_DISPLAY_MAX 60.0   // W/m²

// Géométries (champ mode de Sim) :
//   carré ; flocon de Koch (profondeur 0 à 4) ;
//   carré crénelé n (témoin non fractal : même aire, même périmètre et même
//     taille de détail que le flocon de profondeur n, n = 1 à 4) ;
//   île de Minkowski n (autre fractale, murs alignés sur la grille, n = 1 à 3).
enum { GEOM_CARRE = 0, GEOM_KOCH = 1, GEOM_CRENEAU = 2, GEOM_MINKOWSKI = 3 };

// Profondeur de Koch par défaut (« fractal ») et profondeurs maximales
#define KOCH_DEPTH_REF 2
#define KOCH_DEPTH_MAX 4
#define CRENEAU_MAX 4
#define MINKOWSKI_MAX 3
#define POLY_MAX 4096           // ≥ 4 · 8^MINKOWSKI_MAX sommets

typedef struct { int x, y; } Point;
typedef struct { double x, y; } Vec2;

typedef struct {
    double *T, *Tstar;      // température de chaque case (°C), étape intermédiaire
    double *lambda;         // conductivité (W/(m·K)) : air brassé ou béton
    double *cap;            // capacité de chaque case : ρc·Δx² (J/(m·K))
    double *src;            // puissance injectée dans chaque case (W/m)
    // Conductances des faces (W/(m·K)) : kx[i] entre i et i+1, ky[i] entre
    // i et i+GRID_SIZE
    double *kx, *ky;
    double *g_ext;          // case de mur : conductance vers l'extérieur
    double *work;           // tableaux de l'algorithme de Thomas (4 · GRID_SIZE)
    uint8_t *zones;
    uint32_t *pixels;       // image de la grille (fenêtre)

    // Contour intérieur de la pièce (m) et position de la source
    Vec2 *poly;
    int n_poly;
    double area, perimeter; // du polygone (m², m)
    double src_x, src_y;    // centre de la source (centre de gravité de la pièce)
    int n_src;              // cases de la source

    double q_int;           // W/m
    double r_si, r_se;      // m²·K/W
    int raw_surfaces;       // 1 : sans corrections de surface (test uniquement)
    double dt;              // s

    long iteration;
    double sim_time;        // s

    int mode;               // GEOM_CARRE, GEOM_KOCH, GEOM_CRENEAU, GEOM_MINKOWSKI
    int depth;              // profondeur (sauf carré)
    double wall;            // épaisseur du mur (m)
    double room_area;       // aire intérieure imposée (m²)

    FILE *csv;              // enregistrement en cours (fenêtre)

    // Statistiques
    double stat_T_air;      // température moyenne de l'air
    double stat_T_air_max;  // température maximale de l'air (source)
    double stat_T_wall;     // température moyenne des murs
    double stat_p_mur;      // puissance de l'air vers les murs (W/m)
    double stat_p_ext;      // puissance des murs vers l'extérieur (W/m)

    // Géométrie courante
    int n_int, n_wall;      // pixels d'air et de mur
    int n_faces_int_wall;   // faces air/mur
    int n_faces_wall_ext;   // faces mur/extérieur
    // Longueurs vraies des surfaces (m), mesurées sur la grille : somme des
    // faces divisées par le facteur d'escalier. len_int ≈ périmètre du
    // polygone ; len_ext : surface extérieure du mur, que le mur raccourcit
    // en bouchant les creux.
    double len_int, len_ext;
} Sim;

// Géométrie
void init_geometry(Sim *s);
void geometry_from_polygon(Sim *s);
void surface_geometry(const Sim *s, double x, double y, int outer, double *f, double *delta);

// Physique
void compute_conductances(Sim *s);
void step_physics(Sim *s);
int solve_steady(Sim *s, double tol, int max_iter);
void compute_stats(Sim *s);

// Entrées-sorties (io.c)
int sim_alloc(Sim *s);
void sim_free(Sim *s);
void sim_reset(Sim *s);
void geom_name(const Sim *s, char *buf, size_t n);
int parse_geometry(const char *g, int *mode, int *depth);
void csv_header(FILE *f, const Sim *s);
void csv_row(FILE *f, const Sim *s);

// Rendu (render.c)
enum { VIEW_TEMP = 0, VIEW_FLUX = 1 };
double wall_flux(const Sim *s, int i);
SDL_Color temp_color(double T);
SDL_Color flux_color(double q);
void fill_pixels(Sim *s, int view);

// Fenêtre interactive (ui.c), taille en points
#define WIN_W 1400
#define WIN_H 780
int run_gui(void);
int capture_gui(const char *path, int depth, int view, double hours);
int film_gui(int depth, int frames);

#endif
