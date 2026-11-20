#include "thermal_diffusion.h"
#include <math.h>

/*
 * Couleurs de la fenêtre. Deux vues :
 *   - température : de T_EXT (bleu) à T_DISPLAY_MAX (rouge), air et murs ;
 *   - flux : intensité du flux de chaleur dans les murs (W/m²), de 0 (noir) à
 *     FLUX_DISPLAY_MAX (jaune pâle) ; l'air est gris.
 */

typedef struct { double t; uint8_t r, g, b; } Stop;

static const Stop TEMP_MAP[] = {
    { 0.00,  34,  51, 107 },
    { 0.20,  59, 111, 182 },
    { 0.45,  95, 191, 192 },
    { 0.65, 242, 227, 148 },
    { 0.85, 242, 156,  80 },
    { 1.00, 216,  66,  58 },
};

static const Stop FLUX_MAP[] = {
    { 0.00,   0,   0,   4 },
    { 0.20,  59,  15, 112 },
    { 0.40, 140,  41, 129 },
    { 0.60, 222,  73, 104 },
    { 0.80, 254, 159, 109 },
    { 1.00, 252, 253, 191 },
};

static SDL_Color colormap(const Stop *m, int n, double t) {
    if (!(t > 0.0)) t = 0.0;  // NaN compris
    if (t > 1.0) t = 1.0;
    int k = 1;
    while (k < n - 1 && t > m[k].t) k++;
    double u = (t - m[k - 1].t) / (m[k].t - m[k - 1].t);
    return (SDL_Color){
        (uint8_t)(m[k - 1].r + u * (m[k].r - m[k - 1].r)),
        (uint8_t)(m[k - 1].g + u * (m[k].g - m[k - 1].g)),
        (uint8_t)(m[k - 1].b + u * (m[k].b - m[k - 1].b)),
        255
    };
}

SDL_Color temp_color(double T) {
    return colormap(TEMP_MAP, 6, (T - T_EXT) / (T_DISPLAY_MAX - T_EXT));
}

SDL_Color flux_color(double q) {
    return colormap(FLUX_MAP, 6, q / FLUX_DISPLAY_MAX);
}

// Intensité du flux de chaleur (W/m²) au centre d'une case de mur : moyenne
// des flux des deux faces dans chaque direction, divisée par la longueur Δx
// d'une face.
double wall_flux(const Sim *s, int i) {
    const int G = GRID_SIZE;
    const double *T = s->T;
    double qx = (s->kx[i - 1] * (T[i - 1] - T[i]) + s->kx[i] * (T[i] - T[i + 1])) / (2.0 * DX);
    double qy = (s->ky[i - G] * (T[i - G] - T[i]) + s->ky[i] * (T[i] - T[i + G])) / (2.0 * DX);
    return hypot(qx, qy);
}

static uint32_t pack(SDL_Color c) {
    return (255u << 24) | ((uint32_t)c.r << 16) | ((uint32_t)c.g << 8) | c.b;
}

// Vue température : les murs sont assombris (même teinte, 70 % de la
// luminosité) pour que la forme de la pièce se lise même quand l'air et les
// murs ont la même température.
void fill_pixels(Sim *s, int view) {
    const uint32_t ext = 0xFF101218, air_flux = 0xFF2A2E38;
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        uint8_t z = s->zones[i];
        if (z == ZONE_EXT) s->pixels[i] = ext;
        else if (view == VIEW_TEMP) {
            SDL_Color c = temp_color(s->T[i]);
            if (z == ZONE_WALL) {
                c.r = (uint8_t)(c.r * 7 / 10); c.g = (uint8_t)(c.g * 7 / 10);
                c.b = (uint8_t)(c.b * 7 / 10);
            }
            s->pixels[i] = pack(c);
        }
        else if (z == ZONE_INT) s->pixels[i] = air_flux;
        else s->pixels[i] = pack(flux_color(wall_flux(s, i)));
    }
}
