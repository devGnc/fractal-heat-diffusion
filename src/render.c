#include "thermal_diffusion.h"
#include <math.h>

// Carte de couleurs de T_EXT (bleu sombre) à T_DISPLAY_MAX (jaune). L'air est
// éclairci pour le distinguer des murs ; l'extérieur est noir.
static uint32_t get_color(double T, uint8_t zone) {
    if (zone == ZONE_EXT) return 0xFF101010;

    double norm = (T - T_EXT) / (T_DISPLAY_MAX - T_EXT);
    if (norm < 0.0) norm = 0.0;
    if (norm > 1.0) norm = 1.0;

    double v = sqrt(norm);
    double r = v * 255, g = v * v * 200, b = (1.0 - v) * 120;
    if (zone == ZONE_INT) {
        r = r + 0.3 * (255 - r);
        g = g + 0.3 * (255 - g);
        b = b + 0.3 * (255 - b);
    }
    return (255u << 24) | ((uint32_t)r << 16) | ((uint32_t)g << 8) | (uint32_t)b;
}

void render(SDL_Renderer *ren, Sim *s) {
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        s->pixels[i] = get_color(s->T[i], s->zones[i]);
    }
    SDL_UpdateTexture(s->tex, NULL, s->pixels, GRID_SIZE * sizeof(uint32_t));
    SDL_RenderCopy(ren, s->tex, NULL, NULL);
}
