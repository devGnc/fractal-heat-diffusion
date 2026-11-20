#include "thermal_diffusion.h"
#include <math.h>

uint32_t get_color(float T, uint8_t zone) {
    if (zone == ZONE_WALL) return 0xFF6B6B6B;
    
    if (zone == ZONE_EXT) return 0xFF101010;
    
    float norm = (T - AMBIENT_TEMP) / (SOURCE_TEMP - AMBIENT_TEMP);
    if (norm < 0.0f) norm = 0.0f;
    if (norm > 1.0f) norm = 1.0f;
    
    float v = sqrtf(norm);
    uint8_t r = (uint8_t)(v * 255);
    uint8_t g = (uint8_t)(v * v * 200);
    uint8_t b = (uint8_t)((1.0f - v) * 120);
    
    return (255u << 24) | (r << 16) | (g << 8) | b;
}

void render(SDL_Renderer *ren, Sim *s) {
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        s->pixels[i] = get_color(s->T[i], s->zones[i]);
    }
    SDL_UpdateTexture(s->tex, NULL, s->pixels, GRID_SIZE * sizeof(uint32_t));
    SDL_RenderCopy(ren, s->tex, NULL, NULL);
}