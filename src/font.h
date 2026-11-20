#ifndef FONT_H
#define FONT_H

#include <SDL2/SDL.h>

// Police bitmap intégrée (voir font.c). Texte en UTF-8.
int font_init(SDL_Renderer *ren);
void font_free(void);
int text_draw(SDL_Renderer *ren, int x, int y, int scale, SDL_Color col, const char *s);
int text_width(int scale, const char *s);
int text_height(int scale);

#endif
