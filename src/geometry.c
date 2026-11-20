#include "thermal_diffusion.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

void flood_fill_safe(uint8_t *zones, int sx, int sy) {
    int n = GRID_SIZE * GRID_SIZE;
    Point *stack = malloc(n * sizeof(Point));
    if (!stack) return;
    int top = 0;
    
    if (sx < 0 || sx >= GRID_SIZE || sy < 0 || sy >= GRID_SIZE) {
        free(stack);
        return;
    }
    if (zones[sy * GRID_SIZE + sx] != ZONE_EXT) {
        free(stack);
        return;
    }
    
    stack[top++] = (Point){sx, sy};
    zones[sy * GRID_SIZE + sx] = ZONE_INT;
    
    while (top > 0) {
        Point p = stack[--top];
        const int nb[4][2] = {{1,0}, {-1,0}, {0,1}, {0,-1}};
        
        for (int k = 0; k < 4; k++) {
            int nx = p.x + nb[k][0], ny = p.y + nb[k][1];
            if (nx >= 0 && nx < GRID_SIZE && ny >= 0 && ny < GRID_SIZE) {
                int idx = ny * GRID_SIZE + nx;
                if (zones[idx] == ZONE_EXT) {
                    zones[idx] = ZONE_INT;
                    if (top < n - 1) stack[top++] = (Point){nx, ny};
                }
            }
        }
    }
    free(stack);
}

void draw_line(uint8_t *zones, Point p1, Point p2, uint8_t val) {
    int dx = abs(p2.x - p1.x), dy = abs(p2.y - p1.y);
    int sx = (p1.x < p2.x) ? 1 : -1, sy = (p1.y < p2.y) ? 1 : -1;
    int err = dx - dy;
    
    while (1) {
        if (p1.x >= 0 && p1.x < GRID_SIZE && p1.y >= 0 && p1.y < GRID_SIZE)
            zones[p1.y * GRID_SIZE + p1.x] = val;
        if (p1.x == p2.x && p1.y == p2.y) break;
        
        int e2 = 2 * err;
        if (e2 > -dy) { err -= dy; p1.x += sx; }
        if (e2 < dx)  { err += dx; p1.y += sy; }
    }
}

void koch(uint8_t *zones, Point a, Point b, int depth) {
    if (depth == 0) {
        draw_line(zones, a, b, ZONE_WALL);
        return;
    }
    
    Point p3 = { a.x + (b.x - a.x) / 3, a.y + (b.y - a.y) / 3 };
    Point p4 = { a.x + 2 * (b.x - a.x) / 3, a.y + 2 * (b.y - a.y) / 3 };
    
    float dx = p4.x - p3.x, dy = p4.y - p3.y;
    float ang = -M_PI / 3.0f;
    Point tip = {
        p3.x + (int)(dx * cos(ang) - dy * sin(ang)),
        p3.y + (int)(dx * sin(ang) + dy * cos(ang))
    };
    
    koch(zones, a, p3, depth - 1);
    koch(zones, p3, tip, depth - 1);
    koch(zones, tip, p4, depth - 1);
    koch(zones, p4, b, depth - 1);
}

int count_interior_pixels(uint8_t *zones) {
    int count = 0;
    for (int i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
        if (zones[i] == ZONE_INT) count++;
    }
    return count;
}

void init_geometry(Sim *s) {
    int N = GRID_SIZE * GRID_SIZE;
    int cx = GRID_SIZE / 2, cy = GRID_SIZE / 2;
    
    for (int i = 0; i < N; i++) {
        s->zones[i] = ZONE_EXT;
        s->T[i] = AMBIENT_TEMP;
        s->alpha[i] = 0.0f;
    }
    
    if (s->mode == 0) {
        int target_area = s->count_int_fractal;
        
        if (target_area == 0) {
            float r_frac = GRID_SIZE / 3.5f;
            float side_frac = r_frac * sqrtf(3.0f);
            float area_triangle = (sqrtf(3.0f) / 4.0f) * side_frac * side_frac;
            float area_koch = (8.0f / 5.0f) * area_triangle;
            target_area = (int)area_koch;
        }
        
        int side_square = (int)sqrtf((float)target_area);
        int r = side_square / 2;
        
        Point p1 = { cx - r, cy - r };
        Point p2 = { cx + r, cy - r };
        Point p3 = { cx + r, cy + r };
        Point p4 = { cx - r, cy + r };
        
        draw_line(s->zones, p1, p2, ZONE_WALL);
        draw_line(s->zones, p2, p3, ZONE_WALL);
        draw_line(s->zones, p3, p4, ZONE_WALL);
        draw_line(s->zones, p4, p1, ZONE_WALL);
        
        flood_fill_safe(s->zones, cx, cy);
        s->count_int_square = count_interior_pixels(s->zones);
        
        printf("[CARRÉ] Aire intérieure = %d pixels (target = %d, écart = %.2f%%)\n",
               s->count_int_square, target_area,
               100.0f * (s->count_int_square - target_area) / (float)target_area);
    }
    else {
        float r_frac = GRID_SIZE / 3.5f;
        float side_frac = r_frac * sqrtf(3.0f);
        int fractal_depth = 2;
        
        Point p1 = { cx - (int)(side_frac / 2.0f), cy + (int)(r_frac / 2.0f) };
        Point p2 = { cx + (int)(side_frac / 2.0f), cy + (int)(r_frac / 2.0f) };
        Point p3 = { cx, cy - (int)r_frac };
        
        koch(s->zones, p1, p3, fractal_depth);
        koch(s->zones, p3, p2, fractal_depth);
        koch(s->zones, p2, p1, fractal_depth);
        
        flood_fill_safe(s->zones, cx, cy);
        s->count_int_fractal = count_interior_pixels(s->zones);
        
        printf("[FRACTALE] Aire intérieure = %d pixels\n", s->count_int_fractal);
    }
    
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            int idx = y * GRID_SIZE + x;
            
            if (s->zones[idx] == ZONE_INT) {
                s->alpha[idx] = DIFF_AIR;
                s->T[idx] = AMBIENT_TEMP;
            }
            else if (s->zones[idx] == ZONE_WALL) {
                s->alpha[idx] = DIFF_WALL;
                s->T[idx] = WALL_TEMP_INIT;
            }
            else {
                s->alpha[idx] = 0.0f;
                s->T[idx] = AMBIENT_TEMP;
            }
        }
    }
    
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            if (in_source(x, y)) s->T[y * GRID_SIZE + x] = SOURCE_TEMP;
        }
    }
}