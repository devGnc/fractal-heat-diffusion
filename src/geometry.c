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

static Point round_point(Vec2 v) {
    return (Point){ (int)lround(v.x), (int)lround(v.y) };
}

// Calcul en flottants, arrondi uniquement au tracé : pas d'erreur cumulée
// d'une profondeur à l'autre.
void koch(uint8_t *zones, Vec2 a, Vec2 b, int depth) {
    if (depth == 0) {
        draw_line(zones, round_point(a), round_point(b), ZONE_WALL);
        return;
    }
    
    double dx = (b.x - a.x) / 3.0, dy = (b.y - a.y) / 3.0;
    Vec2 p3 = { a.x + dx, a.y + dy };
    Vec2 p4 = { a.x + 2.0 * dx, a.y + 2.0 * dy };
    
    double ang = -M_PI / 3.0;
    Vec2 tip = {
        p3.x + dx * cos(ang) - dy * sin(ang),
        p3.y + dx * sin(ang) + dy * cos(ang)
    };
    
    koch(zones, a, p3, depth - 1);
    koch(zones, p3, tip, depth - 1);
    koch(zones, tip, p4, depth - 1);
    koch(zones, p4, b, depth - 1);
}

// Aire du flocon de profondeur n rapportée à celle du triangle initial :
// A_n / A_0 = 8/5 - 3/5 · (4/9)^n
static double koch_area_factor(int n) {
    return 8.0 / 5.0 - 3.0 / 5.0 * pow(4.0 / 9.0, n);
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
        s->alpha[i] = 0.0;
    }
    
    // Rayon du cercle circonscrit au triangle du flocon de référence
    double r_ref = GRID_SIZE / 3.5;
    
    if (s->mode == 0) {
        int target_area = s->count_int_fractal;
        
        if (target_area == 0) {
            double side_ref = r_ref * sqrt(3.0);
            double area_triangle = (sqrt(3.0) / 4.0) * side_ref * side_ref;
            target_area = (int)(koch_area_factor(KOCH_DEPTH_REF) * area_triangle);
        }
        
        // Murs en cx ± r : l'intérieur fait (2r - 1)² pixels. On choisit r pour
        // que cette aire soit la plus proche possible de la cible.
        int r = (int)lround((sqrt((double)target_area) + 1.0) / 2.0);
        
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
        
        printf("[CARRÉ] Aire intérieure = %d pixels (cible = %d, écart = %+.2f%%)\n",
               s->count_int_square, target_area,
               100.0 * (s->count_int_square - target_area) / (double)target_area);
    }
    else {
        // Même aire théorique pour toutes les profondeurs : on agrandit le
        // triangle de départ quand la profondeur est plus faible.
        double r = r_ref * sqrt(koch_area_factor(KOCH_DEPTH_REF) / koch_area_factor(s->depth));
        double side = r * sqrt(3.0);
        
        // Triangle centré sur son centre de gravité (cx, cy)
        Vec2 p1 = { cx - side / 2.0, cy + r / 2.0 };
        Vec2 p2 = { cx + side / 2.0, cy + r / 2.0 };
        Vec2 p3 = { cx, cy - r };
        
        koch(s->zones, p1, p3, s->depth);
        koch(s->zones, p3, p2, s->depth);
        koch(s->zones, p2, p1, s->depth);
        
        flood_fill_safe(s->zones, cx, cy);
        int area = count_interior_pixels(s->zones);
        if (s->depth == KOCH_DEPTH_REF) s->count_int_fractal = area;
        
        printf("[KOCH %d] Aire intérieure = %d pixels\n", s->depth, area);
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
                s->alpha[idx] = 0.0;
                s->T[idx] = AMBIENT_TEMP;
            }
        }
    }
    
    for (int y = 0; y < GRID_SIZE; y++) {
        for (int x = 0; x < GRID_SIZE; x++) {
            if (in_source(x, y)) s->T[y * GRID_SIZE + x] = SOURCE_TEMP;
        }
    }
    
    compute_conductances(s);
}