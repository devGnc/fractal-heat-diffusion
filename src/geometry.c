#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * La pièce est décrite par son contour intérieur, un polygone en mètres
 * (flocon de Koch ou carré). Sur la grille :
 *   - air : cases dont le centre est à l'intérieur du polygone ;
 *   - mur : cases extérieures à moins de WALL_THICKNESS du polygone (le mur
 *     est l'épaisseur du contour, vers l'extérieur) ;
 *   - extérieur : le reste.
 * Le mur suit le contour vers l'extérieur : il comble les creux plus étroits
 * que son épaisseur, ce qui limite physiquement la profondeur utile du flocon.
 */

// Aire du flocon de profondeur n rapportée à celle du triangle initial :
// A_n / A_0 = 8/5 − 3/5 · (4/9)^n
static double koch_area_factor(int n) {
    return 8.0 / 5.0 - 3.0 / 5.0 * pow(4.0 / 9.0, n);
}

// Ajoute les sommets du côté [a, b] de profondeur n, sauf b (premier sommet
// du côté suivant).
static void koch_side(Vec2 a, Vec2 b, int depth, Vec2 *out, int *n) {
    if (depth == 0) {
        out[(*n)++] = a;
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
    koch_side(a, p3, depth - 1, out, n);
    koch_side(p3, tip, depth - 1, out, n);
    koch_side(tip, p4, depth - 1, out, n);
    koch_side(p4, b, depth - 1, out, n);
}

// Point dans le polygone (règle pair-impair)
static int inside_polygon(const Vec2 *p, int n, double x, double y) {
    int c = 0;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        if ((p[i].y > y) != (p[j].y > y) &&
            x < (p[j].x - p[i].x) * (y - p[i].y) / (p[j].y - p[i].y) + p[i].x)
            c = !c;
    }
    return c;
}

// Distance du point (x, y) au contour du polygone
static double distance_polygon(const Vec2 *p, int n, double x, double y) {
    double best = INFINITY;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double ex = p[i].x - p[j].x, ey = p[i].y - p[j].y;
        double t = ((x - p[j].x) * ex + (y - p[j].y) * ey) / (ex * ex + ey * ey);
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        double dx = x - (p[j].x + t * ex), dy = y - (p[j].y + t * ey);
        double d2 = dx * dx + dy * dy;
        if (d2 < best) best = d2;
    }
    return sqrt(best);
}

/*
 * Corrections de surface, au milieu (x, y) d'une face air/mur (outer = 0) ou
 * mur/extérieur (outer = 1). La surface intérieure est le contour ; la surface
 * extérieure, la courbe parallèle au contour à la distance WALL_THICKNESS.
 *
 * f : facteur d'escalier |n_x| + |n_y| de la normale n à la surface. Sur la
 *   grille, une surface inclinée d'un angle θ est une suite de marches : pour
 *   une longueur L de mur, elle compte L·(|cos θ| + |sin θ|) de faces. Une
 *   résistance posée sur chaque face surestime donc l'échange de ce facteur
 *   (jusqu'à √2 à 45°), et ce biais ne diminue pas en raffinant la grille.
 *   physics.c multiplie les résistances de surface de la face par f.
 *
 * delta : épaisseur de mur manquante entre la face et la vraie surface (m),
 *   négative si la grille met trop de mur. Les cases sont entières : sans
 *   correction, un mur de 20 cm sur des cases de 3,5 cm n'en fait que 17,5.
 *   physics.c ajoute la résistance delta/λ à la face.
 *
 * Normale : celle du segment le plus proche (surface intérieure), ou la
 * direction depuis le point le plus proche du contour (surface extérieure).
 */
void surface_geometry(const Sim *s, double x, double y, int outer, double *f, double *delta) {
    const Vec2 *p = s->poly;
    double best = INFINITY, nx = 0.0, ny = 0.0;
    for (int i = 0, j = s->n_poly - 1; i < s->n_poly; j = i++) {
        double ex = p[i].x - p[j].x, ey = p[i].y - p[j].y;
        double t = ((x - p[j].x) * ex + (y - p[j].y) * ey) / (ex * ex + ey * ey);
        t = t < 0.0 ? 0.0 : (t > 1.0 ? 1.0 : t);
        double qx = p[j].x + t * ex, qy = p[j].y + t * ey;
        double d2 = (x - qx) * (x - qx) + (y - qy) * (y - qy);
        if (d2 < best) {
            best = d2;
            if (outer) { nx = x - qx; ny = y - qy; }
            else { nx = -ey; ny = ex; }
        }
    }
    double len = hypot(nx, ny);
    *f = len > 0.0 ? (fabs(nx) + fabs(ny)) / len : 1.0;
    double d = sqrt(best);
    if (outer) *delta = WALL_THICKNESS - d;
    else *delta = inside_polygon(p, s->n_poly, x, y) ? -d : d;
}

// Remplit zones, matériaux et températures initiales à partir de s->poly,
// puis calcule les conductances. Utilisé aussi par le test de l'anneau.
void geometry_from_polygon(Sim *s) {
    const int G = GRID_SIZE;
    const Vec2 *p = s->poly;
    const int n = s->n_poly;

    s->area = s->perimeter = 0.0;
    double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        s->area += p[j].x * p[i].y - p[i].x * p[j].y;
        s->perimeter += hypot(p[i].x - p[j].x, p[i].y - p[j].y);
        xmin = fmin(xmin, p[i].x); xmax = fmax(xmax, p[i].x);
        ymin = fmin(ymin, p[i].y); ymax = fmax(ymax, p[i].y);
    }
    s->area = fabs(s->area) / 2.0;

    const double e = WALL_THICKNESS;
    if (xmin - e < DX || ymin - e < DX || xmax + e > (G - 1) * DX || ymax + e > (G - 1) * DX) {
        fprintf(stderr, "Erreur : la pièce et ses murs dépassent de la grille\n");
        exit(1);
    }

    for (int y = 0; y < G; y++) {
        for (int x = 0; x < G; x++) {
            int i = y * G + x;
            double cx = (x + 0.5) * DX, cy = (y + 0.5) * DX;
            uint8_t z = ZONE_EXT;
            if (cx >= xmin - e && cx <= xmax + e && cy >= ymin - e && cy <= ymax + e) {
                if (inside_polygon(p, n, cx, cy)) z = ZONE_INT;
                else if (distance_polygon(p, n, cx, cy) <= e) z = ZONE_WALL;
            }
            s->zones[i] = z;
            s->lambda[i] = z == ZONE_WALL ? LAMBDA_WALL : 0.0;
            s->cap[i] = z == ZONE_WALL ? RHOC_WALL * DX * DX : 0.0;
            s->T[i] = s->Tnext[i] = z == ZONE_EXT ? T_EXT : T_INIT;
        }
    }

    compute_conductances(s);
    s->T_air = T_INIT;
    s->C_air = RHOC_AIR * s->n_int * DX * DX;
}

void init_geometry(Sim *s) {
    const double c = GRID_SIZE * DX / 2.0;  // centre du domaine
    s->n_poly = 0;

    if (s->mode == 0) {
        double h = sqrt(ROOM_AREA) / 2.0;
        s->poly[0] = (Vec2){ c - h, c - h };
        s->poly[1] = (Vec2){ c + h, c - h };
        s->poly[2] = (Vec2){ c + h, c + h };
        s->poly[3] = (Vec2){ c - h, c + h };
        s->n_poly = 4;
    } else {
        // Même aire ROOM_AREA à toutes les profondeurs : on règle le rayon r
        // du cercle circonscrit au triangle de départ.
        double r = sqrt(ROOM_AREA / (koch_area_factor(s->depth) * 3.0 * sqrt(3.0) / 4.0));
        double side = r * sqrt(3.0);
        // Triangle centré sur son centre de gravité
        Vec2 p1 = { c - side / 2.0, c + r / 2.0 };
        Vec2 p2 = { c + side / 2.0, c + r / 2.0 };
        Vec2 p3 = { c, c - r };
        koch_side(p1, p3, s->depth, s->poly, &s->n_poly);
        koch_side(p3, p2, s->depth, s->poly, &s->n_poly);
        koch_side(p2, p1, s->depth, s->poly, &s->n_poly);
    }

    // Centrage du cadre englobant (le centre de gravité du triangle n'est pas
    // au milieu de son cadre) : la pièce occupe mieux le domaine.
    double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY;
    for (int k = 0; k < s->n_poly; k++) {
        xmin = fmin(xmin, s->poly[k].x); xmax = fmax(xmax, s->poly[k].x);
        ymin = fmin(ymin, s->poly[k].y); ymax = fmax(ymax, s->poly[k].y);
    }
    double sx = c - (xmin + xmax) / 2.0, sy = c - (ymin + ymax) / 2.0;
    for (int k = 0; k < s->n_poly; k++) {
        s->poly[k].x += sx;
        s->poly[k].y += sy;
    }

    geometry_from_polygon(s);

    char name[16];
    if (s->mode) snprintf(name, sizeof name, "KOCH %d", s->depth);
    else snprintf(name, sizeof name, "CARRÉ");
    printf("[%s] aire %.2f m² (grille : %.2f m²), périmètre intérieur %.2f m, "
           "%d px de mur\n",
           name, s->area, s->n_int * DX * DX, s->perimeter, s->n_wall);
}
