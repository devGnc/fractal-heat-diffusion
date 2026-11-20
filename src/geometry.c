#include "thermal_diffusion.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * La pièce est décrite par son contour intérieur, un polygone en mètres
 * (flocon de Koch ou carré). Sur la grille :
 *   - air : cases dont le centre est à l'intérieur du polygone ;
 *   - mur : cases extérieures à moins de s->wall du polygone (le mur
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
 * extérieure, la courbe parallèle au contour à la distance s->wall.
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
    if (outer) *delta = s->wall - d;
    else *delta = inside_polygon(p, s->n_poly, x, y) ? -d : d;
}

// Remplit zones, matériaux, source et températures initiales à partir de
// s->poly, puis calcule les conductances. Utilisé aussi par le test de
// l'anneau. La source est un disque de rayon SOURCE_R centré sur le centre de
// gravité de la pièce.
void geometry_from_polygon(Sim *s) {
    const int G = GRID_SIZE;
    const Vec2 *p = s->poly;
    const int n = s->n_poly;

    // Aire et centre de gravité (formule du polygone), périmètre, cadre
    double a2 = 0.0, gx = 0.0, gy = 0.0;
    s->perimeter = 0.0;
    double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY;
    for (int i = 0, j = n - 1; i < n; j = i++) {
        double cr = p[j].x * p[i].y - p[i].x * p[j].y;
        a2 += cr;
        gx += (p[j].x + p[i].x) * cr;
        gy += (p[j].y + p[i].y) * cr;
        s->perimeter += hypot(p[i].x - p[j].x, p[i].y - p[j].y);
        xmin = fmin(xmin, p[i].x); xmax = fmax(xmax, p[i].x);
        ymin = fmin(ymin, p[i].y); ymax = fmax(ymax, p[i].y);
    }
    s->area = fabs(a2) / 2.0;
    s->src_x = gx / (3.0 * a2);
    s->src_y = gy / (3.0 * a2);

    const double e = s->wall;
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
            s->lambda[i] = z == ZONE_WALL ? LAMBDA_WALL : (z == ZONE_INT ? LAMBDA_AIR : 0.0);
            s->cap[i] = z == ZONE_WALL ? RHOC_WALL * DX * DX
                      : (z == ZONE_INT ? RHOC_AIR * DX * DX : 0.0);
            s->T[i] = s->Tstar[i] = z == ZONE_EXT ? T_EXT : T_INIT;
            s->src[i] = (z == ZONE_INT && hypot(cx - s->src_x, cy - s->src_y) <= SOURCE_R)
                        ? 1.0 : 0.0;
        }
    }

    // Puissance répartie également entre les cases de la source
    s->n_src = 0;
    for (int i = 0; i < G * G; i++) s->n_src += s->src[i] > 0.0;
    if (s->n_src == 0) {
        fprintf(stderr, "Erreur : la source ne contient aucune case d'air\n");
        exit(1);
    }
    for (int i = 0; i < G * G; i++) s->src[i] *= s->q_int / s->n_src;

    compute_conductances(s);
}

/*
 * Carré crénelé n : témoin non fractal du flocon de profondeur n. Même aire,
 * même périmètre P_n, et des créneaux de la taille des plus petits détails du
 * flocon : sur chaque côté de longueur a, m = 2k + 1 intervalles de largeur
 * w = a/m (créneau, creux, créneau...), avec k = 3^(n−1) créneaux de
 * profondeur d tournés vers l'extérieur, et un creux à chaque extrémité. Aire
 * et périmètre :
 *   A = a² + 4·k·w·d,   P = 4a + 8·k·d
 * En éliminant d = (P − 4a)/(8k) : a²·(1 − 2/m) + a·P/(2m) − A = 0.
 */
static void creneau_polygon(Sim *s, double c) {
    const int n = s->depth;
    const double A = s->room_area;
    // Périmètre du flocon de même aire et de profondeur n
    double r = sqrt(A / (koch_area_factor(n) * 3.0 * sqrt(3.0) / 4.0));
    double P = 3.0 * r * sqrt(3.0) * pow(4.0 / 3.0, n);
    int k = 1;
    for (int i = 1; i < n; i++) k *= 3;
    int m = 2 * k + 1;
    double qa = 1.0 - 2.0 / m, qb = P / (2.0 * m);
    double a = (-qb + sqrt(qb * qb + 4.0 * qa * A)) / (2.0 * qa);
    double w = a / m, d = (P - 4.0 * a) / (8.0 * k);
    if (d <= 0.0) {
        fprintf(stderr, "Erreur : carré crénelé %d impossible\n", n);
        exit(1);
    }
    // Côtés parcourus dans le sens trigonométrique ; normale extérieure
    // (u_y, −u_x) pour la direction u.
    double h = a / 2.0;
    Vec2 corner[4] = { { c - h, c - h }, { c + h, c - h }, { c + h, c + h }, { c - h, c + h } };
    for (int side = 0; side < 4; side++) {
        Vec2 p = corner[side], q = corner[(side + 1) % 4];
        double ux = (q.x - p.x) / a, uy = (q.y - p.y) / a, nx = uy, ny = -ux;
        s->poly[s->n_poly++] = p;
        for (int t = 0; t < k; t++) {
            double x0 = (2 * t + 1) * w, x1 = x0 + w;
            s->poly[s->n_poly++] = (Vec2){ p.x + x0 * ux, p.y + x0 * uy };
            s->poly[s->n_poly++] = (Vec2){ p.x + x0 * ux + d * nx, p.y + x0 * uy + d * ny };
            s->poly[s->n_poly++] = (Vec2){ p.x + x1 * ux + d * nx, p.y + x1 * uy + d * ny };
            s->poly[s->n_poly++] = (Vec2){ p.x + x1 * ux, p.y + x1 * uy };
        }
    }
}

/*
 * Île de Minkowski : chaque segment [A, B] est remplacé par 8 segments de
 * longueur |AB|/4, avec u = (B − A)/4 et v = u tourné de −90° :
 *   u, v, u, −v, −v, u, v, u
 * Une bosse vers l'extérieur, une vers l'intérieur : l'aire du carré de départ
 * est conservée, le périmètre double à chaque profondeur et les détails
 * mesurent a/4^n. Tous les murs restent alignés sur la grille.
 * Ajoute les sommets du côté [a, b] sauf b.
 */
static void minkowski_side(Vec2 a, Vec2 b, int depth, Vec2 *out, int *n) {
    if (depth == 0) {
        out[(*n)++] = a;
        return;
    }
    double ux = (b.x - a.x) / 4.0, uy = (b.y - a.y) / 4.0, vx = uy, vy = -ux;
    const int steps[8][2] = { { 1, 0 }, { 0, 1 }, { 1, 0 }, { 0, -1 },
                              { 0, -1 }, { 1, 0 }, { 0, 1 }, { 1, 0 } };
    Vec2 p = a;
    for (int k = 0; k < 8; k++) {
        Vec2 q = { p.x + steps[k][0] * ux + steps[k][1] * vx,
                   p.y + steps[k][0] * uy + steps[k][1] * vy };
        if (k == 7) q = b;  // évite l'accumulation d'arrondis
        minkowski_side(p, q, depth - 1, out, n);
        p = q;
    }
}

void init_geometry(Sim *s) {
    const double c = GRID_SIZE * DX / 2.0;  // centre du domaine
    const double A = s->room_area;
    s->n_poly = 0;

    if (s->mode == GEOM_CARRE || s->mode == GEOM_MINKOWSKI) {
        double h = sqrt(A) / 2.0;
        Vec2 q[4] = { { c - h, c - h }, { c + h, c - h }, { c + h, c + h }, { c - h, c + h } };
        int depth = s->mode == GEOM_MINKOWSKI ? s->depth : 0;
        for (int k = 0; k < 4; k++) minkowski_side(q[k], q[(k + 1) % 4], depth, s->poly, &s->n_poly);
    } else if (s->mode == GEOM_CRENEAU) {
        creneau_polygon(s, c);
    } else {
        // Même aire A à toutes les profondeurs : on règle le rayon r du
        // cercle circonscrit au triangle de départ.
        double r = sqrt(A / (koch_area_factor(s->depth) * 3.0 * sqrt(3.0) / 4.0));
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
    geom_name(s, name, sizeof name);
    printf("[%s] aire %.2f m² (grille : %.2f m²), périmètre intérieur %.2f m, "
           "surface extérieure %.2f m, mur %.3g m, %d px de mur\n",
           name, s->area, s->n_int * DX * DX, s->perimeter, s->len_ext, s->wall, s->n_wall);
}
