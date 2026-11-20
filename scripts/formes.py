"""Contours des pièces, mêmes constructions que src/geometry.c (en mètres,
centrés sur l'origine), et taille de leurs plus petits détails.

    koch(n, A)       flocon de Koch de profondeur n et d'aire A
    creneau(n, A)    carré crénelé, témoin non fractal du flocon de profondeur n
    minkowski(n, A)  île de Minkowski de profondeur n (carré pour n = 0)
    carre(A)
    detail(famille, n, A)  plus petit détail (m) : segment du flocon, côté
                           d'un créneau, segment de l'île de Minkowski
"""
import math


def koch_area_factor(n):
    return 8 / 5 - 3 / 5 * (4 / 9) ** n


def _koch_side(a, b, depth, out):
    if depth == 0:
        out.append(a)
        return
    dx, dy = (b[0] - a[0]) / 3, (b[1] - a[1]) / 3
    p3 = (a[0] + dx, a[1] + dy)
    p4 = (a[0] + 2 * dx, a[1] + 2 * dy)
    ang = -math.pi / 3
    tip = (p3[0] + dx * math.cos(ang) - dy * math.sin(ang),
           p3[1] + dx * math.sin(ang) + dy * math.cos(ang))
    for u, v in ((a, p3), (p3, tip), (tip, p4), (p4, b)):
        _koch_side(u, v, depth - 1, out)


def koch_triangle_side(n, A=36.0):
    """Côté s_n du triangle de départ du flocon de profondeur n et d'aire A."""
    r = math.sqrt(A / (koch_area_factor(n) * 3 * math.sqrt(3) / 4))
    return r * math.sqrt(3)


def koch(n, A=36.0):
    s = koch_triangle_side(n, A)
    r = s / math.sqrt(3)
    p1, p2, p3 = (-s / 2, r / 2), (s / 2, r / 2), (0, -r)
    pts = []
    for u, v in ((p1, p3), (p3, p2), (p2, p1)):
        _koch_side(u, v, n, pts)
    return pts


def koch_perimeter(n, A=36.0):
    return 3 * koch_triangle_side(n, A) * (4 / 3) ** n


def creneau_params(n, A=36.0):
    """Côté a, largeur w et profondeur d des créneaux, nombre k par côté."""
    P = koch_perimeter(n, A)
    k = 3 ** (n - 1)
    m = 2 * k + 1
    qa, qb = 1 - 2 / m, P / (2 * m)
    a = (-qb + math.sqrt(qb * qb + 4 * qa * A)) / (2 * qa)
    return a, a / m, (P - 4 * a) / (8 * k), k


def creneau(n, A=36.0):
    a, w, d, k = creneau_params(n, A)
    h = a / 2
    corners = [(-h, -h), (h, -h), (h, h), (-h, h)]
    pts = []
    for s in range(4):
        p, q = corners[s], corners[(s + 1) % 4]
        ux, uy = (q[0] - p[0]) / a, (q[1] - p[1]) / a
        nx, ny = uy, -ux
        pts.append(p)
        for t in range(k):
            x0, x1 = (2 * t + 1) * w, (2 * t + 2) * w
            pts += [(p[0] + x0 * ux, p[1] + x0 * uy),
                    (p[0] + x0 * ux + d * nx, p[1] + x0 * uy + d * ny),
                    (p[0] + x1 * ux + d * nx, p[1] + x1 * uy + d * ny),
                    (p[0] + x1 * ux, p[1] + x1 * uy)]
    return pts


def _minkowski_side(a, b, depth, out):
    if depth == 0:
        out.append(a)
        return
    ux, uy = (b[0] - a[0]) / 4, (b[1] - a[1]) / 4
    vx, vy = uy, -ux
    p = a
    for k, (su, sv) in enumerate(((1, 0), (0, 1), (1, 0), (0, -1),
                                  (0, -1), (1, 0), (0, 1), (1, 0))):
        q = b if k == 7 else (p[0] + su * ux + sv * vx, p[1] + su * uy + sv * vy)
        _minkowski_side(p, q, depth - 1, out)
        p = q


def minkowski(n, A=36.0):
    h = math.sqrt(A) / 2
    corners = [(-h, -h), (h, -h), (h, h), (-h, h)]
    pts = []
    for s in range(4):
        _minkowski_side(corners[s], corners[(s + 1) % 4], n, pts)
    return pts


def carre(A=36.0):
    return minkowski(0, A)


def detail(famille, n, A=36.0):
    if famille == "koch":
        return koch_triangle_side(n, A) / 3 ** n
    if famille == "minkowski":
        return math.sqrt(A) / 4 ** n
    if famille == "creneau":
        return creneau_params(n, A)[1]
    return math.sqrt(A)


def perimetre(pts):
    return sum(math.dist(pts[i - 1], pts[i]) for i in range(len(pts)))
