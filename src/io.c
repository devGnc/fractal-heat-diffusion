#include "thermal_diffusion.h"
#include <stdlib.h>
#include <string.h>

#define CSV_HEADER "Iter,Time_h,T_air,T_air_max,T_mur,P_int,P_mur,P_ext\n"

int sim_alloc(Sim *s) {
    int N = GRID_SIZE * GRID_SIZE;
    s->T = calloc(N, sizeof(double));
    s->Tstar = calloc(N, sizeof(double));
    s->lambda = calloc(N, sizeof(double));
    s->cap = calloc(N, sizeof(double));
    s->src = calloc(N, sizeof(double));
    s->kx = calloc(N, sizeof(double));
    s->ky = calloc(N, sizeof(double));
    s->g_ext = calloc(N, sizeof(double));
    s->work = calloc(4 * GRID_SIZE, sizeof(double));
    s->zones = calloc(N, 1);
    s->pixels = malloc(N * sizeof(uint32_t));
    s->poly = malloc(POLY_MAX * sizeof(Vec2));
    s->q_int = Q_INT;
    s->wall = WALL_THICKNESS;
    s->room_area = ROOM_AREA;
    s->r_si = R_SI;
    s->r_se = R_SE;
    s->dt = DT;
    return s->T && s->Tstar && s->lambda && s->cap && s->src && s->kx && s->ky &&
           s->g_ext && s->work && s->zones && s->pixels && s->poly;
}

void sim_free(Sim *s) {
    free(s->T); free(s->Tstar); free(s->lambda); free(s->cap); free(s->src);
    free(s->kx); free(s->ky); free(s->g_ext); free(s->work);
    free(s->zones); free(s->pixels); free(s->poly);
}

// (Ré)initialise la géométrie du mode courant et remet le temps à zéro.
void sim_reset(Sim *s) {
    s->iteration = 0;
    s->sim_time = 0.0;
    init_geometry(s);
    compute_stats(s);
}

void geom_name(const Sim *s, char *buf, size_t n) {
    switch (s->mode) {
    case GEOM_CARRE: snprintf(buf, n, "carre"); break;
    case GEOM_CRENEAU: snprintf(buf, n, "creneau%d", s->depth); break;
    case GEOM_MINKOWSKI: snprintf(buf, n, "minkowski%d", s->depth); break;
    default: snprintf(buf, n, "koch%d", s->depth);
    }
}

// Lit un nom de géométrie : carre, koch0 … koch4, fractal (= koch2),
// creneau1 … creneau4, minkowski1 … minkowski3. Renvoie 0 si inconnu.
int parse_geometry(const char *g, int *mode, int *depth) {
    static const struct { const char *prefix; int mode, min, max; } kinds[] = {
        { "koch", GEOM_KOCH, 0, KOCH_DEPTH_MAX },
        { "creneau", GEOM_CRENEAU, 1, CRENEAU_MAX },
        { "minkowski", GEOM_MINKOWSKI, 1, MINKOWSKI_MAX },
    };
    *depth = 0;
    if (strcmp(g, "carre") == 0) { *mode = GEOM_CARRE; return 1; }
    if (strcmp(g, "fractal") == 0) { *mode = GEOM_KOCH; *depth = KOCH_DEPTH_REF; return 1; }
    for (size_t k = 0; k < sizeof kinds / sizeof kinds[0]; k++) {
        size_t l = strlen(kinds[k].prefix);
        if (strncmp(g, kinds[k].prefix, l) == 0 && g[l] >= '0' + kinds[k].min
            && g[l] <= '0' + kinds[k].max && g[l + 1] == '\0') {
            *mode = kinds[k].mode;
            *depth = g[l] - '0';
            return 1;
        }
    }
    return 0;
}

// En-tête : une ligne de commentaire décrivant la géométrie (lue par les
// scripts avec pandas.read_csv(..., comment="#")), puis les noms de colonnes.
void csv_header(FILE *f, const Sim *s) {
    char name[16];
    geom_name(s, name, sizeof name);
    fprintf(f, "# geometrie=%s profondeur=%d aire_m2=%.4f perimetre_m=%.4f "
               "surface_int_m=%.4f surface_ext_m=%.4f "
               "px_air=%d px_mur=%d faces_air_mur=%d faces_mur_ext=%d "
               "dx=%g dt=%g epaisseur=%g lambda_air=%g t_ext=%g q_int=%g\n",
            name, s->mode != GEOM_CARRE ? s->depth : -1, s->area, s->perimeter,
            s->len_int, s->len_ext,
            s->n_int, s->n_wall, s->n_faces_int_wall, s->n_faces_wall_ext,
            (double)DX, s->dt, s->wall, LAMBDA_AIR, T_EXT, s->q_int);
    fprintf(f, CSV_HEADER);
}

void csv_row(FILE *f, const Sim *s) {
    fprintf(f, "%ld,%.4f,%.5f,%.5f,%.5f,%.3f,%.4f,%.4f\n",
            s->iteration, s->sim_time / 3600.0, s->stat_T_air, s->stat_T_air_max,
            s->stat_T_wall, s->q_int, s->stat_p_mur, s->stat_p_ext);
}
