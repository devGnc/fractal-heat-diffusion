#include "thermal_diffusion.h"
#include <stdlib.h>

#define CSV_HEADER "Iter,Time_h,T_air,T_mur,P_int,P_mur,P_ext\n"

int sim_alloc(Sim *s) {
    int N = GRID_SIZE * GRID_SIZE;
    s->T = calloc(N, sizeof(double));
    s->Tnext = calloc(N, sizeof(double));
    s->lambda = calloc(N, sizeof(double));
    s->cap = calloc(N, sizeof(double));
    s->kx = calloc(N, sizeof(double));
    s->ky = calloc(N, sizeof(double));
    s->g_air = calloc(N, sizeof(double));
    s->g_ext = calloc(N, sizeof(double));
    s->zones = calloc(N, 1);
    s->pixels = malloc(N * sizeof(uint32_t));
    s->poly = malloc(POLY_MAX * sizeof(Vec2));
    s->q_int = Q_INT;
    s->r_si = R_SI;
    s->r_se = R_SE;
    return s->T && s->Tnext && s->lambda && s->cap && s->kx && s->ky &&
           s->g_air && s->g_ext && s->zones && s->pixels && s->poly;
}

void sim_free(Sim *s) {
    free(s->T); free(s->Tnext); free(s->lambda); free(s->cap);
    free(s->kx); free(s->ky); free(s->g_air); free(s->g_ext);
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
    if (s->mode == 0) snprintf(buf, n, "carre");
    else snprintf(buf, n, "koch%d", s->depth);
}

// En-tête : une ligne de commentaire décrivant la géométrie (lue par les
// scripts avec pandas.read_csv(..., comment="#")), puis les noms de colonnes.
void csv_header(FILE *f, const Sim *s) {
    char name[16];
    geom_name(s, name, sizeof name);
    fprintf(f, "# geometrie=%s profondeur=%d aire_m2=%.4f perimetre_m=%.4f "
               "px_air=%d px_mur=%d faces_air_mur=%d faces_mur_ext=%d "
               "dx=%g dt=%g epaisseur=%g t_ext=%g q_int=%g\n",
            name, s->mode ? s->depth : -1, s->area, s->perimeter,
            s->n_int, s->n_wall, s->n_faces_int_wall, s->n_faces_wall_ext,
            (double)DX, (double)DT, (double)WALL_THICKNESS, T_EXT, s->q_int);
    fprintf(f, CSV_HEADER);
}

void csv_row(FILE *f, const Sim *s) {
    fprintf(f, "%ld,%.4f,%.5f,%.5f,%.3f,%.4f,%.4f\n",
            s->iteration, s->sim_time / 3600.0, s->T_air, s->stat_T_wall,
            s->q_int, s->stat_p_mur, s->stat_p_ext);
}
