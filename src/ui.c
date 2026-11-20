#include "thermal_diffusion.h"
#include "font.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/*
 * Fenêtre interactive : le carré (référence) et un flocon de Koch, de même
 * aire et avec les mêmes apports de chaleur, simulés côte à côte. Deux vues
 * (température, flux de chaleur dans les murs), chiffres clés sous chaque
 * pièce, légende, courbe de la température de l'air et commandes.
 *
 * Mise en page en unités logiques (points) ; sur un écran Retina, tout est
 * dessiné en pixels physiques (facteur ui.dpi) pour rester net.
 */

#define WIN_W 1400
#define WIN_H 700
#define MARGIN 20
#define PANEL 460                       // côté d'une pièce à l'écran
#define PANEL_Y 96
#define SIDE_X (2 * PANEL + 3 * MARGIN) // colonne de droite
#define SIDE_W (WIN_W - SIDE_X - MARGIN)
#define HIST_MAX (24 * 20)              // 20 jours, un point par heure

// Itérations par image : 30 s simulées par itération
static const int SPEEDS[] = { 5, 10, 20, 40, 80, 160, 320 };
#define N_SPEEDS (int)(sizeof SPEEDS / sizeof SPEEDS[0])

static const SDL_Color C_BG     = {  18,  20,  26, 255 };
static const SDL_Color C_PANEL  = {  28,  31,  40, 255 };
static const SDL_Color C_GRID   = {  55,  60,  74, 255 };
static const SDL_Color C_TEXT   = { 228, 231, 238, 255 };
static const SDL_Color C_DIM    = { 140, 147, 162, 255 };
static const SDL_Color C_REF    = { 150, 190, 255, 255 };  // carré
static const SDL_Color C_KOCH   = { 255, 172,  92, 255 };  // flocon
static const SDL_Color C_GOOD   = { 125, 220, 145, 255 };
static const SDL_Color C_BAD    = { 240, 112, 100, 255 };

typedef struct {
    Sim sim;
    SDL_Texture *tex;
    float hist[HIST_MAX];   // T_air, une valeur par heure simulée
    int n_hist;
} Room;

static struct {
    SDL_Renderer *ren;
    double dpi;             // pixels physiques par point
} ui;

/* ------------------------------------------------------------------------ */
/* Dessin de base                                                           */
/* ------------------------------------------------------------------------ */

static int L(double v) { return (int)lround(v * ui.dpi); }

static void set_color(SDL_Color c) {
    SDL_SetRenderDrawColor(ui.ren, c.r, c.g, c.b, c.a);
}

static void fill_rect(double x, double y, double w, double h, SDL_Color c) {
    set_color(c);
    SDL_Rect r = { L(x), L(y), L(x + w) - L(x), L(y + h) - L(y) };
    SDL_RenderFillRect(ui.ren, &r);
}

static void line(double x1, double y1, double x2, double y2, SDL_Color c) {
    set_color(c);
    SDL_RenderDrawLine(ui.ren, L(x1), L(y1), L(x2), L(y2));
}

// Taille de police : nombre de points par pixel de la police (5 × 9)
static int font_px(double size) {
    int px = (int)lround(size * ui.dpi);
    return px < 1 ? 1 : px;
}

static double text(double x, double y, double size, SDL_Color c, const char *s) {
    return text_draw(ui.ren, L(x), L(y), font_px(size), c, s) / ui.dpi;
}

static double text_w(double size, const char *s) {
    return text_width(font_px(size), s) / ui.dpi;
}

static void text_right(double x, double y, double size, SDL_Color c, const char *s) {
    text(x - text_w(size, s), y, size, c, s);
}

// Nombres à la française : virgule décimale, vrai signe moins
static void fr(char *buf, size_t n) {
    char tmp[256];
    size_t k = 0;
    for (const char *p = buf; *p && k + 4 < sizeof tmp; p++) {
        if (*p == '.') tmp[k++] = ',';
        else if (*p == '-' && (p == buf || p[-1] == ' ' || p[-1] == '(')) {
            tmp[k++] = (char)0xE2; tmp[k++] = (char)0x88; tmp[k++] = (char)0x92;
        }
        else tmp[k++] = *p;
    }
    tmp[k] = '\0';
    snprintf(buf, n, "%s", tmp);
}

#define FMT(buf, ...) do { snprintf(buf, sizeof buf, __VA_ARGS__); fr(buf, sizeof buf); } while (0)

/* ------------------------------------------------------------------------ */
/* Simulation                                                               */
/* ------------------------------------------------------------------------ */

static void room_reset(Room *r) {
    sim_reset(&r->sim);
    r->n_hist = 0;
    r->hist[r->n_hist++] = (float)r->sim.T_air;
}

static void room_step(Room *r, int n_iter) {
    for (int k = 0; k < n_iter; k++) {
        step_physics(&r->sim);
        if (r->sim.iteration % STATS_EVERY == 0) {
            compute_stats(&r->sim);
            if (r->n_hist < HIST_MAX) r->hist[r->n_hist++] = (float)r->sim.T_air;
            if (r->sim.recording) {
                csv_row(r->sim.csv, &r->sim);
                fflush(r->sim.csv);
            }
        }
    }
    compute_stats(&r->sim);
}

static double room_G(const Room *r) {
    double dT = r->sim.T_air - T_EXT;
    return dT > 0.05 ? r->sim.stat_p_ext / dT : NAN;
}

static int room_steady(const Room *r) {
    return fabs(r->sim.stat_p_ext - r->sim.q_int) < 1e-3 * r->sim.q_int;
}

static void stop_recording(Room *r) {
    if (!r->sim.recording) return;
    fclose(r->sim.csv);
    r->sim.csv = NULL;
    r->sim.recording = 0;
}

static void start_recording(Room *r) {
    char g[16], name[64];
    geom_name(&r->sim, g, sizeof g);
    snprintf(name, sizeof name, "res_%s.csv", g);
    r->sim.csv = fopen(name, "w");
    if (!r->sim.csv) return;
    csv_header(r->sim.csv, &r->sim);
    csv_row(r->sim.csv, &r->sim);
    r->sim.recording = 1;
    printf("Enregistrement vers %s\n", name);
}

/* ------------------------------------------------------------------------ */
/* Éléments de l'interface                                                  */
/* ------------------------------------------------------------------------ */

static void draw_room(Room *r, double x, int view, const char *label, SDL_Color accent) {
    // Titre au-dessus de la pièce
    fill_rect(x, PANEL_Y - 24, 6, 14, accent);
    text(x + 14, PANEL_Y - 25, 1.75, C_TEXT, label);

    // Image de la grille
    fill_pixels(&r->sim, view);
    SDL_UpdateTexture(r->tex, NULL, r->sim.pixels, GRID_SIZE * sizeof(uint32_t));
    SDL_Rect dst = { L(x), L(PANEL_Y), L(x + PANEL) - L(x), L(PANEL_Y + PANEL) - L(PANEL_Y) };
    SDL_RenderCopy(ui.ren, r->tex, NULL, &dst);

    // Barre d'échelle de 1 m
    double m = PANEL / (GRID_SIZE * DX);  // points par mètre
    double bx = x + 14, by = PANEL_Y + PANEL - 16;
    line(bx, by, bx + m, by, C_TEXT);
    line(bx, by - 4, bx, by + 1, C_TEXT);
    line(bx + m, by - 4, bx + m, by + 1, C_TEXT);
    text(bx + m + 6, by - 7, 1.25, C_TEXT, "1 m");
}

static void draw_room_info(const Room *r, const Room *ref, double x) {
    const Sim *s = &r->sim;
    char buf[128];
    double y = PANEL_Y + PANEL + 14;

    text(x, y, 1.5, C_DIM, "Température de l'air");
    FMT(buf, "%.1f °C", s->T_air);
    text(x, y + 18, 4.0, temp_color(s->T_air), buf);

    // G = P/(T_air − T_ext) n'est la conductance de la pièce qu'en régime
    // établi : pendant la montée en température, la valeur est provisoire.
    double cx = x + 200;
    double G = room_G(r);
    int steady = room_steady(r) && room_steady(ref);
    if (isnan(G)) snprintf(buf, sizeof buf, "G = ...");
    else if (steady) FMT(buf, "G = %.1f W/(m·K)", G);
    else FMT(buf, "G provisoire : %.1f W/(m·K)", G);
    text(cx, y - 2, 1.5, steady ? C_TEXT : C_DIM, buf);
    FMT(buf, "Longueur de mur : %.1f m", s->perimeter);
    text(cx, y + 16, 1.5, C_DIM, buf);
    FMT(buf, "Chaleur évacuée : %.0f W/m", s->stat_p_ext);
    text(cx, y + 34, 1.5, C_DIM, buf);
    FMT(buf, "Par mètre de mur : %.1f W/m²", s->stat_p_ext / s->perimeter);
    text(cx, y + 52, 1.5, C_DIM, buf);

    y += 74;
    if (r == ref) {
        text(x, y, 1.5, C_DIM, "Référence : la forme la plus compacte");
        text(x, y + 20, 1.5, C_DIM, "(périmètre minimal pour 4 côtés droits)");
        return;
    }
    double dT = s->T_air - ref->sim.T_air;
    FMT(buf, "Écart avec le carré : %+.1f °C", dT);
    text(x, y, 1.5, dT <= 0 ? C_GOOD : C_BAD, buf);
    double Gr = room_G(ref);
    if (!steady) {
        text(x, y + 20, 1.5, C_DIM, "Gain par degré : attendre le régime établi");
    } else if (!isnan(G) && !isnan(Gr)) {
        double gain = 100.0 * (G / Gr - 1.0);
        if (gain >= 0) FMT(buf, "Évacue %.0f %% de chaleur en plus par degré", gain);
        else FMT(buf, "Évacue %.0f %% de chaleur en moins par degré", -gain);
        text(x, y + 20, 1.5, gain >= 0 ? C_GOOD : C_BAD, buf);
    }
}

static void draw_legend(double x, double y, int view) {
    const double w = SIDE_W, h = 16;
    text(x, y, 1.75, C_TEXT, view == VIEW_TEMP ? "Vue : température" : "Vue : flux de chaleur");
    y += 26;

    // Dégradé, colonne de pixels physiques par colonne
    int x0 = L(x), x1 = L(x + w);
    for (int px = x0; px < x1; px++) {
        double t = (px - x0) / (double)(x1 - x0 - 1);
        SDL_Color c = view == VIEW_TEMP ? temp_color(T_EXT + t * (T_DISPLAY_MAX - T_EXT))
                                        : flux_color(t * FLUX_DISPLAY_MAX);
        set_color(c);
        SDL_RenderDrawLine(ui.ren, px, L(y), px, L(y + h) - 1);
    }

    // Graduations
    char buf[32];
    int n = view == VIEW_TEMP ? 5 : 6;
    for (int k = 0; k <= n; k++) {
        double tx = x + w * k / n;
        line(tx, y + h, tx, y + h + 4, C_DIM);
        if (view == VIEW_TEMP) FMT(buf, "%.0f", T_EXT + k * (T_DISPLAY_MAX - T_EXT) / n);
        else FMT(buf, "%.0f", k * FLUX_DISPLAY_MAX / n);
        double tw = text_w(1.25, buf);
        double lx = k == 0 ? tx : (k == n ? tx - tw : tx - tw / 2);
        text(lx, y + h + 7, 1.25, C_DIM, buf);
    }
    text_right(x + w, y + h + 22, 1.25, C_DIM, view == VIEW_TEMP ? "°C" : "W/m²");

    y += h + 42;
    const char *lines_t[] = {
        "L'air est brassé : une seule température.",
        "Le dégradé dans les murs montre la chaleur",
        "qui traverse le béton vers l'extérieur.",
        "Touche F : voir où la chaleur sort.",
    };
    const char *lines_f[] = {
        "Intensité du flux de chaleur dans les murs.",
        "Clair : le mur évacue beaucoup (pointes) ;",
        "sombre : il évacue peu (creux étroits).",
        "Chaque mètre de mur du flocon évacue moins,",
        "mais le flocon a beaucoup plus de mur.",
    };
    int n_lines = view == VIEW_TEMP ? 4 : 5;
    for (int k = 0; k < n_lines; k++)
        text(x, y + 18 * k, 1.5, C_TEXT, view == VIEW_TEMP ? lines_t[k] : lines_f[k]);
}

static void draw_graph(double x, double y, const Room *ref, const Room *koch) {
    const double w = SIDE_W, h = 190;
    char buf[32];
    text(x, y, 1.75, C_TEXT, "Température de l'air (°C)");
    y += 24;
    fill_rect(x, y, w, h, C_PANEL);

    const double gx = x + 34, gy = y + 10, gw = w - 44, gh = h - 50;
    const double tmin = T_EXT, tmax = T_DISPLAY_MAX;

    // Durée affichée : au moins 2 jours, arrondie au jour supérieur
    int n = ref->n_hist > koch->n_hist ? ref->n_hist : koch->n_hist;
    int days = (int)ceil((n - 1) / 24.0);
    if (days < 2) days = 2;
    int step = days > 10 ? 2 : 1;

    for (int k = 0; k <= 5; k++) {
        double ty = gy + gh - gh * k / 5.0;
        line(gx, ty, gx + gw, ty, C_GRID);
        FMT(buf, "%.0f", tmin + k * (tmax - tmin) / 5.0);
        text_right(gx - 6, ty - 5, 1.25, C_DIM, buf);
    }
    for (int d = 0; d <= days; d += step) {
        double tx = gx + gw * d / days;
        line(tx, gy + gh, tx, gy + gh + 4, C_DIM);
        FMT(buf, "%d", d);
        text(tx - text_w(1.25, buf) / 2, gy + gh + 7, 1.25, C_DIM, buf);
    }
    const char *xl = "temps (jours)";
    text(gx + (gw - text_w(1.25, xl)) / 2, gy + gh + 23, 1.25, C_DIM, xl);

    const Room *rooms[2] = { ref, koch };
    const SDL_Color cols[2] = { C_REF, C_KOCH };
    for (int r = 0; r < 2; r++) {
        SDL_Point pts[HIST_MAX];
        int m = rooms[r]->n_hist;
        for (int k = 0; k < m; k++) {
            double t = fmin(fmax(rooms[r]->hist[k], tmin), tmax);
            pts[k].x = L(gx + gw * (k / 24.0) / days);
            pts[k].y = L(gy + gh - gh * (t - tmin) / (tmax - tmin));
        }
        set_color(cols[r]);
        for (int dy = 0; dy < (ui.dpi > 1.5 ? 2 : 1); dy++) {
            for (int k = 0; k < m; k++) pts[k].y += dy;
            SDL_RenderDrawLines(ui.ren, pts, m);
        }
    }

    // Légende
    double ly = y + h + 8;
    fill_rect(x, ly + 5, 18, 3, C_REF);
    text(x + 24, ly, 1.5, C_TEXT, "Carré");
    FMT(buf, "Koch, profondeur %d", koch->sim.depth);
    fill_rect(x + 110, ly + 5, 18, 3, C_KOCH);
    text(x + 134, ly, 1.5, C_TEXT, buf);

}

static void draw_controls(double x, double y, int speed, int paused, int recording) {
    char buf[64];
    text(x, y, 1.75, C_TEXT, "Commandes");
    y += 26;
    double mins = SPEEDS[speed] * DT / 60.0;
    char rate[32];
    if (mins < 60) FMT(rate, "%.0f min", mins);
    else FMT(rate, "%.1f h", mins / 60.0);
    FMT(buf, "vitesse (%s par image)", rate);
    const char *keys[] = { "0 à 4", "F", "Espace", "+ / -", "R", "S", "Q" };
    const char *what[] = { "profondeur du flocon", "température / flux",
                           paused ? "reprendre" : "pause", buf, "recommencer",
                           recording ? "arrêter l'enregistrement" : "enregistrer en CSV",
                           "quitter" };
    for (int k = 0; k < 7; k++) {
        text(x, y + 19 * k, 1.5, C_KOCH, keys[k]);
        text(x + 80, y + 19 * k, 1.5, C_TEXT, what[k]);
    }
}

static void draw_header(const Room *ref, const Room *koch, int paused, int recording) {
    char buf[128];
    text(MARGIN, 14, 2.5, C_TEXT, "Même aire, même chaleur produite : quelle forme évacue le mieux ?");
    FMT(buf, "Pièces de %.0f m², murs en béton de %.0f cm, %.0f W/m d'apports internes, "
             "extérieur à %.0f °C", ROOM_AREA, 100 * WALL_THICKNESS, Q_INT, T_EXT);
    text(MARGIN, 44, 1.5, C_DIM, buf);

    double t = ref->sim.sim_time;
    int d = (int)(t / 86400.0), h = (int)fmod(t / 3600.0, 24.0);
    snprintf(buf, sizeof buf, "t = %d j %02d h", d, h);
    text_right(WIN_W - MARGIN, 12, 2.5, C_TEXT, buf);
    const char *status = paused ? "PAUSE" : (room_steady(ref) && room_steady(koch)) ?
                         "régime établi" : "en cours";
    SDL_Color sc = paused ? C_BAD : (room_steady(ref) && room_steady(koch)) ? C_GOOD : C_DIM;
    text_right(WIN_W - MARGIN, 44, 1.5, sc, status);
    if (recording) {
        double w = text_w(1.5, status);
        text_right(WIN_W - MARGIN - w - 16, 44, 1.5, C_BAD, "REC");
    }
}

/* ------------------------------------------------------------------------ */
/* Boucle principale                                                        */
/* ------------------------------------------------------------------------ */

typedef struct {
    Room rooms[2];          // 0 : carré, 1 : flocon
    int view, speed, paused;
} App;

static int app_init(App *a, SDL_Renderer *ren, int depth) {
    ui.ren = ren;
    SDL_SetRenderDrawBlendMode(ren, SDL_BLENDMODE_BLEND);
    if (!font_init(ren)) {
        fprintf(stderr, "Police : %s\n", SDL_GetError());
        return 0;
    }
    SDL_SetHint(SDL_HINT_RENDER_SCALE_QUALITY, "1");  // pièces lissées
    for (int k = 0; k < 2; k++) {
        if (!sim_alloc(&a->rooms[k].sim)) {
            fprintf(stderr, "Allocation mémoire impossible\n");
            return 0;
        }
        a->rooms[k].tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
                                            SDL_TEXTUREACCESS_STREAMING, GRID_SIZE, GRID_SIZE);
        if (!a->rooms[k].tex) return 0;
    }
    a->rooms[0].sim.mode = 0;
    a->rooms[1].sim.mode = 1;
    a->rooms[1].sim.depth = depth;
    room_reset(&a->rooms[0]);
    room_reset(&a->rooms[1]);
    a->view = VIEW_TEMP;
    a->speed = 2;
    a->paused = 0;
    return 1;
}

static void app_free(App *a) {
    for (int k = 0; k < 2; k++) {
        stop_recording(&a->rooms[k]);
        if (a->rooms[k].tex) SDL_DestroyTexture(a->rooms[k].tex);
        sim_free(&a->rooms[k].sim);
    }
    font_free();
}

static void draw_frame(App *a) {
    Room *ref = &a->rooms[0], *koch = &a->rooms[1];
    set_color(C_BG);
    SDL_RenderClear(ui.ren);

    char label[64];
    snprintf(label, sizeof label, "Flocon de Koch, profondeur %d", koch->sim.depth);
    draw_header(ref, koch, a->paused, ref->sim.recording);
    draw_room(ref, MARGIN, a->view, "Carré", C_REF);
    draw_room(koch, 2 * MARGIN + PANEL, a->view, label, C_KOCH);
    draw_room_info(ref, ref, MARGIN);
    draw_room_info(koch, ref, 2 * MARGIN + PANEL);
    draw_legend(SIDE_X, PANEL_Y - 25, a->view);
    draw_graph(SIDE_X, PANEL_Y + 168, ref, koch);
    draw_controls(SIDE_X, PANEL_Y + 428, a->speed, a->paused, ref->sim.recording);
}

int run_gui(void) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init : %s\n", SDL_GetError());
        return 1;
    }
    SDL_Window *win = SDL_CreateWindow("Diffusion thermique et géométries fractales",
                                       SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       WIN_W, WIN_H, SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_Renderer *ren = win ? SDL_CreateRenderer(win, -1, SDL_RENDERER_ACCELERATED |
                                                          SDL_RENDERER_PRESENTVSYNC) : NULL;
    if (!ren) {
        fprintf(stderr, "Fenêtre SDL : %s\n", SDL_GetError());
        return 1;
    }
    int ow, oh;
    SDL_GetRendererOutputSize(ren, &ow, &oh);
    ui.dpi = ow / (double)WIN_W;

    static App app;
    if (!app_init(&app, ren, KOCH_DEPTH_REF + 1)) return 1;
    Room *ref = &app.rooms[0], *koch = &app.rooms[1];
    check_stability(&koch->sim);
    printf("Fenêtre : carré et flocon de Koch côte à côte (Δx = %.1f cm, Δt = %.0f s).\n"
           "Touches : 0-4 profondeur, F vue, Espace pause, +/- vitesse, R, S, Q.\n",
           100.0 * DX, (double)DT);

    int running = 1;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
            if (e.type != SDL_KEYDOWN) continue;
            SDL_Keycode k = e.key.keysym.sym;
            if (k == SDLK_q || k == SDLK_ESCAPE) running = 0;
            else if (k == SDLK_f) app.view = !app.view;
            else if (k == SDLK_SPACE) app.paused = !app.paused;
            else if (k == SDLK_PLUS || k == SDLK_EQUALS || k == SDLK_KP_PLUS) {
                if (app.speed < N_SPEEDS - 1) app.speed++;
            }
            else if (k == SDLK_MINUS || k == SDLK_KP_MINUS) {
                if (app.speed > 0) app.speed--;
            }
            else if (k == SDLK_r || (k >= SDLK_0 && k <= SDLK_0 + KOCH_DEPTH_MAX)) {
                // Les deux pièces repartent ensemble : les courbes restent comparables
                stop_recording(ref);
                stop_recording(koch);
                if (k != SDLK_r) koch->sim.depth = k - SDLK_0;
                room_reset(ref);
                room_reset(koch);
            }
            else if (k == SDLK_s) {
                if (ref->sim.recording) {
                    stop_recording(ref);
                    stop_recording(koch);
                    printf("Enregistrement arrêté.\n");
                } else {
                    start_recording(ref);
                    start_recording(koch);
                }
            }
        }

        if (!app.paused) {
            room_step(ref, SPEEDS[app.speed]);
            room_step(koch, SPEEDS[app.speed]);
        }
        draw_frame(&app);
        SDL_RenderPresent(ren);
    }

    app_free(&app);
    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}

// Image de la fenêtre sans l'ouvrir (rendu logiciel, échelle 2 comme un écran
// Retina), après `hours` heures simulées : vérification et illustrations.
int capture_gui(const char *path, int depth, int view, double hours) {
    SDL_Surface *surf = SDL_CreateRGBSurfaceWithFormat(0, 2 * WIN_W, 2 * WIN_H, 32,
                                                       SDL_PIXELFORMAT_ARGB8888);
    SDL_Renderer *ren = surf ? SDL_CreateSoftwareRenderer(surf) : NULL;
    if (!ren) {
        fprintf(stderr, "Rendu logiciel : %s\n", SDL_GetError());
        return 1;
    }
    ui.dpi = 2.0;
    static App app;
    if (!app_init(&app, ren, depth)) return 1;
    app.view = view;
    long n = (long)lround(hours * 3600.0 / DT);
    room_step(&app.rooms[0], (int)n);
    room_step(&app.rooms[1], (int)n);
    draw_frame(&app);
    SDL_RenderPresent(ren);
    int ok = SDL_SaveBMP(surf, path) == 0;
    if (!ok) fprintf(stderr, "Enregistrement de %s : %s\n", path, SDL_GetError());
    else printf("Image enregistrée : %s (t = %.0f h)\n", path, hours);
    app_free(&app);
    SDL_DestroyRenderer(ren);
    SDL_FreeSurface(surf);
    return !ok;
}
