#include "thermal_diffusion.h"
#include "font.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

/*
 * Fenêtre interactive : le carré et un flocon de Koch, de même aire et avec la
 * même source de chaleur, simulés côte à côte. L'image d'abord : deux grandes
 * pièces, la température moyenne de l'air sous chacune, une barre de couleurs
 * et le temps. L'aide (touche H) liste les commandes.
 *
 * Vitesse automatique : le temps simulé par image croît avec le temps écoulé
 * (2 % de t par image). On voit d'abord la chaleur quitter la source et
 * remplir la pièce en quelques minutes, puis les murs se réchauffer en
 * quelques jours. Le pas de temps est plafonné à 30 s (au-delà, l'air,
 * de capacité très faible, n'est plus amorti correctement) : une image peut
 * enchaîner plusieurs pas.
 *
 * Mise en page en points ; sur un écran Retina, tout est dessiné en pixels
 * physiques (facteur ui.dpi) pour rester net.
 */

#define MARGIN 20
#define PANEL 670
#define STRIP_Y (MARGIN + PANEL)
#define DT_MAX 30.0                     // s : pas de temps maximal
#define REC_EVERY 600.0                 // s simulées entre deux lignes de CSV

// Vitesses manuelles : temps simulé par image (s)
static const double SPEEDS[] = { 1, 2, 5, 10, 30, 60, 150, 300, 600, 1200 };
#define N_SPEEDS (int)(sizeof SPEEDS / sizeof SPEEDS[0])

static const SDL_Color C_BG    = {  14,  16,  21, 255 };
static const SDL_Color C_TEXT  = { 228, 231, 238, 255 };
static const SDL_Color C_DIM   = { 125, 132, 148, 255 };
static const SDL_Color C_REF   = { 150, 190, 255, 255 };
static const SDL_Color C_KOCH  = { 255, 172,  92, 255 };
static const SDL_Color C_REC   = { 240,  90,  80, 255 };

typedef struct {
    Sim sim;
    SDL_Texture *tex;
} Room;

typedef struct {
    Room rooms[2];          // 0 : carré, 1 : flocon
    int view, iso, help, paused;
    int speed;              // -1 : automatique, sinon indice dans SPEEDS
    char rec_dir[128];      // dossier de la session d'enregistrement ("" : aucune)
    double next_rec;        // prochain instant d'enregistrement (s)
} App;

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

// Nombres à la française : virgule décimale
static void fr(char *buf) {
    for (char *p = buf; *p; p++) if (*p == '.') *p = ',';
}

static void format_time(char *buf, size_t n, double t) {
    int d = (int)(t / 86400.0), h = (int)fmod(t / 3600.0, 24.0), m = (int)fmod(t / 60.0, 60.0);
    if (d > 0) snprintf(buf, n, "%d j %02d h", d, h);
    else if (h > 0) snprintf(buf, n, "%d h %02d min", h, m);
    else snprintf(buf, n, "%d min %02d s", m, (int)fmod(t, 60.0));
}

/* ------------------------------------------------------------------------ */
/* Simulation et enregistrement                                             */
/* ------------------------------------------------------------------------ */

// Temps simulé par image
static double frame_span(const App *a) {
    if (a->speed >= 0) return SPEEDS[a->speed];
    double span = 0.02 * a->rooms[0].sim.sim_time;
    return span < 1.0 ? 1.0 : (span > 1200.0 ? 1200.0 : span);
}

static void stop_recording(App *a) {
    if (!a->rec_dir[0]) return;
    for (int k = 0; k < 2; k++) {
        if (a->rooms[k].sim.csv) fclose(a->rooms[k].sim.csv);
        a->rooms[k].sim.csv = NULL;
    }
    printf("Enregistrement terminé : %s\n", a->rec_dir);
    a->rec_dir[0] = '\0';
}

// Une session = un dossier enregistrements/AAAA-MM-JJ_HHhMMmSS contenant un
// CSV par pièce (même format que le mode sans fenêtre, lu par les scripts).
static void start_recording(App *a) {
    time_t now = time(NULL);
    char stamp[32];
    strftime(stamp, sizeof stamp, "%Y-%m-%d_%Hh%Mm%S", localtime(&now));
    mkdir("enregistrements", 0755);
    snprintf(a->rec_dir, sizeof a->rec_dir, "enregistrements/%s", stamp);
    if (mkdir(a->rec_dir, 0755) != 0) {
        perror(a->rec_dir);
        a->rec_dir[0] = '\0';
        return;
    }
    for (int k = 0; k < 2; k++) {
        Sim *s = &a->rooms[k].sim;
        char g[16], path[192];
        geom_name(s, g, sizeof g);
        snprintf(path, sizeof path, "%s/%s.csv", a->rec_dir, g);
        s->csv = fopen(path, "w");
        if (!s->csv) { perror(path); continue; }
        compute_stats(s);
        s->dt = DT_MAX;  // pas maximal (la vitesse le fait varier en dessous)
        csv_header(s->csv, s);
        csv_row(s->csv, s);
    }
    a->next_rec = a->rooms[0].sim.sim_time + REC_EVERY;
    printf("Enregistrement dans %s/ (une ligne toutes les %.0f min simulées)\n",
           a->rec_dir, REC_EVERY / 60.0);
}

static void reset_rooms(App *a) {
    stop_recording(a);
    for (int k = 0; k < 2; k++) {
        sim_reset(&a->rooms[k].sim);
        a->rooms[k].sim.dt = DT_MAX;
    }
}

// Avance les deux pièces de `span` secondes, par pas d'au plus DT_MAX
static void advance(App *a, double span) {
    int n = (int)ceil(span / DT_MAX - 1e-9);
    double dt = span / n;
    for (int s = 0; s < n; s++) {
        for (int k = 0; k < 2; k++) {
            a->rooms[k].sim.dt = dt;
            step_physics(&a->rooms[k].sim);
        }
        if (a->rec_dir[0] && a->rooms[0].sim.sim_time >= a->next_rec - 1e-6) {
            for (int k = 0; k < 2; k++) {
                Sim *sim = &a->rooms[k].sim;
                compute_stats(sim);
                if (sim->csv) { csv_row(sim->csv, sim); fflush(sim->csv); }
            }
            a->next_rec += REC_EVERY;
        }
    }
    for (int k = 0; k < 2; k++) compute_stats(&a->rooms[k].sim);
}

/* ------------------------------------------------------------------------ */
/* Dessin                                                                   */
/* ------------------------------------------------------------------------ */

static void draw_room(App *a, int k, double x) {
    Room *r = &a->rooms[k];
    Sim *s = &r->sim;
    fill_pixels(s, a->view);

    // Isothermes tous les 1 °C : on assombrit les cases où le degré entier
    // change par rapport au voisin de droite ou du dessous.
    if (a->iso && a->view == VIEW_TEMP) {
        const int G = GRID_SIZE;
        for (int i = 0; i < G * (G - 1); i++) {
            if (s->zones[i] == ZONE_EXT || (i % G) == G - 1) continue;
            int li = (int)floor(s->T[i] - T_EXT);
            int lr = s->zones[i + 1] != ZONE_EXT ? (int)floor(s->T[i + 1] - T_EXT) : li;
            int ld = s->zones[i + G] != ZONE_EXT ? (int)floor(s->T[i + G] - T_EXT) : li;
            if (li != lr || li != ld) {
                uint32_t p = s->pixels[i];
                uint32_t rr = ((p >> 16) & 0xFF) * 3 / 5, gg = ((p >> 8) & 0xFF) * 3 / 5,
                         bb = (p & 0xFF) * 3 / 5;
                s->pixels[i] = 0xFF000000u | (rr << 16) | (gg << 8) | bb;
            }
        }
    }

    SDL_UpdateTexture(r->tex, NULL, s->pixels, GRID_SIZE * sizeof(uint32_t));
    SDL_Rect dst = { L(x), L(MARGIN), L(x + PANEL) - L(x), L(MARGIN + PANEL) - L(MARGIN) };
    SDL_RenderCopy(ui.ren, r->tex, NULL, &dst);

    // Sous la pièce : son nom et la température moyenne de l'air
    char name[48], temp[32];
    if (s->mode == GEOM_KOCH) snprintf(name, sizeof name, "Flocon de Koch %d", s->depth);
    else snprintf(name, sizeof name, "Carré");
    snprintf(temp, sizeof temp, "%.1f °C", s->stat_T_air);
    fr(temp);
    double y = STRIP_Y + 12;
    fill_rect(x, y + 3, 5, 14, k ? C_KOCH : C_REF);
    double w = text(x + 13, y + 2, 2.0, C_TEXT, name);
    text(x + 13 + w + 18, y - 3, 3.0, temp_color(s->stat_T_air), temp);
}

static void draw_strip(const App *a) {
    char buf[64];
    double y = STRIP_Y + 54;

    // Barre de couleurs au centre
    const double w = 420, h = 10, x0 = (WIN_W - w) / 2;
    int p0 = L(x0), p1 = L(x0 + w);
    for (int px = p0; px < p1; px++) {
        double t = (px - p0) / (double)(p1 - p0 - 1);
        SDL_Color c = a->view == VIEW_TEMP ? temp_color(T_EXT + t * (T_DISPLAY_MAX - T_EXT))
                                           : flux_color(t * FLUX_DISPLAY_MAX);
        set_color(c);
        SDL_RenderDrawLine(ui.ren, px, L(y), px, L(y + h) - 1);
    }
    const char *lo = a->view == VIEW_TEMP ? "20 °C" : "0";
    char hi[40];
    if (a->view == VIEW_TEMP) snprintf(hi, sizeof hi, "%.0f °C", T_DISPLAY_MAX);
    else snprintf(hi, sizeof hi, "%.0f W/m² dans les murs", FLUX_DISPLAY_MAX);
    text(x0 - text_w(1.5, lo) - 10, y - 2, 1.5, C_DIM, lo);
    text(x0 + w + 10, y - 2, 1.5, C_DIM, hi);

    // Temps à gauche, aide à droite
    format_time(buf, sizeof buf, a->rooms[0].sim.sim_time);
    double tx = MARGIN + text(MARGIN, y - 4, 2.0, C_TEXT, buf);
    if (a->paused) text(tx + 14, y - 1, 1.5, C_DIM, "pause");
    else if (a->rec_dir[0]) text(tx + 14, y - 1, 1.5, C_REC, "REC");
    const char *hint = "H : aide";
    text(WIN_W - MARGIN - text_w(1.5, hint), y - 1, 1.5, C_DIM, hint);
}

static void draw_help(const App *a) {
    const double w = 520, h = 300, x = (WIN_W - w) / 2, y = (STRIP_Y - h) / 2;
    fill_rect(x, y, w, h, (SDL_Color){ 10, 12, 16, 235 });
    text(x + 24, y + 20, 2.0, C_TEXT, "Commandes");
    char speed[64];
    if (a->speed < 0) snprintf(speed, sizeof speed, "vitesse (auto)");
    else snprintf(speed, sizeof speed, "vitesse (%g s par image)", SPEEDS[a->speed]);
    const char *keys[] = { "Espace", "+ / -", "A", "0 à 4", "F", "I", "S", "R", "Q" };
    const char *what[] = { "pause", speed, "vitesse automatique", "profondeur du flocon",
                           "température / flux dans les murs", "isothermes (1 °C)",
                           a->rec_dir[0] ? "arrêter l'enregistrement" : "enregistrer (CSV)",
                           "recommencer", "quitter" };
    for (int k = 0; k < 9; k++) {
        text(x + 24, y + 56 + 24 * k, 1.5, C_KOCH, keys[k]);
        text(x + 124, y + 56 + 24 * k, 1.5, C_TEXT, what[k]);
    }
}

static void draw_frame(App *a) {
    set_color(C_BG);
    SDL_RenderClear(ui.ren);
    draw_room(a, 0, MARGIN);
    draw_room(a, 1, 2 * MARGIN + PANEL);
    draw_strip(a);
    if (a->help) draw_help(a);
}

/* ------------------------------------------------------------------------ */
/* Initialisation et boucle                                                 */
/* ------------------------------------------------------------------------ */

static int app_init(App *a, SDL_Renderer *ren, int depth) {
    memset(a, 0, sizeof *a);
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
    a->view = VIEW_TEMP;
    a->iso = 1;
    a->speed = -1;
    reset_rooms(a);
    return 1;
}

static void app_free(App *a) {
    stop_recording(a);
    for (int k = 0; k < 2; k++) {
        if (a->rooms[k].tex) SDL_DestroyTexture(a->rooms[k].tex);
        sim_free(&a->rooms[k].sim);
    }
    font_free();
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
    printf("Carré et flocon de Koch côte à côte. H : aide, S : enregistrer, Q : quitter.\n");

    int running = 1;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = 0;
            if (e.type != SDL_KEYDOWN) continue;
            SDL_Keycode k = e.key.keysym.sym;
            if (k == SDLK_q || k == SDLK_ESCAPE) running = 0;
            else if (k == SDLK_h) app.help = !app.help;
            else if (k == SDLK_f) app.view = !app.view;
            else if (k == SDLK_i) app.iso = !app.iso;
            else if (k == SDLK_SPACE) app.paused = !app.paused;
            else if (k == SDLK_a) app.speed = -1;
            else if (k == SDLK_PLUS || k == SDLK_EQUALS || k == SDLK_KP_PLUS ||
                     k == SDLK_MINUS || k == SDLK_KP_MINUS) {
                int up = !(k == SDLK_MINUS || k == SDLK_KP_MINUS);
                if (app.speed < 0) {  // depuis l'automatique : vitesse la plus proche
                    double span = frame_span(&app);
                    app.speed = 0;
                    while (app.speed < N_SPEEDS - 1 && SPEEDS[app.speed] < span) app.speed++;
                }
                if (up && app.speed < N_SPEEDS - 1) app.speed++;
                if (!up && app.speed > 0) app.speed--;
            }
            else if (k == SDLK_r || (k >= SDLK_0 && k <= SDLK_0 + KOCH_DEPTH_MAX)) {
                if (k != SDLK_r) app.rooms[1].sim.depth = k - SDLK_0;
                reset_rooms(&app);
            }
            else if (k == SDLK_s) {
                if (app.rec_dir[0]) stop_recording(&app);
                else start_recording(&app);
            }
        }

        if (!app.paused) advance(&app, frame_span(&app));
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
    if (hours > 0) advance(&app, hours * 3600.0);
    draw_frame(&app);
    SDL_RenderPresent(ren);
    int ok = SDL_SaveBMP(surf, path) == 0;
    if (!ok) fprintf(stderr, "Enregistrement de %s : %s\n", path, SDL_GetError());
    else printf("Image enregistrée : %s (t = %g h)\n", path, hours);
    app_free(&app);
    SDL_DestroyRenderer(ren);
    SDL_FreeSurface(surf);
    return !ok;
}

// Film de la fenêtre sans l'ouvrir : `frames` images brutes (BGRA, WIN_W × WIN_H,
// échelle 1) écrites sur la sortie standard, à envoyer à ffmpeg (make animation).
// Le temps simulé par image croît de 5 % de t (au moins 5 s) : quelques
// minutes pour voir la chaleur quitter la source, puis une semaine en
// quelques secondes de film.
int film_gui(int depth, int frames) {
    SDL_Surface *surf = SDL_CreateRGBSurfaceWithFormat(0, WIN_W, WIN_H, 32,
                                                       SDL_PIXELFORMAT_ARGB8888);
    SDL_Renderer *ren = surf ? SDL_CreateSoftwareRenderer(surf) : NULL;
    if (!ren) {
        fprintf(stderr, "Rendu logiciel : %s\n", SDL_GetError());
        return 1;
    }
    // Les images partent sur la sortie standard d'origine ; les messages de la
    // simulation (printf) sont renvoyés vers l'erreur standard.
    FILE *out = fdopen(dup(STDOUT_FILENO), "wb");
    if (!out || dup2(STDERR_FILENO, STDOUT_FILENO) < 0) {
        perror("Sortie standard");
        return 1;
    }
    ui.dpi = 1.0;
    static App app;
    if (!app_init(&app, ren, depth)) return 1;
    int ok = 1;
    for (int k = 0; k < frames && ok; k++) {
        if (k > 0) {
            double span = 0.05 * app.rooms[0].sim.sim_time;
            advance(&app, span < 5.0 ? 5.0 : span);
        }
        draw_frame(&app);
        SDL_RenderPresent(ren);
        for (int y = 0; y < WIN_H && ok; y++)
            ok = fwrite((char *)surf->pixels + (size_t)y * surf->pitch, 4, WIN_W, out)
                 == (size_t)WIN_W;
    }
    if (fclose(out) != 0) ok = 0;
    fprintf(stderr, "%s : t = %.1f h\n", ok ? "Film terminé" : "Écriture interrompue",
            app.rooms[0].sim.sim_time / 3600.0);
    app_free(&app);
    SDL_DestroyRenderer(ren);
    SDL_FreeSurface(surf);
    return !ok;
}
