/* Referee (Section 16, fourth pass): candidate junctions in d = 3 (spatial slice Z^2) containing p = (bd q, s-1), counted by a method that differs
 * from the deficiency-driven searches of the authors (next_order_count.py, order_counts.py) and of the third referee (referee-section16c/):
 * the fronts are placed first, and the plaquettes are then read off the faces of W = A u B.
 *
 * Usage: ./junc_faces FRONTS_DUMP n R k1 [k2 ...]
 *   FRONTS_DUMP: the fronts through the origin of Z^2 ("k x y dir x y dir ..." per line, e.g. from referee-section16c/fronts_redelmeier.c);
 *   n: the number of plaquettes |J|_C; R: the placement radius (R >= n suffices, see below); k1 k2 ...: the sizes of the further front-pairs.
 *
 * A candidate junction (P, Phi) at slice s: P a non-empty set of n plaquettes of Z^2, Phi = {(bd q, 0)} u {(G_i, l_i)} with fronts G_i (connected,
 * all degrees >= 2) of the given sizes and layers l_i in {0 = s-1, 1 = s}, such that (Definition mg:def:junction, Lemma mg:lem:junction):
 *   (a) distinct front-pairs; fronts of one layer vertex-disjoint;
 *   (b) n_J(e) = [e in A] + [e in B] + #{c in P: e in c} >= 2 for every e in E(J) = A u B u E(P)  (A, B: the fronts of layers 0, 1);
 *   (c) every front shares a vertex with E(P) (touched by a marker; no supports);
 *   (d) H_s connected: elements = link-components of P (clusters) and front-pairs, adjacent when they share a vertex.
 * Two consequences, used only to organise the search (every candidate found is checked against (a)-(d) from scratch):
 *   (F1) In Z^2 every link lies in exactly two plaquettes. By (b) a link outside W := A u B lies in 0 or 2 plaquettes of P, so P is a union of
 *        faces of W (classes of plaquettes joined across links not in W), and of bounded faces as P is finite.
 *   (F2) A cluster C has a non-empty mod-2 boundary, which lies in W by (b): C shares a link with some front. So, by (d), the graph on the fronts
 *        in which two fronts are adjacent when they share a vertex or are joined through a chain of clusters (consecutive ones sharing a vertex) is
 *        connected; such a chain has at most n plaquettes and a vertex-connected set of m unit squares spans at most m in each coordinate. Hence in
 *        a breadth-first order from bd q every further front has a vertex within L-infinity distance n of a vertex of an earlier front.
 * The search places the further fronts in every order of the sizes, each with a vertex within distance R of the vertices placed so far (each
 * front once per step: u is the smallest of its vertices in that neighbourhood), in either layer, with (a); then computes the bounded faces of W,
 * takes every union of them with n plaquettes as P, checks (a)-(d), and deduplicates.
 * Output: the number of candidate junctions, by the layers of the further fronts (relative to p: 'same' = layer s-1, 'other' = layer s).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define G 192
#define OFF 96
#define MAXK 16
#define MAXF 6
static inline int vid(int x, int y) { return (x + OFF) * G + (y + OFF); }
static inline int eid(int x, int y, int dir) { return vid(x, y) * 2 + dir; }

typedef struct { int k, nv; int ex[MAXK], ey[MAXK], ed[MAXK]; int vx[2 * MAXK], vy[2 * MAXK]; } Tmpl;
static Tmpl *T; static int nT = 0, capT = 0;
static int *bysize[64], nby[64], capby[64];
static int want[64];                          /* the sizes to load from the dump */

static int n_plaq, R, nsz, sizes[MAXF], order[MAXF];
static unsigned char *inA, *inB;              /* edge flags */
static int *occ[2];                            /* vertex occupancy per layer */
static int *nbc;                               /* number of placed vertices within distance R */
static int pl_t[MAXF], pl_x[MAXF], pl_y[MAXF], pl_l[MAXF], npl = 0;
static long long nconf = 0, nfaces = 0;

/* ---------- templates ---------- */
static int tmpl_check(Tmpl *t) {                /* connected, all degrees >= 2, contains the origin, distinct edges */
    int deg[2 * MAXK] = {0}, par[2 * MAXK];
    for (int i = 0; i < t->k; i++) for (int j = i + 1; j < t->k; j++)
        if (t->ex[i] == t->ex[j] && t->ey[i] == t->ey[j] && t->ed[i] == t->ed[j]) return 0;
    int org = 0;
    for (int v = 0; v < t->nv; v++) { par[v] = v; if (t->vx[v] == 0 && t->vy[v] == 0) org = 1; }
    if (!org) return 0;
    for (int i = 0; i < t->k; i++) {
        int x2 = t->ex[i] + (t->ed[i] == 0), y2 = t->ey[i] + (t->ed[i] == 1), a = -1, b = -1;
        for (int v = 0; v < t->nv; v++) { if (t->vx[v] == t->ex[i] && t->vy[v] == t->ey[i]) a = v; if (t->vx[v] == x2 && t->vy[v] == y2) b = v; }
        deg[a]++; deg[b]++;
        while (par[a] != a) a = par[a]; while (par[b] != b) b = par[b]; par[a] = b;
    }
    int roots = 0;
    for (int v = 0; v < t->nv; v++) { if (deg[v] < 2) return 0; if (par[v] == v) roots++; }
    return roots == 1;
}
static void load(const char *path) {
    FILE *f = fopen(path, "r"); if (!f) { perror(path); exit(1); }
    memset(nby, 0, sizeof nby); memset(capby, 0, sizeof capby);
    int k;
    while (fscanf(f, "%d", &k) == 1) {
        if (k < 1 || k > MAXK) { fprintf(stderr, "front size %d not supported\n", k); exit(1); }
        if (!want[k]) { int a, b, c; for (int i = 0; i < k; i++) if (fscanf(f, "%d %d %d", &a, &b, &c) != 3) { fprintf(stderr, "bad dump\n"); exit(1); } continue; }
        if (nT == capT) { capT = capT ? 2 * capT : 4096; T = realloc(T, (size_t)capT * sizeof(Tmpl)); if (!T) { fprintf(stderr, "out of memory\n"); exit(1); } }
        Tmpl *t = &T[nT]; t->k = k; t->nv = 0;
        for (int i = 0; i < k; i++) {
            if (fscanf(f, "%d %d %d", &t->ex[i], &t->ey[i], &t->ed[i]) != 3) { fprintf(stderr, "bad dump\n"); exit(1); }
            int xs[2] = {t->ex[i], t->ex[i] + (t->ed[i] == 0)}, ys[2] = {t->ey[i], t->ey[i] + (t->ed[i] == 1)};
            for (int s = 0; s < 2; s++) {
                int found = 0;
                for (int v = 0; v < t->nv; v++) if (t->vx[v] == xs[s] && t->vy[v] == ys[s]) found = 1;
                if (!found) { t->vx[t->nv] = xs[s]; t->vy[t->nv] = ys[s]; t->nv++; }
            }
        }
        if (!tmpl_check(t)) { fprintf(stderr, "template %d is not a front through the origin\n", nT); exit(1); }
        if (nby[k] == capby[k]) { capby[k] = capby[k] ? 2 * capby[k] : 1024; bysize[k] = realloc(bysize[k], (size_t)capby[k] * sizeof(int)); }
        bysize[k][nby[k]++] = nT; nT++;
    }
    fclose(f);
}

/* ---------- deduplication: 128-bit fingerprints ---------- */
#define HS (1 << 22)
static uint64_t *H1, *H2; static unsigned char *Hu; static long long nsol = 0;
static uint64_t mix(uint64_t h) { h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33; return h; }
static int insert(uint64_t a, uint64_t b) {
    uint64_t i = a & (HS - 1);
    while (Hu[i]) { if (H1[i] == a && H2[i] == b) return 0; i = (i + 1) & (HS - 1); }
    Hu[i] = 1; H1[i] = a; H2[i] = b; nsol++;
    if (nsol > HS / 2) { fprintf(stderr, "hash table full\n"); exit(1); }
    return 1;
}
static int cmp64(const void *a, const void *b) { uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b; return x < y ? -1 : x > y; }
static int cmpint(const void *a, const void *b) { return *(const int *)a - *(const int *)b; }

/* ---------- statistics by layer pattern ---------- */
#define MAXPAT 256
static char patname[MAXPAT][128]; static long long patcnt[MAXPAT]; static int npat = 0;
static void addpat(const char *s) {
    for (int i = 0; i < npat; i++) if (!strcmp(patname[i], s)) { patcnt[i]++; return; }
    strcpy(patname[npat], s); patcnt[npat++] = 1;
}

/* ---------- evaluation of a full placement ---------- */
static int *sqlab; static int *pi_cnt; static unsigned char *vP;
static int fx[MAXF + 1][2 * MAXK], fy[MAXF + 1][2 * MAXK], fnv[MAXF + 1];   /* vertices of bd q and the further fronts */

static int find(int *par, int a) { while (par[a] != a) { par[a] = par[par[a]]; a = par[a]; } return a; }

static void check_P(int *Psq, int np, int x0, int y0, int x1, int y1) {
    /* (b): n_J(e) >= 2 on E(J) */
    int Pe[64], nPe = 0;
    for (int i = 0; i < np; i++) {
        int x = Psq[i] / G - OFF, y = Psq[i] % G - OFF;
        int es[4] = {eid(x, y, 0), eid(x, y + 1, 0), eid(x, y, 1), eid(x + 1, y, 1)};
        for (int j = 0; j < 4; j++) { if (pi_cnt[es[j]] == 0) Pe[nPe++] = es[j]; pi_cnt[es[j]]++; }
    }
    int ok = 1;
    for (int i = 0; i < nPe && ok; i++) if (inA[Pe[i]] + inB[Pe[i]] + pi_cnt[Pe[i]] < 2) ok = 0;
    for (int f = 0; f <= npl && ok; f++) {                     /* the links of the fronts (bd q is front 0) */
        const Tmpl *t = (f == 0) ? NULL : &T[pl_t[f - 1]];
        int k = (f == 0) ? 4 : t->k;
        for (int j = 0; j < k && ok; j++) {
            int e;
            if (f == 0) { int qe[4] = {eid(0, 0, 0), eid(0, 1, 0), eid(0, 0, 1), eid(1, 0, 1)}; e = qe[j]; }
            else e = eid(t->ex[j] + pl_x[f - 1], t->ey[j] + pl_y[f - 1], t->ed[j]);
            if (inA[e] + inB[e] + pi_cnt[e] < 2) ok = 0;
        }
    }
    /* (c): every front touches E(P) */
    for (int i = 0; i < np; i++) { int x = Psq[i] / G - OFF, y = Psq[i] % G - OFF; vP[vid(x, y)] = vP[vid(x + 1, y)] = vP[vid(x, y + 1)] = vP[vid(x + 1, y + 1)] = 1; }
    for (int f = 0; f <= npl && ok; f++) {
        int touch = 0;
        for (int v = 0; v < fnv[f]; v++) if (vP[vid(fx[f][v], fy[f][v])]) touch = 1;
        if (!touch) ok = 0;
    }
    /* (d): connectivity of clusters and fronts through common vertices */
    if (ok) {
        int par[64 + MAXF + 1], ne = np + npl + 1;
        for (int i = 0; i < ne; i++) par[i] = i;
        for (int i = 0; i < np; i++) for (int j = i + 1; j < np; j++) {      /* clusters: plaquettes sharing a link */
            int xi = Psq[i] / G - OFF, yi = Psq[i] % G - OFF, xj = Psq[j] / G - OFF, yj = Psq[j] % G - OFF;
            if ((xi == xj && abs(yi - yj) == 1) || (yi == yj && abs(xi - xj) == 1)) par[find(par, i)] = find(par, j);
        }
        /* an element per vertex: plaquettes (their clusters) and fronts */
        static int owner[G * G]; static int touched[4096]; int nt = 0;
        for (int i = 0; i < np; i++) {
            int x = Psq[i] / G - OFF, y = Psq[i] % G - OFF, vs[4] = {vid(x, y), vid(x + 1, y), vid(x, y + 1), vid(x + 1, y + 1)};
            for (int j = 0; j < 4; j++) { int v = vs[j]; if (owner[v] == 0) { owner[v] = i + 1; touched[nt++] = v; } else par[find(par, owner[v] - 1)] = find(par, i); }
        }
        for (int f = 0; f <= npl; f++) for (int v0 = 0; v0 < fnv[f]; v0++) {
            int v = vid(fx[f][v0], fy[f][v0]), el = np + f;
            if (owner[v] == 0) { owner[v] = el + 1; touched[nt++] = v; } else par[find(par, owner[v] - 1)] = find(par, el);
        }
        for (int i = 0; i < nt; i++) owner[touched[i]] = 0;
        int r0 = find(par, 0);
        for (int i = 1; i < ne; i++) if (find(par, i) != r0) { ok = 0; break; }
    }
    for (int i = 0; i < np; i++) { int x = Psq[i] / G - OFF, y = Psq[i] % G - OFF; vP[vid(x, y)] = vP[vid(x + 1, y)] = vP[vid(x, y + 1)] = vP[vid(x + 1, y + 1)] = 0; }
    for (int i = 0; i < nPe; i++) pi_cnt[Pe[i]] = 0;
    if (!ok) return;
    /* a canonical fingerprint: the further fronts (sorted edge lists with layer) and the plaquettes */
    uint64_t fh[MAXF];
    for (int f = 0; f < npl; f++) {
        const Tmpl *t = &T[pl_t[f]]; int es[MAXK];
        for (int j = 0; j < t->k; j++) es[j] = eid(t->ex[j] + pl_x[f], t->ey[j] + pl_y[f], t->ed[j]);
        qsort(es, t->k, sizeof(int), cmpint);
        uint64_t h = 1469598103934665603ULL ^ (uint64_t)pl_l[f];
        for (int j = 0; j < t->k; j++) h = mix(h ^ (uint64_t)es[j]) + 0x9e3779b97f4a7c15ULL;
        fh[f] = h;
    }
    qsort(fh, npl, sizeof(uint64_t), cmp64);
    int ps[64]; memcpy(ps, Psq, np * sizeof(int)); qsort(ps, np, sizeof(int), cmpint);
    uint64_t a = 0x12345678ULL, b = 0x9abcdef0ULL;
    for (int f = 0; f < npl; f++) { a = mix(a ^ fh[f]); b = mix(b + 3 * fh[f] + 7); }
    for (int i = 0; i < np; i++) { a = mix(a ^ (uint64_t)(ps[i] + 1000003)); b = mix(b + 5 * (uint64_t)ps[i] + 11); }
    if (insert(a, b)) {
        /* layer pattern: sorted (size, same/other) */
        char s[128] = ""; int idx[MAXF];
        for (int f = 0; f < npl; f++) idx[f] = T[pl_t[f]].k * 2 + pl_l[f];
        qsort(idx, npl, sizeof(int), cmpint);
        for (int f = 0; f < npl; f++) { char w[24]; sprintf(w, "%d%s ", idx[f] / 2, (idx[f] % 2) ? "o" : "s"); strcat(s, w); }
        addpat(s);
    }
}

static int faceP[256], nfaceP, faceArea[256], faceStart[256], faceSq[8192];
static void choose_faces(int fi, int area, int *Psq, int np, int x0, int y0, int x1, int y1) {
    if (area == n_plaq) { check_P(Psq, np, x0, y0, x1, y1); return; }
    if (fi == nfaceP) return;
    choose_faces(fi + 1, area, Psq, np, x0, y0, x1, y1);
    if (area + faceArea[fi] <= n_plaq) {
        for (int j = 0; j < faceArea[fi]; j++) Psq[np + j] = faceSq[faceStart[fi] + j];
        choose_faces(fi + 1, area + faceArea[fi], Psq, np + faceArea[fi], x0, y0, x1, y1);
    }
}

static void evaluate(void) {
    nconf++;
    /* bd q and the fronts: vertices, the edge flags */
    int qe[4] = {eid(0, 0, 0), eid(0, 1, 0), eid(0, 0, 1), eid(1, 0, 1)};
    int sizeA = 4, sizeB = 0, inter = 0;
    for (int j = 0; j < 4; j++) inA[qe[j]] = 1;
    for (int f = 0; f < npl; f++) {
        const Tmpl *t = &T[pl_t[f]];
        for (int j = 0; j < t->k; j++) {
            int e = eid(t->ex[j] + pl_x[f], t->ey[j] + pl_y[f], t->ed[j]);
            if (pl_l[f] == 0) inA[e] = 1; else inB[e] = 1;
        }
        if (pl_l[f] == 0) sizeA += t->k; else sizeB += t->k;
    }
    int x0 = 0, y0 = 0, x1 = 1, y1 = 1;
    fnv[0] = 4; fx[0][0] = 0; fy[0][0] = 0; fx[0][1] = 1; fy[0][1] = 0; fx[0][2] = 0; fy[0][2] = 1; fx[0][3] = 1; fy[0][3] = 1;
    for (int f = 0; f < npl; f++) {
        const Tmpl *t = &T[pl_t[f]]; fnv[f + 1] = t->nv;
        for (int v = 0; v < t->nv; v++) {
            int x = t->vx[v] + pl_x[f], y = t->vy[v] + pl_y[f]; fx[f + 1][v] = x; fy[f + 1][v] = y;
            if (x < x0) x0 = x; if (x > x1) x1 = x; if (y < y0) y0 = y; if (y > y1) y1 = y;
        }
    }
    /* |A n B| (A, B are disjoint unions of their fronts) */
    for (int f = 0; f < npl; f++) if (pl_l[f] == 1) {
        const Tmpl *t = &T[pl_t[f]];
        for (int j = 0; j < t->k; j++) if (inA[eid(t->ex[j] + pl_x[f], t->ey[j] + pl_y[f], t->ed[j])]) inter++;
    }
    if (sizeA + sizeB - 2 * inter <= 4 * n_plaq) {
        /* faces of W in the box of squares [x0-1, x1] x [y0-1, y1]; the ring is unbounded */
        int W = x1 - x0 + 2, H = y1 - y0 + 2, nsq = W * H;
        static int queue[1 << 16];
        for (int i = 0; i < nsq; i++) sqlab[i] = 0;
        int qh = 0, qt = 0;
        for (int i = 0; i < W; i++) for (int j = 0; j < H; j++)
            if (i == 0 || j == 0 || i == W - 1 || j == H - 1) { sqlab[i * H + j] = -1; queue[qt++] = i * H + j; }
        int lab = 0; nfaceP = 0; int nfs = 0;
        for (int pass = 0; pass < 2; pass++) {
            if (pass == 1) {
                for (int s0 = 0; s0 < nsq; s0++) if (sqlab[s0] == 0) {
                    lab++; sqlab[s0] = lab; queue[qt++] = s0; faceStart[nfaceP] = nfs; faceArea[nfaceP] = 0;
                    while (qh < qt) {
                        int s = queue[qh++], i = s / H, j = s % H, x = i + x0 - 1, y = j + y0 - 1;
                        faceSq[nfs++] = vid(x, y); faceArea[nfaceP]++;
                        int nb[4][3] = {{i, j - 1, eid(x, y, 0)}, {i, j + 1, eid(x, y + 1, 0)}, {i - 1, j, eid(x, y, 1)}, {i + 1, j, eid(x + 1, y, 1)}};
                        for (int m = 0; m < 4; m++) {
                            int ii = nb[m][0], jj = nb[m][1], e = nb[m][2];
                            if (ii < 0 || jj < 0 || ii >= W || jj >= H) continue;
                            if (inA[e] || inB[e]) continue;
                            if (sqlab[ii * H + jj] == 0) { sqlab[ii * H + jj] = lab; queue[qt++] = ii * H + jj; }
                        }
                    }
                    if (faceArea[nfaceP] <= n_plaq) nfaceP++; else nfs = faceStart[nfaceP];
                }
                break;
            }
            while (qh < qt) {
                int s = queue[qh++], i = s / H, j = s % H, x = i + x0 - 1, y = j + y0 - 1;
                int nb[4][3] = {{i, j - 1, eid(x, y, 0)}, {i, j + 1, eid(x, y + 1, 0)}, {i - 1, j, eid(x, y, 1)}, {i + 1, j, eid(x + 1, y, 1)}};
                for (int m = 0; m < 4; m++) {
                    int ii = nb[m][0], jj = nb[m][1], e = nb[m][2];
                    if (ii < 0 || jj < 0 || ii >= W || jj >= H) continue;
                    if (inA[e] || inB[e]) continue;
                    if (sqlab[ii * H + jj] == 0) { sqlab[ii * H + jj] = -1; queue[qt++] = ii * H + jj; }
                }
            }
            qh = qt = 0;
        }
        int tot = 0; for (int i = 0; i < nfaceP; i++) tot += faceArea[i];
        if (tot >= n_plaq) { nfaces++; int Psq[64]; choose_faces(0, 0, Psq, 0, x0, y0, x1, y1); }
    }
    for (int j = 0; j < 4; j++) inA[qe[j]] = 0;
    for (int f = 0; f < npl; f++) {
        const Tmpl *t = &T[pl_t[f]];
        for (int j = 0; j < t->k; j++) { int e = eid(t->ex[j] + pl_x[f], t->ey[j] + pl_y[f], t->ed[j]); inA[e] = 0; inB[e] = 0; }
    }
}

/* ---------- placement ---------- */
static int bx0, by0, bx1, by1;     /* bounding box of the placed vertices */
static void nb_add(int x, int y, int d) {
    for (int i = x - R; i <= x + R; i++) for (int j = y - R; j <= y + R; j++) nbc[vid(i, j)] += d;
}
static void place(int level) {
    if (level == nsz) { evaluate(); return; }
    int k = order[level];
    int sx0 = bx0, sy0 = by0, sx1 = bx1, sy1 = by1;
    for (int ux = sx0 - R; ux <= sx1 + R; ux++) for (int uy = sy0 - R; uy <= sy1 + R; uy++) {
        int u = vid(ux, uy);
        if (!nbc[u]) continue;
        for (int ti = 0; ti < nby[k]; ti++) {
            const Tmpl *t = &T[bysize[k][ti]];
            int canon = 1;
            for (int v = 0; v < t->nv; v++) { int w = vid(t->vx[v] + ux, t->vy[v] + uy); if (nbc[w] && w < u) { canon = 0; break; } }
            if (!canon) continue;
            for (int l = 0; l < 2; l++) {
                int free_ = 1;
                for (int v = 0; v < t->nv; v++) if (occ[l][vid(t->vx[v] + ux, t->vy[v] + uy)]) { free_ = 0; break; }
                if (!free_) continue;
                pl_t[npl] = bysize[k][ti]; pl_x[npl] = ux; pl_y[npl] = uy; pl_l[npl] = l; npl++;
                for (int v = 0; v < t->nv; v++) occ[l][vid(t->vx[v] + ux, t->vy[v] + uy)]++;
                int last = (level + 1 == nsz);
                int ox0 = bx0, oy0 = by0, ox1 = bx1, oy1 = by1;
                if (!last) for (int v = 0; v < t->nv; v++) {
                    int x = t->vx[v] + ux, y = t->vy[v] + uy; nb_add(x, y, 1);
                    if (x < bx0) bx0 = x; if (x > bx1) bx1 = x; if (y < by0) by0 = y; if (y > by1) by1 = y;
                }
                place(level + 1);
                if (!last) for (int v = 0; v < t->nv; v++) nb_add(t->vx[v] + ux, t->vy[v] + uy, -1);
                bx0 = ox0; by0 = oy0; bx1 = ox1; by1 = oy1;
                for (int v = 0; v < t->nv; v++) occ[l][vid(t->vx[v] + ux, t->vy[v] + uy)]--;
                npl--;
            }
        }
    }
}

static void permute(int i, int *used, int *buf) {     /* every distinct order of the multiset of sizes */
    if (i == nsz) { memcpy(order, buf, sizeof(int) * nsz); place(0); return; }
    for (int j = 0; j < nsz; j++) {
        if (used[j]) continue;
        int dup = 0; for (int m = 0; m < j; m++) if (!used[m] && sizes[m] == sizes[j]) dup = 1;
        if (dup) continue;
        used[j] = 1; buf[i] = sizes[j]; permute(i + 1, used, buf); used[j] = 0;
    }
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: %s DUMP n R k1 [k2 ...]\n", argv[0]); return 1; }
    n_plaq = atoi(argv[2]); R = atoi(argv[3]); nsz = argc - 4;
    if (nsz > MAXF) { fprintf(stderr, "too many fronts\n"); return 1; }
    for (int i = 0; i < nsz; i++) { sizes[i] = atoi(argv[4 + i]); want[sizes[i]] = 1; }
    load(argv[1]);
    for (int i = 0; i < nsz; i++) fprintf(stderr, "%s", "");
    for (int i = 0; i < nsz; i++) for (int j = i + 1; j < nsz; j++) if (sizes[j] < sizes[i]) { int t = sizes[i]; sizes[i] = sizes[j]; sizes[j] = t; }
    inA = calloc(G * G * 2, 1); inB = calloc(G * G * 2, 1); occ[0] = calloc(G * G, sizeof(int)); occ[1] = calloc(G * G, sizeof(int));
    nbc = calloc(G * G, sizeof(int)); sqlab = calloc(1 << 16, sizeof(int)); pi_cnt = calloc(G * G * 2, sizeof(int)); vP = calloc(G * G, 1);
    H1 = malloc(sizeof(uint64_t) * HS); H2 = malloc(sizeof(uint64_t) * HS); Hu = calloc(HS, 1);
    /* bd q in layer 0 */
    int qv[4][2] = {{0, 0}, {1, 0}, {0, 1}, {1, 1}};
    for (int i = 0; i < 4; i++) { occ[0][vid(qv[i][0], qv[i][1])]++; nb_add(qv[i][0], qv[i][1], 1); }
    bx0 = 0; by0 = 0; bx1 = 1; by1 = 1;
    int used[MAXF] = {0}, buf[MAXF];
    permute(0, used, buf);
    printf("n=%d R=%d sizes=", n_plaq, R);
    for (int i = 0; i < nsz; i++) printf("%d%s", sizes[i], i + 1 < nsz ? "," : "");
    printf(": %lld junctions  [", nsol);
    for (int i = 0; i < npat; i++) printf("%s%s: %lld", i ? "; " : "", patname[i], patcnt[i]);
    printf("]  (%lld placements, %lld with faces)\n", nconf, nfaces);
    return 0;
}
