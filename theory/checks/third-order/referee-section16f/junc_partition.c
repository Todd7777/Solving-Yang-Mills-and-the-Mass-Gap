/* Referee (Section 16, sixth pass): candidate junctions containing p = (bd q, s-1), counted from the definitions by a partition search.
 * Usage: junc_partition FRONTS_DUMP D NMAX SIZES [noprune]      (SIZES: sizes of the further fronts, e.g. 4,10 or 4,4,4; FRONTS_DUMP from fronts_enum)
 *
 * A candidate junction is J = (P, Phi): P a finite set of plaquettes of Z^D, Phi a set of front-pairs (G, layer), layer 0 = s-1, 1 = s, containing
 * p0 = (bd q, 0), q the plaquette at the origin in the plane (0,1), with further fronts of exactly the sizes SIZES, and (Definition mg:def:junction,
 * Lemma mg:lem:junction):
 *   (1) every G a front (connected, all degrees >= 2): the translates of the classes of the dump;
 *   (2) distinct fronts of one layer vertex-disjoint;
 *   (3) n_J(e) = [e in E_0] + [e in E_1] + #{c in P: e in c} >= 2 for every link e of E(J);
 *   (4) every front shares a vertex with a plaquette of P (so P is not empty);
 *   (5) H_s connected: nodes = link-components of P and the front-pairs, adjacency = a common vertex (front-front adjacency included).
 * Counted by n = |P| <= NMAX.
 *
 * Search: a partition of the set of junctions extending the current state X (a subset of J), each junction reached exactly once (a set of
 * fingerprints checks this).  If X has a link of coverage one, let l be the smallest: J contains a further element through l, and branch i
 * means "c_i is the first candidate through l that J contains", so c_1..c_{i-1} are excluded below.  If X has no such link, either J = X, or
 * J contains an element sharing a vertex with X (the elements of J are connected by common vertices, by (5)), and branch i means "c_i is the
 * first such candidate in J".  Exclusions are kept in a LIFO hash set.  This is neither the authors' search (all elements through l decided at
 * once, exclusion of plaquettes only, memo of closed states) nor the third referee's (one element per step, memo of all states).
 * Prune (optional, a necessary condition): every link of coverage one needs a further element through it, and the further elements have at
 * most 4 (NMAX - |P|) + (sum of the remaining front sizes) links.  It is applied to the state and to each candidate before branching (a candidate
 * that fails cannot lie in any J extending X, so it need not be excluded either).  With "noprune" it is off. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef unsigned long long u64;

static int D, NMAX, R, W, NPL, PRUNE = 1;
static long stride[4], NVERT, NLINK, NPLAQ, ORIGIN;
static int pi_[3], pj_[3];
typedef struct { int k, nv; long loff[16]; int ldir[16]; long voff[40]; } FClass;
static FClass *cls; static int ncls;
static int rem[64], maxsize;
typedef struct { int c; long off; } TL;
static TL *thruL[64][3]; static int nthruL[64][3];
static TL *thruV[64]; static int nthruV[64];
static unsigned char *cov, *occ[2], *inP, *exP, *vmark, *candmark;
static long Pl[64]; static int nP;
typedef struct { int c; long t; int lay; } Fr;
static Fr Fs[32]; static int nF;
static long long found[64], nodes = 0, boundary_hits = 0, dups = 0;

/* ---------- LIFO hash sets (linear probing; deletions in reverse order of insertion just clear the slot) ---------- */
typedef struct { u64 *tab; long mask, n; } HSet;
static void hs_init(HSet *h, int bits) { h->tab = calloc(1L << bits, 8); h->mask = (1L << bits) - 1; h->n = 0; }
static u64 mix(u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static int hs_has(HSet *h, u64 key) { long i = mix(key) & h->mask; while (h->tab[i]) { if (h->tab[i] == key) return 1; i = (i + 1) & h->mask; } return 0; }
static int hs_add(HSet *h, u64 key) {
    long i = mix(key) & h->mask;
    while (h->tab[i]) { if (h->tab[i] == key) return 0; i = (i + 1) & h->mask; }
    h->tab[i] = key; h->n++;
    if (2 * h->n > h->mask) { fprintf(stderr, "hash set full\n"); exit(2); }
    return 1;
}
static void hs_del(HSet *h, u64 key) { long i = mix(key) & h->mask; while (h->tab[i] != key) i = (i + 1) & h->mask; h->tab[i] = 0; h->n--; }
static HSet EXF, CAND, FOUND;

/* ---------- lattice ---------- */
static int coord(long v, int i) { return (int)((v / stride[i]) % W) - R; }
static int inbox(long v, int margin) { for (int i = 0; i < D; i++) { int c = coord(v, i); if (c < -R + margin || c > R - margin) return 0; } return 1; }
static long plinks(long p, long *out) {
    long x = p / NPL; int pl = (int)(p % NPL), i = pi_[pl], j = pj_[pl];
    out[0] = x * D + i; out[1] = (x + stride[i]) * D + j; out[2] = (x + stride[j]) * D + i; out[3] = x * D + j; return 4;
}
static void pverts(long p, long *out) {
    long x = p / NPL; int pl = (int)(p % NPL), i = pi_[pl], j = pj_[pl];
    out[0] = x; out[1] = x + stride[i]; out[2] = x + stride[j]; out[3] = x + stride[i] + stride[j];
}
static int plane_of(int i, int j) { for (int t = 0; t < NPL; t++) if (pi_[t] == i && pj_[t] == j) return t; return -1; }
static int plaqs_through(long l, long *out) {
    long x = l / D; int k = (int)(l % D), m = 0;
    for (int a = 0; a < D; a++) {
        if (a == k) continue;
        int i = k < a ? k : a, j = k < a ? a : k, pl = plane_of(i, j);
        out[m++] = x * NPL + pl; out[m++] = (x - stride[a]) * NPL + pl;
    }
    return m;
}
static int plaqs_at(long v, long *out) {
    int m = 0;
    for (int pl = 0; pl < NPL; pl++) {
        int i = pi_[pl], j = pj_[pl];
        out[m++] = v * NPL + pl; out[m++] = (v - stride[i]) * NPL + pl; out[m++] = (v - stride[j]) * NPL + pl;
        out[m++] = (v - stride[i] - stride[j]) * NPL + pl;
    }
    return m;
}
static u64 fkey(int c, long t, int lay) { return (((u64)(c + 1)) << 34) | ((u64)t << 1) | (u64)lay; }
static int front_fits(int c, long t) {         /* all vertices well inside the box */
    for (int a = 0; a < cls[c].nv; a++) if (!inbox(t + cls[c].voff[a], 2)) return 0;
    return 1;
}

/* ---------- state ---------- */
static int capacity(void) { int s = 4 * (NMAX - nP); for (int k = 0; k <= maxsize; k++) if (rem[k] > 0) s += k * rem[k]; return s; }
static int front_ok_layer(int c, long t, int lay) { for (int a = 0; a < cls[c].nv; a++) if (occ[lay][t + cls[c].voff[a]]) return 0; return 1; }
static void add_plaq(long p) { long L[4]; plinks(p, L); for (int a = 0; a < 4; a++) cov[L[a]]++; inP[p] = 1; Pl[nP++] = p; }
static void del_plaq(long p) { long L[4]; plinks(p, L); for (int a = 0; a < 4; a++) cov[L[a]]--; inP[p] = 0; nP--; }
static void add_front(int c, long t, int lay, int counted) {
    for (int a = 0; a < cls[c].k; a++) cov[(t + cls[c].loff[a]) * D + cls[c].ldir[a]]++;
    for (int a = 0; a < cls[c].nv; a++) occ[lay][t + cls[c].voff[a]] = 1;
    Fs[nF].c = c; Fs[nF].t = t; Fs[nF].lay = lay; nF++;
    if (counted) rem[cls[c].k]--;
}
static void del_front(int counted) {
    nF--; int c = Fs[nF].c; long t = Fs[nF].t; int lay = Fs[nF].lay;
    for (int a = 0; a < cls[c].k; a++) cov[(t + cls[c].loff[a]) * D + cls[c].ldir[a]]--;
    for (int a = 0; a < cls[c].nv; a++) occ[lay][t + cls[c].voff[a]] = 0;
    if (counted) rem[cls[c].k]++;
}
/* deficient links of X: smallest one, and how many */
static int deficient(long *lmin) {
    int nd = 0; *lmin = -1;
    long L[4];
    for (int a = 0; a < nP; a++) { plinks(Pl[a], L); for (int b = 0; b < 4; b++) if (cov[L[b]] == 1) { nd++; if (*lmin < 0 || L[b] < *lmin) *lmin = L[b]; } }
    for (int a = 0; a < nF; a++) {
        FClass *C = &cls[Fs[a].c];
        for (int b = 0; b < C->k; b++) { long l = (Fs[a].t + C->loff[b]) * D + C->ldir[b]; if (cov[l] == 1) { nd++; if (*lmin < 0 || l < *lmin) *lmin = l; } }
    }
    /* a link of coverage one lies in exactly one element, so it was counted once */
    return nd;
}
/* deficient count after adding a plaquette / a front (prune of candidates) */
static int def_after_plaq(int nd, long p) { long L[4]; plinks(p, L); for (int a = 0; a < 4; a++) { if (cov[L[a]] == 1) nd--; else if (cov[L[a]] == 0) nd++; } return nd; }
static int def_after_front(int nd, int c, long t) {
    for (int a = 0; a < cls[c].k; a++) { long l = (t + cls[c].loff[a]) * D + cls[c].ldir[a]; if (cov[l] == 1) nd--; else if (cov[l] == 0) nd++; }
    return nd;
}

/* ---------- final check of (2)-(5) from scratch; (1) by construction, (3) re-checked ---------- */
static int final_check(void) {
    if (nP == 0) return 0;
    for (int k = 0; k <= maxsize; k++) if (rem[k] != 0) return 0;
    long L[4], V[4];
    /* (3) */
    for (int a = 0; a < nP; a++) { plinks(Pl[a], L); for (int b = 0; b < 4; b++) if (cov[L[b]] < 2) return 0; }
    for (int a = 0; a < nF; a++) { FClass *C = &cls[Fs[a].c]; for (int b = 0; b < C->k; b++) if (cov[(Fs[a].t + C->loff[b]) * D + C->ldir[b]] < 2) return 0; }
    /* (2) */
    for (int a = 0; a < nF; a++) for (int b = a + 1; b < nF; b++) {
        if (Fs[a].lay != Fs[b].lay) continue;
        if (Fs[a].c == Fs[b].c && Fs[a].t == Fs[b].t) return 0;
        FClass *A = &cls[Fs[a].c], *B = &cls[Fs[b].c];
        for (int x = 0; x < A->nv; x++) for (int y = 0; y < B->nv; y++) if (Fs[a].t + A->voff[x] == Fs[b].t + B->voff[y]) return 0;
    }
    /* (4) */
    for (int a = 0; a < nP; a++) { pverts(Pl[a], V); for (int b = 0; b < 4; b++) vmark[V[b]] = 1; }
    int ok = 1;
    for (int a = 0; a < nF && ok; a++) { FClass *C = &cls[Fs[a].c]; int t = 0; for (int b = 0; b < C->nv; b++) if (vmark[Fs[a].t + C->voff[b]]) { t = 1; break; } if (!t) ok = 0; }
    for (int a = 0; a < nP; a++) { pverts(Pl[a], V); for (int b = 0; b < 4; b++) vmark[V[b]] = 0; }
    if (!ok) return 0;
    /* (5): union-find over plaquettes (common link) then nodes (common vertex) */
    int nn = nP + nF, par[96];
    for (int a = 0; a < nn; a++) par[a] = a;
#define FIND(x) ({ int _r = (x); while (par[_r] != _r) _r = par[_r]; _r; })
    for (int a = 0; a < nP; a++) for (int b = a + 1; b < nP; b++) {
        long L2[4]; plinks(Pl[a], L); plinks(Pl[b], L2); int sh = 0;
        for (int x = 0; x < 4; x++) for (int y = 0; y < 4; y++) if (L[x] == L2[y]) sh = 1;
        if (sh) { int ra = FIND(a), rb = FIND(b); par[ra] = rb; }
    }
    /* vertex sets of nodes: a marker's vertices are its plaquettes' vertices, so adjacency of nodes = some plaquette/front pair sharing a vertex */
    for (int a = 0; a < nn; a++) for (int b = a + 1; b < nn; b++) {
        long Va[40], Vb[40]; int na, nb;
        if (a < nP) { pverts(Pl[a], Va); na = 4; } else { FClass *C = &cls[Fs[a - nP].c]; na = C->nv; for (int x = 0; x < na; x++) Va[x] = Fs[a - nP].t + C->voff[x]; }
        if (b < nP) { pverts(Pl[b], Vb); nb = 4; } else { FClass *C = &cls[Fs[b - nP].c]; nb = C->nv; for (int x = 0; x < nb; x++) Vb[x] = Fs[b - nP].t + C->voff[x]; }
        int sh = 0;
        for (int x = 0; x < na && !sh; x++) for (int y = 0; y < nb; y++) if (Va[x] == Vb[y]) { sh = 1; break; }
        if (sh) { int ra = FIND(a), rb = FIND(b); par[ra] = rb; }
    }
    int r0 = FIND(0);
    for (int a = 1; a < nn; a++) if (FIND(a) != r0) return 0;
    return 1;
}
static void record(void) {
    /* fingerprint: sorted plaquettes and fronts, hashed */
    long ps[64]; int np = nP; memcpy(ps, Pl, sizeof(long) * nP);
    for (int a = 0; a < np; a++) for (int b = a + 1; b < np; b++) if (ps[b] < ps[a]) { long t = ps[a]; ps[a] = ps[b]; ps[b] = t; }
    u64 fs[32]; int nf = nF;
    for (int a = 0; a < nf; a++) fs[a] = fkey(Fs[a].c, Fs[a].t, Fs[a].lay);
    for (int a = 0; a < nf; a++) for (int b = a + 1; b < nf; b++) if (fs[b] < fs[a]) { u64 t = fs[a]; fs[a] = fs[b]; fs[b] = t; }
    u64 h = 1469598103934665603ULL;
    for (int a = 0; a < np; a++) h = mix(h ^ (u64)ps[a]) + 0x9e3779b97f4a7c15ULL;
    h = mix(h ^ 0xabcdefULL);
    for (int a = 0; a < nf; a++) h = mix(h ^ fs[a]) + 0x9e3779b97f4a7c15ULL;
    if (h == 0) h = 1;
    if (!hs_add(&FOUND, h)) { dups++; return; }
    found[nP]++;
}

/* candidates */
typedef struct { int kind; long p; int c; long t; int lay; } Cand;

static void rec(void);
static void try_cands(Cand *cs, int nc) {
    /* branch i: add cs[i]; then exclude it */
    for (int i = 0; i < nc; i++) {
        Cand *x = &cs[i];
        if (x->kind == 0) { add_plaq(x->p); rec(); del_plaq(x->p); exP[x->p] = 1; }
        else { add_front(x->c, x->t, x->lay, 1); rec(); del_front(1); hs_add(&EXF, fkey(x->c, x->t, x->lay)); }
    }
    for (int i = nc - 1; i >= 0; i--) {
        Cand *x = &cs[i];
        if (x->kind == 0) exP[x->p] = 0; else hs_del(&EXF, fkey(x->c, x->t, x->lay));
    }
}

static void rec(void) {
    nodes++;
    long lmin; int nd = deficient(&lmin);
    int cap = capacity();
    if (PRUNE && nd > cap) return;
    int ncap = 1 << 12; Cand *cs = malloc(sizeof(Cand) * ncap); int nc = 0;
#define PUSH(...) do { if (nc == ncap) { ncap *= 2; cs = realloc(cs, sizeof(Cand) * ncap); } cs[nc++] = (Cand){__VA_ARGS__}; } while (0)
    if (nd == 0) {
        if (final_check()) record();
        /* elements sharing a vertex with X */
        long Vx[1024]; int nvx = 0; long V[4];
        for (int a = 0; a < nP; a++) { pverts(Pl[a], V); for (int b = 0; b < 4; b++) if (!vmark[V[b]]) { vmark[V[b]] = 1; Vx[nvx++] = V[b]; } }
        for (int a = 0; a < nF; a++) { FClass *C = &cls[Fs[a].c]; for (int b = 0; b < C->nv; b++) { long v = Fs[a].t + C->voff[b]; if (!vmark[v]) { vmark[v] = 1; Vx[nvx++] = v; } } }
        for (int a = 0; a < nvx; a++) vmark[Vx[a]] = 0;
        if (nP < NMAX) {
            long Q[64];
            for (int a = 0; a < nvx; a++) {
                int m = plaqs_at(Vx[a], Q);
                for (int b = 0; b < m; b++) {
                    long p = Q[b];
                    if (inP[p] || exP[p] || candmark[p]) continue;
                    long PV[4]; pverts(p, PV); int inb = 1; for (int z = 0; z < 4; z++) if (!inbox(PV[z], 2)) inb = 0;
                    if (!inb) { boundary_hits++; continue; }
                    if (PRUNE) { nP++; int c2 = capacity(); nP--; if (def_after_plaq(nd, p) > c2 - 0) { candmark[p] = 1; PUSH(0, p, 0, 0, 0); cs[nc - 1].kind = -1; continue; } }
                    candmark[p] = 1; PUSH(0, p, 0, 0, 0);
                }
            }
        }
        for (int k = 0; k <= maxsize; k++) {
            if (rem[k] <= 0) continue;
            for (int a = 0; a < nvx; a++) for (int b = 0; b < nthruV[k]; b++) {
                int c = thruV[k][b].c; long t = Vx[a] - thruV[k][b].off;
                for (int lay = 0; lay < 2; lay++) {
                    u64 key = fkey(c, t, lay);
                    if (hs_has(&CAND, key)) continue;
                    hs_add(&CAND, key);
                    if (!front_fits(c, t)) { boundary_hits++; PUSH(-2, 0, c, t, lay); continue; }
                    if (hs_has(&EXF, key) || !front_ok_layer(c, t, lay)) { PUSH(-2, 0, c, t, lay); continue; }
                    if (PRUNE && def_after_front(nd, c, t) > cap - k) { PUSH(-2, 0, c, t, lay); continue; }
                    PUSH(1, 0, c, t, lay);
                }
            }
        }
    } else {
        if (nP < NMAX) {
            long Q[8]; int m = plaqs_through(lmin, Q);
            for (int b = 0; b < m; b++) {
                long p = Q[b];
                if (inP[p] || exP[p]) continue;
                long PV[4]; pverts(p, PV); int inb = 1; for (int z = 0; z < 4; z++) if (!inbox(PV[z], 2)) inb = 0;
                if (!inb) { boundary_hits++; continue; }
                if (PRUNE) { nP++; int c2 = capacity(); nP--; if (def_after_plaq(nd, p) > c2) continue; }
                PUSH(0, p, 0, 0, 0);
            }
        }
        long x = lmin / D; int dir = (int)(lmin % D);
        for (int k = 0; k <= maxsize; k++) {
            if (rem[k] <= 0) continue;
            for (int b = 0; b < nthruL[k][dir]; b++) {
                int c = thruL[k][dir][b].c; long t = x - thruL[k][dir][b].off;
                if (!front_fits(c, t)) { boundary_hits++; continue; }
                for (int lay = 0; lay < 2; lay++) {
                    if (!front_ok_layer(c, t, lay)) continue;
                    if (hs_has(&EXF, fkey(c, t, lay))) continue;
                    if (PRUNE && def_after_front(nd, c, t) > cap - k) continue;
                    PUSH(1, 0, c, t, lay);
                }
            }
        }
    }
    /* clean the candidate marks of the closed-state enumeration, keep only real candidates */
    if (nd == 0) {
        for (int i = nc - 1; i >= 0; i--) {
            if (cs[i].kind == 0 || cs[i].kind == -1) candmark[cs[i].p] = 0;
            else hs_del(&CAND, fkey(cs[i].c, cs[i].t, cs[i].lay));
        }
        int m = 0;
        for (int i = 0; i < nc; i++) if (cs[i].kind == 0 || cs[i].kind == 1) cs[m++] = cs[i];
        nc = m;
    }
    try_cands(cs, nc);
    free(cs);
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: junc_partition DUMP D NMAX SIZES [noprune]\n"); return 1; }
    D = atoi(argv[2]); NMAX = atoi(argv[3]);
    int sizes[16], ns = 0; maxsize = 4;
    char sizestr[64]; strncpy(sizestr, argv[4], 63); sizestr[63] = 0;
    for (char *s = strtok(argv[4], ","); s; s = strtok(NULL, ",")) { sizes[ns] = atoi(s); rem[sizes[ns]]++; if (sizes[ns] > maxsize) maxsize = sizes[ns]; ns++; }
    if (argc > 5 && !strcmp(argv[5], "noprune")) PRUNE = 0;
    int sumk = 0; for (int a = 0; a < ns; a++) sumk += sizes[a];
    R = NMAX + sumk + 4; W = 2 * R + 1;
    stride[0] = 1; for (int i = 1; i <= D; i++) stride[i] = stride[i - 1] * W;
    NVERT = stride[D]; NLINK = NVERT * D; NPL = D * (D - 1) / 2; NPLAQ = NVERT * NPL;
    { int t = 0; for (int i = 0; i < D; i++) for (int j = i + 1; j < D; j++) { pi_[t] = i; pj_[t] = j; t++; } }
    ORIGIN = 0; for (int i = 0; i < D; i++) ORIGIN += R * stride[i];
    cov = calloc(NLINK, 1); occ[0] = calloc(NVERT, 1); occ[1] = calloc(NVERT, 1); inP = calloc(NPLAQ, 1); exP = calloc(NPLAQ, 1);
    vmark = calloc(NVERT, 1); candmark = calloc(NPLAQ, 1);
    hs_init(&EXF, 23); hs_init(&CAND, 23); hs_init(&FOUND, 22);
    /* load the classes */
    FILE *fh = fopen(argv[1], "r"); if (!fh) { perror("dump"); return 1; }
    int cap = 1024; cls = malloc(sizeof(FClass) * cap); ncls = 0;
    int k;
    while (fscanf(fh, "%d", &k) == 1) {
        FClass C; C.k = k; C.nv = 0;
        int use = (k <= maxsize);
        for (int a = 0; a < k; a++) {
            int xs[3] = {0, 0, 0}, d;
            for (int i = 0; i < D; i++) if (fscanf(fh, "%d", &xs[i]) != 1) return 3;
            if (fscanf(fh, "%d", &d) != 1) return 3;
            long off = 0; for (int i = 0; i < D; i++) off += xs[i] * stride[i];
            C.loff[a] = off; C.ldir[a] = d;
            long ends[2] = {off, off + stride[d]};
            for (int e = 0; e < 2; e++) { int seenv = 0; for (int b = 0; b < C.nv; b++) if (C.voff[b] == ends[e]) seenv = 1; if (!seenv) C.voff[C.nv++] = ends[e]; }
        }
        if (!use) continue;
        if (ncls == cap) { cap *= 2; cls = realloc(cls, sizeof(FClass) * cap); }
        cls[ncls++] = C;
    }
    fclose(fh);
    for (int kk = 0; kk <= maxsize; kk++) { for (int d = 0; d < D; d++) { thruL[kk][d] = malloc(sizeof(TL) * 1); nthruL[kk][d] = 0; } thruV[kk] = malloc(sizeof(TL)); nthruV[kk] = 0; }
    for (int c = 0; c < ncls; c++) {
        int kk = cls[c].k;
        for (int a = 0; a < kk; a++) { int d = cls[c].ldir[a]; thruL[kk][d] = realloc(thruL[kk][d], sizeof(TL) * (nthruL[kk][d] + 1)); thruL[kk][d][nthruL[kk][d]++] = (TL){c, cls[c].loff[a]}; }
        for (int a = 0; a < cls[c].nv; a++) { thruV[kk] = realloc(thruV[kk], sizeof(TL) * (nthruV[kk] + 1)); thruV[kk][nthruV[kk]++] = (TL){c, cls[c].voff[a]}; }
    }
    /* p0 = (bd q, 0): the 4-cycle in the plane (0,1) with its smallest vertex at the origin */
    int c0 = -1;
    for (int c = 0; c < ncls && c0 < 0; c++) {
        if (cls[c].k != 4) continue;
        long want[4] = {0 * D + 0, (stride[0]) * D + 1, (stride[1]) * D + 0, 0 * D + 1};
        int okc = 1;
        for (int a = 0; a < 4; a++) { int f = 0; for (int b = 0; b < 4; b++) if (cls[c].loff[b] * D + cls[c].ldir[b] == want[a]) f = 1; if (!f) okc = 0; }
        if (okc) c0 = c;
    }
    if (c0 < 0) { fprintf(stderr, "no plaquette class\n"); return 4; }
    printf("D=%d NMAX=%d sizes=%s prune=%d: %d classes loaded (k<=%d); N^v by size:", D, NMAX, sizestr, PRUNE, ncls, maxsize);
    for (int kk = 4; kk <= maxsize; kk++) if (nthruV[kk]) printf(" %d:%d", kk, nthruV[kk]);
    printf("\n"); fflush(stdout);
    add_front(c0, ORIGIN, 0, 0);
    rec();
    printf("  junctions by n:");
    for (int n = 1; n <= NMAX; n++) printf(" %d:%lld", n, found[n]);
    printf("\n  nodes %lld, duplicates %lld, candidates cut by the box %lld\n", nodes, dups, boundary_hits);
    return 0;
}
