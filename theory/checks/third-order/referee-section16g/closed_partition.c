/* Referee (Section 16, seventh pass): closed polymers through a plaquette of Z^d by a partition search.
 * Usage: closed_partition d NMAX            (build: cc -O2 -o closed_partition closed_partition.c)
 *
 * Counts the sets Y of plaquettes of Z^d with p0 in Y, |Y| <= NMAX, link-connected (two plaquettes adjacent when they share a link: the polymers of
 * Lemma mg:lem:theta7), in which every link of E(Y) lies in at least two plaquettes of Y.
 * Search: a state is (X, Excl) with p0 in X; the targets below it are the closed link-connected Y with X in Y and Y disjoint from Excl.
 *  - If X has links of coverage one, take the one with the fewest candidates (plaquettes through it, not in X, not excluded); Y contains a
 *    candidate; branch i = "c_i is the first candidate (in index order) in Y": add c_i, exclude c_1..c_{i-1}.
 *  - If X is closed, record it (Y = X); otherwise Y contains a plaquette sharing a link with X (link-connected), not excluded; branch i = "c_i is
 *    the first such plaquette in Y".
 * Exclusions are never cleared, so every Y is reached exactly once (a hash set of the sorted plaquettes counts duplicates, which must be 0).
 * Prune: every link of coverage one needs a further plaquette through it, and a plaquette has four links.
 * This differs from the authors' closed_sets.c (smallest link, all subsets at once, exclusions cleared at closed sets, memo of closed sets). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int d, NMAX, NO;
static long R, W, stride[8], ORIG;
static size_t NVERT, NLINK, NPLAQ;
static int oi[32], oj[32], orient[8][8];
static unsigned char *cov; static char *inY, *exY, *mark;
static long Yl[64]; static int nY;
static long *exst; static long nex, capex;
static long long cnt_n[64], dups, nodes;

typedef struct { uint64_t *t; long mask, n; } HS;
static HS H;
static uint64_t mix(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
static int hs_add(uint64_t k) {
    if (2 * (H.n + 1) > H.mask) {
        long nm = 2 * (H.mask + 1) - 1; uint64_t *nt = calloc(nm + 1, 8);
        for (long i = 0; i <= H.mask; i++) if (H.t[i]) { long j = mix(H.t[i]) & nm; while (nt[j]) j = (j + 1) & nm; nt[j] = H.t[i]; }
        free(H.t); H.t = nt; H.mask = nm;
    }
    long i = mix(k) & H.mask;
    while (H.t[i]) { if (H.t[i] == k) return 0; i = (i + 1) & H.mask; }
    H.t[i] = k; H.n++; return 1;
}
static inline void plinks(long P, long *L) {
    long v = P / NO; int o = P % NO, i = oi[o], j = oj[o];
    L[0] = v * d + i; L[1] = (v + stride[i]) * d + j; L[2] = (v + stride[j]) * d + i; L[3] = v * d + j;
}
static int through(long L, long *out) {
    long v = L / d; int k = L % d, c = 0;
    for (int m = 0; m < d; m++) {
        if (m == k) continue;
        int i = k < m ? k : m, j = k < m ? m : k, o = orient[i][j];
        out[c++] = v * NO + o; out[c++] = (v - stride[m]) * NO + o;
    }
    return c;
}
static void add(long P) { long L[4]; plinks(P, L); for (int a = 0; a < 4; a++) cov[L[a]]++; inY[P] = 1; Yl[nY++] = P; }
static void del(long P) { long L[4]; plinks(P, L); for (int a = 0; a < 4; a++) cov[L[a]]--; inY[P] = 0; nY--; }
static void excl(long P) { if (nex == capex) { capex *= 2; exst = realloc(exst, capex * sizeof(long)); } exY[P] = 1; exst[nex++] = P; }
static int cmpl(const void *a, const void *b) { long x = *(const long *)a, y = *(const long *)b; return (x > y) - (x < y); }

static void record(void) {
    long k[64]; memcpy(k, Yl, nY * sizeof(long)); qsort(k, nY, sizeof(long), cmpl);
    uint64_t h = 0x12345678abcdefULL; for (int a = 0; a < nY; a++) h = mix(h ^ (uint64_t)k[a]) + 0x9e3779b97f4a7c15ULL;
    if (!h) h = 1;
    if (!hs_add(h)) dups++;
    cnt_n[nY]++;
}

static void rec(void) {
    nodes++;
    /* links of coverage one */
    long def[256]; int nd = 0;
    for (int a = 0; a < nY; a++) { long L[4]; plinks(Yl[a], L); for (int b = 0; b < 4; b++) if (cov[L[b]] == 1) def[nd++] = L[b]; }
    long cand[2048]; int nc = 0;
    if (nd == 0) {
        record();
        if (nY >= NMAX) return;
        for (int a = 0; a < nY; a++) {
            long L[4]; plinks(Yl[a], L);
            for (int b = 0; b < 4; b++) { long T[16]; int m = through(L[b], T); for (int x = 0; x < m; x++) if (!inY[T[x]] && !exY[T[x]] && !mark[T[x]]) { mark[T[x]] = 1; cand[nc++] = T[x]; } }
        }
        for (int a = 0; a < nc; a++) mark[cand[a]] = 0;
    } else {
        if (nd > 4 * (NMAX - nY)) return;
        int best = -1, bestc = 99; long bc[16]; int nbc = 0;
        for (int a = 0; a < nd; a++) {
            long T[16]; int m = through(def[a], T), c = 0; long cc[16];
            for (int x = 0; x < m; x++) if (!inY[T[x]] && !exY[T[x]]) cc[c++] = T[x];
            if (c < bestc) { bestc = c; best = a; nbc = c; memcpy(bc, cc, c * sizeof(long)); if (c == 0) return; }
        }
        (void)best;
        for (int a = 0; a < nbc; a++) cand[nc++] = bc[a];
    }
    qsort(cand, nc, sizeof(long), cmpl);
    long nex0 = nex;
    for (int i = 0; i < nc; i++) {
        add(cand[i]); rec(); del(cand[i]);
        excl(cand[i]);
    }
    while (nex > nex0) { nex--; exY[exst[nex]] = 0; }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s d NMAX\n", argv[0]); return 1; }
    d = atoi(argv[1]); NMAX = atoi(argv[2]);
    NO = 0; for (int i = 0; i < d; i++) for (int j = i + 1; j < d; j++) { oi[NO] = i; oj[NO] = j; orient[i][j] = NO; NO++; }
    /* a link-connected set of n plaquettes through p0 spans at most n + 1 values of each coordinate, so with p0 at the centre of a box of side
       2 NMAX + 5 every vertex of a searched set, and every plaquette through one of its links, lies inside */
    R = NMAX + 2; W = 2 * R + 1;
    NVERT = 1; for (int i = 0; i < d; i++) { stride[i] = NVERT; NVERT *= W; }
    NLINK = NVERT * d; NPLAQ = NVERT * NO;
    ORIG = 0; for (int i = 0; i < d; i++) ORIG += R * stride[i];
    cov = calloc(NLINK, 1); inY = calloc(NPLAQ, 1); exY = calloc(NPLAQ, 1); mark = calloc(NPLAQ, 1);
    capex = 1 << 16; exst = malloc(capex * sizeof(long));
    H.mask = (1 << 16) - 1; H.t = calloc(H.mask + 1, 8);
    if (!cov || !inY || !exY || !mark) { fprintf(stderr, "out of memory\n"); return 1; }
    long p0 = ORIG * NO + orient[0][1];
    add(p0); rec();
    printf("d=%d NMAX=%d (partition search, box side %ld): closed link-connected sets through a plaquette by size:", d, NMAX, W);
    for (int n = 1; n <= NMAX; n++) if (cnt_n[n]) printf(" %d:%lld", n, cnt_n[n]);
    printf("\n  duplicates %lld, nodes %lld\n", dups, nodes);
    return 0;
}
