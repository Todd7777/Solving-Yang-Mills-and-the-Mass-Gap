/* S_n(root): the number of sets P of n plaquettes of Z^D such that every plaquette of P is joined to V(root) by a chain of plaquettes of P sharing
 * vertices (the sets of the rooted counting in the proofs of Lemmas mg:lem:twofrontcount and mg:lem:jsums6), for a root front of shape
 *   q     the boundary of a plaquette (4 vertices)
 *   flat  a flat domino boundary (6 vertices)
 *   bent  a bent domino boundary (6 vertices; D >= 3)
 *   skew  the skew hexagon 0, e1, e1+e2, e1+e2+e3, e2+e3, e3 (6 vertices; D >= 3)
 * by Redelmeier's method (each set once). Part IV, Lemma mg:lem:smallcounts; the second count is rooted_rs.c (sharp_counts.py runs both).
 * Usage: rooted_sets D NMAX SHAPE   (build: cc -O2 -o rooted_sets rooted_sets.c) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static int D, NMAX, R, W, NO;
static long stride[8], NVERT, NPLAQ, ORIG;
static int oi[28], oj[28];
static long long SN[32];
static char *reached;
static int nP;
static inline void pverts(long P, long *V) {
    long v = P / NO; int o = P % NO, i = oi[o], j = oj[o];
    V[0] = v; V[1] = v + stride[i]; V[2] = v + stride[j]; V[3] = v + stride[i] + stride[j];
}
static int plaqs_at(long v, long *out) {
    int c = 0;
    for (int o = 0; o < NO; o++) for (int si = 0; si < 2; si++) for (int sj = 0; sj < 2; sj++)
        out[c++] = (v - si * stride[oi[o]] - sj * stride[oj[o]]) * NO + o;
    return c;
}
static void redel(long *untried, int nunt) {
    long *mine = malloc(nunt * sizeof(long));
    memcpy(mine, untried, nunt * sizeof(long));
    long *nu = malloc((nunt + 64) * sizeof(long));
    for (int a = 0; a < nunt; a++) {
        long c = mine[a];
        nP++; SN[nP]++;
        if (nP < NMAX) {
            long nb[64]; int nn = 0; long V[4]; pverts(c, V);
            for (int x = 0; x < 4; x++) { long ps[64]; int k = plaqs_at(V[x], ps); for (int y = 0; y < k; y++) if (!reached[ps[y]]) { reached[ps[y]] = 1; nb[nn++] = ps[y]; } }
            int nnu = 0;
            for (int b = a + 1; b < nunt; b++) nu[nnu++] = mine[b];
            for (int b = 0; b < nn; b++) nu[nnu++] = nb[b];
            redel(nu, nnu);
            for (int b = 0; b < nn; b++) reached[nb[b]] = 0;
        }
        nP--;
    }
    free(mine); free(nu);
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s D NMAX q|flat|bent|skew\n", argv[0]); return 1; }
    D = atoi(argv[1]); NMAX = atoi(argv[2]); const char *shape = argv[3];
    NO = 0; for (int i = 0; i < D; i++) for (int j = i + 1; j < D; j++) { oi[NO] = i; oj[NO] = j; NO++; }
    R = 2 * NMAX + 6; W = 2 * R + 1;
    NVERT = 1; for (int i = 0; i < D; i++) { stride[i] = NVERT; NVERT *= W; }
    NPLAQ = NVERT * NO;
    ORIG = 0; for (int i = 0; i < D; i++) ORIG += (long)R * stride[i];
    reached = calloc(NPLAQ, 1);
    long RV[8]; int nr = 0; long o = ORIG, e1 = stride[0], e2 = stride[1], e3 = D >= 3 ? stride[2] : 0;
    if (!strcmp(shape, "q"))         { long v[] = {o, o + e1, o + e2, o + e1 + e2}; nr = 4; memcpy(RV, v, sizeof v); }
    else if (!strcmp(shape, "flat")) { long v[] = {o, o + e1, o + 2 * e1, o + e2, o + e1 + e2, o + 2 * e1 + e2}; nr = 6; memcpy(RV, v, sizeof v); }
    else if (!strcmp(shape, "bent") && D >= 3) { long v[] = {o, o + e1, o + e2, o + e1 + e2, o + e3, o + e1 + e3}; nr = 6; memcpy(RV, v, sizeof v); }
    else if (!strcmp(shape, "skew") && D >= 3) { long v[] = {o, o + e1, o + e1 + e2, o + e1 + e2 + e3, o + e2 + e3, o + e3}; nr = 6; memcpy(RV, v, sizeof v); }
    else { fprintf(stderr, "unknown shape\n"); return 1; }
    long start[512]; int ns = 0;
    for (int x = 0; x < nr; x++) { long ps[64]; int k = plaqs_at(RV[x], ps); for (int y = 0; y < k; y++) if (!reached[ps[y]]) { reached[ps[y]] = 1; start[ns++] = ps[y]; } }
    printf("D=%d, root %s (%d vertices, %d plaquettes with a vertex in it):", D, shape, nr, ns);
    nP = 0; redel(start, ns);
    for (int n = 1; n <= NMAX; n++) printf(" %lld", SN[n]);
    printf("\n");
    return 0;
}
