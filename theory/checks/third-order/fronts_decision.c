/* Section 16 of Part IV: exact N^v_k (fronts through a vertex of Z^D), the decision search of fronts_decision.py in C.
 * Usage: fronts_decision D K            (build: cc -O2 -o fronts_decision fronts_decision.c)
 *
 * A front is a connected set of edges of Z^D in which every vertex has degree >= 2. The search produces every front through the origin with at most K
 * edges exactly once: it decides all edges at one vertex at a time (each undecided incident edge goes into the front or is excluded for good), first at
 * the origin (at least two edges), then always at an undecided vertex of the current set, chosen by a rule that depends only on the current state
 * (degree one first, then fewest undecided incident edges, then the smallest index). For a given front the sequence of decisions is forced, so it is
 * reached by exactly one path; the set stays connected, and when no undecided vertex of positive degree remains every vertex has all its edges decided
 * and degree >= 2. Prune: each vertex of degree one needs another edge, and one edge serves at most two of them.
 * It prints N^v_k (fronts of k edges through the origin) and N_l(k) (fronts of k edges containing the edge from the origin in direction 0).
 * Python version: fronts_decision.py (same rule up to the tie-break; the counts agree for D = 2, K <= 12 and D = 3, K <= 12). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int D, K, W;
static long NV, stride[8];
static char *inS, *inX, *dec;
static unsigned char *deg;
static long *active; static int nact;      /* vertices of positive degree, in order of first appearance */
static char *isact;
static int nS;
static long long Nv[64], Nl[64];
static long origin, link0;

static inline long lnk(long v, int i) { return v * D + i; }
/* the 2D edges at v: (v, i) to v + e_i and (v - e_i, i) to v */
static inline int incident(long v, long *L, long *other) {
    int c = 0;
    for (int i = 0; i < D; i++) {
        L[c] = lnk(v, i); other[c] = v + stride[i]; c++;
        L[c] = lnk(v - stride[i], i); other[c] = v - stride[i]; c++;
    }
    return c;
}

static void rec(long v) {
    long L[16], O[16], U[16], UO[16]; int nu = 0;
    int m = incident(v, L, O);
    for (int j = 0; j < m; j++) if (!inS[L[j]] && !inX[L[j]]) { U[nu] = L[j]; UO[nu] = O[j]; nu++; }
    for (int mask = 0; mask < (1 << nu); mask++) {
        int r = __builtin_popcount(mask);
        if (nS + r > K) continue;
        if (deg[v] + r < 2) continue;
        /* apply */
        int nact0 = nact;
        for (int j = 0; j < nu; j++) {
            if (mask >> j & 1) {
                inS[U[j]] = 1;
                long w = UO[j];
                deg[v]++; deg[w]++;
                if (!isact[w]) { isact[w] = 1; active[nact++] = w; }
            } else inX[U[j]] = 1;
        }
        if (!isact[v]) { isact[v] = 1; active[nact++] = v; }
        dec[v] = 1; nS += r;
        /* next vertex */
        long best = -1; int bk1 = 0, bk2 = 0, ones = 0;
        for (int a = 0; a < nact; a++) {
            long w = active[a];
            if (dec[w]) continue;
            if (deg[w] == 1) ones++;
            int k1 = deg[w] == 1 ? 0 : 1, k2 = 0;
            long L2[16], O2[16]; int m2 = incident(w, L2, O2);
            for (int j = 0; j < m2; j++) if (!inS[L2[j]] && !inX[L2[j]]) k2++;
            if (best < 0 || k1 < bk1 || (k1 == bk1 && (k2 < bk2 || (k2 == bk2 && w < best)))) { best = w; bk1 = k1; bk2 = k2; }
        }
        if (best < 0) {
            Nv[nS]++;
            if (inS[link0]) Nl[nS]++;
        } else if (nS + (ones + 1) / 2 <= K) rec(best);
        /* undo */
        nS -= r; dec[v] = 0;
        for (int j = 0; j < nu; j++) {
            if (mask >> j & 1) { inS[U[j]] = 0; deg[v]--; deg[UO[j]]--; }
            else inX[U[j]] = 0;
        }
        while (nact > nact0) { nact--; isact[active[nact]] = 0; }
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s D K\n", argv[0]); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]);
    if (D < 2 || D > 6 || K < 4 || K > 60) { fprintf(stderr, "bad arguments\n"); return 1; }
    W = 2 * K + 3;                                   /* coordinates -K-1 .. K+1: every vertex of a front through the origin is within distance K */
    NV = 1; for (int i = 0; i < D; i++) { stride[i] = NV; NV *= W; }
    inS = calloc(NV * D, 1); inX = calloc(NV * D, 1); dec = calloc(NV, 1); deg = calloc(NV, 1); isact = calloc(NV, 1);
    active = malloc(sizeof(long) * (K + 4) * 4);
    if (!inS || !inX || !dec || !deg || !isact || !active) { fprintf(stderr, "out of memory\n"); return 1; }
    origin = 0; for (int i = 0; i < D; i++) origin += (long)(K + 1) * stride[i];
    link0 = lnk(origin, 0);
    nact = 0; nS = 0;
    rec(origin);
    printf("D=%d, K=%d: N^v_k:", D, K);
    for (int k = 1; k <= K; k++) if (Nv[k]) printf(" %d:%lld", k, Nv[k]);
    printf("\nN_l(k) (fronts containing a given edge):");
    for (int k = 1; k <= K; k++) if (Nl[k]) printf(" %d:%lld", k, Nl[k]);
    printf("\n");
    return 0;
}
