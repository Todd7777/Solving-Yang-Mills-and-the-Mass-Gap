/* Referee (Section 16, sixth pass): fronts of Z^D up to translation, by Redelmeier's algorithm on the edges.
 * Usage: fronts_enum D K [dumpfile]
 *
 * A front is a finite connected set of links in which every vertex has degree >= 2 (Lemma mg:lem:junction / Definition mg:def:fronts, as used
 * by the junction searches).  Each translation class is generated exactly once, rooted at its lexicographically smallest vertex O: Redelmeier's
 * algorithm enumerates every connected edge set with an edge at O and all vertices >= O (lex), each exactly once (untried set, edges excluded
 * once popped).  Prunes (both necessary): a vertex of degree one whose incident links are all in the set, excluded or not allowed is dead;
 * and #(vertices of degree one) <= 2 (K - size).  Output: for each k, the number of classes and N^v_k = sum over classes of |V(front)|, the
 * number of fronts with k edges through a given vertex.  With a dump file, every class is written as "k  x y z dir  x y z dir ...".
 * This is not the authors' decision search (fronts_decision.py) nor the third referee's code (fronts_redelmeier.c, rooted at a vertex). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int D, K, W, OFF;          /* box [-OFF, OFF]^D, W = 2*OFF+1 */
static long NV;                    /* number of vertices in the box */
static unsigned char *inA, *seen, *inU, *allowed;
static int *deg;
static int asize = 0, nverts = 0, ndeg1 = 0;
static int animal[64];
static long long ncls[64], nv[64];
static FILE *dump = NULL;
static int stride[4];

static int coord(long v, int i) { return (int)((v / stride[i]) % W) - OFF; }
static long vadd(long v, int i, int s) { return v + s * stride[i]; }
static long e_of(long v, int i) { return v * D + i; }
static long e_tail(long e) { return e / D; }
static int e_dir(long e) { return (int)(e % D); }
static long e_head(long e) { return vadd(e / D, (int)(e % D), 1); }

static int lex_ge_origin(long v) {
    for (int i = 0; i < D; i++) { int c = coord(v, i); if (c > 0) return 1; if (c < 0) return 0; }
    return 1;
}
static int inbox(long v) {
    for (int i = 0; i < D; i++) { int c = coord(v, i); if (c < -OFF + 1 || c > OFF - 1) return 0; }
    return 1;
}

/* incident edges of v: up to 2D */
static int incident(long v, long *out) {
    int m = 0;
    for (int i = 0; i < D; i++) {
        out[m++] = e_of(v, i);
        out[m++] = e_of(vadd(v, i, -1), i);
    }
    return m;
}

static void add_edge(long e) {
    long a = e_tail(e), b = e_head(e);
    inA[e] = 1; animal[asize++] = (int)e;
    long ends[2] = {a, b};
    for (int t = 0; t < 2; t++) {
        long v = ends[t];
        if (deg[v] == 0) nverts++;
        if (deg[v] == 1) ndeg1--;
        deg[v]++;
        if (deg[v] == 1) ndeg1++;
    }
}
static void remove_edge(long e) {
    long a = e_tail(e), b = e_head(e);
    inA[e] = 0; asize--;
    long ends[2] = {a, b};
    for (int t = 0; t < 2; t++) {
        long v = ends[t];
        if (deg[v] == 1) ndeg1--;
        deg[v]--;
        if (deg[v] == 1) ndeg1++;
        if (deg[v] == 0) nverts--;
    }
}

static int dead(void) {
    /* a vertex of degree one with no incident edge in the current untried set */
    long inc[8];
    for (int t = 0; t < asize; t++) {
        long e = animal[t];
        long ends[2] = {e_tail(e), e_head(e)};
        for (int s = 0; s < 2; s++) {
            long v = ends[s];
            if (deg[v] != 1) continue;
            int m = incident(v, inc), ok = 0;
            for (int j = 0; j < m; j++) if (inU[inc[j]]) { ok = 1; break; }
            if (!ok) return 1;
        }
    }
    return 0;
}

static void record(void) {
    ncls[asize]++; nv[asize] += nverts;
    if (dump) {
        fprintf(dump, "%d", asize);
        for (int t = 0; t < asize; t++) {
            long e = animal[t];
            fprintf(dump, " ");
            for (int i = 0; i < D; i++) fprintf(dump, " %d", coord(e_tail(e), i));
            fprintf(dump, " %d", e_dir(e));
        }
        fprintf(dump, "\n");
    }
}

static void rec(long *U, int nU) {
    long *U2 = malloc(sizeof(long) * (nU + 4 * D * 2 + 8));
    while (nU > 0) {
        long e = U[--nU];
        inU[e] = 0;                               /* popped: excluded from now on at this level */
        add_edge(e);
        if (ndeg1 == 0) record();
        if (asize < K && ndeg1 <= 2 * (K - asize)) {
            int n2 = 0;
            for (int j = 0; j < nU; j++) U2[n2++] = U[j];
            int start = n2;
            long ends[2] = {e_tail(e), e_head(e)}, inc[8];
            for (int s = 0; s < 2; s++) {
                int m = incident(ends[s], inc);
                for (int j = 0; j < m; j++) {
                    long f = inc[j];
                    if (!allowed[f] || seen[f] || inA[f]) continue;
                    seen[f] = 1; inU[f] = 1; U2[n2++] = f;
                }
            }
            if (!dead()) rec(U2, n2);
            for (int j = 0; j < start; j++) inU[U2[j]] = 1;   /* restore the inherited untried flags (rec popped them) */
            for (int j = start; j < n2; j++) { seen[U2[j]] = 0; inU[U2[j]] = 0; }
        }
        remove_edge(e);
    }
    free(U2);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: fronts_enum D K [dumpfile]\n"); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]);
    if (argc > 3) dump = fopen(argv[3], "w");
    OFF = K + 2; W = 2 * OFF + 1;
    stride[0] = 1; for (int i = 1; i <= D; i++) stride[i] = stride[i - 1] * W;
    NV = stride[D];
    long NE = NV * D;
    inA = calloc(NE, 1); seen = calloc(NE, 1); inU = calloc(NE, 1); allowed = calloc(NE, 1); deg = calloc(NV, sizeof(int));
    for (long v = 0; v < NV; v++) {
        if (!inbox(v) || !lex_ge_origin(v)) continue;
        for (int i = 0; i < D; i++) if (inbox(vadd(v, i, 1))) allowed[e_of(v, i)] = 1;
    }
    long O = 0; for (int i = 0; i < D; i++) O += OFF * stride[i];
    long inc[8]; int m = incident(O, inc);
    long U[8]; int nU = 0;
    for (int j = 0; j < m; j++) if (allowed[inc[j]]) { U[nU++] = inc[j]; seen[inc[j]] = 1; inU[inc[j]] = 1; }
    rec(U, nU);
    printf("D=%d K=%d\n k  classes  N^v_k\n", D, K);
    for (int k = 1; k <= K; k++) if (ncls[k]) printf("%2d %8lld %9lld\n", k, ncls[k], nv[k]);
    if (dump) fclose(dump);
    return 0;
}
