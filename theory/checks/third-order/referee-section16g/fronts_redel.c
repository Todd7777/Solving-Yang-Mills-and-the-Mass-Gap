/* Referee (Section 16, seventh pass): N^v_k, the fronts of Z^D with k edges through a vertex, by Redelmeier's enumeration of edge sets.
 * Usage: fronts_redel D K [dumpfile]  (build: cc -O2 -o fronts_redel fronts_redel.c)
 * With a dump file, every class of fronts is written as one line "k  x.. dir  x.. dir ..." (coordinates of the tail of each edge, then its
 * direction), the format read by the sixth referee's junc_partition.c.
 *
 * A front is a finite non-empty connected set of edges (two edges adjacent when they share an endpoint) in which every vertex has degree >= 2.
 * Every translation class of connected edge sets is enumerated once, with its smallest edge (in the order index = vertex index * D + direction,
 * the vertex index being lexicographic, hence translation invariant) at the origin: Redelmeier's algorithm on the line graph, the allowed edges
 * being those greater than the root.  N^v_k is the sum of |V(phi)| over the classes that are fronts (each class has |V| translates through a
 * given vertex); N_l(k), the fronts through a given edge of direction 0, is the sum of the numbers of edges of direction 0.
 * Prunes (necessary conditions, so no front is lost): (a) a vertex of degree one needs a further edge, one edge serves at most two such vertices, so
 * #deg-1 vertices <= 2 (K - size); (b) in Redelmeier's scheme every edge adjacent to the current set has been reached, so a further edge at a vertex
 * of the set can only come from the current untried list: a vertex of degree one with no incident edge in that list is dead.
 * This is neither the authors' decision search (fronts_decision.c) nor the sixth referee's fronts_enum.c (rooted at the smallest vertex). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int D, K, W;
static size_t NV, NE; static long stride[8], ORIG;
static char *reached, *inS, *inU;
static unsigned char *deg;
static long rootE;
static long long cls[64], Nv[64], Nl[64], nodes;
static FILE *dump = NULL;
static long Sl[64]; static int nS;

static inline long eend(long e, int s) { long v = e / D; return s ? v + stride[e % D] : v; }
static int inc(long v, long *E) { int c = 0; for (int i = 0; i < D; i++) { E[c++] = v * D + i; E[c++] = (v - stride[i]) * D + i; } return c; }
static int adj(long e, long *out) {   /* edges sharing an endpoint with e */
    int c = 0; long E[16];
    for (int s = 0; s < 2; s++) { int m = inc(eend(e, s), E); for (int a = 0; a < m; a++) if (E[a] != e) out[c++] = E[a]; }
    return c;
}

static void rec(long *U, int nu) {
    for (int a = 0; a < nu; a++) {
        long e = U[a];
        /* add e */
        Sl[nS++] = e; inS[e] = 1;
        long v0 = eend(e, 0), v1 = eend(e, 1);
        deg[v0]++; deg[v1]++;
        nodes++;
        /* statistics of the vertices of S */
        int ones = 0, nv = 0; long V[130];
        for (int b = 0; b < nS; b++) for (int s = 0; s < 2; s++) {
            long w = eend(Sl[b], s); int seen = 0;
            for (int c = 0; c < nv; c++) if (V[c] == w) { seen = 1; break; }
            if (!seen) { V[nv++] = w; if (deg[w] == 1) ones++; }
        }
        if (ones == 0) {
            cls[nS]++; Nv[nS] += nv;
            for (int b = 0; b < nS; b++) if (Sl[b] % D == 0) Nl[nS]++;
            if (dump) {
                fprintf(dump, "%d", nS);
                for (int b = 0; b < nS; b++) {
                    long v = Sl[b] / D; fprintf(dump, " ");
                    for (int i = 0; i < D; i++) fprintf(dump, " %ld", (long)((v / stride[i]) % W) - (K + 2));
                    fprintf(dump, " %ld", Sl[b] % D);
                }
                fprintf(dump, "\n");
            }
        }
        if (nS < K && ones <= 2 * (K - nS)) {
            /* the child's untried list: U[a+1..] and the new neighbours of e (not reached, greater than the root) */
            long *nu2 = malloc(sizeof(long) * (nu + 4 * D + 4)); int m = 0;
            for (int b = a + 1; b < nu; b++) nu2[m++] = U[b];
            long N[32]; int nn = adj(e, N); int m0 = m;
            for (int b = 0; b < nn; b++) if (N[b] > rootE && !reached[N[b]]) { reached[N[b]] = 1; nu2[m++] = N[b]; }
            /* dead-vertex prune */
            for (int b = 0; b < m; b++) inU[nu2[b]] = 1;
            int dead = 0;
            for (int c = 0; c < nv && !dead; c++) {
                if (deg[V[c]] != 1) continue;
                long E[16]; int k = inc(V[c], E), ok = 0;
                for (int x = 0; x < k; x++) if (inU[E[x]]) { ok = 1; break; }
                if (!ok) dead = 1;
            }
            for (int b = 0; b < m; b++) inU[nu2[b]] = 0;
            if (!dead) rec(nu2, m);
            for (int b = m0; b < m; b++) reached[nu2[b]] = 0;
            free(nu2);
        }
        /* remove e */
        deg[v0]--; deg[v1]--; inS[e] = 0; nS--;
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s D K\n", argv[0]); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]);
    if (argc > 3) dump = fopen(argv[3], "w");
    W = 2 * K + 5; NV = 1; for (int i = 0; i < D; i++) { stride[i] = NV; NV *= W; }
    NE = NV * D; ORIG = 0; for (int i = 0; i < D; i++) ORIG += (long)(K + 2) * stride[i];
    reached = calloc(NE, 1); inS = calloc(NE, 1); inU = calloc(NE, 1); deg = calloc(NV, 1);
    for (int i = 0; i < D; i++) {
        rootE = ORIG * D + i;
        reached[rootE] = 1;
        long U[1] = {rootE};
        rec(U, 1);
        reached[rootE] = 0;
        memset(reached, 0, NE);
    }
    if (dump) fclose(dump);
    printf("D=%d K=%d (Redelmeier, smallest edge at the origin): %lld nodes\n", D, K, nodes);
    printf("classes of fronts by k:"); for (int k = 1; k <= K; k++) if (cls[k]) printf(" %d:%lld", k, cls[k]); printf("\n");
    printf("N^v_k:"); for (int k = 1; k <= K; k++) if (Nv[k]) printf(" %d:%lld", k, Nv[k]); printf("\n");
    printf("N_l(k):"); for (int k = 1; k <= K; k++) if (Nl[k]) printf(" %d:%lld", k, Nl[k]); printf("\n");
    return 0;
}
