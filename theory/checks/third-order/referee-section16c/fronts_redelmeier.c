/* Referee (Section 16, third pass): N^v_k, the number of fronts with k edges through a vertex of Z^D
 * (a front: a non-empty connected set of edges in which every vertex has degree >= 2; Section 7 of Part IV).
 *
 * Algorithm (independent of the authors' decision search fronts_decision.py and of ref_fronts.py): Redelmeier's
 * enumeration of connected edge sets ("bond animals") containing a root edge at the origin, each animal exactly once.
 * A front through the origin contains at least two of the 2D edges at the origin; it is counted under the smallest
 * of them (the smaller origin edges are forbidden for that root). Prune: a set S with d1 vertices of degree one has
 * no front superset with fewer than |S| + ceil(d1/2) edges (each further edge raises at most two degrees).
 *
 * Usage: ./fronts_redelmeier D K [dumpfile kmin]   (dump all fronts with kmin <= k <= K, one per line:
 *        "k  x_1 .. x_D dir  x_1 .. x_D dir ...", the edge from x to x + e_dir)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int D, K, S, NV, NE;
static int off;                 /* coordinate offset */
static int *deg; static unsigned char *inS, *seen;
static long long cnt[64];
static int d1 = 0;              /* number of vertices of degree exactly one in the current set */
static int cur[64], ncur = 0;
static FILE *dump = NULL; static int kmin_dump = 1 << 30;
static int stride[4];

static int vid(const int *x) { int v = 0; for (int i = 0; i < D; i++) v += (x[i] + off) * stride[i]; return v; }
static void vcoord(int v, int *x) { for (int i = 0; i < D; i++) { x[i] = (v / stride[i]) % S - off; } }

static int nbrs(int e, int *out) {        /* edges sharing an endpoint with e = (v, i) */
    int v = e / D, i = e % D, n = 0;
    int w = v + stride[i];
    for (int j = 0; j < D; j++) {
        int a;
        a = v * D + j;                   if (a != e) out[n++] = a;       /* (v, j) */
        a = (v - stride[j]) * D + j;     out[n++] = a;                   /* (v - e_j, j) */
        a = w * D + j;                   out[n++] = a;                   /* (w, j) */
        a = (w - stride[j]) * D + j;     if (a != e) out[n++] = a;       /* (w - e_j, j) */
    }
    return n;
}
static void add(int e) {
    int v = e / D, w = v + stride[e % D];
    int ends[2] = {v, w};
    for (int t = 0; t < 2; t++) { int u = ends[t]; deg[u]++; if (deg[u] == 1) d1++; else if (deg[u] == 2) d1--; }
    inS[e] = 1; cur[ncur++] = e;
}
static void rem(int e) {
    int v = e / D, w = v + stride[e % D];
    int ends[2] = {v, w};
    for (int t = 0; t < 2; t++) { int u = ends[t]; deg[u]--; if (deg[u] == 0) d1--; else if (deg[u] == 1) d1++; }
    inS[e] = 0; ncur--;
}
static void emit(void) {
    fprintf(dump, "%d", ncur);
    for (int t = 0; t < ncur; t++) {
        int e = cur[t], x[4]; vcoord(e / D, x);
        fprintf(dump, " ");
        for (int i = 0; i < D; i++) fprintf(dump, " %d", x[i]);
        fprintf(dump, " %d", e % D);
    }
    fprintf(dump, "\n");
}
static void rec(int *untried, int nu) {
    int nb[32];
    while (nu > 0) {
        int e = untried[--nu];
        add(e);
        if (d1 == 0) { cnt[ncur]++; if (dump && ncur >= kmin_dump) emit(); }
        if (ncur < K && ncur + (d1 + 1) / 2 <= K) {
            int *nu2 = malloc(sizeof(int) * (nu + 32));
            memcpy(nu2, untried, sizeof(int) * nu);
            int m = nu, k = nbrs(e, nb), added[32], na = 0;
            for (int t = 0; t < k; t++) if (!seen[nb[t]]) { seen[nb[t]] = 1; nu2[m++] = nb[t]; added[na++] = nb[t]; }
            rec(nu2, m);
            for (int t = 0; t < na; t++) seen[added[t]] = 0;
            free(nu2);
        }
        rem(e);
    }
}
int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s D K [dumpfile kmin]\n", argv[0]); return 1; }
    D = atoi(argv[1]); K = atoi(argv[2]);
    if (argc >= 5) { dump = fopen(argv[3], "w"); kmin_dump = atoi(argv[4]); }
    off = K + 2; S = 2 * off + 1;
    stride[0] = 1; for (int i = 1; i < D; i++) stride[i] = stride[i - 1] * S;
    NV = stride[D - 1] * S; NE = NV * D;
    deg = calloc(NV, sizeof(int)); inS = calloc(NE, 1); seen = calloc(NE, 1);
    int o[4] = {0, 0, 0, 0}, ov = vid(o);
    int orig[8], no = 0;                       /* the 2D edges at the origin, in a fixed order */
    for (int j = 0; j < D; j++) { orig[no++] = ov * D + j; orig[no++] = (ov - stride[j]) * D + j; }
    for (int r = 0; r < no; r++) {
        memset(seen, 0, NE);
        for (int t = 0; t < r; t++) seen[orig[t]] = 1;    /* forbidden: smaller origin edges */
        seen[orig[r]] = 1;
        int u[1] = {orig[r]};
        rec(u, 1);
    }
    printf("D=%d K=%d: N^v_k =", D, K);
    for (int k = 1; k <= K; k++) if (cnt[k]) printf(" %d:%lld", k, cnt[k]);
    printf("\n");
    if (dump) fclose(dump);
    return 0;
}
