/* A second count of the rooted sets of rooted_sets.c, by a different method (reverse search; Avis-Fukuda): a valid set P (every plaquette joined to the
 * root's vertex set V0 by a chain of plaquettes of P sharing vertices) has the parent P \ {x}, x the largest plaquette of P (by its code) whose removal
 * leaves a valid set; the children of a valid set Q are the sets Q u {y}, y a plaquette with a vertex in V0 or in a plaquette of Q, such that no
 * plaquette larger than y is removable from Q u {y}. Every valid set is generated exactly once, from the empty set.
 * Plaquettes are coded as (x_1..x_D, i, j) with coordinates offset by R; validity is checked by a breadth-first search.
 * Usage: rooted_rs D NMAX q|flat|bent|skew   (build: cc -O2 -o rooted_rs rooted_rs.c) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct { int x[4]; int i, j; } Pl;
static int D, NMAX, OFF;
static int V0[8][4], nV0;
static long long CNT[32];
static long code(const Pl *p) { long c = 0; for (int a = 0; a < D; a++) c = c * 64 + (p->x[a] + OFF); return (c * 8 + p->i) * 8 + p->j; }
static void verts(const Pl *p, int V[4][4]) {
    for (int s = 0; s < 4; s++) { for (int a = 0; a < D; a++) V[s][a] = p->x[a]; if (s & 1) V[s][p->i]++; if (s & 2) V[s][p->j]++; }
}
static int veq(const int *a, const int *b) { for (int c = 0; c < D; c++) if (a[c] != b[c]) return 0; return 1; }
static int share(const Pl *p, const Pl *q) {
    int A[4][4], B[4][4]; verts(p, A); verts(q, B);
    for (int s = 0; s < 4; s++) for (int t = 0; t < 4; t++) if (veq(A[s], B[t])) return 1;
    return 0;
}
static int touches_root(const Pl *p) {
    int A[4][4]; verts(p, A);
    for (int s = 0; s < 4; s++) for (int t = 0; t < nV0; t++) if (veq(A[s], V0[t])) return 1;
    return 0;
}
/* valid: every element of S (size m, skipping index skip) joined to the root */
static int valid(Pl *S, int m, int skip) {
    int seen[32] = {0}, queue[32], h = 0, t = 0;
    for (int a = 0; a < m; a++) if (a != skip && touches_root(&S[a])) { seen[a] = 1; queue[t++] = a; }
    while (h < t) { int u = queue[h++]; for (int a = 0; a < m; a++) if (a != skip && !seen[a] && share(&S[u], &S[a])) { seen[a] = 1; queue[t++] = a; } }
    for (int a = 0; a < m; a++) if (a != skip && !seen[a]) return 0;
    return 1;
}
static Pl S[32];
static void plaqs_at(const int *v, Pl *out, int *n) {
    for (int i = 0; i < D; i++) for (int j = i + 1; j < D; j++) for (int si = 0; si < 2; si++) for (int sj = 0; sj < 2; sj++) {
        Pl p; for (int a = 0; a < D; a++) p.x[a] = v[a]; p.x[i] -= si; p.x[j] -= sj; p.i = i; p.j = j; out[(*n)++] = p;
    }
}
static void rec(int m) {
    CNT[m]++;
    if (m == NMAX) return;
    /* candidates y: plaquettes with a vertex in V0 or in a plaquette of S, not in S; distinct */
    static Pl cand[32][2048]; int nc = 0;
    Pl *C = cand[m];
    int W[64][4], nw = 0;
    for (int t = 0; t < nV0; t++) { memcpy(W[nw++], V0[t], sizeof W[0]); }
    for (int a = 0; a < m; a++) { int A[4][4]; verts(&S[a], A); for (int s = 0; s < 4; s++) memcpy(W[nw++], A[s], sizeof W[0]); }
    for (int w = 0; w < nw; w++) {
        Pl tmp[64]; int nt = 0; plaqs_at(W[w], tmp, &nt);
        for (int u = 0; u < nt; u++) {
            long cu = code(&tmp[u]); int dup = 0;
            for (int a = 0; a < m && !dup; a++) if (code(&S[a]) == cu) dup = 1;
            for (int a = 0; a < nc && !dup; a++) if (code(&C[a]) == cu) dup = 1;
            if (!dup) C[nc++] = tmp[u];
        }
    }
    for (int u = 0; u < nc; u++) {
        S[m] = C[u]; long cy = code(&C[u]);
        int ok = 1;
        for (int a = 0; a < m && ok; a++) if (code(&S[a]) > cy && valid(S, m + 1, a)) ok = 0;   /* a larger removable plaquette: not the parent */
        if (ok) rec(m + 1);
    }
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s D NMAX q|flat|bent|skew\n", argv[0]); return 1; }
    D = atoi(argv[1]); NMAX = atoi(argv[2]); OFF = 2 * NMAX + 8;
    const char *sh = argv[3];
    int o[4] = {0, 0, 0, 0};
    #define SET(k, a, b, c) do { V0[k][0] = a; V0[k][1] = b; if (D > 2) V0[k][2] = c; } while (0)
    if (!strcmp(sh, "q")) { nV0 = 4; SET(0,0,0,0); SET(1,1,0,0); SET(2,0,1,0); SET(3,1,1,0); }
    else if (!strcmp(sh, "flat")) { nV0 = 6; SET(0,0,0,0); SET(1,1,0,0); SET(2,2,0,0); SET(3,0,1,0); SET(4,1,1,0); SET(5,2,1,0); }
    else if (!strcmp(sh, "bent") && D >= 3) { nV0 = 6; SET(0,0,0,0); SET(1,1,0,0); SET(2,0,1,0); SET(3,1,1,0); SET(4,0,0,1); SET(5,1,0,1); }
    else if (!strcmp(sh, "skew") && D >= 3) { nV0 = 6; SET(0,0,0,0); SET(1,1,0,0); SET(2,1,1,0); SET(3,1,1,1); SET(4,0,1,1); SET(5,0,0,1); }
    else { fprintf(stderr, "bad shape\n"); return 1; }
    (void)o;
    rec(0);
    printf("reverse search, D=%d, root %s:", D, sh);
    for (int n = 1; n <= NMAX; n++) printf(" %lld", CNT[n]);
    printf("\n");
    return 0;
}
