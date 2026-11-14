/* admissible_split.c: the fifth referee's referee-section16e/admissible_redelmeier.c (Lemma mg:lem:uppertails, the admissible sets), unchanged except for a
 * split into independent first-level branches, for counts that take hours (Remark mg:rem:su5(iv): 11 plaquettes in Z^4).
 * Usage: ./admissible_split d A|Q nmax [prune [k J]]     (build: cc -O2 -o admissible_split admissible_split.c)
 * With k J, only the first-level branches of Redelmeier's search with index i = k mod J are run: the branch of the root plaquette roots[i] has the untried
 * set roots[0..i), and the roots after it stay marked, exactly as in the unsplit search; so the counts of the J pieces add up to the unsplit counts (checked:
 * A and Q in Z^4 to 9 plaquettes, J = 2, 3). run_admissible_pieces.sh runs the pieces resumably.
 */
/* Referee check (Part IV, Lemma mg:lem:uppertails, Remark mg:rem:upperexact): the admissible sets by a second, independent method.
 *
 * Usage: ./admissible_redelmeier d A|Q nmax [prune [k J]]   (split copy: only the first-level branches with index = k mod J)
 *
 * Redelmeier's algorithm enumerates, exactly once each, the sets Y of plaquettes of Z^d (outside C) such that Y + {o} is connected in
 * the graph where plaquettes sharing a link are adjacent and a virtual root o is adjacent to every plaquette of P_S \ C.  These are
 * exactly the sets satisfying conditions (1) and (3) of the definition of "admissible" (Lemma mg:lem:uppertails): Y and C disjoint,
 * every component of Y contains a plaquette of P_S \ C.  Multi-component sets are generated directly.  At every node condition (2),
 * "every link of E(Y) \ S lies in at least two plaquettes of Y", is tested by the number D of links of E(Y) \ S covered exactly once.
 * Pruning: a superset needs at least ceil(D/4) further plaquettes (each covers at most four links); with a fourth argument "2" also at least
 * |I| further plaquettes, for any set I of deficient links no two of which lie in a common plaquette (a greedy such I; each of its links needs
 * its own new plaquette).  Plaquettes of C are pre-marked.
 * The geometry (p0, Q, S, C) is set up as in admissible_partition.py: p0 at the origin in directions (1,2), Q the unit cube at the
 * origin in directions (0,1,2).  The box is large enough that no set of nmax plaquettes reaches its boundary (checked).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXD 4
static int D_, NMAX, SIDE, OFF, NPAIR;
static int pairidx[MAXD][MAXD], pairmu[6], pairnu[6];
static long NV, NPL, NLK;
static long stride[MAXD];
static unsigned char *mark, *isS, *inC;
static unsigned char *cnt;
static long counts[64], counts2[64], nodes = 0;
static int Dcur = 0;
static long cur[64];
static int boundary_hit = 0;
static int PRUNE2 = 0;
static long Ycur_links_count = 0;

static long vid(const int *x) { long v = 0; for (int i = 0; i < D_; i++) v += (long)(x[i] + OFF) * stride[i]; return v; }
static void vcoord(long v, int *x) { for (int i = 0; i < D_; i++) { x[i] = (int)(v % SIDE) - OFF; v /= SIDE; } }
static long plid(const int *x, int m, int n) { return vid(x) * NPAIR + pairidx[m][n]; }
static long lkid(const int *x, int m) { return vid(x) * D_ + m; }
static int inbox(const int *x) { for (int i = 0; i < D_; i++) if (x[i] + OFF < 1 || x[i] + OFF > SIDE - 3) return 0; return 1; }

static void plinks(long p, long *L) {
    int x[MAXD], y[MAXD]; long v = p / NPAIR; int k = (int)(p % NPAIR), m = pairmu[k], n = pairnu[k];
    vcoord(v, x);
    L[0] = lkid(x, m);
    memcpy(y, x, sizeof x); y[m]++; L[1] = lkid(y, n);
    memcpy(y, x, sizeof x); y[n]++; L[2] = lkid(y, m);
    L[3] = lkid(x, n);
}
/* plaquettes containing link l */
static int lplaqs(long l, long *P) {
    int x[MAXD], y[MAXD], c = 0; long v = l / D_; int m = (int)(l % D_);
    vcoord(v, x);
    for (int n = 0; n < D_; n++) if (n != m) {
        int a = m < n ? m : n, b = m < n ? n : m;
        P[c++] = plid(x, a, b);
        memcpy(y, x, sizeof x); y[n]--; P[c++] = plid(y, a, b);
    }
    return c;
}
static int neighbours(long p, long *Nb) {
    long L[4], P[2 * MAXD]; int c = 0;
    plinks(p, L);
    for (int i = 0; i < 4; i++) {
        int k = lplaqs(L[i], P);
        for (int j = 0; j < k; j++) if (P[j] != p) {
            int dup = 0; for (int t = 0; t < c; t++) if (Nb[t] == P[j]) { dup = 1; break; }
            if (!dup) Nb[c++] = P[j];
        }
    }
    return c;
}
/* two links (coordinates x, y; directions m, n) lie in a common plaquette */
static int common_plaquette_xy(const int *x, int m, const int *y, int n) {
    if (m != n) {      /* perpendicular: they must share a vertex */
        for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) {
            int same = 1;
            for (int i = 0; i < D_; i++) if (x[i] + (a && i == m) != y[i] + (b && i == n)) { same = 0; break; }
            if (same) return 1;
        }
        return 0;
    }
    /* parallel: opposite edges of a plaquette, y = x +- e_l with l != m */
    int diff = 0, dl = -1;
    for (int i = 0; i < D_; i++) if (x[i] != y[i]) { diff++; dl = i; }
    return diff == 1 && dl != m && abs(x[dl] - y[dl]) == 1;
}
/* greedy lower bound: deficient links (outside S, covered once) pairwise in no common plaquette */
static int indep_bound(int size) {
    int ni = 0; int Ix[64][MAXD]; int Im[64];
    for (int i = 0; i < size; i++) { long L[4]; plinks(cur[i], L);
        for (int j = 0; j < 4; j++) if (!isS[L[j]] && cnt[L[j]] == 1) {
            int x[MAXD]; vcoord(L[j] / D_, x); int m = (int)(L[j] % D_);
            int ok = 1;
            for (int k = 0; k < ni; k++) if (common_plaquette_xy(x, m, Ix[k], Im[k])) { ok = 0; break; }
            if (ok) { memcpy(Ix[ni], x, sizeof x); Im[ni] = m; ni++; }
        } }
    return ni;
}
static void addp(long p) {
    long L[4]; plinks(p, L);
    for (int i = 0; i < 4; i++) { unsigned char c = ++cnt[L[i]]; if (!isS[L[i]]) { if (c == 1) Dcur++; else if (c == 2) Dcur--; } }
}
static void remp(long p) {
    long L[4]; plinks(p, L);
    for (int i = 0; i < 4; i++) { unsigned char c = --cnt[L[i]]; if (!isS[L[i]]) { if (c == 1) Dcur++; else if (c == 0) Dcur--; } }
}
/* number of link-connected components of cur[0..n) */
static int ncomp(int n) {
    int par[64]; for (int i = 0; i < n; i++) par[i] = i;
    long Li[64][4]; for (int i = 0; i < n; i++) plinks(cur[i], Li[i]);
    for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) {
        int share = 0; for (int a = 0; a < 4 && !share; a++) for (int b = 0; b < 4; b++) if (Li[i][a] == Li[j][b]) { share = 1; break; }
        if (share) { int ri = i, rj = j; while (par[ri] != ri) ri = par[ri]; while (par[rj] != rj) rj = par[rj]; if (ri != rj) par[ri] = rj; }
    }
    int c = 0; for (int i = 0; i < n; i++) if (par[i] == i) c++;
    return c;
}
static void rec(long *untried, int nu, int size) {
    long *nw = malloc(sizeof(long) * (nu + 64 * 24));
    while (nu > 0) {
        long p = untried[--nu];
        cur[size] = p; addp(p); nodes++;
        int s1 = size + 1;
        if (Dcur == 0) { counts[s1]++; if (ncomp(s1) == 2) counts2[s1]++; }
        if (s1 < NMAX && s1 + (Dcur + 3) / 4 <= NMAX && (!PRUNE2 || s1 + indep_bound(s1) <= NMAX)) {
            memcpy(nw, untried, sizeof(long) * nu);
            int nn = nu, start = nu;
            long Nb[64]; int k = neighbours(p, Nb);
            for (int j = 0; j < k; j++) if (!mark[Nb[j]]) {
                int x[MAXD]; vcoord(Nb[j] / NPAIR, x);
                if (!inbox(x)) { boundary_hit = 1; continue; }
                mark[Nb[j]] = 1; nw[nn++] = Nb[j];
            }
            rec(nw, nn, s1);
            for (int j = start; j < nn; j++) mark[nw[j]] = 0;
        }
        remp(p);
    }
    free(nw);
}
/* one iteration of rec's loop at size 0: take p = untried[nu-1], with untried[0..nu-1) left for its subtree */
static void rec_first(long *untried, int nu) {
    long *nw = malloc(sizeof(long) * (nu + 64 * 24));
    long p = untried[--nu]; int size = 0;
    cur[size] = p; addp(p); nodes++;
    int s1 = size + 1;
    if (Dcur == 0) { counts[s1]++; if (ncomp(s1) == 2) counts2[s1]++; }
    if (s1 < NMAX && s1 + (Dcur + 3) / 4 <= NMAX && (!PRUNE2 || s1 + indep_bound(s1) <= NMAX)) {
        memcpy(nw, untried, sizeof(long) * nu);
        int nn = nu, start = nu;
        long Nb[64]; int k = neighbours(p, Nb);
        for (int j = 0; j < k; j++) if (!mark[Nb[j]]) {
            int x[MAXD]; vcoord(Nb[j] / NPAIR, x);
            if (!inbox(x)) { boundary_hit = 1; continue; }
            mark[Nb[j]] = 1; nw[nn++] = Nb[j];
        }
        rec(nw, nn, s1);
        for (int j = start; j < nn; j++) mark[nw[j]] = 0;
    }
    remp(p);
    free(nw);
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s d A|Q nmax\n", argv[0]); return 1; }
    D_ = atoi(argv[1]); char kind = argv[2][0]; NMAX = atoi(argv[3]); if (argc > 4) PRUNE2 = (atoi(argv[4]) == 2);
    OFF = NMAX + 3; SIDE = 2 * OFF + 3;
    NPAIR = 0; for (int m = 0; m < D_; m++) for (int n = m + 1; n < D_; n++) { pairidx[m][n] = NPAIR; pairmu[NPAIR] = m; pairnu[NPAIR] = n; NPAIR++; }
    NV = 1; for (int i = 0; i < D_; i++) { stride[i] = NV; NV *= SIDE; }
    NPL = NV * NPAIR; NLK = NV * D_;
    mark = calloc(NPL, 1); inC = calloc(NPL, 1); isS = calloc(NLK, 1); cnt = calloc(NLK, 1);
    int O[MAXD] = {0, 0, 0, 0};
    long C[6]; int nC = 0;
    if (kind == 'A') C[nC++] = plid(O, 1, 2);
    else {
        int pr[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
        for (int t = 0; t < 3; t++) { int x[MAXD] = {0, 0, 0, 0}; C[nC++] = plid(x, pr[t][0], pr[t][1]); x[pr[t][2]] = 1; C[nC++] = plid(x, pr[t][0], pr[t][1]); }
    }
    for (int i = 0; i < nC; i++) { inC[C[i]] = 1; mark[C[i]] = 1; long L[4]; plinks(C[i], L); for (int j = 0; j < 4; j++) isS[L[j]] = 1; }
    /* P_S \ C, the neighbours of the root */
    long *roots = malloc(sizeof(long) * 1000); int nr = 0, nPS = 0;
    for (int i = 0; i < nC; i++) {
        long L[4]; plinks(C[i], L);
        for (int j = 0; j < 4; j++) {
            long P[2 * MAXD]; int k = lplaqs(L[j], P);
            for (int t = 0; t < k; t++) {
                int dup = 0; for (int s = 0; s < nr; s++) if (roots[s] == P[t]) dup = 1;
                if (!dup && !inC[P[t]]) { roots[nr++] = P[t]; mark[P[t]] = 1; }
            }
        }
    }
    /* count |P_S| for the record */
    nPS = nr + nC;
    /* split (argv 5, 6 = k J): only the first-level branches with index i = k mod J (Redelmeier: the branch of roots[i] has the untried set roots[0..i)) */
    int SPK = argc > 6 ? atoi(argv[5]) : 0, SPJ = argc > 6 ? atoi(argv[6]) : 1;
    if (SPJ == 1) rec(roots, nr, 0);
    else {
        for (int i = nr - 1; i >= 0; i--) {
            if (i % SPJ != SPK) continue;
            long one[1]; one[0] = roots[i];
            /* the branch of roots[i]: untried = roots[0..i) and roots[i] taken first; emulate rec's first iteration */
            long *u = malloc(sizeof(long) * (i + 1)); memcpy(u, roots, sizeof(long) * i); u[i] = roots[i];
            rec_first(u, i + 1);
            free(u);
        }
    }
    printf("d=%d kind=%c nmax=%d prune=%d |C|=%d |P_S|=%d J'=%d nodes=%ld boundary_hit=%d\n", D_, kind, NMAX, PRUNE2 ? 2 : 1, nC, nPS, nr, nodes, boundary_hit);
    printf("admissible sets by size:");
    for (int n = 1; n <= NMAX; n++) if (counts[n]) printf(" %d:%ld", n, counts[n]);
    printf("\nof which with two components:");
    for (int n = 1; n <= NMAX; n++) if (counts2[n]) printf(" %d:%ld", n, counts2[n]);
    printf("\n");
    return 0;
}
