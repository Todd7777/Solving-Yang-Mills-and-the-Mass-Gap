/* Referee (Section 16h): the admissible sets of Lemma mg:lem:uppertails, by a forced-link search -- a different method from the fifth referee's
 * Redelmeier enumeration (admissible_redelmeier.c / admissible_split.c), written from the definition.
 * Usage: ./adm_forced d A|Q nmax          (build: cc -O2 -o adm_forced adm_forced.c)
 *
 * Y is admissible if Y and C are disjoint, every link of E(Y)\S lies in at least two plaquettes of Y, and every component of Y (plaquettes adjacent when
 * they share a link, any link) contains a plaquette of P_S \ C.  Geometry as in the fifth referee's code: (S_A, {p0}) with p0 the plaquette at the
 * origin in directions (1,2); (S_Q, Q) with Q the unit cube at the origin in directions (0,1,2).
 * Phase 1 -- the connected admissible sets.  A state is a link-connected set P containing a root (a plaquette of P_S \ C).  If some link outside S is
 * covered exactly once ("deficient"), take the smallest one, l: every admissible Y containing P contains a further plaquette through l, so branch over
 * the non-empty sets T of the available plaquettes through l (not in P, not in C, not excluded) and exclude the others below.  If no link is deficient,
 * P is a connected admissible set: record it once (a hash table of the sorted plaquette lists) and, with the exclusions lifted, extend it by every single
 * plaquette sharing a link with it.  Every connected admissible Y with |Y| <= nmax is reached from any root it contains: along the forced steps the state
 * stays inside Y (T = the plaquettes of Y through l, and the excluded ones are not in Y), and a closed state P != Y has a plaquette of Y \ P sharing a
 * link with it.  Prunes (necessary conditions): a further plaquette removes at most four deficient links, and a set of deficient links no two of which
 * lie in a common plaquette needs that many further plaquettes.
 * Phase 2 -- the admissible sets with several components: unions of connected admissible sets that share no link (each has >= 5 plaquettes, Lemma
 * mg:lem:five, so with nmax <= 14 at most two).  The output lists the admissible sets by size and those with two components. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXD 4
#define MAXN 16
static int D, NMAX, SIDE, OFF, NPAIR;
static int pidx[MAXD][MAXD], pmu[6], pnu[6];
static long NV, NPL, NLK, stride[MAXD];
static unsigned char *inS, *inC, *inP, *excl, *cnt;
static long nodes = 0, closed_hits = 0;
static int boundary = 0;
static long cur[MAXN]; static int n = 0;
static int Dcur = 0;                         /* deficient links: outside S, covered exactly once */

static long vid(const int *x) { long v = 0; for (int i = 0; i < D; i++) v += (long)(x[i] + OFF) * stride[i]; return v; }
static void vco(long v, int *x) { for (int i = 0; i < D; i++) { x[i] = (int)(v % SIDE) - OFF; v /= SIDE; } }
static int near_edge(const int *x) { for (int i = 0; i < D; i++) if (x[i] + OFF < 2 || x[i] + OFF > SIDE - 3) return 1; return 0; }
static long plid(const int *x, int m, int k) { return vid(x) * NPAIR + pidx[m][k]; }
static long lkid(const int *x, int m) { return vid(x) * D + m; }
static void plinks(long p, long *L) {
    int x[MAXD], y[MAXD]; long v = p / NPAIR; int k = (int)(p % NPAIR), m = pmu[k], q = pnu[k];
    vco(v, x);
    L[0] = lkid(x, m); L[3] = lkid(x, q);
    memcpy(y, x, sizeof x); y[m]++; L[1] = lkid(y, q);
    memcpy(y, x, sizeof x); y[q]++; L[2] = lkid(y, m);
}
static int lplaqs(long l, long *P) {              /* the 2(D-1) plaquettes through the link l */
    int x[MAXD], y[MAXD], c = 0; long v = l / D; int m = (int)(l % D);
    vco(v, x);
    for (int k = 0; k < D; k++) if (k != m) {
        int a = m < k ? m : k, b = m < k ? k : m;
        P[c++] = plid(x, a, b);
        memcpy(y, x, sizeof x); y[k]--; P[c++] = plid(y, a, b);
    }
    return c;
}
static void addp(long p) {
    long L[4]; plinks(p, L); inP[p] = 1; cur[n++] = p;
    int x[MAXD]; vco(p / NPAIR, x); if (near_edge(x)) boundary = 1;
    for (int i = 0; i < 4; i++) { int c = ++cnt[L[i]]; if (!inS[L[i]]) { if (c == 1) Dcur++; else if (c == 2) Dcur--; } }
}
static void remp(void) {
    long p = cur[--n]; long L[4]; plinks(p, L); inP[p] = 0;
    for (int i = 0; i < 4; i++) { int c = --cnt[L[i]]; if (!inS[L[i]]) { if (c == 1) Dcur++; else if (c == 0) Dcur--; } }
}
/* two links share a plaquette */
static int share_plaq(long a, long b) {
    long P[2 * MAXD], Q[2 * MAXD]; int ka = lplaqs(a, P), kb = lplaqs(b, Q);
    for (int i = 0; i < ka; i++) for (int j = 0; j < kb; j++) if (P[i] == Q[j]) return 1;
    return 0;
}
static int indep_deficient(void) {
    long I[64]; int ni = 0;
    for (int i = 0; i < n; i++) { long L[4]; plinks(cur[i], L);
        for (int j = 0; j < 4; j++) if (!inS[L[j]] && cnt[L[j]] == 1) {
            int ok = 1; for (int t = 0; t < ni; t++) if (I[t] == L[j] || share_plaq(I[t], L[j])) { ok = 0; break; }
            if (ok) I[ni++] = L[j];
        } }
    return ni;
}
/* ---------- the closed states: hash table of sorted plaquette lists ---------- */
typedef struct { int k; long p[MAXN]; } Set;
static Set *sets; static long nsets = 0, capsets = 0;
static long *htab; static long hmask;
static uint64_t hkey(const long *s, int k) { uint64_t h = 1469598103934665603ULL ^ (uint64_t)k; for (int i = 0; i < k; i++) { h ^= (uint64_t)s[i]; h *= 1099511628211ULL; h ^= h >> 29; } return h; }
static int cmpl(const void *a, const void *b) { long x = *(const long *)a, y = *(const long *)b; return (x > y) - (x < y); }
static int record_closed(void) {                   /* 1 if new */
    long s[MAXN]; memcpy(s, cur, sizeof(long) * n); qsort(s, n, sizeof(long), cmpl);
    uint64_t h = hkey(s, n); long i = (long)(h & (uint64_t)hmask);
    while (htab[i] >= 0) {
        Set *t = &sets[htab[i]];
        if (t->k == n && !memcmp(t->p, s, sizeof(long) * n)) return 0;
        i = (i + 1) & hmask;
    }
    if (nsets == capsets) { capsets = capsets ? 2 * capsets : 1 << 16; sets = realloc(sets, sizeof(Set) * capsets); }
    sets[nsets].k = n; memcpy(sets[nsets].p, s, sizeof(long) * n); htab[i] = nsets++;
    if (2 * nsets > hmask) { fprintf(stderr, "hash table full\n"); exit(2); }
    return 1;
}
static long conn[MAXN + 1];
static long exstack[1 << 16]; static int exn = 0;
static void rec(void);
static void closed_state(void) {
    closed_hits++;
    if (!record_closed()) return;
    conn[n]++;
    if (n >= NMAX) return;
    /* lift the exclusions, then extend by every plaquette sharing a link with the set */
    int saved = exn; long *sv = malloc(sizeof(long) * (saved + 1)); memcpy(sv, exstack, sizeof(long) * saved);
    for (int i = 0; i < saved; i++) excl[sv[i]] = 0;
    exn = 0;
    long seeds[MAXN * 4 * 2 * MAXD]; int ns = 0;
    for (int i = 0; i < n; i++) { long L[4]; plinks(cur[i], L);
        for (int j = 0; j < 4; j++) { long P[2 * MAXD]; int k = lplaqs(L[j], P);
            for (int t = 0; t < k; t++) { long c = P[t]; if (inP[c] || inC[c]) continue;
                int dup = 0; for (int s = 0; s < ns; s++) if (seeds[s] == c) { dup = 1; break; }
                if (!dup) seeds[ns++] = c; } } }
    for (int s = 0; s < ns; s++) { addp(seeds[s]); rec(); remp(); }
    for (int i = 0; i < saved; i++) excl[sv[i]] = 1;
    memcpy(exstack, sv, sizeof(long) * saved); exn = saved; free(sv);
}
static void rec(void) {
    nodes++;
    if (Dcur == 0) { closed_state(); return; }
    if (Dcur > 4 * (NMAX - n)) return;
    if (indep_deficient() > NMAX - n) return;
    /* the smallest deficient link */
    long l = -1;
    for (int i = 0; i < n; i++) { long L[4]; plinks(cur[i], L); for (int j = 0; j < 4; j++) if (!inS[L[j]] && cnt[L[j]] == 1 && (l < 0 || L[j] < l)) l = L[j]; }
    long P[2 * MAXD], cand[2 * MAXD]; int k = lplaqs(l, P), nc = 0;
    for (int t = 0; t < k; t++) if (!inP[P[t]] && !inC[P[t]] && !excl[P[t]]) cand[nc++] = P[t];
    for (int mask = 1; mask < (1 << nc); mask++) {
        int pc = __builtin_popcount(mask); if (n + pc > NMAX) continue;
        int ex0 = exn;
        for (int t = 0; t < nc; t++) { if (mask >> t & 1) addp(cand[t]); else { excl[cand[t]] = 1; exstack[exn++] = cand[t]; } }
        rec();
        for (int t = nc - 1; t >= 0; t--) if (mask >> t & 1) remp();
        while (exn > ex0) excl[exstack[--exn]] = 0;
    }
}
int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s d A|Q nmax\n", argv[0]); return 1; }
    D = atoi(argv[1]); char kind = argv[2][0]; NMAX = atoi(argv[3]);
    if (NMAX >= 15 || NMAX > MAXN) { fprintf(stderr, "nmax <= 14 (at most two components)\n"); return 1; }
    OFF = NMAX + 3; SIDE = 2 * OFF + 3;
    NPAIR = 0; for (int a = 0; a < D; a++) for (int b = a + 1; b < D; b++) { pidx[a][b] = NPAIR; pmu[NPAIR] = a; pnu[NPAIR] = b; NPAIR++; }
    NV = 1; for (int i = 0; i < D; i++) { stride[i] = NV; NV *= SIDE; }
    NPL = NV * NPAIR; NLK = NV * D;
    inS = calloc(NLK, 1); cnt = calloc(NLK, 1); inC = calloc(NPL, 1); inP = calloc(NPL, 1); excl = calloc(NPL, 1);
    hmask = (1L << 22) - 1; htab = malloc(sizeof(long) * (hmask + 1)); memset(htab, 0xff, sizeof(long) * (hmask + 1));
    int O[MAXD] = {0};
    long C[6]; int nC = 0;
    if (kind == 'A') C[nC++] = plid(O, 1, 2);
    else { int pr[3][3] = {{0, 1, 2}, {0, 2, 1}, {1, 2, 0}};
        for (int t = 0; t < 3; t++) { int x[MAXD] = {0}; C[nC++] = plid(x, pr[t][0], pr[t][1]); x[pr[t][2]] = 1; C[nC++] = plid(x, pr[t][0], pr[t][1]); } }
    for (int i = 0; i < nC; i++) { inC[C[i]] = 1; long L[4]; plinks(C[i], L); for (int j = 0; j < 4; j++) inS[L[j]] = 1; }
    long roots[256]; int nr = 0;
    for (int i = 0; i < nC; i++) { long L[4]; plinks(C[i], L);
        for (int j = 0; j < 4; j++) { long P[2 * MAXD]; int k = lplaqs(L[j], P);
            for (int t = 0; t < k; t++) { if (inC[P[t]]) continue; int dup = 0; for (int s = 0; s < nr; s++) if (roots[s] == P[t]) dup = 1; if (!dup) roots[nr++] = P[t]; } } }
    for (int r = 0; r < nr; r++) { addp(roots[r]); rec(); remp(); }
    /* phase 2: two link-disjoint connected admissible sets */
    long two[MAXN + 1] = {0};
    long (*links)[4 * MAXN] = malloc(sizeof(*links) * nsets);
    for (long i = 0; i < nsets; i++) { for (int a = 0; a < sets[i].k; a++) plinks(sets[i].p[a], &links[i][4 * a]); qsort(links[i], 4 * sets[i].k, sizeof(long), cmpl); }
    for (long i = 0; i < nsets; i++) for (long j = i + 1; j < nsets; j++) {
        int s = sets[i].k + sets[j].k; if (s > NMAX) continue;
        long *a = links[i], *b = links[j]; int na = 4 * sets[i].k, nb = 4 * sets[j].k, x = 0, y = 0, common = 0;
        while (x < na && y < nb) { if (a[x] == b[y]) { common = 1; break; } if (a[x] < b[y]) x++; else y++; }
        if (!common) two[s]++;
    }
    printf("d=%d kind=%c nmax=%d |C|=%d J'=%d nodes=%ld closed_hits=%ld connected_sets=%ld boundary=%d\n", D, kind, NMAX, nC, nr, nodes, closed_hits, nsets, boundary);
    printf("admissible sets by size:"); for (int s = 1; s <= NMAX; s++) if (conn[s] + two[s]) printf(" %d:%ld", s, conn[s] + two[s]);
    printf("\nof which with two components:"); for (int s = 1; s <= NMAX; s++) if (two[s]) printf(" %d:%ld", s, two[s]);
    printf("\n");
    return 0;
}
