/* Section 16 of Part IV: closed polymers through a plaquette of Z^d, closed_sets_ext.py in C.
 * Usage: closed_sets d NMAX          (build: cc -O2 -o closed_sets closed_sets.c)
 *
 * Counts the link-connected sets Y of plaquettes of Z^d containing a fixed plaquette p0, |Y| <= NMAX, in which every link of E(Y) lies in at least
 * two plaquettes of Y (the closed polymers of Lemma mg:lem:theta7; polymers are connected through shared links). The search is that of
 * closed_sets_ext.py: at the smallest link lying in exactly one plaquette, branch over the nonempty sets of further plaquettes through it (the others
 * through it are excluded below); at a closed set, record it (once: the closed sets are stored with their sorted plaquettes) and add each plaquette
 * sharing a link with it, the exclusions cleared. A link-connected closed Y with a closed proper subset X contains a plaquette outside X sharing a link
 * with X, so every Y is reached. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static int d, NMAX, R, W, NO;
static long stride[8], NVERT, NLINK, NPLAQ, ORIG;
static int oi[32], oj[32], orient[8][8];
static unsigned char *cnt; static long tl[4096]; static long ntl;
static char *inY, *exY; static long Yl[64]; static int nY;
static long exst[8192]; static int nex;
static long long found_n[64];

typedef struct { uint64_t h; long *key; } Slot;
typedef struct { Slot *s; long cap, n; } Set;
static Set found;
static uint64_t hkey(const long *k, long len) {
    uint64_t h = 1469598103934665603ULL;
    for (long i = 0; i < len; i++) { h ^= (uint64_t)k[i] + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); h *= 1099511628211ULL; }
    return h ? h : 1;
}
static int set_insert(Set *S, const long *k, long len) {
    if (2 * (S->n + 1) > S->cap) {
        long nc = S->cap ? 2 * S->cap : 1 << 16; Slot *ns = calloc(nc, sizeof(Slot));
        for (long i = 0; i < S->cap; i++) if (S->s[i].h) { long j = S->s[i].h & (nc - 1); while (ns[j].h) j = (j + 1) & (nc - 1); ns[j] = S->s[i]; }
        free(S->s); S->s = ns; S->cap = nc;
    }
    uint64_t h = hkey(k, len); long j = h & (S->cap - 1);
    while (S->s[j].h) {
        if (S->s[j].h == h && S->s[j].key[0] == len && memcmp(S->s[j].key + 1, k, len * sizeof(long)) == 0) return 0;
        j = (j + 1) & (S->cap - 1);
    }
    S->s[j].h = h; S->s[j].key = malloc((len + 1) * sizeof(long)); S->s[j].key[0] = len; memcpy(S->s[j].key + 1, k, len * sizeof(long));
    S->n++; return 1;
}
static inline void plinks(long P, long *L) {
    long v = P / NO; int o = P % NO, i = oi[o], j = oj[o];
    L[0] = v * d + i; L[1] = (v + stride[i]) * d + j; L[2] = (v + stride[j]) * d + i; L[3] = v * d + j;
}
static int plaqs_through(long L, long *out) {
    long v = L / d; int k = L % d, c = 0;
    for (int m = 0; m < d; m++) {
        if (m == k) continue;
        int i = k < m ? k : m, j = k < m ? m : k, o = orient[i][j];
        out[c++] = v * NO + o; out[c++] = (v - stride[m]) * NO + o;
    }
    return c;
}
static inline void inc_link(long L) { if (cnt[L]++ == 0) tl[ntl++] = L; }
static inline void dec_link(long L) { if (--cnt[L] == 0) { ntl--; if (tl[ntl] != L) { fprintf(stderr, "undo order\n"); exit(2); } } }
static void add_p(long P) { long L[4]; plinks(P, L); for (int a = 0; a < 4; a++) inc_link(L[a]); inY[P] = 1; Yl[nY++] = P; }
static void del_p(long P) { long L[4]; plinks(P, L); for (int a = 3; a >= 0; a--) dec_link(L[a]); inY[P] = 0; nY--; }
static int cmpl(const void *a, const void *b) { long x = *(const long *)a, y = *(const long *)b; return (x > y) - (x < y); }
static long *pstamp; static long pstampc;

static void rec(void) {
    long ndef = 0, lmin = -1;
    for (long a = 0; a < ntl; a++) if (cnt[tl[a]] == 1) { ndef++; if (lmin < 0 || tl[a] < lmin) lmin = tl[a]; }
    if (ndef == 0) {
        long key[64]; for (int a = 0; a < nY; a++) key[a] = Yl[a];
        qsort(key, nY, sizeof(long), cmpl);
        if (!set_insert(&found, key, nY)) return;
        found_n[nY]++;
        if (nY >= NMAX) return;
        /* the plaquettes sharing a link with Y, each once */
        long seeds[4096]; int ns = 0; pstampc++;
        for (int a = 0; a < nY; a++) {
            long L[4]; plinks(Yl[a], L);
            for (int c = 0; c < 4; c++) { long th[16]; int k = plaqs_through(L[c], th); for (int x = 0; x < k; x++) if (!inY[th[x]] && pstamp[th[x]] != pstampc) { pstamp[th[x]] = pstampc; seeds[ns++] = th[x]; } }
        }
        int save = nex; long saved[8192]; memcpy(saved, exst, nex * sizeof(long));
        for (int a = 0; a < nex; a++) exY[exst[a]] = 0;
        nex = 0;
        for (int a = 0; a < ns; a++) { add_p(seeds[a]); rec(); del_p(seeds[a]); }
        nex = save; memcpy(exst, saved, save * sizeof(long)); for (int a = 0; a < nex; a++) exY[exst[a]] = 1;
        return;
    }
    if (ndef > 4 * (NMAX - nY)) return;           /* each further plaquette covers at most four of these links */
    long th[16]; int k = plaqs_through(lmin, th);
    long cand[16]; int nc = 0;
    for (int x = 0; x < k; x++) if (!inY[th[x]] && !exY[th[x]]) cand[nc++] = th[x];
    int room = NMAX - nY;
    for (int mask = 1; mask < (1 << nc); mask++) {
        int r = __builtin_popcount(mask);
        if (r > room) continue;
        for (int a = 0; a < nc; a++) if (mask >> a & 1) add_p(cand[a]);
        int nex0 = nex;
        for (int a = 0; a < nc; a++) if (!(mask >> a & 1)) { exY[cand[a]] = 1; exst[nex++] = cand[a]; }
        rec();
        while (nex > nex0) { nex--; exY[exst[nex]] = 0; }
        for (int a = nc - 1; a >= 0; a--) if (mask >> a & 1) del_p(cand[a]);
    }
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s d NMAX\n", argv[0]); return 1; }
    d = atoi(argv[1]); NMAX = atoi(argv[2]);
    NO = 0; for (int i = 0; i < d; i++) for (int j = i + 1; j < d; j++) { oi[NO] = i; oj[NO] = j; orient[i][j] = NO; NO++; }
    R = NMAX + 4; W = 2 * R + 1;
    NVERT = 1; for (int i = 0; i < d; i++) { stride[i] = NVERT; NVERT *= W; }
    NLINK = NVERT * d; NPLAQ = NVERT * NO;
    ORIG = 0; for (int i = 0; i < d; i++) ORIG += (long)R * stride[i];
    cnt = calloc(NLINK, 1); inY = calloc(NPLAQ, 1); exY = calloc(NPLAQ, 1); pstamp = calloc(NPLAQ, sizeof(long));
    if (!cnt || !inY || !exY || !pstamp) { fprintf(stderr, "out of memory\n"); return 1; }
    long p0 = ORIG * NO + orient[0][1];
    add_p(p0);
    rec();
    printf("d=%d, NMAX=%d: closed link-connected sets through a plaquette by size:", d, NMAX);
    for (int n = 1; n <= NMAX; n++) if (found_n[n]) printf(" %d:%lld", n, found_n[n]);
    printf("\n");
    return 0;
}
