"""Independent exact verifier for the Z^4 block-comparison links.

Written from the definitions in the task statement only.

Conventions
-----------
A site x = (x1,x2,x3,x4) of B = {0,1}^4 is encoded as the integer
    idx(x) = 8*x1 + 4*x2 + 2*x3 + x4,
so coordinate x_i is bit (4-i) of idx, and iota(C) = min_{x in C} idx(x).
Subsets of B are 16-bit masks.

Faces F_{i,s} (i in 1..4, s in {0,1}) are numbered f = 2*(i-1) + s.
Opposite face of f is f ^ 1.

All probabilities are handled as exact integers: for p = a/b,
    pi_p(eta) * b^16 = a^|eta| (b-a)^(16-|eta|).
"""
from fractions import Fraction
from collections import defaultdict

NSITES = 16
FULL = (1 << NSITES) - 1


def coords(idx):
    return tuple((idx >> (4 - i)) & 1 for i in range(1, 5))


def idx_of(x):
    return 8 * x[0] + 4 * x[1] + 2 * x[2] + x[3]


FACES = [(i, s) for i in range(1, 5) for s in (0, 1)]  # f -> (i,s)


def face_mask(f):
    i, s = FACES[f]
    m = 0
    for v in range(NSITES):
        if coords(v)[i - 1] == s:
            m |= 1 << v
    return m


FACE_MASK = [face_mask(f) for f in range(8)]


def opposite(f):
    i, s = FACES[f]
    return FACES.index((i, 1 - s))


def mirror_site(f, y):
    """m_F : Fbar -> F, replace the i-th coordinate by s."""
    i, s = FACES[f]
    x = list(coords(y))
    assert x[i - 1] == 1 - s, "mirror applied to a site not on the opposite face"
    x[i - 1] = s
    return idx_of(x)


def mirror_mask(f, ymask):
    out = 0
    v = ymask
    while v:
        low = v & -v
        y = low.bit_length() - 1
        out |= 1 << mirror_site(f, y)
        v ^= low
    return out


# local indexing of a face: the sites of face f listed in increasing idx order
FACE_SITES = [[v for v in range(NSITES) if FACE_MASK[f] >> v & 1] for f in range(8)]


def to_local(f, mask):
    """global mask (subset of face f) -> 8-bit local mask."""
    assert mask & ~FACE_MASK[f] == 0
    out = 0
    for k, v in enumerate(FACE_SITES[f]):
        if mask >> v & 1:
            out |= 1 << k
    return out


def from_local(f, lmask):
    out = 0
    for k, v in enumerate(FACE_SITES[f]):
        if lmask >> k & 1:
            out |= 1 << v
    return out


NEIGH = [[v ^ (1 << k) for k in range(4)] for v in range(NSITES)]


def clusters(eta):
    """open clusters of eta (connectivity inside B) as list of masks."""
    seen = 0
    out = []
    for v in range(NSITES):
        if (eta >> v & 1) and not (seen >> v & 1):
            comp = 1 << v
            stack = [v]
            while stack:
                u = stack.pop()
                for w in NEIGH[u]:
                    if (eta >> w & 1) and not (comp >> w & 1):
                        comp |= 1 << w
                        stack.append(w)
            seen |= comp
            out.append(comp)
    return out


def popcount(m):
    return bin(m).count("1")


def phi(C):
    return sum(popcount(C & FACE_MASK[f]) for f in range(8))


def iota(C):
    return (C & -C).bit_length() - 1  # index of lowest set bit = min idx


def designated(eta, minsize=3):
    cl = clusters(eta)
    if not cl:
        return 0
    best = max(cl, key=lambda C: (popcount(C), phi(C), -iota(C)))
    if popcount(best) < minsize:
        return 0
    return best


def enumerate_G(minsize=3):
    """returns dict G -> list of 17 counts (number of eta with G(eta)=G and |eta|=n)."""
    table = defaultdict(lambda: [0] * (NSITES + 1))
    ties = 0
    for eta in range(1 << NSITES):
        G = designated(eta, minsize)
        table[G][popcount(eta)] += 1
    return dict(table)


def count_size_ties(minsize=3):
    """number of eta where the largest cluster (size>=minsize) is not unique in size."""
    n = 0
    for eta in range(1 << NSITES):
        cl = clusters(eta)
        if not cl:
            continue
        sizes = sorted((popcount(C) for C in cl), reverse=True)
        if sizes[0] >= minsize and len(sizes) > 1 and sizes[1] == sizes[0]:
            n += 1
    return n


# ----------------------------------------------------------------------------
# exact probabilities at rational p = a/b
# ----------------------------------------------------------------------------

class Model:
    """Exact quantities for density p = Fraction(a, b)."""

    def __init__(self, Gtable, p):
        p = Fraction(p)
        self.p = p
        a, b = p.numerator, p.denominator
        self.a, self.b = a, b
        self.B16 = b ** NSITES
        apow = [a ** n for n in range(NSITES + 1)]
        cpow = [(b - a) ** n for n in range(NSITES + 1)]
        w = [apow[n] * cpow[NSITES - n] for n in range(NSITES + 1)]
        # Pn[G] = P_p(G(eta)=G) * b^16
        self.Pn = {G: sum(c[n] * w[n] for n in range(NSITES + 1) if c[n]) for G, c in Gtable.items()}
        assert sum(self.Pn.values()) == self.B16
        self.Glist = sorted(self.Pn)
        self._compute_r()

    def _compute_r(self):
        """R[f][lT] = r_f(T) * b^16 for every T subset of face f (lT = local mask)."""
        self.R = []
        for f in range(8):
            fb = opposite(f)
            # distribution of W = m_f(T_fbar(eps')) as local masks of face f
            agg = [0] * 256
            for G, P in self.Pn.items():
                Wg = mirror_mask(f, G & FACE_MASK[fb])
                agg[to_local(f, Wg)] += P
            # zeta: Z[V] = sum_{W subset V} agg[W]
            Z = agg[:]
            for k in range(8):
                bit = 1 << k
                for V in range(256):
                    if V & bit:
                        Z[V] += Z[V ^ bit]
            # r_f(T) = 1 - P(W subset complement(T))
            Rf = [self.B16 - Z[255 ^ lT] for lT in range(256)]
            assert Rf[0] == 0
            self.R.append(Rf)

    def r_exact(self, f, Tmask):
        return Fraction(self.R[f][to_local(f, Tmask)], self.B16)

    def entry_laws(self, f):
        """For entry face f returns (children, laws) where
        children = list of the 7 child faces (sigma bit j <-> children[j]),
        laws[lU] = (Nnum, Zn) for every local nonempty U (1..255):
            Nnum[sigma] = b^128 * Z_U * nu_{f,U}(sigma),   Zn = b^16 * Z_U,
        with Z_U = P_p(T_f(eta) cap U != empty).
        """
        children = [g for g in range(8) if g != f]
        B16 = self.B16
        M = [None] * 256  # M[lS] = 128-vector, sum over G with footprint S on f
        MP = [0] * 256
        Rch = [self.R[g] for g in children]
        loc_cache = {}
        for G in self.Glist:
            S = G & FACE_MASK[f]
            if not S:
                continue
            P = self.Pn[G]
            lS = to_local(f, S)
            vec = [P]
            for j, g in enumerate(children):
                key = (g, G & FACE_MASK[g])
                lT = loc_cache.get(key)
                if lT is None:
                    lT = to_local(g, key[1])
                    loc_cache[key] = lT
                r = Rch[j][lT]
                rc = B16 - r
                vec = [v * rc for v in vec] + [v * r for v in vec]
            if M[lS] is None:
                M[lS] = vec
            else:
                M[lS] = [x + y for x, y in zip(M[lS], vec)]
            MP[lS] += P
        zero = [0] * 128
        Z = [m if m is not None else zero for m in M]
        ZP = MP[:]
        for k in range(8):
            bit = 1 << k
            for V in range(256):
                if V & bit:
                    Z[V] = [x + y for x, y in zip(Z[V], Z[V ^ bit])]
                    ZP[V] += ZP[V ^ bit]
        Tot, TotP = Z[255], ZP[255]
        laws = [None] * 256
        b112 = B16 ** 7
        for lU in range(1, 256):
            C = Z[255 ^ lU]
            Nnum = [x - y for x, y in zip(Tot, C)]
            Zn = TotP - ZP[255 ^ lU]
            assert Zn > 0
            assert sum(Nnum) == Zn * b112, "normalisation self-check failed"
            laws[lU] = (Nnum, Zn)
        return children, laws


# ----------------------------------------------------------------------------
# minimum of a modular function over nontrivial up-sets of {0,1}^k (exact)
# ----------------------------------------------------------------------------

def _maxflow(n, s, t, edges, inf):
    """Dinic max-flow with exact integer capacities. edges: list of (u,v,cap).
    returns (flow, set of nodes reachable from s in final residual graph)."""
    to, cap, adj = [], [], [[] for _ in range(n)]
    for u, v, c in edges:
        adj[u].append(len(to)); to.append(v); cap.append(c)
        adj[v].append(len(to)); to.append(u); cap.append(0)
    flow = 0
    while True:
        level = [-1] * n
        level[s] = 0
        q = [s]
        for u in q:
            for e in adj[u]:
                if cap[e] > 0 and level[to[e]] < 0:
                    level[to[e]] = level[u] + 1
                    q.append(to[e])
        if level[t] < 0:
            reach = {u for u in range(n) if level[u] >= 0}
            # certificate check: recover the flow and verify it is feasible
            net = [0] * n
            for idx, (u, v, c) in enumerate(edges):
                fe = c - cap[2 * idx]
                assert 0 <= fe <= c and fe == cap[2 * idx + 1], "infeasible flow"
                net[u] += fe
                net[v] -= fe
            assert net[s] == flow and net[t] == -flow
            assert all(net[u] == 0 for u in range(n) if u != s and u != t), "conservation"
            # and the residual-reachable set is a cut of capacity == flow
            cutcap = sum(c for (u, v, c) in edges if u in reach and v not in reach)
            assert cutcap == flow, "max-flow/min-cut certificate failed"
            return flow, reach
        it = [0] * n

        def dfs(u, f):
            if u == t:
                return f
            while it[u] < len(adj[u]):
                e = adj[u][it[u]]
                v = to[e]
                if cap[e] > 0 and level[v] == level[u] + 1:
                    d = dfs(v, f if f < cap[e] else cap[e])
                    if d > 0:
                        cap[e] -= d
                        cap[e ^ 1] += d
                        return d
                it[u] += 1
            return 0

        while True:
            f = dfs(s, inf)
            if f == 0:
                break
            flow += f


def min_nontrivial_upset(W, k):
    """W: list of 2^k exact numbers (int), sum(W) == 0 expected.
    Returns (minval, A) where minval = min over up-sets A of {0,1}^k with
    A != empty, A != whole cube of sum_{x in A} W[x], and A is a minimiser (set of ints).

    Nontrivial up-sets are exactly A = {top} u A' with A' an up-set of the
    middle poset {0,1}^k minus {bottom, top}. Minimising a modular function
    over up-sets (closures) is solved exactly by a min s-t cut (Picard 1976).
    """
    n = 1 << k
    top = n - 1
    s, t = n, n + 1
    INF = sum(abs(w) for w in W) + 1  # exceeds the cut {s}, so never in a min cut
    edges = []
    pos = 0
    for x in range(1, top):
        c = -W[x]  # maximise sum of c over closures
        if c > 0:
            edges.append((s, x, c))
            pos += c
        elif c < 0:
            edges.append((x, t, -c))
        for j in range(k):
            y = x | (1 << j)
            if y != x and y != top:
                edges.append((x, y, INF))  # x in A' => y in A'
    flow, reach = _maxflow(n + 2, s, t, edges, INF)
    Aprime = {x for x in reach if 1 <= x < top}
    maxclosure = pos - flow
    assert -sum(W[x] for x in Aprime) == maxclosure, "duality check failed"
    A = Aprime | {top}
    # sanity: A is an up-set
    for x in A:
        for j in range(k):
            assert (x | (1 << j)) in A
    assert 0 not in A
    return W[top] - maxclosure, A


def minimal_elements(A, k):
    return sorted(x for x in A if not any((x & ~(1 << j)) in A for j in range(k) if x >> j & 1))


def bern_weights(q, k):
    """returns (num list over sigma, den) with Bernoulli(q)^k(sigma) = num/den."""
    q = Fraction(q)
    c, d = q.numerator, q.denominator
    return [c ** popcount(s) * (d - c) ** (k - popcount(s)) for s in range(1 << k)], d ** k
