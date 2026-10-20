"""Self-tests for core.py.  Run:  PYTHONUTF8=1 python selftest.py"""
import itertools
import pickle
import random
import time
from fractions import Fraction

from core import (NSITES, FACES, FACE_MASK, opposite, mirror_site, mirror_mask, to_local,
                  from_local, clusters, designated, enumerate_G, popcount, Model,
                  min_nontrivial_upset, bern_weights, coords, idx_of)

T0 = time.time()
rng = random.Random(20261020)


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


# ---------------------------------------------------------------- geometry
def test_geometry():
    # 1. faces: 8 faces of 8 sites; every site lies on exactly 4 faces
    assert all(popcount(m) == 8 for m in FACE_MASK)
    for v in range(NSITES):
        assert sum(FACE_MASK[f] >> v & 1 for f in range(8)) == 4
    # 2. mirror + adjacency across blocks, checked against the actual Z^4 geometry
    for f in range(8):
        i, s = FACES[f]
        e = [0, 0, 0, 0]
        e[i - 1] = 1 if s == 1 else -1  # w = v + e_i if F = F_{i,1}, v - e_i if F = F_{i,0}
        v = (0, 0, 0, 0)
        w = tuple(e)
        fb = opposite(f)
        for x in range(NSITES):
            for y in range(NSITES):
                X = [coords(x)[k] + 2 * v[k] for k in range(4)]
                Y = [coords(y)[k] + 2 * w[k] for k in range(4)]
                adj = sum(abs(X[k] - Y[k]) for k in range(4)) == 1
                pred = (FACE_MASK[fb] >> y & 1) == 1 and mirror_site(f, y) == x
                assert adj == pred, (f, x, y)
        # mirror is a bijection Fbar -> F with inverse m_Fbar
        img = {mirror_site(f, y) for y in range(NSITES) if FACE_MASK[fb] >> y & 1}
        assert sum(1 << u for u in img) == FACE_MASK[f]
        for y in range(NSITES):
            if FACE_MASK[fb] >> y & 1:
                assert mirror_site(fb, mirror_site(f, y)) == y
        for lm in range(256):
            assert to_local(f, from_local(f, lm)) == lm
    log("geometry/mirror/adjacency OK")


# ---------------------------------------------------------------- designated cluster
def connected(mask):
    if mask == 0:
        return False
    cl = clusters(mask)
    return len(cl) == 1


def designated_alt(eta):
    """independent re-implementation with union-find and explicit tuple keys"""
    parent = list(range(NSITES))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    open_ = [v for v in range(NSITES) if eta >> v & 1]
    for u in open_:
        for v in open_:
            x, y = coords(u), coords(v)
            if sum(abs(x[k] - y[k]) for k in range(4)) == 1:
                ru, rv = find(u), find(v)
                if ru != rv:
                    parent[ru] = rv
    comps = {}
    for u in open_:
        comps.setdefault(find(u), []).append(u)
    if not comps:
        return 0
    best = None
    for c in comps.values():
        size = len(c)
        ph = 0
        for u in c:
            x = coords(u)
            for i in range(4):
                for s in (0, 1):
                    if x[i] == s:
                        ph += 1
        io = min(8 * coords(u)[0] + 4 * coords(u)[1] + 2 * coords(u)[2] + coords(u)[3] for u in c)
        key = (size, ph, -io)
        if best is None or key > best[0]:
            best = (key, c)
    if best[0][0] < 3:
        return 0
    return sum(1 << u for u in best[1])


def test_designated(Gt):
    for eta in range(1 << NSITES):
        assert designated(eta) == designated_alt(eta), eta
    log("designated cluster: two implementations agree on all 2^16 states")
    nconn = sum(1 for m in range(1 << NSITES) if popcount(m) >= 3 and connected(m))
    log(f"connected subsets of Q4 with >=3 sites: {nconn}; distinct G incl. empty: {len(Gt)}")
    assert len(Gt) == nconn + 1
    assert sum(sum(c) for c in Gt.values()) == 1 << NSITES
    assert designated((1 << NSITES) - 1) == (1 << NSITES) - 1
    # all-open block has G = B (needed for the root step)


# ---------------------------------------------------------------- r and nu by brute force
def brute_r_all(Gtab_by_eta, p, f):
    """r_f(T) for all 256 T, straight from the definition (sum over eps'), exact Fractions"""
    fb = opposite(f)
    a, b = p.numerator, p.denominator
    agg = {}
    for eps in range(1 << NSITES):
        W = mirror_mask(f, Gtab_by_eta[eps] & FACE_MASK[fb])
        n = popcount(eps)
        agg[W] = agg.get(W, 0) + a ** n * (b - a) ** (NSITES - n)
    out = {}
    for lt in range(256):
        T = from_local(f, lt)
        out[T] = Fraction(sum(w for W, w in agg.items() if W & T), b ** NSITES)
    return out


def brute_nu(Gtab_by_eta, rtab, p, f, U):
    """nu_{f,U}(sigma) straight from the definition: a sum over all 2^16 eta. exact."""
    a, b = p.numerator, p.denominator
    children = [g for g in range(8) if g != f]
    B16 = b ** NSITES
    # integer numerators of r (denominator b^16)
    Rint = {g: {T: int(v * B16) for T, v in rtab[g].items()} for g in children}
    num = [0] * 128
    Z = 0
    cache = {}
    for eta in range(1 << NSITES):
        G = Gtab_by_eta[eta]
        if not (G & FACE_MASK[f] & U):
            continue
        n = popcount(eta)
        w = a ** n * (b - a) ** (NSITES - n)
        Z += w
        vec = cache.get(G)
        if vec is None:
            rs = [Rint[g][G & FACE_MASK[g]] for g in children]
            vec = []
            for sig in range(128):
                pr = 1
                for j in range(7):
                    pr *= rs[j] if sig >> j & 1 else B16 - rs[j]
                vec.append(pr)
            cache[G] = vec
        for sig in range(128):
            num[sig] += w * vec[sig]
    # nu = (num / b^(16*8)) / (Z / b^16)
    return [Fraction(x, Z * B16 ** 7) for x in num]


def test_r_nu(Gt):
    Gtab_by_eta = [designated(eta) for eta in range(1 << NSITES)]
    for p in (Fraction(9, 25), Fraction(3583, 10000)):
        model = Model(Gt, p)
        rtab = {}
        for f in range(8):
            rtab[f] = brute_r_all(Gtab_by_eta, p, f)
            for T, v in rtab[f].items():
                assert model.r_exact(f, T) == v, (f, T)
        log(f"p={p}: all 8 x 256 values r_F(T) match a direct sum over eps' (exact)")
        cases = [(0, 1), (3, 255), (5, 0b10010110), (7, rng.randrange(1, 256)), (6, 0b00010000)]
        if p != Fraction(9, 25):
            cases = cases[:2]
        b128 = model.B16 ** 8
        cache_laws = {}
        for f, lU in cases:
            if f not in cache_laws:
                cache_laws[f] = model.entry_laws(f)
            children, laws = cache_laws[f]
            Nnum, Zn = laws[lU]
            nu_fast = [Fraction(x, b128) / Fraction(Zn, model.B16) for x in Nnum]
            nu_slow = brute_nu(Gtab_by_eta, rtab, p, f, from_local(f, lU))
            assert nu_fast == nu_slow, (f, lU)
            assert sum(nu_fast) == 1
            log(f"p={p}: nu_(F={FACES[f]}, U local={lU:08b}) matches direct sum over eta, all 128 entries exact")


# ---------------------------------------------------------------- up-set minimisation
def all_upsets(k):
    """all up-sets of {0,1}^k, as frozensets (recursive slice construction)"""
    if k == 0:
        return [frozenset(), frozenset({0})]
    prev = all_upsets(k - 1)
    out = []
    hi = 1 << (k - 1)
    for A0 in prev:
        for A1 in prev:
            if A0 <= A1:  # slice with x_k=0 must be inside slice with x_k=1
                out.append(frozenset(A0) | frozenset(x | hi for x in A1))
    return out


def is_upset(A, k):
    return all((x | (1 << j)) in A for x in A for j in range(k))


def test_upsets():
    dedekind = [2, 3, 6, 20, 168, 7581]
    for k in range(6):
        U = all_upsets(k)
        assert len(U) == dedekind[k] == len(set(U)), (k, len(U))
        if k <= 4:
            for A in U:
                assert is_upset(A, k)
            # brute force: count all up-sets among all subsets
            cnt = sum(1 for m in range(1 << (1 << k)) if is_upset({x for x in range(1 << k) if m >> x & 1}, k))
            assert cnt == dedekind[k]
    log("up-set generator reproduces Dedekind numbers M(0..5)")
    for k in (2, 3, 4, 5):
        ups = [A for A in all_upsets(k) if len(A) not in (0, 1 << k)]
        for trial in range(300 if k < 5 else 120):
            W = [rng.randint(-50, 50) for _ in range(1 << k)]
            if trial % 3 == 0:  # small spread -> many ties / nonnegative min
                W = [rng.randint(0, 3) for _ in range(1 << k)]
            W[0] -= sum(W)
            mv, A = min_nontrivial_upset(W, k)
            bf = min(sum(W[x] for x in B) for B in ups)
            assert mv == bf, (k, W, mv, bf)
            assert sum(W[x] for x in A) == mv
    log("min over nontrivial up-sets: max-flow == brute force on 1020 random instances (k=2..5)")
    # domination semantics on k=7 product measures
    q = Fraction(3, 7)
    Wn, Wd = bern_weights(q, 7)
    for trial in range(40):
        rs = [q + Fraction(rng.randint(0, 50), 1000) for _ in range(7)]
        bad = None
        if trial % 2:
            bad = rng.randrange(7)
            rs[bad] = q - Fraction(rng.randint(1, 30), 1000)
        nu = []
        for sig in range(128):
            pr = Fraction(1)
            for j in range(7):
                pr *= rs[j] if sig >> j & 1 else 1 - rs[j]
            nu.append(pr)
        assert sum(nu) == 1
        mu = [Fraction(w, Wd) for w in Wn]
        den = 1
        for x in nu + mu:
            den = den * x.denominator // __import__("math").gcd(den, x.denominator)
        W = [int((nu[s] - mu[s]) * den) for s in range(128)]
        mv, A = min_nontrivial_upset(W, 7)
        mvF = Fraction(mv, den)
        if bad is None:
            assert mvF >= 0
        else:
            assert mvF == rs[bad] - q, (mvF, rs[bad] - q)
            assert A == {s for s in range(128) if s >> bad & 1}
    log("k=7 product-measure domination tests OK (exact minima as predicted)")


if __name__ == "__main__":
    test_geometry()
    Gt = enumerate_G(3)
    pickle.dump(Gt, open("Gtable.pkl", "wb"))
    test_designated(Gt)
    test_upsets()
    test_r_nu(Gt)
    log("ALL SELF-TESTS PASSED")
