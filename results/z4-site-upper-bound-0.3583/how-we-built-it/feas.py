"""Feasibility estimates for the block-coupling upper bound on p_c^site(Z^d). Floating point, NOT certificates.

Block B = prod range(L_a), 2d faces, 2d - 1 children per explored block (Part II of papers/z3-site-bounds).
f(p) = min( min_F r_F(F),  min over entry faces F and test sets U of q*(nu_{F,U}) ), where
q*(nu) = max q with nu >=_st Bernoulli(q)^(2d-1), decided by a max-flow on the Boolean lattice (Strassen).
Test sets U: every single site of the entry face, plus the whole face (in d = 3 the binding case was a single site),
so the f printed here is an upper estimate of the true f (the full check uses every nonempty U).
States: all 2^n with exact weights when n <= 20, else Monte Carlo samples (same uniforms for every p).

usage: python feas.py L1,L2,..,Ld rule p1,p2,... [samples] [seed] [facemode] [test|all]
  A 7th argument switches to the aggregated path: "all" tests every nonempty U, "test" the sets above; states are
  enumerated exactly (in chunks) up to 24 sites. It reproduces the default path's values (2x2x2x2, 3x3x3 MC).
  rule: largest:m      largest open cluster with >= m sites (ties: face weight, then smallest site index)
        core:m         largest open cluster meeting the central core (sites with every coordinate in 1..L_a-2),
                       if it has >= m sites (same ties)
        corepri:m      the cluster containing the first open core site in index order, if it has >= m sites
                       (needs no cluster sizes beyond the minimum)
        coremost:m     the cluster with the most core sites (ties: smallest site index), if it has >= m sites
        coreseed:c     the cluster of the first open core site whose cluster contains >= c core sites
        filters (size-free for a slice-by-slice verifier): write coremost:fh2, corepri:fw6, coreseed:2:fh2, ...
                       sN = >= N sites; fwN = >= N face sites (with multiplicity); fhN = touches >= N faces
  facemode also: pinwheel (centre 2x2 + one arm per edge on every face), belt / beltnc (full / noncorner faces
        perpendicular to the last axis, middle-slice belt on the other faces)
  facemode (optional 6th argument): full (default); center = blocks join only through the central sites of the
        shared face (every other coordinate in 1..L_a-2); noncorner = every face site except the face's corners;
        mixed / mixednc = full faces perpendicular to the last axis (the sweep direction), center / noncorner on
        the other faces
"""
import itertools
import sys
import time

import numpy as np

INF = 1 << 62


def build(L, facemode="full"):
    d = len(L)
    coords = list(itertools.product(*[range(l) for l in L]))  # index order = lexicographic, first coord major
    idx = {c: i for i, c in enumerate(coords)}
    edges = []
    for c in coords:
        for a in range(d):
            if c[a] + 1 < L[a]:
                c2 = list(c)
                c2[a] += 1
                edges.append((idx[c], idx[tuple(c2)]))
    def keep(c, a):
        others = [b for b in range(d) if b != a]
        center = all(1 <= c[b] <= L[b] - 2 for b in others)
        corner = all(c[b] in (0, L[b] - 1) for b in others)
        mode = facemode
        if facemode in ("mixed", "mixednc"):
            mode = "full" if a == d - 1 else ("center" if facemode == "mixed" else "noncorner")
        if facemode in ("belt", "beltnc"):
            if a == d - 1:
                mode = "full" if facemode == "belt" else "noncorner"
            else:
                return 1 <= c[d - 1] <= L[d - 1] - 2  # side faces: the middle slices of the sweep axis
        if facemode == "pinwheel":  # 4x4 faces: centre 2x2 plus one arm site per edge, in face coords (i, j)
            ij = tuple(c[b] for b in others)
            return ij in {(1, 1), (1, 2), (2, 1), (2, 2), (0, 1), (1, 3), (3, 2), (2, 0)}
        return mode == "full" or (mode == "center" and center) or (mode == "noncorner" and not corner)

    faces = []  # order (a,0),(a,L-1) for a = 0..d-1, so the opposite face of f is f ^ 1
    for a in range(d):
        for s in (0, L[a] - 1):
            faces.append([i for i, c in enumerate(coords) if c[a] == s and keep(c, a)])
    mult = np.zeros(len(coords))
    for fs in faces:
        mult[fs] += 1
    core = np.array([all(1 <= c[a] <= L[a] - 2 for a in range(d)) for c in coords])
    return d, len(coords), edges, faces, mult, core


def label(S, edges):
    """S: N x n bool. Returns cluster labels = smallest site index in the cluster (n for closed sites)."""
    n = S.shape[1]
    lab = np.where(S, np.arange(n)[None, :], n)
    changed = True
    while changed:
        changed = False
        for i, j in edges:
            both = S[:, i] & S[:, j]
            m = np.minimum(lab[:, i], lab[:, j])
            upd = both & ((lab[:, i] != m) | (lab[:, j] != m))
            if upd.any():
                changed = True
                lab[upd, i] = m[upd]
                lab[upd, j] = m[upd]
    return lab


def apply_filter(G, filt, faces, mult):
    if not filt:
        return G
    kind, val = (filt[:2], int(filt[2:])) if filt[:2] in ("fw", "fh") else (filt[:1], int(filt[1:]))
    if kind == "s":
        ok = G.sum(axis=1) >= val
    elif kind == "fw":
        ok = (G * mult[None, :]).sum(axis=1) >= val
    else:
        ok = np.stack([G[:, fs].any(axis=1) for fs in faces], axis=1).sum(axis=1) >= val
    return G & ok[:, None]


def designated(S, edges, mult, core, rule):
    N, n = S.shape
    lab = label(S, edges)
    rows = np.repeat(np.arange(N), n) * (n + 1) + lab.ravel()
    size = np.bincount(rows, minlength=N * (n + 1)).reshape(N, n + 1)[:, :n]
    facew = np.bincount(rows, weights=np.tile(mult, N), minlength=N * (n + 1)).reshape(N, n + 1)[:, :n]
    score = size * 1000.0 + facew
    parts = rule.split(":")
    kind = parts[0]
    m = int(parts[1]) if parts[1].isdigit() else 1  # a non-numeric spec means: no size filter here
    if kind in ("coremost", "coreseed"):
        hitc = np.bincount(rows, weights=np.tile(core.astype(float), N), minlength=N * (n + 1)).reshape(N, n + 1)
        if kind == "coremost":
            sc = np.where(size > 0, hitc[:, :n], 0.0)
            pick = np.argmax(sc, axis=1)
            ok = (sc[np.arange(N), pick] > 0) & (size[np.arange(N), pick] >= (m if kind == "coremost" else 1))
            return (lab == pick[:, None]) & ok[:, None]
        coreidx = np.where(core)[0]
        pick = np.full(N, n)
        for j in coreidx[::-1]:  # the first valid core site in index order wins
            lj = lab[:, j]
            valid = S[:, j] & (hitc[np.arange(N), lj] >= m)
            pick = np.where(valid, lj, pick)
        ok = pick < n
        return (lab == pick[:, None]) & ok[:, None]
    if kind == "corepri":
        coreidx = np.where(core)[0]
        openc = S[:, coreidx]
        first = coreidx[np.argmax(openc, axis=1)]
        pick = lab[np.arange(N), first]
        ok = openc.any(axis=1) & (size[np.arange(N), np.minimum(pick, n - 1)] >= m)
        return (lab == pick[:, None]) & ok[:, None]
    if kind == "core":
        hit = np.bincount(rows, weights=np.tile(core.astype(float), N), minlength=N * (n + 1)).reshape(N, n + 1)[:, :n]
        score = np.where(hit > 0, score, -1.0)
    score[:, :][size == 0] = -1.0
    pick = np.argmax(score, axis=1)  # first maximum = smallest label = cluster with the smallest site index
    ok = (score[np.arange(N), pick] >= 0) & (size[np.arange(N), pick] >= m)
    return (lab == pick[:, None]) & ok[:, None]


def footprints(G, faces):
    return np.stack([(G[:, fs].astype(np.int64) << np.arange(len(fs))).sum(axis=1) for fs in faces], axis=1)


def zeta(P):
    Z = P.astype(float).copy()
    m = Z.size.bit_length() - 1
    for b in range(m):
        Zr = Z.reshape(-1, 2, 1 << b)
        Zr[:, 1, :] += Zr[:, 0, :]
    return Z


def maxflow_ok(Bc, Vc, k):
    """Bc, Vc: integer masses on {0,1}^k. True iff all of Bc can flow upward into Vc (domination)."""
    N = 1 << k
    S, T = N, N + 1
    to, cap, head = [], [], [[] for _ in range(N + 2)]

    def add(u, v, c):
        head[u].append(len(to)); to.append(v); cap.append(c)
        head[v].append(len(to)); to.append(u); cap.append(0)

    for x in range(N):
        if Bc[x] > 0:
            add(S, x, Bc[x])
        if Vc[x] > 0:
            add(x, T, Vc[x])
        for i in range(k):
            if not (x >> i) & 1:
                add(x, x | (1 << i), INF)
    need, flow = sum(Bc), 0
    while True:
        level = [-1] * (N + 2)
        level[S] = 0
        q = [S]
        for u in q:
            for e in head[u]:
                if cap[e] > 0 and level[to[e]] < 0:
                    level[to[e]] = level[u] + 1
                    q.append(to[e])
        if level[T] < 0:
            break
        it = [0] * (N + 2)

        def dfs(u, f):
            if u == T:
                return f
            while it[u] < len(head[u]):
                e = head[u][it[u]]
                v = to[e]
                if cap[e] > 0 and level[v] == level[u] + 1:
                    g = dfs(v, min(f, cap[e]))
                    if g:
                        cap[e] -= g
                        cap[e ^ 1] += g
                        return g
                it[u] += 1
            return 0

        while True:
            g = dfs(S, INF)
            if not g:
                break
            flow += g
    return flow >= need - N  # both sides were floored: allow one unit of rounding per state (scale is 2^40)


def qstar(law, k, scale=1 << 40, iters=22):
    pc = np.array([bin(x).count("1") for x in range(1 << k)])
    Vc = [int(v * scale) for v in law]
    marg = [law[[x for x in range(1 << k) if (x >> i) & 1]].sum() for i in range(k)]
    lo, hi = 0.0, min(marg)
    for _ in range(iters):
        q = (lo + hi) / 2
        B = q ** pc * (1 - q) ** (k - pc)
        Bc = [int(b * scale) for b in B]
        if maxflow_ok(Bc, Vc, k):
            lo = q
        else:
            hi = q
    return lo


def f_estimate(Tst, w, faces, d):
    k = 2 * d - 1
    nf = len(faces)
    combos, inv = np.unique(Tst, axis=0, return_inverse=True)
    W = np.bincount(inv.ravel(), weights=w)
    r = []
    for g in range(nf):
        o = g ^ 1
        P = np.bincount(combos[:, o], weights=W, minlength=1 << len(faces[o]))
        Z = zeta(P / P.sum())
        full = (1 << len(faces[g])) - 1
        r.append(1.0 - Z[full ^ np.arange(1 << len(faces[g]))])  # r_g(T) = 1 - P(T_opp subset of complement T)
    root = min(r[g][(1 << len(faces[g])) - 1] for g in range(nf))
    R = np.stack([r[g][combos[:, g]] for g in range(nf)], axis=1)
    pats = np.array([[(x >> i) & 1 for i in range(k)] for x in range(1 << k)], dtype=bool)
    best = (root, "root")
    for f in range(nf):
        ch = [g for g in range(nf) if g != f]
        Pp = np.ones((len(combos), 1 << k))
        for i, g in enumerate(ch):
            Pp *= np.where(pats[None, :, i], R[:, g][:, None], 1 - R[:, g][:, None])
        tests = [(1 << j, "site%d" % j) for j in range(len(faces[f]))] + [((1 << len(faces[f])) - 1, "face")]
        for U, name in tests:
            msk = (combos[:, f] & U) != 0
            wt = W * msk
            if wt.sum() == 0:
                continue
            law = (wt @ Pp) / wt.sum()
            q = qstar(law, k)
            if q < best[0]:
                best = (q, "entry face %d, U=%s" % (f, name))
    return best


def exact_table(n, edges, mult, core, rule, filt, faces, chunk=1 << 20):
    """All 2^n states, in chunks. Returns the distinct designated sets G as site bitmasks, and for each state class
    (G, k open sites) its index into them, k and the number of states."""
    bits = np.arange(n, dtype=np.int64)
    keys, cnts = [], []
    for s in range(0, 1 << n, chunk):
        x = np.arange(s, min(s + chunk, 1 << n), dtype=np.int64)
        S = ((x[:, None] >> bits[None, :]) & 1).astype(bool)
        G = apply_filter(designated(S, edges, mult, core, rule), filt, faces, mult)
        u, c = np.unique((G.astype(np.int64) << bits).sum(axis=1) * (n + 1) + S.sum(axis=1), return_counts=True)
        keys.append(u)
        cnts.append(c)
    uk, inv = np.unique(np.concatenate(keys), return_inverse=True)
    c = np.bincount(inv.ravel(), weights=np.concatenate(cnts))
    Gm, ginv = np.unique(uk // (n + 1), return_inverse=True)
    return Gm, ginv.ravel(), uk % (n + 1), c


def dominates(law, q, k, scale=1 << 40):
    pc = np.array([bin(x).count("1") for x in range(1 << k)])
    B = q ** pc * (1 - q) ** (k - pc)
    return maxflow_ok([int(b * scale) for b in B], [int(v * scale) for v in law], k)


def zeta0(A):
    """Subset sums along axis 0 (index = bitmask), any trailing shape."""
    Z = A.astype(float).copy()
    m = Z.shape[0].bit_length() - 1
    for b in range(m):
        Zr = Z.reshape((-1, 2, 1 << b) + Z.shape[1:])
        Zr[:, 1] += Zr[:, 0]
    return Z


def f_estimate_all(Tst, W, faces, d, umode, chunk=1 << 18):
    """As f_estimate, but every nonempty U when umode == "all". For each entry face the children's law is
    aggregated by the entry footprint t, so that the law given {t meets U} is a subset-sum lookup."""
    k = 2 * d - 1
    nf = len(faces)
    C = len(W)
    r = []
    for g in range(nf):
        P = np.bincount(Tst[:, g ^ 1], weights=W, minlength=1 << len(faces[g ^ 1]))
        Z = zeta(P / P.sum())
        full = (1 << len(faces[g])) - 1
        r.append(1.0 - Z[full ^ np.arange(1 << len(faces[g]))])
    root = min(r[g][(1 << len(faces[g])) - 1] for g in range(nf))
    best = (root, "root")
    nU = 0
    for f in range(nf):
        ch = [g for g in range(nf) if g != f]
        nb = len(faces[f])
        T = Tst[:, f]
        A = np.zeros((1 << nb, 1 << k))
        for s in range(0, C, chunk):
            e = min(s + chunk, C)
            Pp = W[s:e, None].copy()
            for g in ch:
                rg = r[g][Tst[s:e, g]][:, None]
                Pp = np.concatenate([Pp * (1 - rg), Pp * rg], axis=1)  # bit i of the index = child i
            for x in range(1 << k):
                A[:, x] += np.bincount(T[s:e], weights=Pp[:, x], minlength=1 << nb)
        ZA = zeta0(A)
        Atot = ZA[-1]
        full = (1 << nb) - 1
        tests = [1 << j for j in range(nb)] + [full]
        if umode == "all":
            tests += [U for U in range(1, full) if U & (U - 1)]
        for U in tests:
            vec = Atot - ZA[full ^ U]  # mass of the states whose entry footprint meets U, per child pattern
            mass = vec.sum()
            if mass <= 1e-15 * Atot.sum():
                continue
            law = np.maximum(vec, 0) / mass
            nU += 1
            if dominates(law, best[0], k):
                continue
            q = qstar(law, k)
            if q < best[0]:
                best = (q, "entry face %d, U=%s" % (f, format(U, "0%db" % nb)[::-1]))
    return best[0], best[1] + " [%d sets U]" % nU


def main_all(L, rule, ps, samples, seed, facemode, umode):
    parts = rule.split(":")
    filt = next((x for x in parts[1:] if not x.isdigit()), "")
    d, n, edges, faces, mult, core = build(L, facemode)
    exact = n <= 24
    print("L=%s rule=%s faces=%s (%s sites) n=%d %s, U=%s" % (L, rule, facemode, "/".join(str(len(f)) for f in faces),
          n, "exact" if exact else "Monte Carlo, %d samples" % samples, umode))
    t0 = time.time()
    if exact:
        Gm, ginv, kk, cnt = exact_table(n, edges, mult, core, rule, filt, faces)
        Tst = np.stack([sum(((Gm >> j) & 1) << i for i, j in enumerate(fs)) for fs in faces], axis=1)
        print("  %d distinct designated sets (%.0fs)" % (len(Gm), time.time() - t0), flush=True)
    else:
        Us = np.random.default_rng(seed).random((samples, n))
    for p in ps:
        t0 = time.time()
        if exact:
            W = np.bincount(ginv, weights=cnt * p ** kk * (1 - p) ** (n - kk), minlength=len(Gm))
        else:
            G = apply_filter(designated(Us < p, edges, mult, core, rule), filt, faces, mult)
            Tst, inv = np.unique(footprints(G, faces), axis=0, return_inverse=True)
            W = np.bincount(inv.ravel(), minlength=len(Tst)) / samples
        q, where = f_estimate_all(Tst, W, faces, d, umode)
        print("p=%.5f  f(p)~%.5f  f-p=%+.5f  binding: %s  (%.0fs)" % (p, q, q - p, where, time.time() - t0), flush=True)


def main():
    L = [int(x) for x in sys.argv[1].split(",")]
    rule = sys.argv[2]
    ps = [float(x) for x in sys.argv[3].split(",")]
    samples = int(sys.argv[4]) if len(sys.argv) > 4 else 400000
    seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    facemode = sys.argv[6] if len(sys.argv) > 6 else "full"
    if len(sys.argv) > 7:  # "test" or "all": the aggregated path (exact up to 24 sites)
        return main_all(L, rule, ps, samples, seed, facemode, sys.argv[7])
    parts = rule.split(":")
    filt = next((x for x in parts[1:] if not x.isdigit()), "")
    d, n, edges, faces, mult, core = build(L, facemode)
    exact = n <= 20
    print("L=%s rule=%s faces=%s (%d sites each) n=%d %s" % (L, rule, facemode, len(faces[0]), n,
          "exact" if exact else "Monte Carlo, %d samples" % samples))
    if exact:
        allS = ((np.arange(1 << n)[:, None] >> np.arange(n)[None, :]) & 1).astype(bool)
        G = apply_filter(designated(allS, edges, mult, core, rule), filt, faces, mult)
        Tst = footprints(G, faces)
        kopen = allS.sum(axis=1)
    else:
        U = np.random.default_rng(seed).random((samples, n))
    for p in ps:
        t0 = time.time()
        if exact:
            w = p ** kopen * (1 - p) ** (n - kopen)
        else:
            S = U < p
            G = apply_filter(designated(S, edges, mult, core, rule), filt, faces, mult)
            Tst = footprints(G, faces)
            w = np.full(samples, 1.0 / samples)
        q, where = f_estimate(Tst, w, faces, d)
        print("p=%.5f  f(p)~%.5f  f-p=%+.5f  binding: %s  (%.0fs)" % (p, q, q - p, where, time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
