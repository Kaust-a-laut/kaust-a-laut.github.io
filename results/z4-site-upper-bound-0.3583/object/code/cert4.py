"""Certified check of the d = 4 chain: blocks {0,1}^4 of Z^4, exact integer arithmetic (Python standard library).

Block sites are i = 8*x1 + 4*x2 + 2*x3 + x4, so two sites are adjacent iff their indices differ in one bit.
G(eta) = the largest open cluster of eta; ties go to the cluster containing the smallest index; G = {} if it has
fewer than 3 sites. (Every site lies on exactly four faces, so Part II's face-weight tie-break is vacuous here.)
Faces are listed as (a,0),(a,1) for a = 1..4, so the opposite of face f is f ^ 1; the j-th site of a face in index
order is the mirror of the j-th site of the opposite face.

For each link (p, q) the program checks, as in Proposition (finite) of the paper with 7 children:
  (root) r_F(F) >= q for all 8 faces;
  (laws) for each face F and nonempty U on F, nu_{F,U} dominates Bern(q)^7. The certificate is a flow on the
         cover graph of {0,1}^7 that carries all of Bern(q)^7 upward into nu_{F,U}; the flow is validated
         directly (capacities, conservation, saturation), so the max-flow routine that finds it need not be trusted.
It also reports the smallest margin nu(A) - Bern(q)^7(A) over nontrivial up-sets A, as the value of a second
maximum flow, matched by an up-set A attaining it (so the margin printed is exact).

usage: python cert4.py counts OUT.txt              write the state counts (deterministic text)
       python cert4.py chain COUNTS.txt q1 p1 p2 ...  certify the links (p1,q1), (p2,p1), ...; fractions as a/b
"""
import sys
import time
from fractions import Fraction

NS, K, MINSIZE = 16, 7, 3
FACES = [[i for i in range(NS) if (i >> (3 - a)) & 1 == s] for a in range(4) for s in (0, 1)]


def check(cond, what):
    """Every certificate condition goes through here, never through assert (python -O would skip an assert)."""
    if not cond:
        print("CERTIFICATE: FAIL (%s)" % what, flush=True)
        sys.exit(1)


def designated(eta):
    seen, best, bestsize = 0, 0, 0
    for i in range(NS):  # a cluster is first reached at its smallest index, so ">" keeps the smallest on ties
        if (eta >> i) & 1 and not (seen >> i) & 1:
            comp, stack = 1 << i, [i]
            seen |= 1 << i
            while stack:
                u = stack.pop()
                for b in range(4):
                    v = u ^ (1 << b)
                    if (eta >> v) & 1 and not (seen >> v) & 1:
                        seen |= 1 << v
                        comp |= 1 << v
                        stack.append(v)
            size = bin(comp).count("1")
            if size > bestsize:
                best, bestsize = comp, size
    return best if bestsize >= MINSIZE else 0


def footprint(G, face):
    return sum(((G >> s) & 1) << j for j, s in enumerate(face))


def write_counts(path):
    counts = {}
    for eta in range(1 << NS):
        G = designated(eta)
        counts.setdefault(G, [0] * (NS + 1))[bin(eta).count("1")] += 1
    with open(path, "w", newline="\n") as fh:
        for G in sorted(counts):
            fh.write("%d %s\n" % (G, " ".join(map(str, counts[G]))))
    print("%d designated sets, %d states" % (len(counts), sum(map(sum, counts.values()))))


def read_counts(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            v = list(map(int, line.split()))
            rows.append((v[0], v[1:]))
    check(sum(sum(c) for _, c in rows) == 1 << NS, "the census does not cover 2^16 states")
    return rows


def zeta(Z, nbits):
    """In place: Z[s] <- sum of Z[t] over t subset of s (entries are ints or lists of ints)."""
    for b in range(nbits):
        for s in range(1 << nbits):
            if (s >> b) & 1:
                t = s ^ (1 << b)
                Z[s] = Z[s] + Z[t] if isinstance(Z[s], int) else [u + v for u, v in zip(Z[s], Z[t])]
    return Z


class Net:
    def __init__(self, n):
        self.n, self.head, self.to, self.cap = n, [[] for _ in range(n)], [], []

    def add(self, u, v, c):
        self.head[u].append(len(self.to)); self.to.append(v); self.cap.append(c)
        self.head[v].append(len(self.to)); self.to.append(u); self.cap.append(0)

    def maxflow(self, s, t):
        self.orig = list(self.cap)
        flow = 0
        while True:
            level = [-1] * self.n
            level[s] = 0
            queue = [s]
            for u in queue:
                for e in self.head[u]:
                    if self.cap[e] > 0 and level[self.to[e]] < 0:
                        level[self.to[e]] = level[u] + 1
                        queue.append(self.to[e])
            if level[t] < 0:
                return flow
            it = [0] * self.n

            def push(u, f):
                if u == t:
                    return f
                while it[u] < len(self.head[u]):
                    e = self.head[u][it[u]]
                    v = self.to[e]
                    if self.cap[e] > 0 and level[v] == level[u] + 1:
                        g = push(v, min(f, self.cap[e]))
                        if g:
                            self.cap[e] -= g
                            self.cap[e ^ 1] += g
                            return g
                    it[u] += 1
                return 0

            while True:
                g = push(s, self.orig_total + 1)
                if not g:
                    break
                flow += g

    def validate(self, s, t):
        """Check the current flow directly: 0 <= flow <= capacity on every arc, conservation at every inner node.
        Returns the net flow out of s."""
        net = [0] * self.n
        for e in range(0, len(self.to), 2):
            fl = self.orig[e] - self.cap[e]
            check(0 <= fl <= self.orig[e] and self.cap[e ^ 1] == fl, "flow outside an arc's capacity")
            u, v = self.to[e ^ 1], self.to[e]
            net[u] -= fl
            net[v] += fl
        check(all(net[x] == 0 for x in range(self.n) if x not in (s, t)), "flow not conserved")
        return -net[s]


def law_network(mu, nu, forced):
    """Nodes 0..127 = {0,1}^7, source 128, sink 129; supplies mu at x, capacities nu at y, cover arcs upward."""
    total = sum(mu)
    net = Net((1 << K) + 2)
    net.orig_total = total + sum(nu)
    S, T = 1 << K, (1 << K) + 1
    for x in range(1 << K):
        if mu[x]:
            net.add(S, x, mu[x])
        if nu[x]:
            net.add(x, T, nu[x])
        for i in range(K):
            if not (x >> i) & 1:
                net.add(x, x | (1 << i), net.orig_total)
    if forced:  # the top in the source side, the bottom in the sink side: only nontrivial up-sets are cuts
        net.add(S, (1 << K) - 1, net.orig_total)
        net.add(0, T, net.orig_total)
    return net, S, T, total


def coupling_ok(mu, nu):
    net, S, T, total = law_network(mu, nu, False)
    net.maxflow(S, T)
    value = net.validate(S, T)
    saturated = all(net.cap[e] == 0 for e in net.head[S] if e % 2 == 0)  # all of mu carried upward into nu
    return saturated and value == total


def min_margin(mu, nu):
    """min over nontrivial up-sets A of nu(A) - mu(A), with the up-set attaining it."""
    net, S, T, total = law_network(mu, nu, True)
    value = net.maxflow(S, T)
    check(net.validate(S, T) == value, "margin flow value")
    side, queue = {S}, [S]  # residual reachability gives a minimum cut
    for u in queue:
        for e in net.head[u]:
            if net.cap[e] > 0 and net.to[e] not in side:
                side.add(net.to[e])
                queue.append(net.to[e])
    A = [x for x in range(1 << K) if x in side]
    check(all(x | (1 << i) in side for x in A for i in range(K)), "margin cut is not an up-set")
    check((1 << K) - 1 in side and 0 not in side, "margin cut is a trivial up-set")
    gap = sum(nu[x] for x in A) - sum(mu[x] for x in A)
    check(gap == value - total, "margin cut and flow differ")  # flow value = cut capacity: both optimal
    return gap, A


def link(rows, p, q):
    t0 = time.time()
    a, b = p.numerator, p.denominator
    c, e = q.numerator, q.denominator
    B16 = b ** NS
    w = [a ** k * (b - a) ** (NS - k) for k in range(NS + 1)]
    Gs = [G for G, _ in rows]
    W = [sum(n * wk for n, wk in zip(cnt, w)) for _, cnt in rows]
    check(sum(W) == B16, "state weights do not sum to 1")
    fp = [[footprint(G, FACES[f]) for f in range(8)] for G in Gs]
    check(designated((1 << NS) - 1) == (1 << NS) - 1, "G of the all-open block is not B")  # used at the first step
    # R[g][T] = numerator of r_g(T) = P(T meets the mirror of the neighbour's footprint on the opposite face)
    R = []
    for g in range(8):
        P = [0] * 256
        for i in range(len(Gs)):
            P[fp[i][g ^ 1]] += W[i]
        Z = zeta(P, 8)
        R.append([B16 - Z[255 ^ T] for T in range(256)])
    root = min(R[g][255] for g in range(8))
    check(root * e >= c * B16, "root condition r_F(F) >= q fails")
    root_gap = Fraction(root, B16) - q
    mu_unit = [c ** bin(x).count("1") * (e - c) ** (K - bin(x).count("1")) for x in range(1 << K)]
    worst = None
    nlaws = 0
    for f in range(8):
        ch = [g for g in range(8) if g != f]
        A = [[0] * (1 << K) for _ in range(256)]
        for i in range(len(Gs)):
            vals = [W[i]]
            for g in ch:  # bit j of the pattern x = value of the j-th child
                r = R[g][fp[i][g]]
                s = B16 - r
                vals = [v * s for v in vals] + [v * r for v in vals]
            acc = A[fp[i][f]]
            for x in range(1 << K):
                acc[x] += vals[x]
        Wt = [0] * 256
        for i in range(len(Gs)):
            Wt[fp[i][f]] += W[i]
        ZA, ZW = zeta(A, 8), zeta(Wt, 8)
        for U in range(1, 256):
            nuN = [u - v for u, v in zip(ZA[255], ZA[255 ^ U])]  # states whose entry footprint meets U
            mass = ZW[255] - ZW[255 ^ U]
            check(mass > 0 and min(nuN) >= 0 and sum(nuN) == mass * B16 ** K, "law (F=%d, U=%d) is not a law" % (f, U))
            D = sum(nuN)
            mu = [m * D for m in mu_unit]
            nu = [v * e ** K for v in nuN]
            check(sum(mu) == sum(nu), "unequal masses")
            check(coupling_ok(mu, nu), "law (F=%d, U=%d) does not dominate Bern(q)^7" % (f, U))
            gap, upset = min_margin(mu, nu)
            m = Fraction(gap, sum(mu))
            nlaws += 1
            if worst is None or m < worst[0]:
                worst = (m, f, U, len(upset))
    check(nlaws == 8 * 255 and worst[0] >= 0 and root_gap >= 0, "a margin or the root gap is negative")
    print("CERTIFICATE: PASS (p = %s, q = %s): root gap %.4f, %d laws, smallest margin %.6e at face %d, U = %s "
          "(up-set of %d patterns) [%.0fs]" % (p, q, float(root_gap), nlaws, float(worst[0]), worst[1],
                                                format(worst[2], "08b")[::-1], worst[3], time.time() - t0), flush=True)
    return worst


def main():
    if sys.flags.optimize:
        sys.exit("refusing to run under python -O: run plain python")
    if sys.argv[1] == "counts":
        write_counts(sys.argv[2])
        return
    rows = read_counts(sys.argv[2])
    qs = [Fraction(x) for x in sys.argv[3:]]
    for q, p in zip(qs, qs[1:]):
        check(p < q, "the chain must decrease")
        link(rows, p, q)


if __name__ == "__main__":
    main()
