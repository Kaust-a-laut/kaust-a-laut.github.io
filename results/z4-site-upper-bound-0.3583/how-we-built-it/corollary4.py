"""Site bounds for d = 5..9 from Gomes-Pereira-Sanchis (Lemma 4.2, Proposition 3.4), seeded with
p_c(Z^3) <= 4427/10000 and a direct bound p_c(Z^4) <= P4 (the d = 4 chain). Same rounding and interval checks as
corollary.py (from the Z^3 upper-bound work; needs mpmath), whose functions are reused.

usage: python corollary4.py P4        (P4 as a fraction, e.g. 3583/10000)
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from corollary import GPS, lemma42, sol  # noqa: E402

JL26 = {5: "0.2839", 6: "0.2088", 7: "0.1676", 8: "0.1409", 9: "0.1221"}


def table(p2, p3, p4, dmax=9):
    p = {2: (p2, "input"), 3: (p3, "input"), 4: (p4, "input (d = 4 chain)")}
    for d in range(5, dmax + 1):
        cands = [(sol(p[d - 1][0], F(2 * (d - 1), 2 * (d - 1) - 1)), "Prop 3.4 from d=%d" % (d - 1))]
        for k in range(2, d):
            if d % k == 0:
                cands.append((lemma42(p[k][0], k, d), "Lemma 4.2, k=%d" % k))
        p[d] = min(cands)
    return p


def table_unrounded(inputs, dmax=9):
    """The same recursion with no intermediate rounding (200-bit floating point; for comparison only).
    inputs: {d: value} for the dimensions given as inputs."""
    from mpmath import mp, mpf
    mp.prec = 200
    v = {d: mpf(x.numerator) / x.denominator for d, x in inputs.items()}
    for d in range(max(inputs) + 1, dmax + 1):
        a, b = v[d - 1], mpf(2 * (d - 1)) / (2 * (d - 1) - 1)
        lo, hi = mpf(0), mpf(1)
        for _ in range(200):  # Proposition 3.4: the root of p - a [p + (1 - p)^b]
            mid = (lo + hi) / 2
            lo, hi = (lo, mid) if mid - a * (mid + (1 - mid) ** b) > 0 else (mid, hi)
        cands = [hi] + [1 - (1 - v[k]) ** (mpf(k) / d) for k in range(2, d) if d % k == 0]  # Lemma 4.2
        v[d] = min(cands)
    return v


if __name__ == "__main__":
    t = table_unrounded({2: F(68, 100), 3: F(1, 2)})
    print("Without intermediate rounding, from GPS26's inputs p_c(2) <= 0.68, p_c(3) <= 1/2:")
    for d in range(4, 10):
        print("  d=%d  %.6f   GPS26 table: %s" % (d, float(t[d]), GPS[d]))
    t = table_unrounded({2: F(68, 100), 3: F(4427, 10000), 4: F(sys.argv[1])})
    print("Without intermediate rounding, from p_c(3) <= 4427/10000 and p_c(4) <= %s:" % sys.argv[1])
    for d in range(5, 10):
        print("  d=%d  %.6f" % (d, float(t[d])))
    print()
    p4 = F(sys.argv[1])
    for name, p2 in (("p_c(2) <= 0.68 (as in GPS26)", F(68, 100)), ("p_c(2) <= 0.666894 (WO24)", F(666894, 1000000))):
        print("With %s, p_c(3) <= 4427/10000, p_c(4) <= %s:" % (name, p4))
        t = table(p2, F(4427, 10000), p4)
        for d in range(4, 10):
            print("  d=%d  %.4f  (%s)   GPS26 table: %s   JL26: %s" % (d, float(t[d][0]), t[d][1], GPS[d],
                                                                       JL26.get(d, "0.4719")))
