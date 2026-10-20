"""Higher-dimensional site bounds from Gomes-Pereira-Sanchis (J. Appl. Probab., arXiv:2106.10388v2), seeded with
p_c^site(Z^3) <= 4427/10000.

Inputs (theirs): Lemma 4.2: if k | d and (1 - pt) = (1 - p)^(d/k), then theta_d(p) >= theta_k(pt); so
p_c(d) <= 1 - (1 - a_k)^(k/d) whenever theta_k(a_k) > 0 (or a_k > p_c(k)). Proposition 3.4: if
p > p_c(d) [p + (1-p)^(2d/(2d-1))] then theta_{d+1}(p) > 0; with an upper bound a >= p_c(d), every p > Sol(a, b)
works (b = 2d/(2d-1)), since f(p) = p - a[p + (1-p)^b] is increasing with f(Sol) = 0.

Every bound is a 4-decimal rational, rounded UP, and certified with mpmath interval arithmetic:
- Lemma 4.2 route: x = ceil4(1 - (1-a)^(k/d)), checked (1 - x) <= (1 - a)^(k/d) in intervals (so p = x gives
  pt >= a). Then p_c(d) <= x and theta_d(x) > 0 whenever theta_k(a) > 0.
- Prop 3.4 route: x = ceil4(Sol(a, b)), checked f(x) > 0 in intervals. Then p_c(d+1) <= x (strictly below x
  holds too, so theta_{d+1}(x) > 0).
usage: python corollary.py
"""
from fractions import Fraction as F
from math import ceil

from mpmath import iv, mp

iv.prec = 200
mp.prec = 200


def ceil4(x):
    return F(ceil(F(str(x)) * 10000 - F(1, 10 ** 30)), 10000)


def lemma42(a, k, d):
    x = ceil4(1 - (1 - mp.mpf(a.numerator) / a.denominator) ** (mp.mpf(k) / d))
    lhs = 1 - iv.mpf([x.numerator, x.numerator]) / x.denominator
    rhs = (1 - iv.mpf([a.numerator, a.numerator]) / a.denominator) ** (iv.mpf(k) / d)
    assert lhs.b <= rhs.a, "lemma 4.2 rounding check failed"
    return x


def sol(a, b):
    af = mp.mpf(a.numerator) / a.denominator
    lo, hi = mp.mpf(0), mp.mpf(1)
    for _ in range(200):
        mid = (lo + hi) / 2
        if mid - af * (mid + (1 - mid) ** b) > 0:
            hi = mid
        else:
            lo = mid
    x = ceil4(hi)
    X = iv.mpf([x.numerator, x.numerator]) / x.denominator
    A = iv.mpf([a.numerator, a.numerator]) / a.denominator
    B = iv.mpf([b.numerator, b.numerator]) / b.denominator
    fx = X - A * (X + (1 - X) ** B)
    assert fx.a > 0, "Sol rounding check failed"
    return x


def table(p2, p3, dmax=9):
    p = {2: (p2, "input"), 3: (p3, "input")}
    for d in range(4, dmax + 1):
        cands = [(sol(p[d - 1][0], F(2 * (d - 1), 2 * (d - 1) - 1)), "Prop 3.4 from d=%d" % (d - 1))]
        for k in range(2, d):
            if d % k == 0:
                cands.append((lemma42(p[k][0], k, d), "Lemma 4.2, k=%d" % k))
        p[d] = min(cands)
    return p


GPS = {3: "0.5000", 4: "0.4344", 5: "0.4156", 6: "0.2929", 7: "0.2866", 8: "0.2479", 9: "0.2063"}

if __name__ == "__main__":
    print("Reproduction of Gomes-Pereira-Sanchis Table 1 (site) from their inputs p_c(2) <= 0.68, p_c(3) <= 1/2:")
    t0 = table(F(68, 100), F(1, 2))
    for d in range(3, 10):
        print("  d=%d  %s  (%s)   table: %s" % (d, float(t0[d][0]), t0[d][1], GPS[d]))
    print("\nWith p_c(3) <= 4427/10000 (this work), p_c(2) <= 0.68 as in their paper:")
    t1 = table(F(68, 100), F(4427, 10000))
    for d in range(3, 10):
        mark = "  improved" if t1[d][0] < F(GPS[d]) else ""
        print("  d=%d  %.4f  (%s)   table: %s%s" % (d, float(t1[d][0]), t1[d][1], GPS[d], mark))
    print("\nAlso with p_c(2) <= 0.666894 (Wierman and Oberly 2024) in place of 0.68:")
    t2 = table(F(666894, 1000000), F(4427, 10000))
    for d in range(3, 10):
        print("  d=%d  %.4f  (%s)" % (d, float(t2[d][0]), t2[d][1]))
