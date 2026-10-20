"""Design the d = 4 chain (floating point; the certified check is separate). From q, take the smallest p on the
1e-4 grid with f(p) >= q + gap, where f is the all-U estimate for 2x2x2x2 blocks and largest:3.

usage: python design_chain.py q1 gap
"""
import sys

import numpy as np

import feas

L, RULE = [2, 2, 2, 2], "largest:3"
d, n, edges, faces, mult, core = feas.build(L)
Gm, ginv, kk, cnt = feas.exact_table(n, edges, mult, core, RULE, "", faces)
Tst = np.stack([sum(((Gm >> j) & 1) << i for i, j in enumerate(fs)) for fs in faces], axis=1)
cache = {}


def f(m):  # p = m / 10000
    if m not in cache:
        p = m / 10000
        W = np.bincount(ginv, weights=cnt * p ** kk * (1 - p) ** (n - kk), minlength=len(Gm))
        cache[m] = feas.f_estimate_all(Tst, W, faces, d, "all")[0]
    return cache[m]


def main():
    q, gap = float(sys.argv[1]), float(sys.argv[2])
    chain = []
    while True:
        lo, hi = 3000, int(round(q * 10000))  # f(lo) < q, look for the smallest m with f(m) >= q + gap
        if f(hi) < q + gap:
            break
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if f(mid) >= q + gap:
                hi = mid
            else:
                lo = mid
        if chain and hi >= chain[-1][0]:
            break
        chain.append((hi, q, f(hi)))
        print("p=%d/10000  q=%.5f  f(p)~%.6f  gap %.6f" % (hi, q, f(hi), f(hi) - q), flush=True)
        q = hi / 10000


if __name__ == "__main__":
    main()
