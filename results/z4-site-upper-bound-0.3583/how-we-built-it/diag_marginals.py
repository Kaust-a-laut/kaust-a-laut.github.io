"""Print the seven child marginals of one law nu_{F,U} (floating point, a diagnostic only).

usage: python diag_marginals.py L1,..,Ld rule p entry_face site_index [samples]
  exact when the block has <= 16 sites, else Monte Carlo (seed 1)
"""
import sys

import numpy as np

import feas


def main():
    L = [int(x) for x in sys.argv[1].split(",")]
    rule, p, f, j = sys.argv[2], float(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    samples = int(sys.argv[6]) if len(sys.argv) > 6 else 1000000
    d, n, edges, faces, mult, core = feas.build(L)
    if n <= 16:
        S = ((np.arange(1 << n)[:, None] >> np.arange(n)[None, :]) & 1).astype(bool)
        k = S.sum(axis=1)
        w = p ** k * (1 - p) ** (n - k)
    else:
        S = np.random.default_rng(1).random((samples, n)) < p
        w = np.full(samples, 1.0 / samples)
    G = feas.designated(S, edges, mult, core, rule)
    Tst = feas.footprints(G, faces)
    nf = len(faces)
    r = []
    for g in range(nf):
        P = np.bincount(Tst[:, g ^ 1], weights=w, minlength=1 << len(faces[g ^ 1]))
        Z = feas.zeta(P / P.sum())
        full = (1 << len(faces[g])) - 1
        r.append(1.0 - Z[full ^ np.arange(1 << len(faces[g]))])
    msk = (Tst[:, f] >> j) & 1 == 1
    print("L=%s rule=%s p=%.4f entry face %d (axis %d), entry site %d; P(entry site in G) = %.5f"
          % (L, rule, p, f, f // 2, j, (w * msk).sum()))
    print("root r_g(full face):", " ".join("%.4f" % r[g][-1] for g in range(nf)))
    for g in range(nf):
        if g == f:
            continue
        m = (w * msk * r[g][Tst[:, g]]).sum() / (w * msk).sum()
        hit = (w * msk * (Tst[:, g] > 0)).sum() / (w * msk).sum()
        print("  child face %d (axis %d, side %s, %2d sites): marginal %.5f   P(G touches this face | entry) %.5f"
              % (g, g // 2, "low " if g % 2 == 0 else "high", len(faces[g]), m, hit))


if __name__ == "__main__":
    main()
