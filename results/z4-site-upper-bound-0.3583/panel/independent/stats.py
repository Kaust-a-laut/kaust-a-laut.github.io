"""Side statistics on the designated-cluster rule (tie frequency, phi redundancy)."""
from core import *
ties = 0; ties_weighted = {}
phi_redundant = True
for eta in range(1 << NSITES):
    cl = clusters(eta)
    for C in cl:
        if phi(C) != 4 * popcount(C):
            phi_redundant = False
    if len(cl) > 1:
        s = sorted((popcount(C) for C in cl), reverse=True)
        if s[0] >= 3 and s[1] == s[0]:
            ties += 1
            n = popcount(eta)
            ties_weighted[n] = ties_weighted.get(n, 0) + 1
print("phi(C) == 4|C| for every cluster of every state:", phi_redundant)
print("states whose largest cluster (>=3 sites) is tied in size:", ties)
p = 0.36
print("P_p(size tie) at p=0.36: %.3e" % sum(c * p**n * (1-p)**(16-n) for n, c in ties_weighted.items()))
print("tie counts by |eta|:", dict(sorted(ties_weighted.items())))
