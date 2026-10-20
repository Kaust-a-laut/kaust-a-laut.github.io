# Independent re-check

`independent/` is a second implementation of the whole verification, written during an AI review from the
definitions alone (block, faces, mirror map, designated cluster, r_F, nu_{F,U}, the pass condition) and without
access to object/code or to the margins of the certified runs. It shares no code with object/code/cert4.py. For
each law it computes min over nontrivial up-sets A of nu(A) - Bernoulli(q)^7(A) exactly, as a minimum-weight closure
problem solved by a minimum cut, checks the flow, the cut and their equality, and decides domination from the sign of
the minimum. This is the dual form of the coupling flow in object/code, so the independence is in the code and the
computation of the laws, not in the underlying max-flow method. Run it with plain python, not python -O: its flow and
cut checks are asserts.

Results (exact decisions; decimals for reporting):
- the census: 37246 designated sets, equal to 1 + the number of connected subsets of the 4-cube with at least 3 sites;
- all nine links PASS, with the same smallest margins as outputs/chain_full.txt (results_links.txt);
- at p = 9/25 the largest passing q on the grid 1/10000 is 3621/10000, and the exact threshold lies in
  [0.3621709283, 0.3621709284] (results_sharp.txt);
- every link's q is one or two grid steps below the largest passing q (results_slack.txt), and the chain cannot be
  continued below 0.3583 on that grid (results_below.txt; the fixed point lies between 0.3581 and 0.3582).

The reviewer did not check the two links of the alternative start (from 4428/10000).

Run from `independent/` (Python 3.8+, standard library; each script writes its own outputs and a cache file
Gtable.pkl next to itself):

    python selftest.py        # about two minutes; last line ALL SELF-TESTS PASSED
    python run.py links       # the nine links, about seven minutes
    python run.py sharp       # the threshold at p = 9/25
