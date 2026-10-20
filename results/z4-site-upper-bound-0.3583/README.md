# p_c^site(Z^4) <= 0.3583: an explicit upper bound for site percolation on the four-dimensional lattice

## What

In words: the critical probability of site percolation on the hypercubic lattice Z^4 is at most 3583/10000. The
best published rigorous upper bound we found was 0.4344 (Gomes, Pereira and Sanchis, Table 1); Jiang and Lang give
0.4719, and the bound p_c^site(Z^3) <= 0.4427 gives 0.4164 through Gomes, Pereira and Sanchis, Proposition 3.4.

The proof is Section 8 of Part II of the paper *Explicit rigorous bounds for the site percolation threshold of the
simple cubic lattice* (Caleb Ashton). It runs the block argument of Part II (a growth process in the style of van den
Berg and Ermakov, with Strassen couplings) on Z^4 with the blocks {0,1}^4. A block's designated cluster is its largest
open cluster, ties going to the cluster that contains the smallest index 8x1+4x2+2x3+x4, if it has at least 3 sites.
Each step of the exploration must dominate Bernoulli(q)^7 (seven children); this reduces to root conditions on the 8
faces and to 8 x 255 = 2040 laws on {0,1}^7. That cube has 2,414,682,040,998 up-sets, too many to list, so each
domination is certified by an explicit flow on the cover graph of {0,1}^7 carrying Bernoulli(q)^7 upward into the law
(a monotone coupling), and the flow is checked directly. All 2^16 block states are enumerated, and everything is exact
integer arithmetic. A chain of nine such certificates goes from theta_4(0.4345) > 0 (the bound 0.4344) down to
p = 0.3583. An alternative start, theta_4(0.4428) > 0 from p_c^site(Z^3) <= 0.4427 (Part II of the same paper; Z^3
is a subgraph of Z^4), adds two links that join the chain at its second link.

## Check it yourself

    sh verify-it-yourself/run_checks.sh             # about a minute; last line PASS
    FULL=1 sh verify-it-yourself/run_checks.sh      # about ten minutes; all 9 + 2 links

Python 3.8 or later, standard library only. The default run rebuilds the census of all 2^16 states, compares it byte
for byte with object/counts_b2222_s3.txt and certifies the last link; FULL=1 certifies every link. Each link prints
`CERTIFICATE: PASS (p = ..., q = ...)` with its smallest margin; compare with outputs/chain_full.txt and
outputs/chain_from_z3.txt. The check also rejects: outputs/sharpness.txt shows (p, q) = (0.36, 0.3621) accepted and
(0.36, 0.3622) rejected, against a floating-point threshold of 0.36217.

## Files

- object/code/cert4.py: the certificate (enumeration, exact laws, flow certificates validated directly, margins as
  minimum cuts matched by flows).
- object/counts_b2222_s3.txt: the census; one line per designated set (a 16-bit mask), followed by the number of
  states with 0, 1, ..., 16 open sites that have it. 37246 lines, 65536 states.
- object/bound.json: the claim, the chain and its margins.
- outputs/: logs of the certified runs (chain_full.txt, chain_from_z3.txt), the accept/reject test (sharpness.txt),
  the floating-point design of the links (design_chain.txt) and the d = 5..9 consequences (corollary4_output.txt).
- panel/: an independent re-check in separate code, written from the definitions alone during an AI review; it
  decides domination from the exact minimum over up-sets, computed as a minimum cut (the dual of the coupling flow),
  and confirms all nine links with the same margins (panel/README.md).
- how-we-built-it/: the floating-point exploration that chose the block and the rule (feas.py, design_chain.py,
  diag_marginals.py; RESULTS.md is the log of every round, including a 4x4x4 study for Z^3 that was set aside) and
  the d = 5..9 script (corollary4.py with corollary.py; needs mpmath).

## Compared to

| source | cell | value |
|---|---|---|
| P. A. Gomes, A. Pereira and R. Sanchis, J. Appl. Probab. (2026), doi:10.1017/jpr.2025.10041 (arXiv:2106.10388v2), Table 1 | site percolation, d = 4 | 0.4344 |
| J. Jiang and S. Lang, arXiv:2609.29806 (2026), Table 2 | site percolation, d = 4 | 0.4719 |

Monte Carlo estimates put p_c^site(Z^4) near 0.197; that is an estimate, not a bound.

## Not claimed

- Not the value of p_c, and not a bound below 0.3583: this block and rule have a fixed point near 0.3581.
- No new bound in dimensions 5 to 9: through Gomes, Pereira and Sanchis this gives 0.3446 (d = 5) and 0.1990
  (d = 8), above Jiang and Lang's 0.2839 and 0.1409.
- Nothing about bond percolation.

## Credits

- J. van den Berg and A. Ermakov, A new lower bound for the critical probability of site percolation on the square
  lattice, Random Structures Algorithms 8 (1996) 199-212 (the growth-process and Strassen-coupling method).
- V. Strassen, The existence of probability measures with given marginals, Ann. Math. Statist. 36 (1965) 423-439.
- P. A. Gomes, A. Pereira and R. Sanchis, Upper bounds for critical probabilities in Bernoulli percolation models,
  J. Appl. Probab. (2026) (the base 0.4344, the table compared against, and Lemma 4.2 and Proposition 3.4).
- J. Jiang and S. Lang, Coexistence of infinite clusters for percolation and Ising model on Z^d, arXiv:2609.29806.
- Found by Caleb Ashton (GitHub: Kaust-a-laut), with Claude Code (Anthropic) running claude-opus-5-5.

Licence: source code (`*.py`, `*.sh`) under the GNU Affero General Public License v3.0 only; everything else under
Creative Commons Attribution 4.0 International (the licence files are in the record that holds this package).
