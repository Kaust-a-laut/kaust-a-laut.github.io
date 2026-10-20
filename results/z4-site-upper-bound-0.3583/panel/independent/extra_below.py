"""Extra: does the chain continue below 0.3583?  largest passing grid q for p = 3582, 3581, 3580 (exact)."""
from fractions import Fraction
from run import Laws, load_G, largest_grid_q, log
import os
fh = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_below.txt"), "w")
Gt = load_G()
for pk in (3582, 3581, 3580):
    L = Laws(Gt, Fraction(pk, 10000))
    k = largest_grid_q(L, pk, pk + 10)
    r = L.check(Fraction(pk + 1, 10000))
    log(f"p={pk}/10000: largest passing q = {k}/10000; link ({pk},{pk+1}) {'PASS' if r['pass'] else 'FAIL'}, min margin {float(r['worst'][0]):.4g}", fh=fh)
