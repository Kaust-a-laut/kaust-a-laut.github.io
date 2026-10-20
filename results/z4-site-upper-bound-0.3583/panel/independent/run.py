"""Main verification run.

  PYTHONUTF8=1 python run.py links     -> step 1 + step 2 (nine links), writes results_links.{txt,json}
  PYTHONUTF8=1 python run.py sharp     -> step 3 (p = 9/25), writes results_sharp.{txt,json}
  PYTHONUTF8=1 python run.py slack     -> extra: largest grid q for each link's p, writes results_slack.{txt,json}

Every PASS/FAIL decision below is made on exact integers.
"""
import json
import os
import pickle
import sys
import time
from fractions import Fraction

from core import (FACES, FACE_MASK, NSITES, Model, enumerate_G, from_local, coords,
                  min_nontrivial_upset, minimal_elements, bern_weights)

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()

LINKS = [(3934, 4345), (3745, 3934), (3658, 3745), (3618, 3658), (3599, 3618),
         (3590, 3599), (3586, 3590), (3584, 3586), (3583, 3584)]


def log(*a, fh=None):
    s = " ".join(str(x) for x in a)
    print(f"[{time.time() - T0:8.1f}s] {s}", flush=True)
    if fh:
        fh.write(s + "\n")
        fh.flush()


def face_name(f):
    i, s = FACES[f]
    return f"F_{{{i},{s}}} (x{i}={s})"


def sites_of(mask):
    return [coords(v) for v in range(NSITES) if mask >> v & 1]


def load_G():
    path = os.path.join(HERE, "Gtable.pkl")
    if os.path.exists(path):
        return pickle.load(open(path, "rb"))
    Gt = enumerate_G(3)
    pickle.dump(Gt, open(path, "wb"))
    return Gt


class Laws:
    """all 2040 entry laws at a fixed p, as exact integers."""

    def __init__(self, Gt, p):
        t = time.time()
        self.model = Model(Gt, p)
        self.faces = []
        for f in range(8):
            children, laws = self.model.entry_laws(f)
            self.faces.append((children, laws))
        self.build_time = time.time() - t

    def root_margin(self, q):
        """exact min_F r_F(F) - q"""
        m = self.model
        return min(Fraction(m.R[f][255], m.B16) for f in range(8)) - Fraction(q)

    def law_margin(self, f, lU, q, Wn=None, Wd=None):
        """exact min over nontrivial up-sets A of nu_{f,U}(A) - Bern(q)^7(A); returns (Fraction, A)"""
        if Wn is None:
            Wn, Wd = bern_weights(q, 7)
        children, laws = self.faces[f]
        Nnum, Zn = laws[lU]
        b112 = self.model.B16 ** 7
        K = Zn * b112
        W = [Nnum[s] * Wd - K * Wn[s] for s in range(128)]
        assert sum(W) == 0
        mv, A = min_nontrivial_upset(W, 7)
        # nu(A)-mu(A) = W(A) / (b^112 * d^7 * Zn)
        return Fraction(mv, K * Wd), A

    def check(self, q, order=None, stop_on_fail=False):
        """returns dict with exact results for the link (p, q)"""
        q = Fraction(q)
        Wn, Wd = bern_weights(q, 7)
        res = {"root_margin": self.root_margin(q)}
        worst = None
        nfail = 0
        per = {}
        items = order or [(f, lU) for f in range(8) for lU in range(1, 256)]
        for f, lU in items:
            mv, A = self.law_margin(f, lU, q, Wn, Wd)
            per[(f, lU)] = mv
            if mv < 0:
                nfail += 1
            if worst is None or mv < worst[0]:
                worst = (mv, f, lU, A)
            if stop_on_fail and mv < 0:
                break
        res.update(worst=worst, nfail=nfail, per=per, nlaws=len(per))
        res["pass"] = res["root_margin"] >= 0 and nfail == 0 and len(per) == 2040
        if stop_on_fail and (res["root_margin"] < 0 or nfail):
            res["pass"] = False
        return res


def describe_A(A, children):
    mins = minimal_elements(A, 7)
    out = []
    for s in mins:
        out.append("{" + ", ".join(f"x{FACES[children[j]][0]}={FACES[children[j]][1]}"
                                   for j in range(7) if s >> j & 1) + "}")
    return out


def run_links():
    fh = open(os.path.join(HERE, "results_links.txt"), "w", encoding="utf-8")
    t = time.time()
    Gt = load_G()
    log(f"STEP 1: distinct designated sets G (incl. empty) = {len(Gt)}", fh=fh)
    rows = []
    for (pk, qk) in LINKS:
        p, q = Fraction(pk, 10000), Fraction(qk, 10000)
        t1 = time.time()
        L = Laws(Gt, p)
        t2 = time.time()
        r = L.check(q)
        t3 = time.time()
        mv, f, lU, A = r["worst"]
        children = L.faces[f][0]
        Umask = from_local(f, lU)
        # per-face minima (to show the face asymmetry)
        perface = [min(v for (ff, uu), v in r["per"].items() if ff == f2) for f2 in range(8)]
        rootvals = [Fraction(L.model.R[f2][255], L.model.B16) for f2 in range(8)]
        row = {
            "p": f"{pk}/10000", "q": f"{qk}/10000",
            "pass": r["pass"],
            "root_margin": float(r["root_margin"]),
            "root_values": [float(x) for x in rootvals],
            "min_margin": float(mv),
            "min_margin_exact_num_digits": len(str(mv.numerator)),
            "min_face": face_name(f),
            "min_U_sites": sites_of(Umask),
            "min_A_minimal_elements": describe_A(A, children),
            "laws_checked": r["nlaws"], "laws_failing": r["nfail"],
            "all_minimisers": [(face_name(ff), sites_of(from_local(ff, uu)))
                               for (ff, uu), v in sorted(r["per"].items()) if v == mv],
            "per_face_min_margin": [float(x) for x in perface],
            "time_laws_s": round(t2 - t1, 1), "time_flows_s": round(t3 - t2, 1),
        }
        rows.append(row)
        log(f"LINK p={pk}/10000 q={qk}/10000: {'PASS' if row['pass'] else 'FAIL'}; "
            f"root min r_F(F)-q = {row['root_margin']:.6f}; min margin = {row['min_margin']:.4g} "
            f"at {row['min_face']}, U = {row['min_U_sites']}; A minimal elements = {row['min_A_minimal_elements']}; "
            f"failing laws = {row['laws_failing']}/{row['laws_checked']}; "
            f"all (F,U) attaining the exact minimum = {row['all_minimisers']}; "
            f"per-face minima = {[f'{x:.4g}' for x in row['per_face_min_margin']]}; "
            f"r_F(F) = {[f'{x:.6f}' for x in row['root_values']]}; "
            f"time laws {row['time_laws_s']}s flows {row['time_flows_s']}s", fh=fh)
        del L
    log(f"total time {time.time() - t:.1f}s", fh=fh)
    json.dump({"distinct_G": len(Gt), "links": rows}, open(os.path.join(HERE, "results_links.json"), "w"), indent=1)


def largest_grid_q(L, lo, hi, fh=None):
    """largest k in [lo, hi] with link (p, k/10000) passing; requires pass(lo), not pass(hi)."""
    order = None
    while L.check(Fraction(hi, 10000), stop_on_fail=True)["pass"]:
        lo, hi = hi, min(10000, 2 * hi - lo)
    while not L.check(Fraction(lo, 10000), stop_on_fail=True)["pass"]:
        lo, hi = max(0, 2 * lo - hi), lo
    log(f"   bracket: q={lo}/10000 passes, q={hi}/10000 fails", fh=fh)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        r = L.check(Fraction(mid, 10000), order=order, stop_on_fail=True)
        # check the worst laws first next time
        if r["worst"] is not None:
            f, lU = r["worst"][1], r["worst"][2]
            base = [(ff, uu) for ff in range(8) for uu in range(1, 256)]
            base.remove((f, lU))
            order = [(f, lU)] + base
        log(f"   q={mid}/10000 -> {'pass' if r['pass'] else 'fail'}", fh=fh)
        if r["pass"]:
            lo = mid
        else:
            hi = mid
    return lo


def run_sharp():
    fh = open(os.path.join(HERE, "results_sharp.txt"), "w", encoding="utf-8")
    Gt = load_G()
    p = Fraction(9, 25)
    t = time.time()
    L = Laws(Gt, p)
    log(f"laws at p=9/25 built in {L.build_time:.1f}s", fh=fh)
    log(f"root: min_F r_F(F) = {float(L.root_margin(0)):.6f}", fh=fh)
    k = largest_grid_q(L, 3600, 4000, fh=fh)
    # full exact confirmation at k and k+1 (all 2040 laws, no early stop)
    r1 = L.check(Fraction(k, 10000))
    r2 = L.check(Fraction(k + 1, 10000))
    assert r1["pass"] and not r2["pass"]
    out = {"p": "9/25", "largest_q_on_grid": f"{k}/10000"}
    for name, r, kk in (("at_q", r1, k), ("at_q_plus", r2, k + 1)):
        mv, f, lU, A = r["worst"]
        out[name] = {"q": f"{kk}/10000", "pass": r["pass"], "root_margin": float(r["root_margin"]),
                     "min_margin": float(mv), "face": face_name(f),
                     "U_sites": sites_of(from_local(f, lU)),
                     "A_minimal": describe_A(A, L.faces[f][0]), "laws_failing": r["nfail"]}
        log(f"q={kk}/10000: {'PASS' if r['pass'] else 'FAIL'}; min margin {float(mv):.4g} at {face_name(f)}, "
            f"U={sites_of(from_local(f, lU))}; failing laws {r['nfail']}/2040; A minimal = {describe_A(A, L.faces[f][0])}", fh=fh)
    # finer: bisection on exact rationals to locate the real threshold f(9/25)
    lo, hi = Fraction(k, 10000), Fraction(k + 1, 10000)
    order = [(r2["worst"][1], r2["worst"][2])] + [(ff, uu) for ff in range(8) for uu in range(1, 256)
                                                    if (ff, uu) != (r2["worst"][1], r2["worst"][2])]
    for _ in range(20):
        mid = (lo + hi) / 2
        r = L.check(mid, order=order, stop_on_fail=True)
        if r["pass"]:
            lo = mid
        else:
            hi = mid
    out["threshold_bracket"] = [float(lo), float(hi)]
    log(f"exact bisection: threshold q*(9/25) in [{float(lo):.10f}, {float(hi):.10f}]", fh=fh)
    log(f"step 3 time {time.time() - t:.1f}s", fh=fh)
    json.dump(out, open(os.path.join(HERE, "results_sharp.json"), "w"), indent=1)


def run_slack():
    fh = open(os.path.join(HERE, "results_slack.txt"), "w", encoding="utf-8")
    Gt = load_G()
    out = []
    for (pk, qk) in LINKS:
        L = Laws(Gt, Fraction(pk, 10000))
        k = largest_grid_q(L, qk, 5000)
        out.append({"p": pk, "q_used": qk, "largest_q": k})
        log(f"p={pk}/10000: largest passing q on grid = {k}/10000 (link uses {qk}/10000)", fh=fh)
        del L
    json.dump(out, open(os.path.join(HERE, "results_slack.json"), "w"), indent=1)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "links"
    {"links": run_links, "sharp": run_sharp, "slack": run_slack}[what]()
