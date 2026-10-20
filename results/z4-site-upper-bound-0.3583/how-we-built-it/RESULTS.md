# Feasibility log (floating point; NOT certificates). Tool: feas.py (single-site + whole-face U only).

## Calibration, 2026-10-20: 3x3x3, largest:6, Monte Carlo 400k samples, seed 1
- f(0.4427) ~ 0.44247 (known: certified at q = 0.4428), f(0.45) ~ 0.46307 (paper Fig.: 0.46353),
  f(0.4634) ~ 0.50115 (known ~0.501). Agreement within ~5e-4 (MC noise); binding U = single corner site, as in the
  paper. (First run gave f = 0 everywhere: a rounding bug in the flow check, fixed before these numbers.)

## Z^4, 2x2x2x2 blocks, largest:4, EXACT enumeration (2^16 states)
- f(0.4164) ~ 0.49478 (+0.078), f(0.36) ~ 0.35675 (-0.003), f(0.30) ~ 0.22018 (-0.080).
- So the iteration's fixed point is a little above 0.36 for this rule: a chain started above the known d = 4 bound
  (0.4344 GPS26, or our derived 0.4164) would descend to about 0.36-0.37, if the full check (all U) agrees.
- Context: GPS26 0.4344, Jiang-Lang 0.4719, our derived 0.4164, numerical p_c(Z^4) ~ 0.197.
- Not yet tried: other m, all nonempty U, 2x2x2x3 blocks. Via GPS26 Prop. 3.4, a d = 4 bound near 0.36 gives only
  about 0.35 in d = 5, above Jiang-Lang's 0.2839: no d = 5 impact.

## Z^3, 4x4x4 blocks, gate round 1 (2026-10-20; MC 200k samples, seed 1; log gate_444.log)
Fixed points by linear interpolation of f(p) - p (uncertainty about +-0.002):
- V0 largest:8, full faces: ~0.4116 (reproduces the earlier ~0.410 estimate)
- V1 core:8 (largest cluster meeting the core {1,2}^3), full faces: ~0.4154  -> gate 1 passes, just
- V2 corepri:8 (first open core site's cluster, no sizes), full faces: ~0.426
- V3 core:8, centre 2x2 faces only: f(0.43) ~ 0.28 -> fails completely (fixed point > 0.5)
- V4 corepri:8, centre 2x2 faces only: f(0.43) ~ 0.27 -> fails completely
Conclusion: centre-only joining is dead; a verifier must handle (near-)full faces. Round 2 (the author's ideas):
core-site-count rules (coremost, coreseed:2/3) and richer faces (noncorner 12, mixednc = full sweep-end faces +
noncorner sides, mixed = full sweep-end faces + centre sides).

## Z^3, 4x4x4, round 2 (log gate_444_round2.log; same settings)
- R1 coremost:8 (most core sites; saturating size >= 8), full faces: ~0.419  <- best verifier-friendly rule so far
- R2 coreseed:2, full: ~0.431; R3 coreseed:3, full: > 0.45 (too strict)
- R4 core:8, noncorner faces (12): ~0.425; R5 core:8, mixednc (full sweep-end + noncorner sides): ~0.422
- R6 core:8, mixed (full sweep-end + centre sides): fails (> 0.5)
Note (the author's catch, 2026-10-20): "size >= m" is cheap in a verifier as a saturating counter min(size, m); only
"largest cluster" (exact sizes to 64) is expensive. Round 3 tests the free face-based filters (fhN = touches >= N
faces, fwN = >= N face sites) and the author's 8-site faces: pinwheel (centre 2x2 + one arm per edge) and belt
(full or noncorner sweep-end faces, middle-slice belt on the sides).

## Z^3, 4x4x4, round 3 (log gate_444_round3.log; same settings)
- F1 coremost:fh2 ~0.4256; F2 coremost:fh3 ~0.4227; F3 coremost:fw6 ~0.4229 (free filters; size filter R1 ~0.419)
- F4 corepri:fh2 ~0.439; F5 corepri:fw6 ~0.432
- P1 coremost:8 pinwheel (8/face), B1 belt, B2 beltnc: f(0.44) ~ 0.4115 -> fixed point > 0.45: fail.
  B1 and B2 agree to 5 decimals because the binding constraint is a side child's marginal (side faces are the same
  8-site belt in both); the sweep-end faces never bind.

## Summary of the 4x4x4 gate (fixed points; current certified bound 0.4427)
| verifier-friendly? | rule / faces | fixed point |
|---|---|---|
| no (exact sizes) | largest:8 / full | 0.4116 |
| no (exact sizes) | core:8 / full | 0.4154 |
| yes (core count + saturating size) | coremost:8 / full | 0.419 |
| yes (face filter, no size) | coremost:fh3 or fw6 / full | 0.423 |
| restricted faces | core:8 / full sweep-end + 12-site sides | 0.422 |
| restricted faces | core:8 / 12-site faces | 0.425 |
| restricted faces | any 8-site layout (pinwheel, belt) | > 0.45 (fails) |
| restricted faces | 4-site centre on any face | > 0.5 (fails) |
Reading: the prize for an exact 4x4x4 verifier is 0.4427 -> about 0.419-0.426; no verifier-friendly variant meets
the 0.415 gate. Face information cannot drop below ~12 sites per face, so a verifier must carry near-full side-face
footprints through the sweep: the hard part.

## Z^4, 2x2x2x2, full check of every entry set (2026-10-20; logs z4_2222_allU.log, z4_2222_filters.log)
New aggregated path (`feas.py ... full all`): every nonempty U (8 faces x 255 = 2040 laws), exact over all 2^16
states. Checks: in "test" mode it reproduces the earlier exact values (0.49478, 0.35675, 0.22018); on 3x3x3 MC 400k
with every U it reproduces the d = 3 calibration (f(0.4427) ~ 0.44247, f(0.45) ~ 0.46307; calib_333_allU.log).
Every U agrees with the single-site tests: a single entry site binds throughout.
| rule | f(0.35) | f(0.36) | f(0.37) | fixed point |
|---|---|---|---|---|
| largest:1 | 0.2999 | 0.3192 | 0.3387 | > 0.38 |
| largest:2 | 0.3246 | 0.3441 | 0.3639 | ~0.376 |
| largest:3 | 0.3407 | 0.3622 | 0.3838 | ~0.358  <- best |
| largest:4 | 0.3327 | 0.3568 | 0.3811 | ~0.362 |
| largest:5 | 0.2963 | 0.3215 | 0.3473 | > 0.38 |
| largest:6..10 | | | | far worse |
| largest:fh7 (touches >= 7 of 8 faces) | 0.3311 | 0.3553 | 0.3797 | ~0.363 |
| largest:fh8, fw24 | | | | far worse |
Reading: with 16 sites any rule is checkable by plain enumeration, so "largest" is fine. Best: largest:3, fixed
point ~0.358 against 0.4344 (GPS26, published) and 0.4164 (our Corollary). Not a d = 5 gain (Jiang-Lang 0.2839).
Certification differences from Part II: 7 children, so up-sets of {0,1}^7 (Dedekind number ~2.4e12) cannot be
listed; each domination must be certified by an explicit flow, checked in exact integer arithmetic (all inputs are
rational when p is). Part II's text is written for Z^3 (five children, 3x3x3) and would need restating for Z^d boxes.

## Z^4, 2x2x2x3 (logs z4_2223_allU.log, z4_2223_mc_high.log; diagnostic diag_marginals.py)
- Exact, every U (25080 laws, 7.89M distinct designated sets, ~6.5 min per p), largest:3: f(0.33) ~ 0.1984,
  f(0.34) ~ 0.2205 (run stopped after these two: the grid 0.33-0.37 cannot contain the fixed point).
- MC 400k (single-site tests), largest:3: f(0.40) ~ 0.3741, f(0.42) ~ 0.4326 -> fixed point ~0.414. Loses to 2^4.
- Cause, measured (the author's prediction): the coarse comparison is isotropic, so q* <= the weakest child marginal.
  Entering a 12-site side face at a corner, the far end of the length-3 axis (8 sites) has marginal 0.198 at p = 0.33
  (= f(0.33)) and 0.268 at p = 0.36, against 0.40-0.60 for the other six children; its root value r is 0.60-0.67
  against 0.77-0.84.
- Contrast, 2x2x2x2 at p = 0.36: child marginals 0.401-0.528 but q* = 0.362, so there the binding constraint is
  joint (correlations between children), not a marginal: the up-set / flow check matters in d = 4.
- The tie-break by site index breaks the block's symmetry (2^4 root values 0.6905-0.6969 by face), so a certificate
  checks all 8 faces x 255 sets U, as Part II does for its six faces.
