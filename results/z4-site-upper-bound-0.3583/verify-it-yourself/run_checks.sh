#!/bin/sh
# Run from the package folder: sh verify-it-yourself/run_checks.sh   (last line: PASS or FAIL)
# Default, about a minute (Python 3.8+, standard library only): rebuild the census of all 2^16 block states, compare it
# byte for byte with object/counts_b2222_s3.txt, and certify the last link of the chain, (p, q) = (3583, 3584)/10000.
# FULL=1, about ten minutes: certify all nine links from q1 = 4345/10000, and the two links of the alternative start
# from q = 4428/10000. The certificate runs as python -I (isolated: PYTHON* variables such as PYTHONOPTIMIZE are
# ignored), and cert4.py itself refuses to run under -O; every check in it is explicit, never an assert.
PY=${PYTHON:-python3}
command -v "$PY" >/dev/null 2>&1 || PY=python
"$PY" -c "import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)" 2>/dev/null || PY=python
W=$(mktemp -d)
"$PY" -I object/code/cert4.py counts "$W/counts.txt" || { echo "FAIL (census)"; exit 1; }
cmp "$W/counts.txt" object/counts_b2222_s3.txt || { echo "FAIL (census differs from object/counts_b2222_s3.txt)"; exit 1; }
if [ "${FULL:-0}" = "1" ]; then
  "$PY" -I object/code/cert4.py chain "$W/counts.txt" 4345/10000 3934/10000 3745/10000 3658/10000 3618/10000 \
    3599/10000 3590/10000 3586/10000 3584/10000 3583/10000 | tee "$W/main.txt"
  "$PY" -I object/code/cert4.py chain "$W/counts.txt" 4428/10000 3972/10000 3934/10000 | tee "$W/alt.txt"
  n=$(cat "$W/main.txt" "$W/alt.txt" | grep -c "CERTIFICATE: PASS")
  rm -rf "$W"
  if [ "$n" -ne 11 ]; then echo "FAIL (only $n of 11 links passed)"; exit 1; fi
  echo "PASS (full: census reproduced byte for byte, 9 of 9 links and 2 of 2 alternative links certified)"
else
  "$PY" -I object/code/cert4.py chain "$W/counts.txt" 3584/10000 3583/10000 | tee "$W/last.txt"
  n=$(grep -c "CERTIFICATE: PASS" "$W/last.txt")
  rm -rf "$W"
  if [ "$n" -ne 1 ]; then echo "FAIL (last link)"; exit 1; fi
  echo "PASS (census reproduced byte for byte, last link certified; FULL=1 reruns the whole chain)"
fi
