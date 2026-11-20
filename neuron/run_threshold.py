"""Find the MRG activation threshold for every COMSOL export (tutorial method). Usage: python run_threshold.py [export_folder]"""
import os
import re
import sys
import glob
import numpy as np
from neuron import h, load_mechanisms

HERE = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
PROJECT = os.path.dirname(HERE)
EXPORT_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJECT, "comsol_exports")
RESULTS_DIR = os.path.join(PROJECT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

POLARITY = -1
DELAY = 1.0
PW = 0.1
SIGMA = 0.2
R_DOMAIN = 25.0

if os.getcwd().replace("\\", "/") != HERE:
    load_mechanisms(HERE)
h.load_file("stdrun.hoc")
h.load_file(HERE + "/mrg_axon.hoc")
n_nodes, n_sec = int(h.axonnodes), int(h.axontotal)

h("""
create electrode
objref exIClmp, apc, Vext
electrode exIClmp = new IClamp(0.5)
node[axonnodes-1] apc = new APCount(0.5)
Vext = new Vector(axontotal)

proc advance() {  local j, k
	j = 0
	for k=0,axonnodes-1  { node[k] e_extracellular(0.5) = exIClmp.i*Vext.x[j]  j += 1 }
	for k=0,paranodes1-1 { MYSA[k] e_extracellular(0.5) = exIClmp.i*Vext.x[j]  j += 1 }
	for k=0,paranodes2-1 { FLUT[k] e_extracellular(0.5) = exIClmp.i*Vext.x[j]  j += 1 }
	for k=0,axoninter-1  { STIN[k] e_extracellular(0.5) = exIClmp.i*Vext.x[j]  j += 1 }
	fadvance()
}
""")
setattr(h.exIClmp, "del", DELAY)
h.exIClmp.dur = PW


def skip(fname, why):
    print(f"  SKIPPED {os.path.basename(fname)}: {why}")


def load_comsol(fname, dist):
    d = np.loadtxt([ln.replace(",", " ") for ln in open(fname) if not ln.startswith("%")])
    if d.shape != (n_sec, 4):
        return skip(fname, f"expected {n_sec} rows x 4 columns, got {d.shape}")
    xyz, V = d[:, :3], d[:, 3]
    if not np.allclose(xyz, np.loadtxt(f"{HERE}/coords/axon_xyz_d{dist:g}mm.txt"), atol=1e-4):
        return skip(fname, f"coordinates do not match coords/axon_xyz_d{dist:g}mm.txt "
                           "(wrong coordinate file, unit m vs mm, or row order changed?)")
    if np.isnan(V).any():
        return skip(fname, "NaN in V -> some points are outside the COMSOL domain (R_domain too small?)")
    if np.all(V == 0):
        return skip(fname, "V is 0 everywhere -> no solution exported (study not computed, or wrong dataset?)")
    return V * 1e3


def trial(amp):
    h.exIClmp.amp = POLARITY * amp
    h.run()
    return h.apc.n > 0


def threshold(eps=0.01, lo=1e-3, hi=100.0, start=0.01):
    strength, lb, ub = start, lo, hi
    while lb == lo or ub == hi:
        if trial(strength):
            ub, strength = strength, strength / 2
        else:
            lb, strength = strength, strength * 2
        if lb > ub or strength > ub:
            return np.nan
    while abs(ub - lb) > eps * ub:
        strength = (ub + lb) / 2
        if trial(strength):
            ub = strength
        else:
            lb = strength
    return ub


files = sorted(glob.glob(os.path.join(EXPORT_DIR, "*_d*mm.txt")))
if not files:
    sys.exit(f"No files like {{shape}}_d{{distance}}mm.txt in {EXPORT_DIR}")

rows = []
for f in files:
    shape, dist = re.match(r"(.+)_d([\d.]+)mm\.txt$", os.path.basename(f)).groups()
    dist = float(dist)
    Ve = load_comsol(f, dist)
    if Ve is None:
        continue
    h.Vext.from_python(Ve)
    if shape.startswith("point"):
        analytic = 1e-3 / (4 * np.pi * SIGMA) * (1 / (dist * 1e-3) - 1 / (R_DOMAIN * 1e-3)) * 1e3
        print(f"  [check] {shape} d={dist} mm: Ve above electrode COMSOL {Ve[n_nodes // 2]:.1f} mV,"
              f" analytic {analytic:.1f} mV")
    Ith = threshold() * 1e3
    rows.append((shape, dist, Ith, Ith * PW))
    print(f"{shape:12s} d = {dist:5.2f} mm   threshold = {Ith:8.1f} uA   charge = {Ith * PW:7.2f} nC")

out_csv = os.path.join(RESULTS_DIR, "thresholds.csv")
with open(out_csv, "w") as fh:
    fh.write("shape,distance_mm,threshold_uA,charge_nC\n")
    for r in rows:
        fh.write("%s,%g,%.3f,%.4f\n" % r)
print("saved", out_csv)
