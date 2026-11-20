"""
COMSOL Ve  ->  MRG axon (NEURON)  ->  activation threshold
Same method as the BIOMEDE 417 tutorial (extracellular_stimulation_MRG_axon.hoc):
  e_extracellular(section) = I_stim(t) * Ve(section)
where Ve is the COMSOL potential for a 1 mA electrode current (I_0 = 1 mA in COMSOL).

Input: COMSOL "From file" exports in ../comsol_exports/, named {shape}_d{distance}mm.txt
       (evaluated at coords/axon_xyz_d{distance}mm.txt -> 221 rows: x y z V)

Run:
    python run_threshold.py                 # uses ../comsol_exports
    python run_threshold.py some/folder
"""
import os
import re
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from neuron import h, load_mechanisms

HERE = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")   # hoc needs '/'
PROJECT = os.path.dirname(HERE)
EXPORT_DIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(PROJECT, "comsol_exports")
RESULTS_DIR = os.path.join(PROJECT, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

# ---------------- settings ----------------
POLARITY = -1      # -1 cathodic, +1 anodic
DELAY = 1.0        # ms
PW = 0.1           # ms, pulse width
SIGMA = 0.2        # S/m      (only for the point-source sanity check)
R_DOMAIN = 25.0    # mm       (only for the point-source sanity check)

# ---------------- 1. build the MRG axon ----------------
if os.getcwd().replace("\\", "/") != HERE:          # NEURON auto-loads compiled mechanisms from the cwd
    load_mechanisms(HERE)                           # compiled AXNODE.mod (Windows dll or Mac/Linux folder)
h.load_file("stdrun.hoc")
h.load_file(HERE + "/mrg_axon.hoc")
n_nodes, n_sec = int(h.axonnodes), int(h.axontotal)

# ---------------- 2. Extracellular stimulation code (copied from BME 417 tutorial) ----------------
h("""
create electrode                 // dummy compartment that holds the IClamp
objref exIClmp, apc, Vext
electrode exIClmp = new IClamp(0.5)
node[axonnodes-1] apc = new APCount(0.5)   // AP counter at the last node
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
setattr(h.exIClmp, "del", DELAY)                    # 'del' is a Python keyword
h.exIClmp.dur = PW


def load_comsol(fname, dist):
    """Read a COMSOL 'From file' export -> Ve [mV per mA] at every compartment.
    Same as the tutorial's XYZVmat.scanf(...) + Vext.x[i] = XYZVmat.x[i][3]*1e3,
    plus a check that the rows match our coordinate file (catches unit / order mistakes)."""
    d = np.loadtxt(fname, comments="%")             # columns: x y z V   (header lines start with %)
    if d.shape != (n_sec, 4):
        sys.exit(f"{fname}: expected {n_sec} rows x 4 columns, got {d.shape}")
    xyz, V = d[:, :3], d[:, 3]
    ref = np.loadtxt(f"{HERE}/coords/axon_xyz_d{dist:g}mm.txt")
    if not np.allclose(xyz, ref, atol=1e-4):   # 0.1 um tolerance
        sys.exit(f"{fname}: coordinates do not match coords/axon_xyz_d{dist:g}mm.txt "
                 f"(wrong coordinate file, unit m vs mm, or row order changed?)")
    if np.isnan(V).any():
        sys.exit(f"{fname}: NaN in V -> some points are outside the COMSOL domain (R_domain too small?)")
    return V * 1e3      # V -> mV


def trial(amp):
    """Run one simulation with |I| = amp [mA]. Returns True if an AP reaches the last node."""
    h.exIClmp.amp = POLARITY * amp
    h.run()
    return h.apc.n > 0


def threshold(eps=0.01, lo=1e-3, hi=100.0, start=0.01):
    """Binary search (bracket, then bisect). Returns threshold in mA.
    start is low (tutorial uses 1 mA): close to the electrode, >= 1 mA BLOCKS conduction."""
    strength, lb, ub = start, lo, hi
    while lb == lo or ub == hi:                     # bracket
        if trial(strength):
            ub, strength = strength, strength / 2
        else:
            lb, strength = strength, strength * 2
        if lb > ub or strength > ub:
            return np.nan
    while abs(ub - lb) > eps * ub:                  # bisect
        strength = (ub + lb) / 2
        if trial(strength):
            ub = strength
        else:
            lb = strength
    return ub


# ---------------- 3. loop over all exported files ----------------
files = sorted(glob.glob(os.path.join(EXPORT_DIR, "*_d*mm.txt")))
if not files:
    sys.exit(f"No files like {{shape}}_d{{distance}}mm.txt in {EXPORT_DIR}")

rows = []                                           # (shape, distance mm, threshold uA, charge nC)
for f in files:
    shape, dist = re.match(r"(.+)_d([\d.]+)mm\.txt$", os.path.basename(f)).groups()
    dist = float(dist)
    Ve = load_comsol(f, dist)
    h.Vext.from_python(Ve)
    if shape.startswith("point"):                   # sanity check: grounded sphere, I/(4*pi*sigma)*(1/r - 1/R)
        analytic = 1e-3 / (4 * np.pi * SIGMA) * (1 / (dist * 1e-3) - 1 / (R_DOMAIN * 1e-3)) * 1e3
        print(f"  [check] {shape} d={dist} mm: Ve above electrode COMSOL {Ve[n_nodes // 2]:.1f} mV,"
              f" analytic {analytic:.1f} mV")
    Ith = threshold() * 1e3                         # uA
    rows.append((shape, dist, Ith, Ith * PW))       # uA * ms = nC
    print(f"{shape:10s} d = {dist:5.2f} mm   threshold = {Ith:8.1f} uA   charge = {Ith * PW:7.2f} nC")

# save table
out_csv = os.path.join(RESULTS_DIR, "thresholds.csv")
with open(out_csv, "w") as fh:
    fh.write("shape,distance_mm,threshold_uA,charge_nC\n")
    for r in rows:
        fh.write("%s,%g,%.3f,%.4f\n" % r)
print("saved", out_csv)

# plot: threshold vs distance, one line per shape
fig, ax = plt.subplots(figsize=(5.5, 4))
for shape in sorted({r[0] for r in rows}):
    d_, I_ = zip(*sorted((r[1], r[2]) for r in rows if r[0] == shape))
    ax.plot(d_, I_, "o-", label=shape)
ax.set(xscale="log", yscale="log", xlabel="electrode-to-axon distance (mm)",
       ylabel="threshold current (uA)",
       title=f"MRG {h.fiberD:g} um, {PW} ms {'cathodic' if POLARITY < 0 else 'anodic'} pulse")
plain = plt.FuncFormatter(lambda v, _: f"{v:g}")     # plain numbers instead of 10^x labels
ax.xaxis.set_minor_formatter(plt.NullFormatter())
ax.set_xticks(sorted({r[1] for r in rows}))          # one tick per tested distance
ax.xaxis.set_major_formatter(plain)
ax.yaxis.set_major_formatter(plain)
ax.yaxis.set_minor_formatter(plain)
ax.legend(); ax.grid(alpha=.3, which="both")
fig.tight_layout()
fig.savefig(os.path.join(RESULTS_DIR, "threshold_vs_distance.png"), dpi=150)
