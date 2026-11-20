"""Make all figures in ../results from the COMSOL exports and thresholds.csv. Run after run_threshold.py."""
import os
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.lines import Line2D
from figure_style import STYLE, INK, INK2, AXON, MM, WIDTH_DOUBLE, WIDTH_SINGLE, marker_kw, save

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
EXPORT_DIR = os.path.join(PROJECT, "comsol_exports")
RESULTS_DIR = os.path.join(PROJECT, "results")
D = 0.05
N_NODES = 21
MID = N_NODES // 2

R_DISC = 56.42
R_IN, R_OUT = 30.0, 63.90
R_CYL, H_CYL = 25.39, 50.0
R_CONE, H_CONE = 46.58, 50.0
T_FLAT = 6.0


def load(shape, dist):
    f = os.path.join(EXPORT_DIR, f"{shape}_d{dist:g}mm.txt")
    if not os.path.exists(f):
        return None
    d = np.loadtxt([ln.replace(",", " ") for ln in open(f) if not ln.startswith("%")])
    return None if np.all(d[:, 3] == 0) else (d[:, 0], d[:, 3])


def peak_af(V):
    return -(V[MID - 1] - 2 * V[MID] + V[MID + 1])


def legend_above(ax, shapes, ncol=3, line=True):
    handles = [Line2D([], [], ls=STYLE[s]["ls"] if line else "none", lw=1.0, label=STYLE[s]["label"],
                      **marker_kw(s, ms=3.5)) for s in shapes]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=ncol,
              handlelength=1.8 if line else 1.0, columnspacing=1.2, borderaxespad=0)


def log_ticks(axis):
    axis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}"))
    axis.set_minor_formatter(plt.FuncFormatter(lambda v, _: f"{v:g}" if f"{v:g}"[0] in "25" else ""))


x_all = np.loadtxt(os.path.join(HERE, "coords", f"axon_xyz_d{D:g}mm.txt"))[:, 0]
SPACING = x_all[MID + 1] - x_all[MID]
data = {s: load(s, D) for s in STYLE}
shapes = [s for s in STYLE if data[s] is not None]
rows = [(r["shape"], float(r["distance_mm"]), float(r["threshold_uA"]))
        for r in csv.DictReader(open(os.path.join(RESULTS_DIR, "thresholds.csv")))]
out = lambda name: os.path.join(RESULTS_DIR, name)

# axon_model
fig, ax = plt.subplots(figsize=(WIDTH_DOUBLE, 32 * MM))
ax.plot([-11.9, 11.9], [0, 0], color=AXON, lw=3, solid_capstyle="round", zorder=1)
ax.plot(x_all[:N_NODES], np.zeros(N_NODES), "o", color=INK, ms=3, zorder=2)
ax.plot(0, -1.25, marker="^", color=INK2, ms=5)
ax.text(0, -1.9, "electrode", ha="center", va="top", color=INK2)
ax.annotate("", (0, -0.25), (0, -1.0), arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.6, shrinkA=0, shrinkB=0))
ax.text(0.25, -0.62, "d", color=INK2, va="center", style="italic")
ax.text(-11.9, 0.55, f"MRG axon, 10 µm fiber: 21 nodes of Ranvier, {SPACING:.2f} mm apart", color=INK)
ax.plot([8.9, 10.9], [-1.6, -1.6], color=INK, lw=1, solid_capstyle="butt")
ax.text(9.9, -1.85, "2 mm", ha="center", va="top")
ax.set(xlim=(-12.2, 12.2), ylim=(-2.5, 1.1))
ax.axis("off")
save(fig, out("axon_model"))

# internode_compartments
fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 34 * MM))
length = {"node": 1.0, "MYSA": 3.0, "FLUT": 46.0, "STIN": (SPACING * 1e3 - 1 - 2 * 3 - 2 * 46) / 6}
index = {"node": (0, N_NODES), "MYSA": (N_NODES, 61), "FLUT": (61, 101), "STIN": (101, 221)}
for row, (name, (i0, i1)) in enumerate(index.items()):
    xs = x_all[i0:i1] * 1e3
    xs = xs[(xs > -1) & (xs < SPACING * 1e3 + 1)]
    for xc in xs:
        ax.add_patch(Rectangle((xc - length[name] / 2, row - 0.22), length[name], 0.44, color=AXON, lw=0))
    ax.plot(xs, np.full(len(xs), row), "o", color=INK, ms=2)
ax.set(xlim=(-30, SPACING * 1e3 + 30), ylim=(3.6, -0.6), yticks=range(4), xticks=[0, 500, 1000],
       yticklabels=["Node (1 µm)", "MYSA (2 × 3 µm)", "FLUT (2 × 46 µm)", f"STIN (6 × {length['STIN']:.0f} µm)"],
       xlabel="Distance from node (µm)")
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
save(fig, out("internode_compartments"))

# electrode_geometry
geom = {
    "disc": [Rectangle((-R_DISC, -T_FLAT), 2 * R_DISC, T_FLAT)],
    "ring": [Rectangle((-R_OUT, -T_FLAT), R_OUT - R_IN, T_FLAT), Rectangle((R_IN, -T_FLAT), R_OUT - R_IN, T_FLAT)],
    "cylinder": [Rectangle((-R_CYL, -H_CYL), 2 * R_CYL, H_CYL)],
    "cone": [Polygon([(0, 0), (R_CONE, -H_CONE), (-R_CONE, -H_CONE)])],
}
dims = {
    "disc": f"r = {R_DISC:.0f} µm",
    "ring": f"r = {R_IN:.0f}–{R_OUT:.0f} µm",
    "cylinder": f"r = {R_CYL:.0f} µm, h = {H_CYL:.0f} µm",
    "cone": f"r = {R_CONE:.0f} µm, h = {H_CONE:.0f} µm",
}
d_um = D * 1e3
fig, axes = plt.subplots(1, 4, figsize=(WIDTH_DOUBLE, 34 * MM), gridspec_kw=dict(wspace=0.08))
for i, (shape, ax) in enumerate(zip(geom, axes)):
    ax.add_patch(Rectangle((-135, d_um - 5), 270, 10, color=AXON, lw=0))
    ax.plot(0, d_um, "o", color=INK, ms=2.5)
    for p in geom[shape]:
        p.set(facecolor=STYLE[shape]["color"], edgecolor=INK, lw=0.5)
        ax.add_patch(p)
    ax.text(0, -57, STYLE[shape]["label"], ha="center", va="top", color=INK, fontsize=7)
    ax.text(0, -79, dims[shape], ha="center", va="top", color=INK2, fontsize=6)
    if i == 0:
        ax.annotate("", (-95, d_um - 5), (-95, 0),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=0.6, shrinkA=0, shrinkB=0))
        ax.text(-90, d_um / 2 - 2, f"d = {d_um:.0f} µm", color=INK2, va="center", fontsize=6)
        ax.plot([-135, -85], [-62, -62], color=INK, lw=1, solid_capstyle="butt")
        ax.text(-110, -67, "50 µm", ha="center", va="top", fontsize=6)
    ax.set(xlim=(-140, 140), ylim=(-98, 62), aspect="equal")
    ax.axis("off")
save(fig, out("electrode_geometry"))

# Ve_along_axon
fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 62 * MM))
for s in shapes:
    x, V = data[s]
    order = np.argsort(x)
    z = 1 if s == "point_source" else 2
    ax.plot(x[order], V[order], color=STYLE[s]["color"], ls=STYLE[s]["ls"], zorder=z)
    ax.plot(x[:N_NODES], V[:N_NODES], ls="none", zorder=z + 1, **marker_kw(s))
ax.set(xlim=(-1.3, 1.3), ylim=(0, None), xticks=[-1, 0, 1], xlabel="Position along axon (mm)",
       ylabel=r"$V_{\mathrm{e}}$ (V mA$^{-1}$)")
legend_above(ax, shapes)
save(fig, out("Ve_along_axon"))

# activating_function
fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 62 * MM))
for s in shapes:
    x, V = data[s]
    Vn = V[:N_NODES]
    af = -(Vn[:-2] - 2 * Vn[1:-1] + Vn[2:])
    ax.plot(x[1:N_NODES - 1], af, ls=STYLE[s]["ls"], lw=0.8, zorder=1 if s == "point_source" else 2, **marker_kw(s))
ax.axhline(0, color=INK2, lw=0.5, zorder=0)
ax.set(xlim=(-3.6, 3.6), yticks=[-5, 0, 5, 10, 15], xlabel="Position along axon (mm)",
       ylabel=r"$-\Delta^{2} V_{\mathrm{e}}$ at nodes (V mA$^{-1}$)")
legend_above(ax, shapes)
save(fig, out("activating_function"))

# threshold_vs_distance
fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 62 * MM))
done = [s for s in STYLE if any(r[0] == s for r in rows)]
for s in done:
    d_, I_ = zip(*sorted((r[1] * 1e3, r[2]) for r in rows if r[0] == s))
    ax.plot(d_, I_, ls=STYLE[s]["ls"], zorder=1 if s == "point_source" else 2, **marker_kw(s))
dists = sorted({r[1] * 1e3 for r in rows})
ax.set(xscale="log", yscale="log", xlabel="Electrode–axon distance (µm)", ylabel="Threshold current (µA)")
ax.set_xticks(dists, [f"{v:g}" for v in dists])
ax.xaxis.set_minor_formatter(plt.NullFormatter())
log_ticks(ax.yaxis)
legend_above(ax, done)
save(fig, out("threshold_vs_distance"))

# threshold_vs_activating_function
pts = []
for s, dist, I in rows:
    V = load(s, dist)[1]
    pts.append((s, peak_af(V), I))
af_all = np.array([p[1] for p in pts])
I_all = np.array([p[2] for p in pts])
slope, icpt = np.polyfit(np.log10(af_all), np.log10(I_all), 1)
r2 = np.corrcoef(np.log10(af_all), np.log10(I_all))[0, 1] ** 2
fig, ax = plt.subplots(figsize=(WIDTH_SINGLE, 66 * MM))
xf = np.logspace(np.log10(af_all.min() / 1.3), np.log10(af_all.max() * 1.3), 50)
ax.plot(xf, 10 ** icpt * xf ** slope, color=INK2, lw=0.6, ls=":", zorder=0)
for s in done:
    ax.plot([p[1] for p in pts if p[0] == s], [p[2] for p in pts if p[0] == s], ls="none", zorder=2,
            **marker_kw(s, ms=4))
ax.text(0.97, 0.97, f"fit: slope {slope:.2f}, $R^2$ = {r2:.4f}".replace("-", "−"), transform=ax.transAxes,
        ha="right", va="top", fontsize=6, color=INK2)
ax.set(xscale="log", yscale="log", ylabel="Threshold current (µA)",
       xlabel=r"Peak activating function, $-\Delta^{2} V_{\mathrm{e}}$ (V mA$^{-1}$)")
log_ticks(ax.xaxis)
log_ticks(ax.yaxis)
legend_above(ax, done, line=False)
save(fig, out("threshold_vs_activating_function"))
