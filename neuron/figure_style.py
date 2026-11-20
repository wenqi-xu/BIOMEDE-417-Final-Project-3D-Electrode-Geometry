"""Shared figure style: one fixed color and marker per electrode, publication settings."""
import matplotlib.pyplot as plt

STYLE = {
    "disc":         dict(color="#2a78d6", marker="o", ls="-", label="Disc"),
    "cylinder":     dict(color="#eb6834", marker="s", ls="-", label="Cylinder"),
    "cone":         dict(color="#1baf7a", marker="^", ls="-", label="Cone"),
    "ring":         dict(color="#eda100", marker="D", ls="-", label="Ring"),
    "point_source": dict(color="#7a7974", marker="x", ls="--", label="Point source"),
}

INK, INK2, AXON = "#0b0b0b", "#52514e", "#a9a8a2"
MM = 1 / 25.4
WIDTH_DOUBLE = 180 * MM
WIDTH_SINGLE = 85 * MM

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "mathtext.default": "regular",
    "axes.linewidth": 0.6, "axes.edgecolor": INK, "axes.labelcolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.labelpad": 2,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "xtick.color": INK, "ytick.color": INK,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "legend.frameon": False,
    "pdf.fonttype": 42, "savefig.dpi": 600, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
})


def marker_kw(shape, ms=3.2):
    st = STYLE[shape]
    edge = st["color"] if st["marker"] == "x" else "white"
    return dict(marker=st["marker"], color=st["color"], ms=ms, mec=edge, mew=0.9 if st["marker"] == "x" else 0.4)


def save(fig, path_no_ext):
    for ext in ("png", "pdf"):
        fig.savefig(f"{path_no_ext}.{ext}")
    plt.close(fig)
    print("saved", path_no_ext)
