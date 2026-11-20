"""Write the axon coordinate files (mm) for COMSOL "From file" export: coords/axon_xyz_d{distance}mm.txt."""
import os
import numpy as np
from neuron import h, load_mechanisms

HERE = os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")
DISTANCES_MM = [0.05, 0.1, 0.2, 0.5]

if os.getcwd().replace("\\", "/") != HERE:
    load_mechanisms(HERE)
h.load_file("stdrun.hoc")
h.load_file(HERE + "/mrg_axon.hoc")

secs = list(h.node) + list(h.MYSA) + list(h.FLUT) + list(h.STIN)
mid = h.node[int(h.axonnodes) // 2](0.5)
pos = np.array([h.distance(h.node[0](0.5), s(0.5)) for s in secs])
x_mm = (pos - h.distance(h.node[0](0.5), mid)) * 1e-3

os.makedirs(HERE + "/coords", exist_ok=True)
for d in DISTANCES_MM:
    xyz = np.column_stack([x_mm, np.zeros_like(x_mm), np.full_like(x_mm, d)])
    fname = f"{HERE}/coords/axon_xyz_d{d}mm.txt"
    np.savetxt(fname, xyz, fmt="%.9g", delimiter="\t")
    print(fname, xyz.shape)
