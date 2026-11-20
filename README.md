# 3D Electrode Geometry Effects on Peripheral Nerve Activation

BIOMEDE 417 group project. COMSOL gives the extracellular potential for each electrode shape, then NEURON (MRG axon) finds the threshold current.

How it works

make_coords.py in the neuron folder writes the axon compartment positions to neuron/coords/axon_xyz_d{distance}mm.txt. In COMSOL, go to Results > Export > Data, set "Points to evaluate in" to From file, and pick the coordinate file for that distance. Save the export to comsol_exports/{shape}_d{distance}mm.txt, for example disc_d0.5mm.txt. Then run_threshold.py reads all the exports, finds the threshold for each one, and saves thresholds.csv and threshold_vs_distance.png to the results folder.

If your COMSOL model uses m instead of mm, put the raw exports in comsol_exports/raw_m_units and run python convert_to_mm.py ../comsol_exports/raw_m_units/* from the neuron folder. It writes the mm version to comsol_exports.

COMSOL settings

Same for every geometry. Length unit in mm (the coordinate files are in mm). Applied current 1 mA. Tissue conductivity 0.2 S/m, R_domain at least 25 mm. Exposed area 0.01 mm^2. The electrode axis is the z-axis, the point of the electrode closest to the axon is at the origin, and the rest of the electrode goes down into -z.

The axon runs along x from -11.5 to 11.5 mm, at y = 0 and z = distance, with the middle node right above the electrode.

Setup

Need NEURON 9 or newer, numpy and matplotlib. Compile the MRG mechanism once in the neuron folder: on Windows run mknrndll there (it makes nrnmech.dll), on Mac run nrnivmodl (refer to lecture tutorial). Then:

    cd neuron
    python make_coords.py      (only needed if the axon or distances change)
    python run_threshold.py    (reads ../comsol_exports/)

Mesh check

To compare meshes, put the mesh number in the shape name, for example disc-mesh1_d0.05mm.txt, disc-mesh2_d0.05mm.txt, and so on. Put these files in their own folder (like mesh_test) and run the script on that folder:

    python run_threshold.py ../mesh_test

Each mesh shows up as its own row in thresholds.csv. The mesh has converged when the threshold changes by less than about 2% between two refinements. Note that this run overwrites the files in the results folder.

Files

neuron/mrg_axon.hoc is the MRG axon model, taken from the course tutorial MRGaxon.hoc. neuron/AXNODE.mod is the node membrane mechanism, also from the course tutorial (originally McIntyre et al. 2002). neuron/run_threshold.py does the COMSOL to NEURON threshold search, e_extracellular = I(t) * V. neuron/make_coords.py writes the coordinate files for COMSOL. neuron/convert_to_mm.py converts exports from m to mm.
