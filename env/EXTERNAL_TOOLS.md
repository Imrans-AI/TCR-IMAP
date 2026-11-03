# External tools, not installable from the environment file

| tool | version used | needed for | where to get it |
|---|---|---|---|
| GROMACS | 2023.2 | running the molecular-dynamics replicates (not needed to reuse their analysis outputs) | https://manual.gromacs.org |
| CHARMM-GUI Solution Builder | web service, used 2026 | rebuilding the four simulation systems from `md/systems/BUILD_MANIFEST.tsv` | https://charmm-gui.org |
| DSSP (mkdssp) | 4.2.2 | secondary-structure assignment in the footprint analysis | conda-forge `dssp`, or https://github.com/PDB-REDO/dssp |
| PyMOL | 2.5+, incremental build | the structural panels of Figures 1, 5 and 6 | https://pymol.org — the rendered panels are included, so PyMOL is needed only to redraw them |

None of these is required by the worked example or by any of the stored-result reproduction
commands in `REPRODUCE.md`.
