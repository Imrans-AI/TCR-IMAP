#!/usr/bin/env python3
"""PyMOL structural panels for Figure 5, p84v2.

Ported from 32b_pymol_panels_p82v1.py. Changes:
  * colours come from results/tables/contact_freq_gdomain182_descriptive_p84v2.npz: 182 G-domain
    columns keyed by (IMGT number, insertion code), so 1061/1061A and 1072/1072A are four columns
  * every residue is selected by IMGT number AND insertion code (PyMOL resi "1061" vs "1061A").
    p82v1 selected "resi 1061", which is the base residue, and coloured it with the insertion
    residue's values (CORRECTION_PLAN_P84V2.md section 1)
  * every selected residue is checked: the selection must hit exactly one residue whose resi string
    equals number+insertion code; inside the 82-position window (84 columns) its residue name must
    equal the residue manifest (branch residue-manifest, manifest.tsv sha256 d41a3ead...). The CDR
    loop selections of the complex panel are checked against the manifest TCR rows the same way.
    Any mismatch stops the script before anything is written.
  * TCR rows resolved in fewer than MIN_ROW_COMPLEXES complexes are excluded from the per-position
    sum: they are single-complex insertion rows (CDR3b 111B, 112C, 112B) whose frequencies are
    1/1 or 1/3, not population estimates
  * PyMOL max_threads 4; CPU ray tracing (cmd.png ray=1)

Panel B intensity = per-position contact frequency summed over the TCR rows resolved in at least
MIN_ROW_COMPLEXES complexes, normalised to its own maximum -- a RELATIVE scale, not panel A's
per-pair scale. Sticks: positions above 0.30 of that maximum.

Run headless:  pymol -cq scripts/55_pymol_fig5_p84v2.py
"""
import sys, os, subprocess, hashlib, io
import numpy as np, pandas as pd
from pymol import cmd

sys.path.insert(0, "src")
from viz.style import CMAP_SEQ, MOLECULE_COLORS, pymol_hex

cmd.set("max_threads", 4)
os.makedirs("figures/pymol", exist_ok=True)
os.makedirs("results/review/figures_p84v2", exist_ok=True)
PDB = "1ao7"
MIN_ROW_COMPLEXES = 5
HOT = 0.30

ART = "results/tables/contact_freq_gdomain182_descriptive_p84v2.npz"
_z = np.load(ART, allow_pickle=True)
num = [int(v) for v in _z["mhc_imgt_number"]]
ins = [str(v) for v in _z["mhc_insertion_code"]]
cols = [f"{n}{i}" for n, i in zip(num, ins)]
assert cols == [str(c) for c in _z["mhc_positions"]] and len(cols) == 182
tcr_rows = [str(t) for t in _z["tcr_position_labels"]]
nres_row = (_z["S"] == 1).any(axis=2).sum(axis=0)
rare = nres_row < MIN_ROW_COMPLEXES
rare_rows = [(t, int(n)) for t, n, r in zip(tcr_rows, nres_row, rare) if r]
assert [t for t, _ in rare_rows] == ["cdr3b:111B", "cdr3b:112C", "cdr3b:112B"], rare_rows
mass_all = _z["F"].sum(0) / _z["F"].sum(0).max()
mass = _z["F"][~rare].sum(0)
mass = mass / mass.max()
reg = (pd.read_csv("results/tables/gdomain182_region_classification_p84v2.csv")
         .assign(imgt_position=lambda d: d.imgt_position.astype(str))
         .set_index("imgt_position")["region"].reindex(cols).values)
assert not pd.isna(reg).any()

MAN_SHA = "d41a3eadc637e00e250d2c5fc8555b555d389c9d779a81e808db719af812b65a"
_raw = subprocess.run(["git", "show", "residue-manifest:results/review/residue_manifest/manifest.tsv"],
                      capture_output=True, check=True).stdout
assert hashlib.sha256(_raw).hexdigest() == MAN_SHA
man = pd.read_csv(io.BytesIO(_raw), sep="\t", dtype={"insertion_code": str}, keep_default_na=False)
man = man[(man.pdb == PDB) & (man.present_in_structure.astype(str) == "True")]
man["resi"] = man.imgt_number.astype(str) + man.insertion_code

idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")
r = idx.loc[PDB]
MHC, A, Bc, PEP = str(r.mhc_chain1), str(r.Achain), str(r.Bchain), str(r.antigen_chain)
man_mhc = man[man.chain_role == "MHC"].set_index("resi")
assert set(man_mhc.chain_id) == {MHC} and len(man_mhc) == 84

VIEW = (-0.740773141, -0.615407944,  0.269310892,
         0.671615958, -0.686662316,  0.278256685,
         0.013684267,  0.386998445,  0.921978831,
         0.000000000,  0.000000000,-165.660995483,
        78.402503967,   5.382022858,   2.056062698,
       130.608398438, 200.713592529, -20.000000000)

C_PEP = pymol_hex(MOLECULE_COLORS["peptide"])
C_TCRA = pymol_hex(MOLECULE_COLORS["TCR alpha"])
C_TCRB = pymol_hex(MOLECULE_COLORS["TCR beta"])
C_MHC = pymol_hex(MOLECULE_COLORS["MHC"])


def residues(sel):
    out = []
    cmd.iterate(f"({sel}) and name CA", "out.append((chain, resi, resn))", space={"out": out})
    return out


def base():
    cmd.reinitialize()
    cmd.set("max_threads", 4)
    cmd.load(f"data/structures/stcrdab/imgt/{PDB}.pdb", "cx")
    cmd.remove("solvent or hydrogens")
    cmd.remove(f"not chain {MHC}+{A}+{Bc}+{PEP}")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.38); cmd.set("specular", 0.12)
    cmd.set("cartoon_fancy_helices", 1); cmd.set("cartoon_smooth_loops", 1)
    cmd.set("stick_radius", 0.16)
    cmd.hide("everything")


# ---------- Panel: complex with CDR loops ----------
base()
cmd.show("cartoon", f"chain {MHC} and resi 1-90+1001-1090")
cmd.color(C_MHC, f"chain {MHC}")
cmd.show("sticks", f"chain {PEP} and not name N+C+O")
cmd.color(C_PEP, f"chain {PEP}")
cmd.show("cartoon", f"chain {A}+{Bc} and resi 1-130")
cmd.color(C_TCRA, f"chain {A} and resi 1-130")
cmd.color(C_TCRB, f"chain {Bc} and resi 1-130")
cmd.set("cartoon_transparency", 0.6, f"chain {A}+{Bc}")
CDR = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}
cdr_check = []
for ch, col, role in ((A, C_TCRA, "TCRa"), (Bc, C_TCRB, "TCRb")):
    for loop, (lo, hi) in CDR.items():
        sel = f"chain {ch} and resi {lo}-{hi}"
        got = {(c, ri, rn) for c, ri, rn in residues(sel)}
        m = man[(man.chain_role == role) & (man.imgt_number >= lo) & (man.imgt_number <= hi)]
        want = {(ch, ri, rn) for ri, rn in zip(m.resi, m.resname)}
        cdr_check.append(dict(chain=ch, loop=loop, selection=sel, n_selected=len(got), n_manifest=len(want),
                              inserted=";".join(sorted(ri for _, ri, _ in got if not ri.isdigit())),
                              match=got == want))
        cmd.set("cartoon_transparency", 0, sel)
        cmd.color(col, sel)
        cmd.set("cartoon_tube_radius", 0.4, sel)
cdr_check = pd.DataFrame(cdr_check)
assert cdr_check.match.all(), cdr_check
cmd.set_view(VIEW)
cmd.turn("x", 180)
cmd.zoom("polymer", 2)
cmd.png("figures/pymol/p84v2_panelA_complex.png", width=1600, height=1700, dpi=300, ray=1)

# ---------- Panel: platform footprint coloured by contact frequency ----------
base()
cmd.show("cartoon", f"chain {MHC} and resi 1-90+1001-1090")
cmd.color(C_MHC, f"chain {MHC}")
cmd.show("sticks", f"chain {PEP} and not name N+C+O")
cmd.color(C_PEP, f"chain {PEP}")

rows, hot, missing, zero_floor = [], [], [], 0
for p, n, m, ma, rg in zip(cols, num, mass, mass_all, reg):
    sel = f"chain {MHC} and resi {p}"           # number + insertion code, e.g. resi 1061A
    got = residues(sel)
    in_window = 50 <= n <= 90 or 1050 <= n <= 1090
    man_resn = man_mhc.resname.get(p, "") if in_window else ""
    if not got:
        missing.append(p)
        rows.append(dict(position=p, selection=sel, n_selected=0, resi="", resn="", in_window=in_window,
                         manifest_resn=man_resn, match=(not in_window) or p not in man_mhc.index,
                         rel_mass=round(float(m), 4), rel_mass_all_rows=round(float(ma), 4), hot=False))
        continue
    ok = len(got) == 1 and got[0][1] == p and (not in_window or got[0][2] == man_resn)
    # sqrt scaling: linear mapping leaves most of the platform near-white because a few
    # positions dominate the summed frequency. Ordering and range are unchanged.
    cmd.set_color(f"cf{p}", list(CMAP_SEQ(float(m) ** 0.5)[:3]))
    cmd.color(f"cf{p}", sel)
    if rg == "beta-sheet floor" and m == 0:
        zero_floor += 1
    if m > HOT:
        cmd.show("sticks", f"{sel} and not name N+C+O"); hot.append(p)
    rows.append(dict(position=p, selection=sel, n_selected=len(got), resi=got[0][1], resn=got[0][2],
                     in_window=in_window, manifest_resn=man_resn, match=ok, rel_mass=round(float(m), 4),
                     rel_mass_all_rows=round(float(ma), 4), hot=m > HOT))
chk = pd.DataFrame(rows)
# every manifest window residue must have been selected exactly once
assert set(chk[chk.in_window & (chk.n_selected == 1)].position) == set(man_mhc.index)
assert chk.match.all(), chk[~chk.match]
cmd.set_view(VIEW)
cmd.zoom(f"chain {MHC} and resi 50-90+1050-1090", 1.2)
cmd.png("figures/pymol/p84v2_panelB_footprint.png", width=1700, height=1400, dpi=300, ray=1)
cmd.save("figures/pymol/p84v2_session.pse")

chk.to_csv("results/review/figures_p84v2/fig5b_selection_check.tsv", sep="\t", index=False)
cdr_check.to_csv("results/review/figures_p84v2/fig5_complex_cdr_selection_check.tsv", sep="\t", index=False)
hot_all = [p for p, ma in zip(cols, mass_all) if ma > HOT]
print(f"\nsource          : {ART}")
print(f"rare TCR rows   : {rare_rows} (resolved in < {MIN_ROW_COMPLEXES} complexes; excluded from the sum)")
print(f"positions       : {len(cols)} requested, {len(cols)-len(missing)} present in {PDB}, "
      f"{len(missing)} absent")
if missing:
    print(f"  absent        : {missing[:20]}{' ...' if len(missing) > 20 else ''}")
print(f"selection check : {int(chk.match.sum())}/{len(chk)} OK; window residues verified against the manifest: "
      f"{int((chk.in_window & (chk.n_selected == 1)).sum())}/84 "
      f"(1061 {chk.set_index('position').resn['1061']}, 1061A {chk.set_index('position').resn['1061A']}, "
      f"1072 {chk.set_index('position').resn['1072']}, 1072A {chk.set_index('position').resn['1072A']})")
print(f"CDR loops       : {int(cdr_check.match.sum())}/{len(cdr_check)} selections equal the manifest; "
      f"inserted residues selected: {';'.join(x for x in cdr_check.inserted if x)}")
print(f"beta-floor at 0 : {zero_floor} positions coloured at the colormap minimum")
print(f"hot (>{HOT} rel, rows >= {MIN_ROW_COMPLEXES}) : {hot}")
print(f"hot (>{HOT} rel, all rows)      : {hot_all}")
print("wrote figures/pymol/p84v2_panelA_complex.png, p84v2_panelB_footprint.png, p84v2_session.pse")
