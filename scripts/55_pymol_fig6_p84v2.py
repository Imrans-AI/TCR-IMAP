#!/usr/bin/env python3
"""PyMOL structural panels for Figure 6, p84v2: observed vs TCR-IMAP-predicted contact sites
for the largest positive and largest negative tcr70 pairwise delta vs the pair-frequency prior,
in the style of SG-TPMI (Gao et al., JCIM 2025) Figure 7.

Ported from 32c_pymol_fig6_p82v1.py. Changes:
  * the two complexes are re-derived, by the same rule, from
    results/tables/percomplex_tcr70_paired_p84v2.tsv (max and min delta); under p82v1 they were
    5yxn and 7n1f, under p84v2 they are whatever the table gives (asserted below and reported)
  * predictions: tcr70 OUT-OF-FOLD pairwise predictions only (results/tables/pairmaps_tcr70_p84v2.npz,
    pred_model; every complex is scored by the fold model that never saw it). No inference here.
  * residue mapping: p82v1 re-keyed residues by IMGT number alone (as 18b did). p84v2 reads the
    residue id of every row/column straight from data/processed/pair_maps_p84v2.npz
    (tcr_resid / mhc_resid, IMGT number + insertion code) and selects it in PyMOL by that string.
    Every displayed residue is checked: exactly one C-alpha in the structure, its residue name
    equals the p84v2 token (tcr_aa / mhc_aa) and the residue manifest (branch residue-manifest,
    manifest.tsv sha256 d41a3ead...) at the same (chain, IMGT number, insertion code).
  * the platform frame uses the 84 window columns (82 + 1061A + 1072A), keyed by residue id

Predicted set rule, unchanged: rank the resolved candidate pairs (S == 1) by TCR-IMAP score and
take the top k, k = that complex's number of observed contacts (Y == 1, 4.5 A minimum heavy-atom
distance) -- so observed and predicted sets have equal size.

Rendering unchanged from 32c (view down onto the alpha1/alpha2 platform from the TCR side, PCA
frame, normal toward the TCR, long axis peptide N -> C; MHC surface grey, peptide sticks yellow,
CDR loops thin translucent tubes; MHC site residues coloured by the contacting TCR chain, ties
split by atom; TCR site residues as C-alpha spheres). PyMOL max_threads 4, CPU ray tracing only.

Run headless:  pymol -cq scripts/55_pymol_fig6_p84v2.py
Writes figures/pymol/p84v2_fig6_<pdb>_{observed,predicted}.png and
figures/pymol/fig6_sites_p84v2.tsv (every displayed pair, for the validation table).
"""
import sys, os, io, subprocess, hashlib
import numpy as np, pandas as pd
from pymol import cmd
from Bio.PDB import PDBParser

sys.path.insert(0, "src")
from viz.style import MOLECULE_COLORS, pymol_hex

cmd.set("max_threads", 4)
os.makedirs("figures/pymol", exist_ok=True)
SEG_RANGE = {"cdr1": (27, 38), "cdr2": (56, 65), "cdr3": (105, 117)}   # loops drawn (32b)
SITES_TSV = "figures/pymol/fig6_sites_p84v2.tsv"

pr = pd.read_csv("results/tables/percomplex_tcr70_paired_p84v2.tsv", sep="\t")
ordr = np.argsort(-pr.delta.values, kind="stable")
TARGETS = [pr.pdb.values[ordr[0]], pr.pdb.values[ordr[-1]]]
assert TARGETS == ["3qfj", "3qeq"], TARGETS

pm = np.load("results/tables/pairmaps_tcr70_p84v2.npz", allow_pickle=True)
lab = np.load("data/processed/pair_maps_p84v2.npz", allow_pickle=True)   # p84v2 label provenance
assert (pm["Y"] == lab["Y"]).all() and (pm["S"] == lab["S"]).all() and (pm["pdb"] == lab["pdb"]).all()
assert bool(pm["covered"].all()) and str(lab["candidate_set"]) == "imgt84_p84v2"
POS = [str(x) for x in lab["col_label"]]
TLAB = [str(x) for x in lab["row_label"]]
assert len(POS) == 84 and len(TLAB) == 75 and "1061A" in POS and "1072A" in POS
pdbs = list(pm["pdb"])
idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")

MAN_SHA = "d41a3eadc637e00e250d2c5fc8555b555d389c9d779a81e808db719af812b65a"
_raw = subprocess.run(["git", "show", "residue-manifest:results/review/residue_manifest/manifest.tsv"],
                      capture_output=True, check=True).stdout
assert hashlib.sha256(_raw).hexdigest() == MAN_SHA
man = pd.read_csv(io.BytesIO(_raw), sep="\t", dtype={"insertion_code": str}, keep_default_na=False)
man = man[man.present_in_structure.astype(str) == "True"]
man["resi"] = man.imgt_number.astype(str) + man.insertion_code
MAN = {(r.pdb, r.chain_id, r.resi): (r.resname, r.aa1) for r in man.itertuples()}

C_PEP, C_MHC = pymol_hex(MOLECULE_COLORS["peptide"]), pymol_hex(MOLECULE_COLORS["MHC"])
C_CH = {"a": pymol_hex(MOLECULE_COLORS["TCR alpha"]), "b": pymol_hex(MOLECULE_COLORS["TCR beta"])}


def rid(res): return f"{res.id[1]}{res.id[2].strip()}"


rows, summary, ORIENT, checks = [], [], {}, []
for pdb in TARGETS:
    i = pdbs.index(pdb)
    Y, S, P = pm["Y"][i].astype(int), pm["S"][i].astype(int), pm["pred_model"][i].astype(float)
    TRES, MRES = [str(x) for x in lab["tcr_resid"][i]], [str(x) for x in lab["mhc_resid"][i]]
    TAA, MAA = [str(x) for x in lab["tcr_aa"][i]], [str(x) for x in lab["mhc_aa"][i]]
    r = idx.loc[pdb]
    MHC, CH = str(r.mhc_chain1), {"a": str(r.Achain), "b": str(r.Bchain)}
    PEP = str(r.antigen_chain)

    # ---- predicted set: top-k resolved pairs, k = observed contacts ----------------------
    ti, mj = np.where(S == 1)
    k = int(Y[ti, mj].sum())
    order = np.argsort(-P[ti, mj], kind="stable")
    score_sorted = P[ti, mj][order]
    tie_at_k = bool(np.isclose(score_sorted[k - 1], score_sorted[k]))
    top = set(zip(ti[order[:k]].tolist(), mj[order[:k]].tolist()))
    obs = {(int(a), int(b)) for a, b in zip(*np.where(Y == 1))}
    assert len(top) == len(obs) == k
    rank = {pair: n + 1 for n, pair in enumerate(zip(ti[order].tolist(), mj[order].tolist()))}

    # ---- residue mapping from the p84v2 label build (number + insertion code) -------------
    model = PDBParser(QUIET=True).get_structure(pdb, f"data/structures/stcrdab/imgt/{pdb}.pdb")[0]
    mres = {rid(x): x for x in model[MHC] if x.id[0] == " "}
    tres = {t: {rid(x): x for x in model[CH[t]] if x.id[0] == " "} for t in "ab"}
    def tcr_res(s):
        seg, _ = TLAB[s].split(":"); t = seg[-1]
        assert TRES[s], (pdb, TLAB[s])
        return t, CH[t], TRES[s]
    def mhc_res(j):
        assert MRES[j], (pdb, POS[j])
        return MRES[j]

    sets = {}
    for name, pairs in (("observed", obs), ("predicted", top)):
        tcr_sites = {}; mhc_votes = {}
        for s, j in sorted(pairs):
            t, ch, tr = tcr_res(s); mr = mhc_res(j)
            tcr_sites[(t, ch, tr, s)] = True
            mhc_votes.setdefault((mr, j), {"a": 0, "b": 0})[t] += 1
            rows.append(dict(pdb=pdb, set=name, tcr_slot=TLAB[s], tcr_chain=ch, tcr_resi=tr,
                             mhc_imgt=POS[j], mhc_resi=mr, observed=int(Y[s, j]),
                             score=round(float(P[s, j]), 6), rank=rank[(s, j)]))
        ties = sorted(m for m, v in mhc_votes.items() if v["a"] == v["b"])
        both = sorted(m for m, v in mhc_votes.items() if v["a"] and v["b"])
        sets[name] = (tcr_sites, mhc_votes, ties, both)
    ov_pairs = len(obs & top)
    o_t, p_t = {x[:3] for x in sets["observed"][0]}, {x[:3] for x in sets["predicted"][0]}
    o_m, p_m = {x[0] for x in sets["observed"][1]}, {x[0] for x in sets["predicted"][1]}
    summary.append(dict(pdb=pdb, k=k, overlap_pairs=ov_pairs, tie_at_k=tie_at_k,
                        tcr_sites_obs=len(o_t), tcr_sites_pred=len(p_t), tcr_sites_shared=len(o_t & p_t),
                        mhc_sites_obs=len(o_m), mhc_sites_pred=len(p_m), mhc_sites_shared=len(o_m & p_m),
                        mhc_both_chain_obs=len(sets["observed"][3]), mhc_tie_obs=len(sets["observed"][2]),
                        mhc_both_chain_pred=len(sets["predicted"][3]), mhc_tie_pred=len(sets["predicted"][2])))

    # ---- structure frame from the platform -------------------------------------------------
    win = np.array([mres[p]["CA"].coord for p in POS if p in mres and "CA" in mres[p]])
    cen = win.mean(0)
    _, _, vt = np.linalg.svd(win - cen)
    e1, n = vt[0], vt[2]
    tcr_ca = np.array([res["CA"].coord for t in "ab" for res in tres[t].values()
                       if "CA" in res and res.id[1] <= 130])
    if np.dot(tcr_ca.mean(0) - cen, n) < 0: n = -n
    pep_ca = [res["CA"].coord for res in model[PEP] if res.id[0] == " " and "CA" in res]
    if np.dot(pep_ca[-1] - pep_ca[0], e1) < 0: e1 = -e1
    e1 = e1 - np.dot(e1, n) * n; e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    R = np.vstack([e1, e2, n])
    a1c = np.array([x["CA"].coord for x in mres.values() if 50 <= x.id[1] <= 90 and "CA" in x]).mean(0)
    a1_up = bool(np.dot(a1c - cen, e2) > 0)

    ORIENT[pdb] = "alpha1 helix " + ("top" if a1_up else "bottom") + "; peptide N-terminus left"
    for name in ("observed", "predicted"):
        tcr_sites, mhc_votes, _, _ = sets[name]
        cmd.reinitialize()
        cmd.set("max_threads", 4)
        cmd.load(f"data/structures/stcrdab/imgt/{pdb}.pdb", "cx")
        cmd.remove("solvent or hydrogens or hetatm")
        cmd.alter_state(1, "cx", "(x,y,z)=fr(x,y,z)",
                        space={"fr": lambda x, y, z: tuple(float(v) for v in R @ (np.array([x, y, z]) - cen))})
        cmd.create("mhc", f"cx and chain {MHC} and resi 1-90+1001-1090")
        cmd.create("pep", f"cx and chain {PEP}")
        loops = " or ".join(f"(chain {CH[t]} and resi {lo}-{hi})" for t in "ab" for lo, hi in SEG_RANGE.values())
        cmd.create("cdr", f"cx and ({loops})")
        cmd.delete("cx")
        cmd.bg_color("white")
        cmd.set("ray_opaque_background", 0); cmd.set("ray_shadows", 0)
        cmd.set("ambient", 0.45); cmd.set("specular", 0.08); cmd.set("depth_cue", 0)
        cmd.set("surface_quality", 1)
        cmd.hide("everything")
        cmd.show("surface", "mhc"); cmd.color(C_MHC, "mhc")
        cmd.show("sticks", "pep and not name N+C+O"); cmd.color(C_PEP, "pep")
        cmd.set("stick_radius", 0.22)
        cmd.show("cartoon", "cdr")
        cmd.set("cartoon_loop_radius", 0.38); cmd.set("cartoon_transparency", 0.25)
        for t in "ab":
            cmd.color(C_CH[t], f"cdr and chain {CH[t]}")
        for (mr, j), v in mhc_votes.items():
            sel = f"mhc and resi {mr}"            # number + insertion code, e.g. resi 1061A
            got = []
            cmd.iterate(f"{sel} and name CA", "got.append(resn)", space={"got": got})
            mres_man = MAN.get((pdb, MHC, mr), ("", ""))
            checks.append(dict(pdb=pdb, set=name, side="MHC", column=POS[j], selection=sel, n_ca=len(got),
                               resn=got[0] if got else "", token=MAA[j], manifest=mres_man[0],
                               ok=len(got) == 1 and got[0] == mres_man[0] and MAA[j] == mres_man[1]))
            t = "a" if v["a"] > v["b"] else "b" if v["b"] > v["a"] else None
            if t is not None:
                cmd.color(C_CH[t], sel)
        for (mr, j), v in mhc_votes.items():   # equal pair counts from both chains: split by atom
            if v["a"] == v["b"]:
                cmd.color(C_CH["a"], f"mhc and resi {mr} and name N+CA+C+O")
                cmd.color(C_CH["b"], f"mhc and resi {mr} and not name N+CA+C+O")
        for (t, ch, tr, s) in tcr_sites:
            sel = f"cdr and chain {ch} and resi {tr} and name CA"
            got = []
            cmd.iterate(sel, "got.append(resn)", space={"got": got})
            tres_man = MAN.get((pdb, ch, tr), ("", ""))
            ok = len(got) == 1 and got[0] == tres_man[0] and TAA[s] == tres_man[1]
            # a TCR site outside the drawn CDR ranges (FR4 118-124) is not in the 'cdr' object
            outside = cmd.count_atoms(sel) == 0 and not any(lo <= int(''.join(c for c in tr if c.isdigit())) <= hi
                                                            for lo, hi in SEG_RANGE.values())
            checks.append(dict(pdb=pdb, set=name, side="TCR", column=TLAB[s], selection=sel, n_ca=len(got),
                               resn=got[0] if got else "", token=TAA[s], manifest=tres_man[0],
                               ok=ok or outside, outside_drawn_loops=outside))
            if got:
                cmd.show("spheres", sel); cmd.color(C_CH[t], sel)
        cmd.set("sphere_scale", 0.72)
        cmd.set("sphere_transparency", 0.0)
        cmd.reset()
        cmd.orient("mhc and resi 50-90+1050-1090")
        # orient() rotates; restore the platform frame (identity rotation) and keep its zoom
        v = list(cmd.get_view()); v[0:9] = [1, 0, 0, 0, 1, 0, 0, 0, 1]
        v[12:15] = [0.0, 0.0, 0.0]
        cmd.set_view(v)
        cmd.zoom("mhc and resi 50-90+1050-1090", 3, complete=1)
        out = f"figures/pymol/p84v2_fig6_{pdb}_{name}.png"
        cmd.png(out, width=1400, height=1150, dpi=300, ray=1)
        print(f"wrote {out}")

chk = pd.DataFrame(checks)
chk.to_csv("results/review/figures_p84v2/fig6_selection_check.tsv", sep="\t", index=False)
assert chk.ok.all(), chk[~chk.ok]
pd.DataFrame(rows).to_csv(SITES_TSV, sep="\t", index=False)
print(f"wrote {SITES_TSV}")
n_out = int(chk["outside_drawn_loops"].eq(True).sum()) if "outside_drawn_loops" in chk else 0
inserted = sorted(set(x for x in chk.selection.str.extract(r"resi (\S+)")[0] if not x.isdigit()))
print(f"targets (max / min delta, percomplex_tcr70_paired_p84v2.tsv): {TARGETS}")
print(f"selection check: {int(chk.ok.sum())}/{len(chk)} displayed residues OK "
      f"(MHC {int((chk.side == 'MHC').sum())}, TCR {int((chk.side == 'TCR').sum())}; "
      f"TCR sites outside the drawn CDR ranges, not drawn: {n_out}); inserted residues displayed: {inserted}")
print(pd.DataFrame(summary).to_string(index=False))
print("orientation:", ORIENT)
