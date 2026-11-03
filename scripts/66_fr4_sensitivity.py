#!/usr/bin/env python3
"""Canonical-CDR-only scoring sensitivity: exclude the FR4 rows (IMGT 118-124) and rescore.

The p84v2 row set keeps 14 J-region framework rows (IMGT 118-124 of both chains). Those rows
carry 0 of the 5,281 contacts at 4.5 A but 204,624 of the 893,928 resolved pairs -- 22.9 % of the
evaluated negatives and none of the positives. Including them lowers the positive base rate and so
depresses average precision for every method. This asks whether any CONCLUSION depends on them.

Nothing is retrained: every number comes from prediction arrays already on disk. The artifact
discovery below was corrected on 2 Nov 2026: an earlier version of this script looked only at
results/tables/pairwise_ours_p84v2.npz (per-complex metrics, no predictions) and wrongly concluded
that the tcr70 pairwise predictions had never been stored. They are in
results/tables/pairmaps_{protocol}_p84v2.npz, which hold pred_model, pred_freq, Y, S, pdb, clu70
and a coverage mask for all five protocols. Per-residue predictions for the 88-complex comparator
subset are in results/tables/_partial/p84v2/31b_perseed/.

Coverage, and the one genuine gap:
  * PAIRWISE, all five frozen protocols -- TCR-IMAP and the pair-frequency prior. (pairmaps_*)
  * PAIRWISE, size-matched control and filtered sweep, all seven thresholds. (scripts 60, 62)
  * PAIRWISE TEIM-Res under tcr70. (teim_res_pairwise_p84v2.npz)
  * PER-RESIDUE on the 88-complex comparator subset under tcr70, per training seed, for TCR-IMAP,
    Struct-TCS and Seq-TCS. (31b_perseed)
  * NOT the predeclared per-residue primary over all 174 complexes: scripts/28b_benchmark_p82v1.py
    returns per-complex METRICS only and never writes a prediction array, so that one comparison
    cannot be rescored without retraining. Reported as outstanding, not estimated.

usage: python scripts/66_fr4_sensitivity.py
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, ".")
from src.eval.metrics import per_complex, cluster_bootstrap_diff

ROOT = Path(__file__).resolve().parent.parent
T = ROOT/"results/tables"
SM = ROOT/"results/review/size_matched_control"
FR = ROOT/"results/review/filtered_rerun"
PS = T/"_partial/p84v2/31b_perseed"
OUT = ROOT/"results/review/fr4_sensitivity"; OUT.mkdir(parents=True, exist_ok=True)
N_BOOT, SEED = 20000, 42
PROTOCOLS = ["random", "allele", "double", "temporal", "tcr70"]
PLAB = {"double": "intersection (legacy)"}


def main():
    z = np.load(ROOT/"data/processed/pair_maps_p84v2.npz", allow_pickle=True)
    num = z["row_imgt_number"].astype(int)
    keep = num <= 117
    md = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t")
    print(f"rows kept {int(keep.sum())}/{len(keep)}; dropped {int((~keep).sum())} FR4 rows "
          f"carrying {int(z['Y'][:, ~keep, :].sum())} contacts and "
          f"{int(z['S'][:, ~keep, :].sum()):,} of {int(z['S'].sum()):,} resolved pairs")

    rows = []

    def add(group, label, Y, S, pa, pb, clu, name_a, name_b):
        for scope, sel in (("all rows (reported)", slice(None)), ("canonical CDRs only", keep)):
            ma = per_complex(Y[:, sel, :], pa[:, sel, :], S[:, sel, :])
            mb = per_complex(Y[:, sel, :], pb[:, sel, :], S[:, sel, :])
            d, lo, hi, p = cluster_bootstrap_diff(ma["auprc"], mb["auprc"], clu,
                                                  n_boot=N_BOOT, seed=SEED)
            rows.append(dict(group=group, comparison=label, scope=scope, n=len(clu),
                             a=round(float(np.nanmean(ma["auprc"])), 4),
                             b=round(float(np.nanmean(mb["auprc"])), 4),
                             diff=round(float(d), 4), lo=round(float(lo), 4),
                             hi=round(float(hi), 4), p=round(float(p), 4)))

    # ---- 1. pairwise, every frozen protocol: TCR-IMAP vs the pair-frequency prior ----
    for proto in PROTOCOLS:
        pm = np.load(T/f"pairmaps_{proto}_p84v2.npz", allow_pickle=True)
        assert [str(x) for x in pm["pdb"]] == list(md.pdb), f"{proto}: identifier order differs"
        assert np.array_equal(pm["clu70"], md.clu70.values), f"{proto}: clu70 differs"
        assert np.array_equal(pm["Y"].astype(np.int8), z["Y"]), f"{proto}: labels differ"
        assert np.array_equal(pm["S"].astype(np.int8), z["S"]), f"{proto}: masks differ"
        idx = np.where(pm["covered"])[0]
        add(f"pairwise, {PLAB.get(proto, proto)}", "TCR-IMAP vs pair-frequency prior",
            pm["Y"].astype(np.float32)[idx], pm["S"].astype(np.float32)[idx],
            pm["pred_model"][idx], pm["pred_freq"][idx], pm["clu70"][idx],
            "TCR-IMAP", "prior")

    # ---- 2. pairwise tcr70: TCR-IMAP vs TEIM-Res ----
    pm = np.load(T/"pairmaps_tcr70_p84v2.npz", allow_pickle=True)
    te = np.load(T/"teim_res_pairwise_p84v2.npz", allow_pickle=True)
    order = {p: i for i, p in enumerate(md.pdb)}
    ti = [order[str(p)] for p in te["pdb"]]
    assert ti == list(range(len(md))), "TEIM-Res order differs from metadata"
    add("pairwise, tcr70", "TCR-IMAP vs TEIM-Res",
        pm["Y"].astype(np.float32), pm["S"].astype(np.float32),
        pm["pred_model"], te["pred"], pm["clu70"], "TCR-IMAP", "TEIM-Res")

    # ---- 3. per-residue on the 88-complex comparator subset, tcr70, seed-averaged ----
    ps = np.load(PS/"4p5A_our_tcr70.npz", allow_pickle=True)
    sub = [order[str(p)] for p in ps["pdb"]]
    assert np.array_equal(ps["clu70"], md.clu70.values[sub]), "perseed clu70 differs from metadata"
    Yr = (z["Y"].sum(-1) > 0).astype(np.float32)[sub]      # per-residue labels, as 28b derives them
    Sr = (z["S"].sum(-1) > 0).astype(np.float32)[sub]
    clu = np.asarray(ps["clu70"])
    def seedavg(m): return np.mean([ps[f"{m}_seed{s}"] for s in (42, 43, 44)], 0)
    for comp in ("Struct-TCS", "Seq-TCS"):
        for scope, sel in (("all rows (reported)", slice(None)), ("canonical CDRs only", keep)):
            ma = per_complex(Yr[:, sel], seedavg("Ours")[:, sel], Sr[:, sel])
            mb = per_complex(Yr[:, sel], seedavg(comp)[:, sel], Sr[:, sel])
            d, lo, hi, p = cluster_bootstrap_diff(ma["auprc"], mb["auprc"], clu,
                                                  n_boot=N_BOOT, seed=SEED)
            rows.append(dict(group="per-residue, tcr70, 88-complex comparator subset",
                             comparison=f"TCR-IMAP vs {comp}", scope=scope, n=len(clu),
                             a=round(float(np.nanmean(ma["auprc"])), 4),
                             b=round(float(np.nanmean(mb["auprc"])), 4),
                             diff=round(float(d), 4), lo=round(float(lo), 4),
                             hi=round(float(hi), 4), p=round(float(p), 4)))

    # ---- 4. the control regimes ----
    for thr in (0.95, 0.85, 0.70):
        for regime, f in (("size-matched random", SM/f"predictions_thr{thr}.npz"),
                          ("similarity-filtered", FR/f"filtered_predictions_thr{thr}.npz")):
            d = np.load(f, allow_pickle=True)
            add(f"sweep {thr}, {regime}", "TCR-IMAP vs pair-frequency prior",
                d["Y"], d["S"], d["pred_model"], d["pred_freq"], md.clu70.values,
                "TCR-IMAP", "prior")

    t = pd.DataFrame(rows)
    t.to_csv(OUT/"fr4_sensitivity.tsv", sep="\t", index=False)
    piv = t.pivot_table(index=["group", "comparison"], columns="scope", values="diff")
    piv["change"] = (piv["canonical CDRs only"] - piv["all rows (reported)"]).round(4)
    piv = piv.reset_index()
    piv.to_csv(OUT/"fr4_change_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 230); pd.set_option("display.max_rows", 200)
    print("\n=== full results ===")
    print(t.to_string(index=False))
    print("\n=== effect of excluding FR4 on each comparison ===")
    print(piv.to_string(index=False))
    print(f"\nlargest absolute change in any tested comparison: {piv['change'].abs().max():.4f}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
