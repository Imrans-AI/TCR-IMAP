"""The MD occupancy analysis as two phases, so that occupancy is fixed before any score is loaded.

Phase A  occupancy_phase(): reads trajectories, computes every replicate's occupancy matrices (all
         window frames and docked frames) and docked-frame bookkeeping. Takes no scores.
Phase B  statistics_phase(): takes the phase-A result, the labels and a dict of score matrices
         {"TCR-IMAP": ..., "pair-frequency prior": ...} and computes every statistic fixed in the
         criteria document ("Analysis definitions - fixed 20 October 2026").
Used by run_occupancy_analysis.py (production, once) and tests/smoke_pipeline_npt.py (restrained NPT
with synthetic scores). Neither phase reads a score file; the caller passes scores in.
"""
import numpy as np
import pandas as pd

import stats as st
from occupancy import run_trajectory

PRIMARY = "TCR-IMAP"
CONTROL = "pair-frequency prior"


def occupancy_phase(systems, reps, xtc_of, residue_maps, chain_segments, t_start_ps, t_end_ps,
                    expected_window_frames, log=print, groups_fn=None):
    """Returns {(pdb, rep): run_trajectory result} after checking frame count and t = 0.

    `groups_fn` is forwarded to run_trajectory: None keeps its p82v1 84 x 82 default, and the p84v2
    analysis passes mdmap_p84v2.groups_from_map_p84v2 for the 75 x 84 axes."""
    out = {}
    kw = {} if groups_fn is None else {"groups_fn": groups_fn}
    kw["docked_identity"] = True   # decision 3b: both docked readings in the same pass
    for pdb in systems:
        for rep in reps:
            r = run_trajectory(pdb, xtc_of(pdb, rep), residue_maps[pdb], chain_segments,
                               t_start_ps, t_end_ps, docked=True, **kw)
            assert r["first_time"] == 0.0, f"{pdb} rep{rep}: first frame t = {r['first_time']} ps, expected 0"
            assert r["n_window"] == expected_window_frames, \
                f"{pdb} rep{rep}: {r['n_window']} window frames, expected {expected_window_frames}"
            out[(pdb, rep)] = r
            log(f"  {pdb} rep{rep}: {r['n_window']} frames, docked {r['n_docked']} "
                f"({100 * r['n_docked'] / r['n_window']:.1f} %), t=0 TCR-pMHC contacts {r['ref_count']}")
    return out


def occupancy_tables(occ, labels, tl, mp):
    """Phase-A tables (no scores): per-pair per-replicate occupancy, and docked-frame bookkeeping."""
    rows, dock = [], []
    for (pdb, rep), r in occ.items():
        Y, S = labels[pdb]
        ti, mj = np.where(S.astype(bool))
        od = r["occ_docked"]
        odi = r.get("occ_docked_identity")
        rows.append(pd.DataFrame(dict(pdb=pdb, rep=rep, tcr_slot=tl[ti], mhc_position=mp[mj],
                                      crystal_contact=Y[ti, mj].astype(int), occupancy=r["occ"][ti, mj],
                                      occupancy_docked_count=od[ti, mj] if od is not None else np.nan,
                                      occupancy_docked_identity=odi[ti, mj] if odi is not None else np.nan)))
        ndi = r.get("n_docked_identity", 0)
        pers = r.get("persisting", np.array([]))
        dock.append(dict(pdb=pdb, rep=rep, n_window_frames=r["n_window"],
                         # reading (a), count-based: the predeclared primary
                         n_docked_frames_count=r["n_docked"],
                         docked_fraction_count=r["n_docked"] / r["n_window"],
                         contributes_count=st.docked_contributes(r["n_docked"]),
                         # reading (b), identity-based: the same t = 0 pairs still present
                         n_docked_frames_identity=ndi,
                         docked_fraction_identity=ndi / r["n_window"],
                         contributes_identity=st.docked_contributes(ndi),
                         readings_agree=bool(st.docked_contributes(r["n_docked"]) == st.docked_contributes(ndi)),
                         ref_count_t0=r["ref_count"], n_ref_pairs_t0=r.get("n_ref_pairs"),
                         persisting_min=int(pers.min()) if len(pers) else np.nan,
                         persisting_median=float(np.median(pers)) if len(pers) else np.nan,
                         persisting_max=int(pers.max()) if len(pers) else np.nan,
                         count_min=int(r["counts"].min()), count_median=float(np.median(r["counts"])),
                         count_max=int(r["counts"].max())))
    return pd.concat(rows, ignore_index=True), pd.DataFrame(dock)


def statistics_phase(occ, labels, scores, systems, reps, n_perm=st.N_PERM, alternative="two-sided",
                     dists=None):
    """Returns dict of DataFrames: check1, rho, topk, docked_rho, strata.

    `dists` maps pdb -> crystal-frame minimum heavy-atom distance matrix, for the distance-stratified
    rho of decision 4. The permutation test is two-sided, as decided 24 Oct 2026."""
    check1, rho_rows, topk_rows, drho_rows, strata_rows = [], [], [], [], []
    sys_occ = {p: st.system_occupancy([occ[(p, r)]["occ"] for r in reps]) for p in systems}
    for p in systems:
        Y, S = labels[p]
        for r in reps:
            check1.append(dict(pdb=p, rep=r, **st.positive_occupancy_summary(occ[(p, r)]["occ"], Y, S)))
        check1.append(dict(pdb=p, rep="system", **st.positive_occupancy_summary(sys_occ[p], Y, S)))

    for name, sc in scores.items():
        pooled = []
        for p in systems:
            Y, S = labels[p]
            n_neg = int(st.neg_mask(Y, S).sum())
            neg = st.neg_mask(Y, S)
            r_sys, p_sys = st.perm_p(sc[p], sys_occ[p], Y, S, n_perm=n_perm, alternative=alternative)
            rho_rows.append(dict(score=name, level="system", pdb=p, rep="mean of 3", n_crystal_negative=n_neg,
                                 **st.occupancy_counts(sys_occ[p][neg]),
                                 rho=r_sys, perm_p_approximate=p_sys))
            for r in reps:
                rho_rows.append(dict(score=name, level="replicate", pdb=p, rep=r, n_crystal_negative=n_neg,
                                     **st.occupancy_counts(occ[(p, r)]["occ"][neg]),
                                     rho=st.rho(sc[p], occ[(p, r)]["occ"], Y, S), perm_p_approximate=np.nan))
            topk_rows.append(dict(score=name, pdb=p, **st.top_k_fraction(sc[p], sys_occ[p], Y, S)))
            pooled.append((sc[p], sys_occ[p], Y, S))
            if dists is not None:
                for row in st.rho_by_stratum(sc[p], sys_occ[p], Y, S, dists[p], n_perm=n_perm,
                                             alternative=alternative):
                    strata_rows.append(dict(score=name, level="system", pdb=p, **row))
        r_pool, p_pool = st.pooled_perm_p(pooled, n_perm=n_perm, alternative=alternative)
        rho_rows.append(dict(score=name, level="pooled", pdb="all four", rep="mean of 3",
                             n_crystal_negative=sum(int(st.neg_mask(y, s).sum()) for _, _, y, s in pooled),
                             **st.occupancy_counts(np.concatenate(
                                 [o[st.neg_mask(y, s)] for _, o, y, s in pooled])),
                             rho=r_pool, perm_p_approximate=p_pool))

        if dists is not None:
            for row in st.pooled_rho_by_stratum([(sc[p], sys_occ[p], *labels[p], dists[p]) for p in systems],
                                                n_perm=n_perm, alternative=alternative):
                strata_rows.append(dict(score=name, level="pooled", pdb="all four", **row))

        # docked-frame secondary: BOTH readings of decision 3, reported whether or not they agree
        for reading, okey, nkey in (("count", "occ_docked", "n_docked"),
                                    ("identity", "occ_docked_identity", "n_docked_identity")):
            dpooled = []
            for p in systems:
                Y, S = labels[p]
                dm, n_use = st.docked_system_occupancy([occ[(p, r)].get(okey) for r in reps],
                                                       [occ[(p, r)].get(nkey, 0) for r in reps])
                neg = st.neg_mask(Y, S)
                if dm is None:
                    # no docked occupancy at all, so no occupancy counts either
                    drho_rows.append(dict(score=name, reading=reading, level="system", pdb=p,
                                          n_replicates_contributing=0, n_crystal_negative=int(neg.sum()),
                                          **st.NO_OCCUPANCY_COUNTS, rho=np.nan, perm_p_approximate=np.nan))
                    continue
                r_d, p_d = st.perm_p(sc[p], dm, Y, S, n_perm=n_perm, alternative=alternative)
                drho_rows.append(dict(score=name, reading=reading, level="system", pdb=p,
                                      n_replicates_contributing=n_use, n_crystal_negative=int(neg.sum()),
                                      **st.occupancy_counts(dm[neg]), rho=r_d, perm_p_approximate=p_d))
                dpooled.append((sc[p], dm, Y, S))
            if dpooled:
                r_dp, p_dp = st.pooled_perm_p(dpooled, n_perm=n_perm, alternative=alternative)
            else:
                r_dp = p_dp = np.nan
            drho_rows.append(dict(score=name, reading=reading, level="pooled",
                                  pdb=f"{len(dpooled)} contributing systems",
                                  n_replicates_contributing=np.nan,
                                  n_crystal_negative=sum(int(st.neg_mask(y, s).sum()) for _, _, y, s in dpooled),
                                  **(st.occupancy_counts(np.concatenate(
                                      [o[st.neg_mask(y, s)] for _, o, y, s in dpooled])) if dpooled
                                     else st.NO_OCCUPANCY_COUNTS),
                                  rho=r_dp, perm_p_approximate=p_dp))

    # the description rule's flag on every correlation table (criteria document, 24 October 2026); no other column changes
    return dict(check1=pd.DataFrame(check1), rho=st.add_interpretable(pd.DataFrame(rho_rows)),
                topk=pd.DataFrame(topk_rows), docked_rho=st.add_interpretable(pd.DataFrame(drho_rows)),
                strata=st.add_interpretable(pd.DataFrame(strata_rows)) if strata_rows else pd.DataFrame(strata_rows))


def characterisation_phase(systems, reps, xtc_of, residue_maps, chain_segments, t_start_ps, t_end_ps,
                           root, sasa_stride=50, hbond_stride=1, log=print, hbond_xtc_of=None):
    """Score-independent simulation-quality and interface measures (criteria document, 24 Oct 2026).

    Reads no score file and takes no score argument, so it cannot depend on a prediction.
    `xtc_of` must give a whole, jump-free trajectory (pbc.whole_nojump_xtc) for the coordinate measures;
    `hbond_xtc_of` (default: xtc_of) gives the trajectory for the periodic-aware hydrogen-bond count.
    Returns (per-replicate summary DataFrame, time-series DataFrame, RMSF DataFrame, H-bond pair DataFrame)."""
    import characterise as ch
    summ, series, rmsf_rows, hb_rows = [], [], [], []
    for pdb in systems:
        for rep in reps:
            c = ch.characterise_trajectory(pdb, xtc_of(pdb, rep), residue_maps[pdb], chain_segments,
                                           t_start_ps, t_end_ps, root, sasa_stride=sasa_stride, log=log)
            hb = ch.hbond_counts(pdb, (hbond_xtc_of or xtc_of)(pdb, rep), residue_maps[pdb], chain_segments,
                                 t_start_ps, t_end_ps, stride=hbond_stride)
            row = dict(pdb=pdb, rep=rep, n_window_frames=c["n_window"],
                       backbone_rmsd_mean=float(c["backbone_rmsd"].mean()),
                       backbone_rmsd_sd=float(c["backbone_rmsd"].std()),
                       backbone_rmsd_final=float(c["backbone_rmsd"][-1]),
                       interface_rmsd_mean=float(c["interface_rmsd"].mean()),
                       interface_rmsd_sd=float(c["interface_rmsd"].std()),
                       interface_rmsd_final=float(c["interface_rmsd"][-1]),
                       buried_sasa_mean=float(c["buried_sasa"].mean()),
                       buried_sasa_sd=float(c["buried_sasa"].std()),
                       buried_sasa_n_frames=len(c["buried_sasa"]), sasa_stride=c["sasa_stride"])
            for k, v in c["per_loop_rmsf"].items():
                row[f"rmsf_{k}"] = v
            for tag, h in hb.items():
                row[f"hbond_{tag}_mean"] = h["mean_per_frame"]
                row[f"hbond_{tag}_sd"] = h["sd_per_frame"]
                for pair, occ_ in h["pair_occupancy"].items():
                    hb_rows.append(dict(pdb=pdb, rep=rep, interface=tag, residue_pair=pair, occupancy=occ_))
            summ.append(row)
            series.append(pd.DataFrame(dict(pdb=pdb, rep=rep, time_ps=c["times"],
                                            backbone_rmsd=c["backbone_rmsd"],
                                            interface_rmsd=c["interface_rmsd"])))
            series.append(pd.DataFrame(dict(pdb=pdb, rep=rep, time_ps=c["sasa_times"],
                                            backbone_rmsd=np.nan, interface_rmsd=np.nan,
                                            buried_sasa=c["buried_sasa"])))
            rmsf_rows.append(pd.DataFrame(dict(pdb=pdb, rep=rep, residue=c["rmsf_labels"],
                                               loop=c["rmsf_loop"], rmsf=c["rmsf_per_residue"])))
    return (pd.DataFrame(summ), pd.concat(series, ignore_index=True),
            pd.concat(rmsf_rows, ignore_index=True),
            pd.DataFrame(hb_rows) if hb_rows else pd.DataFrame(columns=["pdb", "rep", "interface", "residue_pair", "occupancy"]))
