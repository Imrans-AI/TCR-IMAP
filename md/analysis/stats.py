"""Statistics of the MD occupancy analysis, exactly as fixed in md/MD_SYSTEM_SELECTION_CRITERIA.md
("Analysis definitions - fixed 20 October 2026"). Pure functions on arrays; no file I/O, no MD
reading. Every function is tested on synthetic inputs with known answers in tests/test_stats.py.

Conventions
  * a system's matrices are 84 x 82 (TCR slot x MHC position); S = resolved (pair universe),
    Y = crystal contact; crystal-negative = S & ~Y, crystal-positive = S & Y;
  * replicate occupancy: fraction of the 9,001 window frames in contact; system occupancy = mean of
    its replicate occupancies;
  * scores: TCR-IMAP tcr70 out-of-fold pairwise scores (pred_model) or, as the control, the
    within-fold pair-frequency prior (pred_freq), same functions.
"""
import numpy as np
from scipy.stats import spearmanr

MEANINGFUL = 0.25         # meaningful occupancy, fixed
N_PERM = 10_000           # within-system permutations
PERM_SEED = 42
MIN_DOCKED_FRAMES = 1000  # >= 20 ns at 20 ps


def neg_mask(Y, S):
    return S.astype(bool) & ~Y.astype(bool)


def pos_mask(Y, S):
    return S.astype(bool) & Y.astype(bool)


def system_occupancy(rep_occs):
    """Mean of replicate occupancy matrices (list of 84 x 82)."""
    return np.mean(np.stack(rep_occs), axis=0)


def rho(score, occ, Y, S):
    """Spearman rho between score and occupancy among crystal-negative pairs (ties: average ranks)."""
    m = neg_mask(Y, S)
    if undefined_reason(score[m], occ[m]):
        return np.nan
    return float(spearmanr(score[m], occ[m]).correlation)


def pooled_rho(systems):
    """systems: list of (score, occ, Y, S). Spearman over the union of all systems' crystal-negative pairs."""
    s = np.concatenate([sc[neg_mask(y, sm)] for sc, oc, y, sm in systems])
    o = np.concatenate([oc[neg_mask(y, sm)] for sc, oc, y, sm in systems])
    return float(spearmanr(s, o).correlation)


def _rank(x):
    from scipy.stats import rankdata
    return rankdata(x)


def n_distinct(x):
    """Number of distinct values in x (reported beside every correlation, so an undefined one explains itself)."""
    return int(len(np.unique(np.asarray(x))))


def occupancy_counts(occ_vals):
    """The occupancy counts reported beside every correlation (criteria document, 24 October 2026), over exactly
    the values the correlation uses: distinct values, pairs with occupancy > 0, pairs with occupancy >= 0.25."""
    v = np.asarray(occ_vals, dtype=float)
    return dict(n_distinct_occupancy=n_distinct(v), n_nonzero_occupancy=int((v > 0).sum()),
                n_occupancy_ge_025=int((v >= MEANINGFUL).sum()))


NO_OCCUPANCY_COUNTS = dict(n_distinct_occupancy=np.nan, n_nonzero_occupancy=np.nan, n_occupancy_ge_025=np.nan)

# Description rule (criteria document, '"Very few" fixed at 10', 24 October 2026): a correlation resting on fewer
# than 10 non-zero occupancies is described as uninterpretable. It is still reported with rho, p and all counts.
MIN_NONZERO_INTERPRETABLE = 10


def add_interpretable(df):
    """Append the boolean `interpretable` = rho defined and n_nonzero_occupancy >= 10. Changes no other column;
    False wherever there is no rho to interpret (undefined, untested, no docked occupancy)."""
    out = df.copy()
    nz = out["n_nonzero_occupancy"].astype(float)
    out["interpretable"] = out["rho"].notna() & (nz >= MIN_NONZERO_INTERPRETABLE)
    return out


def undefined_reason(score_vals, occ_vals):
    """Why a rank correlation is undefined for this group, or None if it is defined (criteria document,
    "Undefined correlations are reported as undefined", 24 October 2026): fewer than two pairs, all
    occupancy values identical (including all zero), or all scores identical."""
    if len(occ_vals) < 2:
        return "n < 2"
    if n_distinct(occ_vals) < 2:
        return "occupancy constant"
    if n_distinct(score_vals) < 2:
        return "scores constant"
    return None


def _rho_from_ranks(rs, ro):
    rs = rs - rs.mean(); ro = ro - ro.mean()
    return float((rs * ro).sum() / np.sqrt((rs ** 2).sum() * (ro ** 2).sum()))


def perm_p(score, occ, Y, S, n_perm=N_PERM, seed=PERM_SEED, alternative="two-sided"):
    """Within-system permutation p for rho: scores permuted among the system's crystal-negative
    pairs. APPROXIMATE - pairs are not independent. p = (1 + #extreme) / (1 + n_perm)."""
    m = neg_mask(Y, S)
    if undefined_reason(score[m], occ[m]):
        return np.nan, np.nan                      # undefined rho: no test, no p
    rs, ro = _rank(score[m]), _rank(occ[m])
    obs = _rho_from_ranks(rs, ro)
    rng = np.random.default_rng(seed)
    null = np.array([_rho_from_ranks(rng.permutation(rs), ro) for _ in range(n_perm)])
    return obs, _p(obs, null, alternative)


def pooled_perm_p(systems, n_perm=N_PERM, seed=PERM_SEED, alternative="two-sided"):
    """Pooled rho with scores permuted WITHIN each system (system structure preserved). Permuting
    within systems leaves the pooled multiset of scores unchanged, so the pooled score ranks are
    computed once and permuted block by block."""
    parts = [(sc[neg_mask(y, sm)], oc[neg_mask(y, sm)]) for sc, oc, y, sm in systems]
    if undefined_reason(np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])):
        return np.nan, np.nan                      # undefined rho: no test, no p
    rs = _rank(np.concatenate([p[0] for p in parts]))
    ro = _rank(np.concatenate([p[1] for p in parts]))
    obs = _rho_from_ranks(rs, ro)
    edges = np.cumsum([0] + [len(p[0]) for p in parts])
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for b in range(n_perm):
        null[b] = _rho_from_ranks(np.concatenate([rng.permutation(rs[a:z]) for a, z in zip(edges[:-1], edges[1:])]), ro)
    return obs, _p(obs, null, alternative)


def _p(obs, null, alternative):
    if not np.isfinite(obs):
        return np.nan                              # an undefined rho has no p (never the floor)
    if alternative == "two-sided":
        k = (np.abs(null) >= abs(obs) - 1e-12).sum()
    elif alternative == "greater":
        k = (null >= obs - 1e-12).sum()
    else:
        raise ValueError(alternative)
    return float((1 + k) / (1 + len(null)))


def top_k_negatives(score, Y, S):
    """The k highest-scoring crystal-negative pairs, k = number of crystal-positive pairs (Figure 6
    rule: stable descending sort in row-major pair order). Returns (index tuple arrays, k, tie_at_k)."""
    ti, mj = np.where(neg_mask(Y, S))
    k = int(pos_mask(Y, S).sum())
    sc = score[ti, mj]
    order = np.argsort(-sc, kind="stable")
    ss = sc[order]
    tie = bool(len(ss) > k and np.isclose(ss[k - 1], ss[k]))
    sel = order[:k]
    return (ti[sel], mj[sel]), k, tie


def top_k_fraction(score, occ, Y, S, thr=MEANINGFUL):
    """Fraction of the top-k crystal-negative pairs with occupancy >= thr, and the base rate among
    all crystal-negative pairs."""
    sel, k, tie = top_k_negatives(score, Y, S)
    frac = float((occ[sel] >= thr).mean())
    base = float((occ[neg_mask(Y, S)] >= thr).mean())
    return dict(k=k, frac_top_k=frac, n_top_k_ge=int((occ[sel] >= thr).sum()), base_rate=base,
                tie_at_k=tie)


def positive_occupancy_summary(occ, Y, S):
    """Check 1: crystal-positive occupancy of one replicate (or system)."""
    v = occ[pos_mask(Y, S)]
    return dict(n=len(v), mean=float(v.mean()), median=float(np.median(v)), min=float(v.min()),
                n_ge_meaningful=int((v >= MEANINGFUL).sum()))


def docked_contributes(n_docked):
    return n_docked >= MIN_DOCKED_FRAMES


def docked_system_occupancy(rep_docked_occs, rep_n_docked):
    """Docked-frame secondary: mean over the replicates that contribute (>= 1,000 docked frames);
    None if no replicate contributes."""
    use = [o for o, n in zip(rep_docked_occs, rep_n_docked) if docked_contributes(n) and o is not None]
    return (np.mean(np.stack(use), axis=0), len(use)) if use else (None, 0)


# --------------------------------------------------------------------------------------------------
# Distance-stratified rho (criteria document, decision 4, 24 October 2026).
# A crystal-negative pair already close to the interface can show high occupancy through proximity
# alone, so a positive pooled rho could reflect geometry rather than transient contact. Strata are of
# the CRYSTAL-FRAME minimum heavy-atom distance and are fixed; they are not adjusted after the fact.
DIST_STRATA = (("4.5-6 A", 4.5, 6.0), ("6-8 A", 6.0, 8.0), ("> 8 A", 8.0, np.inf))
MIN_STRATUM_N = 20        # fewer crystal-negative pairs than this: report n, no rho and no p


def stratum_masks(Y, S, dist):
    """{stratum name: boolean mask} over crystal-negative resolved pairs, by crystal-frame distance."""
    m = neg_mask(Y, S)
    out = {}
    for name, lo, hi in DIST_STRATA:
        with np.errstate(invalid="ignore"):
            out[name] = m & (dist >= lo) & (dist < hi)
    return out


def rho_by_stratum(score, occ, Y, S, dist, n_perm=N_PERM, seed=PERM_SEED, alternative="two-sided"):
    """rho and its approximate permutation p within each distance stratum."""
    rows = []
    for name, sel in stratum_masks(Y, S, dist).items():
        n = int(sel.sum())
        nd = occupancy_counts(occ[sel])
        if n < MIN_STRATUM_N:
            rows.append(dict(stratum=name, n=n, **nd, rho=np.nan, perm_p_approximate=np.nan,
                             note=f"n < {MIN_STRATUM_N}, not tested"))
            continue
        why = undefined_reason(score[sel], occ[sel])
        if why:
            rows.append(dict(stratum=name, n=n, **nd, rho=np.nan, perm_p_approximate=np.nan,
                             note=f"rho undefined: {why}"))
            continue
        rs, ro = _rank(score[sel]), _rank(occ[sel])
        obs = _rho_from_ranks(rs, ro)
        rng = np.random.default_rng(seed)
        null = np.array([_rho_from_ranks(rng.permutation(rs), ro) for _ in range(n_perm)])
        rows.append(dict(stratum=name, n=n, **nd, rho=obs,
                         perm_p_approximate=_p(obs, null, alternative), note=""))
    return rows


def pooled_rho_by_stratum(systems, n_perm=N_PERM, seed=PERM_SEED, alternative="two-sided"):
    """systems: list of (score, occ, Y, S, dist). Pooled within stratum, permuted within system."""
    rows = []
    for name, lo, hi in DIST_STRATA:
        parts = []
        for sc, oc, y, sm, d in systems:
            with np.errstate(invalid="ignore"):
                sel = neg_mask(y, sm) & (d >= lo) & (d < hi)
            if sel.sum():
                parts.append((sc[sel], oc[sel]))
        n = sum(len(p[0]) for p in parts)
        nd = occupancy_counts(np.concatenate([p[1] for p in parts]) if parts else [])
        if n < MIN_STRATUM_N or not parts:
            rows.append(dict(stratum=name, n=n, **nd, rho=np.nan, perm_p_approximate=np.nan,
                             note=f"n < {MIN_STRATUM_N}, not tested"))
            continue
        why = undefined_reason(np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts]))
        if why:
            rows.append(dict(stratum=name, n=n, **nd, rho=np.nan, perm_p_approximate=np.nan,
                             note=f"rho undefined: {why}"))
            continue
        rs = _rank(np.concatenate([p[0] for p in parts]))
        ro = _rank(np.concatenate([p[1] for p in parts]))
        obs = _rho_from_ranks(rs, ro)
        edges = np.cumsum([0] + [len(p[0]) for p in parts])
        rng = np.random.default_rng(seed)
        null = np.empty(n_perm)
        for b in range(n_perm):
            null[b] = _rho_from_ranks(
                np.concatenate([rng.permutation(rs[a:z]) for a, z in zip(edges[:-1], edges[1:])]), ro)
        rows.append(dict(stratum=name, n=n, **nd, rho=obs,
                         perm_p_approximate=_p(obs, null, alternative), note=""))
    return rows
