#!/usr/bin/env python3
"""Comparator-fidelity sensitivity: Struct-TCS with and without the per-head gate.

Authorised in results/review/COMPARATOR_FIDELITY_RERUN_AUTHORISATION.md, committed before this ran.

Gao et al. describe a plain multi-head graph attention network and never mention gating; the
released GATEncoder instantiates the ungated Multihead_GAT. Our reimplementation adds a per-head
sigmoid gate taken from Multihead_GaAtN, a class present in the release but never instantiated.
That gate is our deviation, and the predeclared primary comparison is measured against the
comparator we built, so it is tested here.

VALIDATION GATE: the gated arm must reproduce the committed Struct-TCS macro-AUPRC 0.6969 and the
paired difference +0.0312 to 4 dp. On failure the ungated arm is NOT interpreted.

tcr70 only; Struct-TCS only; TCR-IMAP read from the stored per-complex arrays, not retrained.

usage: python scripts/68_comparator_gate_sensitivity.py
"""
import os
os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import json, sys, time
from pathlib import Path
import numpy as np, pandas as pd, torch, torch.nn as nn
from sklearn.metrics import average_precision_score

sys.path.insert(0, ".")
from src.data import p84v2 as P84
P84.install_grip_shim()
from src.models.grip import VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN   # noqa: F401
from src.eval.metrics import per_complex, cluster_bootstrap_diff

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT/"results/review/comparator_fidelity"; OUT.mkdir(parents=True, exist_ok=True)
EPOCHS, SEEDS = 300, (42, 43, 44)
N_BOOT, BSEED = 20000, 42
EXPECT_STRUCTCS, EXPECT_DIFF = 0.6969, 0.0312

# --- the comparator, exactly as scripts/21_struct_tcs.py defines it, with the gate switchable ---
class GAT(nn.Module):
    """Multi-head graph attention. gated=True adds the per-head sigmoid gate of our
    implementation; gated=False is the plain form the paper describes and the release uses."""
    def __init__(self, d, heads, gated):
        super().__init__()
        self.h, self.gated = heads, gated
        self.W = nn.ModuleList([nn.Linear(d, d, bias=False) for _ in range(heads)])
        self.a = nn.ModuleList([nn.Linear(2*d, 1, bias=False) for _ in range(heads)])
        if gated:
            self.gate = nn.Linear(d, heads)
        self.w0 = nn.Linear(d*heads, d, bias=False)
    def forward(self, h, adj):
        B, N, D = h.shape
        g = torch.sigmoid(self.gate(h)) if self.gated else None
        outs = []
        for k in range(self.h):
            wh = self.W[k](h)
            e = self.a[k](torch.cat([wh.unsqueeze(2).expand(B, N, N, D),
                                     wh.unsqueeze(1).expand(B, N, N, D)], -1)).squeeze(-1)
            e = torch.nn.functional.leaky_relu(e, 0.2).masked_fill(adj == 0, -1e9)
            o = torch.relu(e.softmax(-1) @ wh)
            outs.append(g[:, :, k:k+1]*o if self.gated else o)
        return self.w0(torch.cat(outs, -1))

class AttnBlock(nn.Module):
    def __init__(self, d, heads, dropout=0.1):
        super().__init__()
        self.att = nn.MultiheadAttention(d, heads, batch_first=True, dropout=dropout)
        self.n1, self.n2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.ff = nn.Sequential(nn.Linear(d, d*2), nn.ReLU(), nn.Linear(d*2, d))
    def forward(self, x):
        a, _ = self.att(x, x, x); x = self.n1(x + a)
        return self.n2(self.ff(x) + x)

class GATBlock(nn.Module):
    def __init__(self, d, heads, gated):
        super().__init__()
        self.gat = GAT(d, heads, gated)
        self.n1, self.n2 = nn.LayerNorm(d), nn.LayerNorm(d)
        self.ff = nn.Sequential(nn.Linear(d, d*2), nn.ReLU(), nn.Linear(d*2, d))
    def forward(self, x, adj):
        t = self.n1(self.gat(x, adj) + x)
        return self.n2(self.ff(t) + t)

class TCS(nn.Module):
    def __init__(self, d=64, pep_heads=4, prot_heads=4, peplayer=2, protlayer=2,
                 caheads=4, mhc_len=84, gated=True):
        super().__init__()
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.tlstm = nn.LSTM(d, d, 2, batch_first=True)
        self.mlstm = nn.LSTM(d, d, 2, batch_first=True)
        self.mhc_blocks = nn.ModuleList([AttnBlock(d, pep_heads) for _ in range(peplayer)])
        self.tcr_blocks = nn.ModuleList([GATBlock(d, prot_heads, gated) for _ in range(protlayer)])
        self.ca = nn.MultiheadAttention(d, caheads, batch_first=True)
        self.ffn = nn.Sequential(nn.Linear(d, d*2), nn.ReLU(), nn.Linear(d*2, d), nn.LayerNorm(d))
        self.final = nn.Linear(d, 1)
    def forward(self, tcr, mhc, adj):
        m = self.mlstm(self.memb(mhc))[0]
        for b in self.mhc_blocks: m = b(m)
        x = self.tlstm(self.emb(tcr))[0]
        for b in self.tcr_blocks: x = b(x, adj)
        inter, _ = self.ca(x, m, m)
        return self.final(self.ffn(inter + x)).squeeze(-1)


def cluster_disjoint_validation(clu, rest, rng, frac=0.15):
    cl = clu[rest]; uc = np.unique(cl); target = frac*len(rest); sel, n = [], 0
    for c in uc[rng.permutation(len(uc))]:
        if n >= target: break
        sel.append(c); n += int((cl == c).sum())
    return rest[np.isin(cl, sel)]


def fit(gated, T, M, A, Y, S, tr, va, te, seed, dev):
    torch.manual_seed(seed)
    net = TCS(mhc_len=M.shape[1], gated=gated).to(dev)
    g = lambda i: tuple(torch.tensor(z[i]).to(dev) for z in (T, M, A, Y, S))
    Xa, Ma, Aa, Ya, Sa = g(tr); Xv, Mv, Av, Yv, Sv = g(va); Xe, Me, Ae, _, _ = g(te)
    pw = ((Ya == 0) & (Sa > 0)).sum()/torch.clamp(((Ya == 1) & (Sa > 0)).sum(), min=1)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-4)
    lf = nn.BCEWithLogitsLoss(reduction="none", pos_weight=pw)
    best, bad, bs = -1, 0, None
    for ep in range(EPOCHS):
        net.train(); opt.zero_grad()
        l = (lf(net(Xa, Ma, Aa), Ya)*Sa).sum()/Sa.sum(); l.backward()
        nn.utils.clip_grad_norm_(net.parameters(), 1.0); opt.step()
        if ep % 5 == 0:
            net.eval()
            with torch.no_grad(): pv = torch.sigmoid(net(Xv, Mv, Av))
            k = Sv.bool()
            sc = average_precision_score(Yv[k].cpu(), pv[k].cpu())
            if sc > best: best, bad, bs = sc, 0, {a: b.clone() for a, b in net.state_dict().items()}
            else:
                bad += 1
                if bad > 10: break
    net.load_state_dict(bs); net.eval()
    with torch.no_grad(): return torch.sigmoid(net(Xe, Me, Ae)).cpu().numpy()


def main():
    t0 = time.time()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"device: {dev}", flush=True)
    d = P84.load()
    Y = (d["Y"].sum(-1) > 0).astype(np.float32)     # per-residue labels, as 28b derives them
    S = (d["S"].sum(-1) > 0).astype(np.float32)
    pdbs = d["pdb"]
    T, M = P84.T_for(pdbs), P84.M_for(pdbs)
    A = np.load(ROOT/"data/processed/tcr_graphs_p84v2.npz")["A"]
    md = pd.read_csv(ROOT/"data/processed/complex_metadata.tsv", sep="\t").set_index("pdb")
    meta = md.loc[pdbs]
    clu = meta.clu70.values
    folds = json.load(open(ROOT/"data/splits_v2/protocols.json"))["tcr70"]
    pdb2i = {p: i for i, p in enumerate(pdbs)}

    stored = np.load(ROOT/"results/tables/percomplex_tcr70_p84v2.npz", allow_pickle=True)
    assert [str(p) for p in stored["pdb"]] == list(pdbs), "stored per-complex order differs"
    ours = stored["Ours_auprc"]

    rows, per = [], {}
    for gated in (True, False):
        P = np.full(Y.shape, np.nan, np.float32); used = np.zeros(len(Y), bool)
        rng = np.random.default_rng(0)                   # 28b's discipline: one per protocol
        for fold in folds:
            te = np.array([pdb2i[p] for p in fold if p in pdb2i])
            if len(te) < 3: continue
            rest = np.setdiff1d(np.arange(len(Y)), te)
            va = cluster_disjoint_validation(clu, rest, rng, 0.15)
            tr = np.setdiff1d(rest, va)
            if len(va) < 3 or len(tr) < 20: continue
            assert not (set(clu[tr]) & set(clu[va]))
            used[te] = True
            P[te] = np.mean([fit(gated, T, M, A, Y, S, tr, va, te, s, dev) for s in SEEDS], 0)
        idx = np.where(used)[0]
        mm = per_complex(Y[idx], P[idx], S[idx])
        per[gated] = mm["auprc"]
        diff, lo, hi, p = cluster_bootstrap_diff(ours[idx], mm["auprc"], clu[idx],
                                                 n_boot=N_BOOT, seed=BSEED)
        rows.append(dict(variant="gated (as reported)" if gated else "ungated (as published)",
                         n=len(idx), struct_tcs=round(float(np.nanmean(mm["auprc"])), 4),
                         tcr_imap=round(float(np.nanmean(ours[idx])), 4),
                         diff=round(float(diff), 4), lo=round(float(lo), 4),
                         hi=round(float(hi), 4), p=round(float(p), 4),
                         significant_at_05="yes" if p < 0.05 else "no"))
        print(f"  {'gated  ' if gated else 'ungated'}: Struct-TCS={rows[-1]['struct_tcs']:.4f}  "
              f"diff={diff:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}  ({time.time()-t0:.0f}s)", flush=True)

    t = pd.DataFrame(rows); t.to_csv(OUT/"gate_sensitivity.tsv", sep="\t", index=False)
    g = t[t.variant.str.startswith("gated")].iloc[0]
    ok = (round(g.struct_tcs, 4) == EXPECT_STRUCTCS) and (round(g["diff"], 4) == EXPECT_DIFF)
    print("\n=== VALIDATION GATE ===")
    print(f"gated arm: Struct-TCS {g.struct_tcs:.4f} (expect {EXPECT_STRUCTCS}), "
          f"diff {g['diff']:+.4f} (expect {EXPECT_DIFF:+.4f})  ->  {'PASS' if ok else 'FAIL'}")
    print("\n" + t.to_string(index=False))
    # paired difference BETWEEN the two comparator variants, same complexes and resamples
    dd, dlo, dhi, dp = cluster_bootstrap_diff(per[False], per[True], clu, n_boot=N_BOOT, seed=BSEED)
    print(f"\nungated minus gated Struct-TCS, paired over the same complexes: "
          f"{dd:+.4f} [{dlo:+.4f}, {dhi:+.4f}] p = {dp:.4f}")
    pd.DataFrame([dict(contrast="ungated minus gated Struct-TCS", diff=round(float(dd), 4),
                       lo=round(float(dlo), 4), hi=round(float(dhi), 4), p=round(float(dp), 4))]
                 ).to_csv(OUT/"gate_paired_contrast.tsv", sep="\t", index=False)
    if not ok:
        (OUT/"GATE_FAILED").write_text("gated arm did not reproduce the committed values\n")
        print("\nGATE FAILED: the ungated arm is NOT interpreted.")
        sys.exit(1)
    print(f"\ntotal {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
