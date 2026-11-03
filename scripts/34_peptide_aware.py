#!/usr/bin/env python3
"""
Does peptide information improve CDR-MHC contact prediction, and is any gain
specific to CDR3?

Variants (tcr70, per-complex macro-AUPRC, cluster bootstrap):
  no_peptide     - current model
  peptide        - + peptide sequence branch, cross-attended
  peptide_shuf   - peptide residues shuffled within each complex (control)
Reported overall and split into germline (CDR1/2) vs CDR3 rows.
"""
import json, sys
import numpy as np, pandas as pd, torch, torch.nn as nn
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
sys.path.insert(0,'.')
from src.models.grip import encode_seq, VOCAB, SEG_ORDER, SEG_LENS, TCR_LEN, GERMLINE_SEGS
from src.eval.metrics import per_complex, cluster_bootstrap_diff

PEP_LEN = 15
SEEDS = (42,43,44); EPOCHS = 400

def get_peptides(meta):
    idx = pd.read_csv("data/interim/stcrdab_index.tsv", sep="\t").drop_duplicates("pdb").set_index("pdb")
    P = PDBParser(QUIET=True); out=[]
    for pdb in meta.pdb:
        r = idx.loc[pdb]
        try:
            m = P.get_structure(pdb, f"data/structures/stcrdab/imgt/{pdb}.pdb")[0]
            ch = m[str(r.antigen_chain)]
            s = "".join(seq1(x.get_resname()) for x in ch if x.id[0]==" ")
        except Exception:
            s = ""
        out.append(s[:PEP_LEN])
    return out

class PepNet(nn.Module):
    def __init__(self, mhc_len, use_pep=True, d=128, heads=8, dropout=0.2):
        super().__init__()
        self.use_pep = use_pep
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.seg = nn.Embedding(len(SEG_ORDER), d)
        self.pos = nn.Parameter(torch.randn(1, TCR_LEN, d)*.02)
        self.memb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.mpos = nn.Parameter(torch.randn(1, mhc_len, d)*.02)
        self.pemb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.ppos = nn.Parameter(torch.randn(1, PEP_LEN, d)*.02)
        self.tenc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), 2)
        self.menc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d*2, dropout, batch_first=True, norm_first=True), 1)
        self.cross_m = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.cross_p = nn.MultiheadAttention(d, heads, dropout=dropout, batch_first=True)
        self.norm = nn.LayerNorm(d)
        self.wt = nn.Linear(d, d, bias=False); self.wm = nn.Linear(d, d, bias=False)
        self.pair_pos = nn.Parameter(torch.zeros(TCR_LEN, mhc_len))
        self.scale = nn.Parameter(torch.tensor(1.0)); self.bias = nn.Parameter(torch.tensor(-4.0))
        sid,o = torch.zeros(TCR_LEN,dtype=torch.long),0
        for i,s in enumerate(SEG_ORDER):
            sid[o:o+SEG_LENS[s]]=i; o+=SEG_LENS[s]
        self.register_buffer("sid", sid)
    def forward(self, t, m, p):
        x = self.tenc(self.emb(t) + self.seg(self.sid)[None] + self.pos)
        mm = self.menc(self.memb(m) + self.mpos)
        a,_ = self.cross_m(x, mm, mm); x = self.norm(x + a)
        if self.use_pep:
            pp = self.pemb(p) + self.ppos
            b,_ = self.cross_p(x, pp, pp); x = self.norm(x + b)
        logits = torch.einsum('bid,bjd->bij', self.wt(x), self.wm(mm))/x.shape[-1]**0.5
        return logits*self.scale + self.pair_pos[None] + self.bias

def run(tag, use_pep, shuffle_pep, Y,S,T,M,Pp,meta,proto,dev,rng):
    pdb2i = {p:i for i,p in enumerate(meta.pdb)}
    Pq = Pp.copy()
    if shuffle_pep:
        for i in range(len(Pq)):
            nz = Pq[i][Pq[i]>0]; rng.shuffle(nz); Pq[i][Pq[i]>0] = nz
    out = np.full(Y.shape, np.nan, np.float32); used=np.zeros(len(Y),bool)
    for fold in proto:
        te = np.array([pdb2i[p] for p in fold if p in pdb2i])
        if len(te)<3: continue
        tr = np.setdiff1d(np.arange(len(Y)), te); used[te]=True
        ps=[]
        for sd in SEEDS:
            torch.manual_seed(sd)
            net = PepNet(M.shape[1], use_pep=use_pep).to(dev)
            Xa,Ma,Pa,Ya,Sa = (torch.tensor(z[tr]).to(dev) for z in (T,M,Pq,Y,S))
            pw = ((Ya==0)&(Sa>0)).sum()/torch.clamp(((Ya==1)&(Sa>0)).sum(),min=1)
            opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=0.01)
            lf = nn.BCEWithLogitsLoss(reduction='none', pos_weight=pw)
            for ep in range(EPOCHS):
                net.train(); opt.zero_grad()
                l=(lf(net(Xa,Ma,Pa),Ya)*Sa).sum()/Sa.sum(); l.backward()
                nn.utils.clip_grad_norm_(net.parameters(),1.0); opt.step()
            net.eval()
            with torch.no_grad():
                ps.append(torch.sigmoid(net(torch.tensor(T[te]).to(dev),
                                            torch.tensor(M[te]).to(dev),
                                            torch.tensor(Pq[te]).to(dev))).cpu().numpy())
        out[te]=np.mean(ps,0)
    return out, np.where(used)[0]

if __name__=="__main__":
    dev='cuda' if torch.cuda.is_available() else 'cpu'
    Z=np.load("data/processed/pair_maps.npz")
    Y,S = Z["Y"].astype(np.float32), Z["S"].astype(np.float32)
    meta=pd.read_parquet("data/processed/pair_meta.parquet")
    md=pd.read_csv("data/processed/complex_metadata.tsv",sep="\t").set_index("pdb")
    meta=meta.join(md[["clu70"]],on="pdb")
    peps = get_peptides(meta)
    print("peptide lengths:", pd.Series([len(x) for x in peps]).describe()[["mean","min","max"]].round(1).to_dict())
    T=np.stack([np.concatenate([encode_seq(r[k],SEG_LENS[k]) for k in SEG_ORDER]) for _,r in meta.iterrows()])
    M=np.stack([encode_seq(s,Y.shape[2]) for s in meta.mhc_contact])
    Pp=np.stack([encode_seq(s,PEP_LEN) for s in peps])
    proto=json.load(open("data/splits_v2/protocols.json"))["tcr70"]
    rng=np.random.default_rng(0)

    o,rows = {}, []
    for tag,(up,sh) in {"no_peptide":(False,False),"peptide":(True,False),"peptide_shuf":(True,True)}.items():
        pred, idx = run(tag,up,sh,Y,S,T,M,Pp,meta,proto,dev,rng)
        o[tag]=(pred,idx)
        mm=per_complex(Y[idx],pred[idx],S[idx])
        print(f"  {tag:13s} macro-AUPRC={np.nanmean(mm['auprc']):.4f}")
    idx=o["peptide"][1]; cl=meta.clu70.values[idx]
    base=per_complex(Y[idx],o["no_peptide"][0][idx],S[idx])
    for tag in ("peptide","peptide_shuf"):
        mm=per_complex(Y[idx],o[tag][0][idx],S[idx])
        d,lo,hi,p=cluster_bootstrap_diff(mm["auprc"],base["auprc"],cl)
        print(f"  {tag} - no_peptide: {d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}")
        rows.append(dict(scope=tag, no_pep=round(float(np.nanmean(base['auprc'])),4),
                         pep=round(float(np.nanmean(mm['auprc'])),4),
                         diff=round(d,4), lo=round(lo,4), hi=round(hi,4), pval=round(p,4)))
    # germline vs CDR3 breakdown
    off, gm, cm = 0, [], []
    for s in SEG_ORDER:
        r=list(range(off,off+SEG_LENS[s])); off+=SEG_LENS[s]
        (gm if s in GERMLINE_SEGS else cm).extend(r)
    print("\nby region:")
    for name, rowsel in (("germline(CDR1/2)",gm), ("CDR3",cm)):
        b=per_complex(Y[idx][:,rowsel],o["no_peptide"][0][idx][:,rowsel],S[idx][:,rowsel])
        pk=per_complex(Y[idx][:,rowsel],o["peptide"][0][idx][:,rowsel],S[idx][:,rowsel])
        d,lo,hi,p=cluster_bootstrap_diff(pk["auprc"],b["auprc"],cl)
        print(f"  {name:18s} no_pep={np.nanmean(b['auprc']):.4f} pep={np.nanmean(pk['auprc']):.4f} "
              f"Δ={d:+.4f} [{lo:+.4f},{hi:+.4f}] p={p:.4f}")
        rows.append(dict(scope=name, no_pep=round(float(np.nanmean(b['auprc'])),4),
                         pep=round(float(np.nanmean(pk['auprc'])),4),
                         diff=round(d,4), lo=round(lo,4), hi=round(hi,4), pval=round(p,4)))
    pd.DataFrame(rows).to_csv("results/tables/peptide_aware.csv", index=False)
    np.savez_compressed("results/tables/peptide_preds.npz",
                        pdb=meta.pdb.values, idx=idx,
                        no_peptide=o["no_peptide"][0], peptide=o["peptide"][0],
                        peptide_shuf=o["peptide_shuf"][0])
    print("\nwrote results/tables/peptide_aware.csv + peptide_preds.npz")
