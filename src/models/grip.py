"""TCR-IMAP (T-cell Receptor Interface Mapping and Prediction): dual encoder with a
germline-selective structural pair bias."""
import numpy as np, torch, torch.nn as nn, torch.nn.functional as F

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i + 1 for i, a in enumerate(AA)}   # 0 = pad
VOCAB = len(AA) + 1

# TCR token layout: [CDR1a, CDR2a, CDR3a, CDR1b, CDR2b, CDR3b]
SEG_LENS = {"cdr1a": 12, "cdr2a": 10, "cdr3a": 20, "cdr1b": 12, "cdr2b": 10, "cdr3b": 20}
SEG_ORDER = ["cdr1a", "cdr2a", "cdr3a", "cdr1b", "cdr2b", "cdr3b"]
SEG_ID = {"cdr1a": 0, "cdr2a": 1, "cdr3a": 2, "cdr1b": 3, "cdr2b": 4, "cdr3b": 5}
TCR_LEN = sum(SEG_LENS.values())          # 84
GERMLINE_SEGS = {"cdr1a", "cdr2a", "cdr1b", "cdr2b"}


def encode_seq(s, L):
    v = np.zeros(L, dtype=np.int64)
    for i, c in enumerate(str(s)[:L]):
        v[i] = AA_IDX.get(c, 0)
    return v


def build_bias_matrix(prior_npz, mhc_len, mode="learned_scale", seed=0):
    """Map priors/contact_prior.npz onto the TCR token layout -> (TCR_LEN, mhc_len)."""
    d = np.load(prior_npz, allow_pickle=True)
    B, rows = d["B"], [str(x) for x in d["row_labels"]]
    out = np.zeros((TCR_LEN, mhc_len), dtype=np.float32)
    if mode == "zero":
        return out
    offs, o = {}, 0
    for s in SEG_ORDER:
        offs[s] = o; o += SEG_LENS[s]
    # row_labels look like "A_CDR1_27" -> chain A/B, loop, IMGT position
    for r, lab in enumerate(rows):
        ch, loop, pos = lab.split("_")
        seg = f"{loop.lower()}{'a' if ch == 'A' else 'b'}"
        if seg not in offs:
            continue
        start = 27 if loop == "CDR1" else 56 if loop == "CDR2" else 105
        k = int(pos) - start
        if 0 <= k < SEG_LENS[seg]:
            out[offs[seg] + k, :B.shape[1]] = B[r, :mhc_len]
    if mode == "permuted":                      # same values, scrambled geometry
        rng = np.random.default_rng(seed)
        flat = out.flatten(); rng.shuffle(flat)
        out = flat.reshape(out.shape)
    elif mode == "uniform":
        out = np.full_like(out, out.mean())
    return out


class CrossAttnBias(nn.Module):
    """Cross-attention with an additive pair bias, gated per head."""
    def __init__(self, d, heads, bias_mat, germline_only=True):
        super().__init__()
        self.h, self.dk = heads, d // heads
        self.q, self.k, self.v = (nn.Linear(d, d) for _ in range(3))
        self.o = nn.Linear(d, d)
        self.gamma = nn.Parameter(torch.full((heads,), 8.0))  # bias must be commensurate with logits
        mask = torch.zeros(TCR_LEN, 1)
        if germline_only:
            o = 0
            for s in SEG_ORDER:
                if s in GERMLINE_SEGS:
                    mask[o:o + SEG_LENS[s]] = 1.0
                o += SEG_LENS[s]
        else:
            mask[:] = 1.0
        self.register_buffer("bias", torch.tensor(bias_mat) * mask)

    def forward(self, x, ctx, ctx_mask=None):
        Bn, T, D = x.shape
        q = self.q(x).view(Bn, T, self.h, self.dk).transpose(1, 2)
        k = self.k(ctx).view(Bn, -1, self.h, self.dk).transpose(1, 2)
        v = self.v(ctx).view(Bn, -1, self.h, self.dk).transpose(1, 2)
        att = q @ k.transpose(-2, -1) / self.dk ** 0.5
        att = att + self.gamma.view(1, -1, 1, 1) * self.bias.unsqueeze(0).unsqueeze(0)
        if ctx_mask is not None:
            att = att.masked_fill(~ctx_mask[:, None, None, :], -1e4)
        out = (att.softmax(-1) @ v).transpose(1, 2).reshape(Bn, T, D)
        return self.o(out)


class TCRIMAP(nn.Module):
    def __init__(self, bias_mat, mhc_len, d=128, heads=8, layers=2, dropout=0.2,
                 germline_only=True, use_mhc=True, use_germline=True):
        super().__init__()
        self.use_mhc, self.use_germline = use_mhc, use_germline
        self.emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.seg = nn.Embedding(len(SEG_ORDER), d)
        self.pos = nn.Parameter(torch.randn(1, TCR_LEN, d) * .02)
        self.mhc_emb = nn.Embedding(VOCAB, d, padding_idx=0)
        self.mhc_pos = nn.Parameter(torch.randn(1, mhc_len, d) * .02)
        enc = nn.TransformerEncoderLayer(d, heads, d * 2, dropout, batch_first=True, norm_first=True)
        self.tcr_enc = nn.TransformerEncoder(enc, layers)
        self.mhc_enc = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(d, heads, d * 2, dropout, batch_first=True, norm_first=True), 1)
        self.cross = CrossAttnBias(d, heads, bias_mat, germline_only)
        self.norm = nn.LayerNorm(d)
        self.head = nn.Sequential(nn.Linear(d * 2, d), nn.GELU(), nn.Dropout(dropout), nn.Linear(d, 1))
        seg_ids, o = torch.zeros(TCR_LEN, dtype=torch.long), 0
        for s in SEG_ORDER:
            seg_ids[o:o + SEG_LENS[s]] = SEG_ID[s]; o += SEG_LENS[s]
        self.register_buffer("seg_ids", seg_ids)
        gm, o = torch.ones(TCR_LEN), 0
        for s in SEG_ORDER:
            if s in GERMLINE_SEGS and not use_germline:
                gm[o:o + SEG_LENS[s]] = 0
            o += SEG_LENS[s]
        self.register_buffer("germ_mask", gm)

    def forward(self, tcr, mhc):
        x = self.emb(tcr) + self.seg(self.seg_ids)[None] + self.pos
        x = x * self.germ_mask[None, :, None]
        x = self.tcr_enc(x, src_key_padding_mask=(tcr == 0) & ~self.germ_mask.bool()[None])
        m = self.mhc_enc(self.mhc_emb(mhc) + self.mhc_pos)
        if self.use_mhc:
            x = self.norm(x + self.cross(x, m))
        pooled = torch.cat([x.mean(1), m.mean(1)], -1)
        return self.head(pooled).squeeze(-1)


GRIP = TCRIMAP  # backward-compatible alias for historical artifacts/checkpoints
