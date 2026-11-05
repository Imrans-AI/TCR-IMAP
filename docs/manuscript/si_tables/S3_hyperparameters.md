# Table S3. Architectures, hyperparameters and parameter counts

**Table S3.** Architecture and training settings of TCR-IMAP (pairwise and per-residue variants), the two reimplemented
per-residue comparators (Struct-TCS, Seq-TCS), the reimplemented pairwise comparator (TEIM-Res) and the docking-angle
regressor (AngleNet), as used for every reported (p84v2) result.

| Item | TCR-IMAP (pairwise) | TCR-IMAP (per-residue) | Struct-TCS | Seq-TCS | TEIM-Res | AngleNet |
|---|---|---|---|---|---|---|
| Task | pairwise contact map (75 × 84) | per-residue (75 rows) | per-residue | per-residue | pairwise contact map | docking-angle regression |
| Source class (script) | `PairNet` (19) | `OursPerRes` (21) | `TCS(mode='struct')` (21) | `TCS(mode='seq')` (21) | `TEIMRes` (36b) | `AngleNet` (16b) |
| p84v2 runner | `53_run_all_p84v2.py` | `53_benchmark_p84v2.py` | `53_benchmark_p84v2.py` | `53_benchmark_p84v2.py` | `53_teim_res_p84v2.py` | `53_docking_angle_p84v2.py` |
| Inputs | TCR (75) + MHC (84) tokens | TCR + MHC tokens | TCR + MHC tokens + TCR residue graph (Cβ–Cβ ≤ 8 Å; Cα for Gly) | TCR + MHC tokens | TCR + MHC tokens | TCR + MHC tokens; FR4 rows (IMGT 118–124) padded |
| Width d | 128 | 96 | 64 | 64 | 64 (embedding); 32 (2D channels) | 96 |
| TCR encoder | 2-layer pre-norm Transformer, FF 256 | 2-layer pre-norm Transformer, FF 192 | 2-layer LSTM → 2 gated multi-head GAT blocks | 2-layer LSTM → 2 self-attention blocks | 2 × Conv1d (k = 3) + BatchNorm + ReLU | 2-layer pre-norm Transformer, FF 192 |
| MHC encoder | 1-layer pre-norm Transformer, FF 256 | 1-layer pre-norm Transformer, FF 192 | 2-layer LSTM → 2 self-attention blocks | 2-layer LSTM → 2 self-attention blocks | 2 × Conv1d (k = 3) + BatchNorm + ReLU | embedding + positional only |
| Attention heads | 8 (encoders and cross-attention) | 6 | 4 (GAT, self- and cross-attention) | 4 | n/a | 6 |
| Fusion | TCR→MHC cross-attention, residual + LayerNorm | same as pairwise | cross-attention (TCR query) + FFN | same as Struct-TCS | element-wise outer product → 1 × 1 conv → 3 dilated 2D residual blocks (dilation 1, 2, 4) | TCR→MHC cross-attention, residual + LayerNorm |
| Output head | bilinear s·⟨W_t h_i, W_m g_j⟩/√d + P_ij + b (P 75 × 84, init 0; s init 1; b init −4) | per row: Linear–GELU–Dropout–Linear | per row: Linear | per row: Linear | 1 × 1 conv | pooled [TCR mean, MHC mean]: Linear–GELU–Dropout–Linear |
| Dropout | 0.2 | 0.25 | 0.1 (self-attention blocks) | 0.1 | none | 0.3 |
| Loss | BCE, positive-class weight = neg/pos (resolved training pairs), masked to resolved pairs | same, per row | same, per row | same, per row | same as pairwise | smooth L1 on the standardized angle |
| Optimizer | AdamW | Adam | Adam | Adam | AdamW | AdamW |
| Learning rate | 1 × 10⁻³ | 1 × 10⁻³ | 1 × 10⁻³ | 1 × 10⁻³ | 1 × 10⁻³ | 1 × 10⁻³ |
| Weight decay | 0.01 | 1 × 10⁻⁴ | 1 × 10⁻⁴ | 1 × 10⁻⁴ | 0.01 | 0.05 |
| Gradient clipping (norm) | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| Batch | full batch | full batch | full batch | full batch | full batch | full batch |
| Epochs | 400 (fixed) | ≤ 300 | ≤ 300 | ≤ 300 | 400 (fixed) | 400 (fixed) |
| Validation / early stopping | none | cluster-disjoint validation (whole clu70 clusters, ≥ 15 % of training complexes); AUPRC every 5 epochs; stop after 10 evaluations without improvement; best weights restored | same as per-residue TCR-IMAP | same | none | none |
| Seeds | 42, 43, 44; predictions averaged | 42, 43, 44; averaged | 42, 43, 44; averaged | 42, 43, 44; averaged | 42, 43, 44; averaged | 42, 43, 44; each a separate run |
| Determinism | reproducible without extra settings | reproducible (deterministic settings change no value) | same | same | deterministic cuDNN/cuBLAS settings required (see Methods) | n/a |
| Trainable parameters | 529,310 | 291,073 | 336,905 | 303,105 | 110,593 | 225,505 |

**Notes.**
- Settings are read from the model and training code; the p84v2 runners execute the same sources with the p84v2 input
  layout (TCR axis 75 rows, segment lengths 7/8/22/6/7/25; MHC axis 84 columns) and change no hyperparameter
  (`results/review/REBUILD_P84V2_REPORT.md` step 3; `CORRECTION_PLAN_P84V2.md` §2).
- Parameter counts were taken by `build_si_tables.py` by instantiating each class on the CPU at the p84v2 layout (no
  forward pass) and counting trainable parameters. The counts for TCR-IMAP (both variants), Struct-TCS and Seq-TCS are
  asserted equal to the values recorded during the build; the TEIM-Res and
  AngleNet counts are recorded here for the first time.
- The pairwise models (TCR-IMAP, TEIM-Res) use no validation subset and no early stopping, so no test complex
  influences their training. The per-residue validation split is described in Methods and evaluated in Figure S4.
- Struct-TCS receives the crystal-structure residue graph of each complex, including test complexes; the other models
  receive amino-acid identities only.
- TEIM-Res: the published model's pretrained peptide autoencoder does not apply to TCR–MHC contacts; all encoders are
  trained from scratch.

**Sources.** `scripts/19_train_pairmap.py` (`PairNet`); `scripts/24b_run_all_p82v1.py` (`fit_predict`: optimizer, loss,
epochs, seeds); `scripts/21_struct_tcs.py` (`OursPerRes`, `TCS`, `GatedGAT`, `AttnBlock`, `GATBlock`);
`scripts/28b_benchmark_p82v1.py` (`fit`, `EPOCHS`, `SEEDS`); `scripts/53_benchmark_p84v2.py`
(`cluster_disjoint_validation`); `scripts/20_build_graphs.py` (graph, `CUT = 8.0`); `scripts/36b_teim_res_p82v1.py`
(`TEIMRes`, training loop); `scripts/53_teim_res_p84v2.py` and `results/review/TEIM_DETERMINISM_REPORT.md`
(determinism); `scripts/16b_docking_angle_p82v1.py` (`AngleNet`, `run_one`); `results/review/VALIDATION_FIX_REPORT.md` §3, §5.
