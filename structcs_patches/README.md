# Struct-TCS / Seq-TCS — patches against the released implementation

Four unified diffs against the released code of Gao *et al.*, *J. Chem. Inf. Model.* **2025**,
*65*, 4674–4686 (DOI 10.1021/acs.jcim.5c00055), repository `gaol00034/pan-pMHC` at commit
`c118be5`, as held in `data/raw/sgtpmi/`.

```
0001-execution-import-GAT-not-GATs.patch          Models/modules.py      execution fix
0002-execution-GATEncoder-argument-arity.patch    Models/modules.py      execution fix
0003-execution-optimizer-never-steps.patch        Scripts/site_train.py  execution fix
0004-model-selection-validation-not-test.patch    Scripts/site_train.py  model-selection change
```

Apply in order from the root of a checkout of that repository:

```bash
for p in 000*.patch; do patch -p1 < "$p"; done
```

All four apply cleanly to `c118be5` and the result parses. Each patch carries a header naming the
defect it addresses and whether it is an execution fix or a change of method.

---

## Read this before using them

**These patches are not the code that produced any number in the paper.** Our comparator is a
reimplementation — `scripts/21_struct_tcs.py`, and `scripts/53_structcs_setting_p84v2.py` for the
comparator-subset setting — written from the published description. We never ran a patched copy of
the released code, because the release cannot be run end to end on our data at all: it ships no
feature builder, no graph builder, no dataset, no splits and no trained weights
(`results/review/COMPARATOR_PROVENANCE.md`).

What the patches are for is narrower and worth having anyway:

* they make the defects **concrete and checkable** rather than a claim in prose — anyone can apply
  them to `c118be5` and see what changed;
* they show **what a minimal repair would look like**, which is the evidence for the statement that
  the released code does not run as published;
* patch 0004 documents the one change that is **method, not repair** — separating model selection
  from test evaluation — which is the single difference most likely to move a number.

**A faithful reimplementation is not constructible** from the paper plus the release. The paper does
not give the graph radius numerically; the release gives no feature builder, no subset memberships
and no splits. These patches do not close that gap and do not claim to. The manuscript states
plainly that this is an architecture comparison under common inputs, not a reproduction of published
performance, and that no number of ours is comparable to any number in that paper.

---

## The gated and ungated variants

The manuscript reports a sensitivity between two forms of the graph attention. They are **not two
separate code files**: the difference is a **configuration flag**.

| | where | what it is |
|---|---|---|
| in the **release** | `Models/GAT.py` | two classes exist: `Multihead_GAT` (ungated, line 63) and `Multihead_GaAtN` (gated, line 34). `GATEncoder` instantiates the **ungated** one; the gated class is never used anywhere in the release |
| in the **paper** | §2.2.3.4 | describes a plain multi-head graph attention. **Gating is never mentioned.** The released ungated code therefore agrees with the paper |
| in **our** code | `scripts/21_struct_tcs.py`, class `GatedGAT` | adds a per-head sigmoid gate, `self.gate = nn.Linear(d, heads)`. **This is our own deviation** from both the paper and the release, taken from the unused `Multihead_GaAtN` |
| the **switch** | `scripts/68_comparator_gate_sensitivity.py`, `GAT(d, heads, gated)` | one boolean. `gated=True` is the reported variant; `gated=False` is the plain form the paper describes and the release uses |

Running `python scripts/68_comparator_gate_sensitivity.py` produces both arms. It carries a
validation gate: the gated arm must reproduce the committed Struct-TCS macro-AUPRC of 0.6969 and the
paired difference of +0.0312 to four decimal places, and if it does not, the ungated arm is not
interpreted. Results are in `results/review/comparator_fidelity/`:

* gated, as reported: **+0.0312**, p = 0.088
* ungated, as published: **+0.0354**, p = 0.036
* the two variants against each other: **−0.0042**, p = 0.441 — indistinguishable

The paper reports the **gated** variant, which is the one named in Methods and in the
predeclaration, and treats the flip in significance as evidence that the comparison is fragile at
this sample size. It is never used to pick a verdict.

---

## What these patches do not fix

Found while writing them, and deliberately left alone, because fixing them would mean rewriting the
module rather than repairing it:

* `Scripts/site_train.py` `fit()` references `os` and `modelname`; neither is imported or defined in
  that module, so saving a checkpoint raises `NameError`. The caller that was meant to define them
  is not in the release.
* `Scripts/site_train.py` has no `__main__`, no data loader and no argument parsing, so there is no
  entry point to run even after patching.
* The release contains no feature builder, no graph builder, no dataset, no splits and no weights.

Patch 0003 fixes the two `NameError`s that sit **inside** `train_epoch`, because the optimizer fix
cannot be exercised without them. The rest are listed here rather than patched.

---

## Provenance

* Defects and deviations: `results/review/COMPARATOR_PROVENANCE.md`; `HANDOFF_1.md` §10a.
* Validation-split change: `results/review/VALIDATION_FIX_REPORT.md`; Methods 2.7.
* Upstream licence and redistribution: the released repository is **not redistributed** in this
  package. These are diffs against it, not copies of it. See `THIRD_PARTY_NOTICES.md`.
