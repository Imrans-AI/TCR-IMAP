// Build the Supporting Information document (figures S1-S9 + tables S1-S10) as a separate .docx.
// Figure captions come from docs/manuscript/FIGURE_CAPTIONS.md (verified by scripts/58_verify_captions.py);
// tables are converted from docs/manuscript/si_tables/*.md, which build_si_tables.py generates from the
// committed tables. No number is recomputed here.
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType, ImageRun, PageBreak,
  Table, TableRow, TableCell, WidthType, ShadingType, BorderStyle, PageOrientation,
  LineRuleType,
} = require("docx");

const REPO = "/home/muhammad-imran/PhD_data/TCRMHC_Module_3_TCR-IMAP";
const OUT = path.join(REPO, "docs/manuscript/TCR-IMAP_supporting_information_v3.docx");
const FONT = "Times New Roman";
const SZ = 22, SZ_SMALL = 20, SZ_TBL = 14; // 11 pt, 10 pt, 7 pt

function rt(text, base = {}) {
  const out = [];
  const re = /\*\*(.+?)\*\*|\*(.+?)\*|`(.+?)`/g;
  let last = 0, m;
  const push = (t, extra = {}) => { if (t !== "") out.push(new TextRun({ text: t, font: FONT, size: SZ, ...base, ...extra })); };
  while ((m = re.exec(text)) !== null) {
    push(text.slice(last, m.index));
    if (m[1] !== undefined) push(m[1], { bold: true });
    else if (m[2] !== undefined) push(m[2], { italics: true });
    else if (m[3] !== undefined) push(m[3], { font: "Consolas" });
    last = re.lastIndex;
  }
  push(text.slice(last));
  return out;
}
const P = (t, o = {}) => new Paragraph({
  children: rt(t), alignment: AlignmentType.JUSTIFIED,
  spacing: { after: 140, line: 300, lineRule: LineRuleType.AUTO }, ...o,
});
const H1 = (t) => new Paragraph({
  children: [new TextRun({ text: t, font: FONT, size: 26, bold: true })],
  spacing: { before: 320, after: 160 }, heading: HeadingLevel.HEADING_1,
});
const H2 = (t) => new Paragraph({
  children: [new TextRun({ text: t, font: FONT, size: 23, bold: true })],
  spacing: { before: 280, after: 120 }, heading: HeadingLevel.HEADING_2,
});

function pngSize(buf) { return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) }; }
function figure(file, label, caption, maxHeightIn) {
  const img = fs.readFileSync(file), dim = pngSize(img);
  const maxW = 6.5 * 96, maxH = maxHeightIn * 96;
  let w = maxW, h = (dim.h / dim.w) * maxW;
  if (h > maxH) { h = maxH; w = (dim.w / dim.h) * maxH; }
  return [
    new Paragraph({
      children: [new ImageRun({ data: img, type: "png", transformation: { width: Math.round(w), height: Math.round(h) } })],
      alignment: AlignmentType.CENTER, spacing: { before: 200, after: 80 },
    }),
    new Paragraph({
      children: [new TextRun({ text: label + ". ", font: FONT, size: SZ_SMALL, bold: true }), ...rt(caption, { size: SZ_SMALL })],
      alignment: AlignmentType.JUSTIFIED, spacing: { after: 240, line: 240, lineRule: LineRuleType.AUTO },
    }),
  ];
}

/* ---- markdown table -> docx table (landscape text width 9.0 in = 12960 dxa) ---- */
const TBL_W = 12960;
function mdTable(rows) {
  const n = Math.max(...rows.map((r) => r.length));
  const colW = Math.floor(TBL_W / n);
  const widths = Array(n).fill(colW);
  widths[n - 1] = TBL_W - colW * (n - 1);
  const mk = (cells, header) => new TableRow({
    tableHeader: !!header,
    children: Array.from({ length: n }, (_, i) => new TableCell({
      width: { size: widths[i], type: WidthType.DXA },
      shading: header ? { type: ShadingType.CLEAR, fill: "EEEEEE" } : undefined,
      margins: { top: 40, bottom: 40, left: 60, right: 60 },
      children: [new Paragraph({
        children: rt(cells[i] === undefined ? "" : cells[i], { size: SZ_TBL, bold: !!header }),
        spacing: { after: 0, line: 200, lineRule: LineRuleType.AUTO },
      })],
    })),
  });
  return new Table({
    columnWidths: widths,
    width: { size: TBL_W, type: WidthType.DXA },
    rows: [mk(rows[0], true), ...rows.slice(1).map((r) => mk(r, false))],
  });
}

function splitRow(line) {
  let t = line.trim();
  if (t.startsWith("|")) t = t.slice(1);
  if (t.endsWith("|")) t = t.slice(0, -1);
  return t.split("|").map((c) => c.trim());
}
const isSep = (l) => /^\s*\|?[\s:|-]+\|[\s:|-]*$/.test(l) && l.includes("-");

function convertTableDoc(file) {
  const lines = fs.readFileSync(file, "utf8").split("\n");
  const out = [];
  let buf = [];
  const flushTable = () => {
    if (buf.length) { out.push(mdTable(buf)); out.push(new Paragraph({ text: "", spacing: { after: 120 } })); buf = []; }
  };
  // Markdown wraps a paragraph across source lines; emit one Word paragraph per markdown
  // paragraph, not per line, so inline bold spanning a line break is not broken in two.
  let para = [];
  let bullet = [];
  const flushBullet = () => {
    if (bullet.length) {
      out.push(new Paragraph({
        children: rt(bullet.join(" "), { size: SZ_SMALL }),
        bullet: { level: 0 }, alignment: AlignmentType.JUSTIFIED,
        spacing: { after: 60, line: 240, lineRule: LineRuleType.AUTO },
      }));
      bullet = [];
    }
  };
  const flushPara = () => {
    flushBullet();
    if (para.length) {
      out.push(P(para.join(" "), { spacing: { after: 120, line: 260, lineRule: LineRuleType.AUTO } }));
      para = [];
    }
  };
  for (let i = 0; i < lines.length; i++) {
    const l = lines[i];
    if (l.trim().startsWith("|")) {
      flushPara();
      if (isSep(l)) continue;
      buf.push(splitRow(l));
      continue;
    }
    flushTable();
    const t = l.trim();
    if (t === "") { flushPara(); continue; }
    if (t.startsWith("# ")) { flushPara(); out.push(H2(t.slice(2))); continue; }
    if (t.startsWith("## ")) { flushPara(); out.push(H2(t.slice(3))); continue; }
    if (t.startsWith("- ") || t.startsWith("* ")) {
      flushPara();
      bullet.push(t.slice(2));
      continue;
    }
    // an indented continuation line belongs to the bullet or paragraph above it
    if (bullet.length && /^\s/.test(l)) { bullet.push(t); continue; }
    flushBullet();
    para.push(t);
  }
  flushPara();
  flushTable();
  return out;
}

/* ================================ front matter + figures (portrait) ========= */
const front = [];
front.push(new Paragraph({
  children: [new TextRun({ text: "Supporting Information", font: FONT, size: 30, bold: true })],
  alignment: AlignmentType.CENTER, spacing: { after: 160 },
}));
front.push(new Paragraph({
  children: [new TextRun({
    text: "TCR-IMAP: Homology-Controlled Evaluation of TCR–MHC Contact Prediction",
    font: FONT, size: 26, bold: true, italics: true,
  })],
  alignment: AlignmentType.CENTER, spacing: { after: 160, line: 300, lineRule: LineRuleType.AUTO },
}));
front.push(new Paragraph({
  children: [new TextRun({ text: "Muhammad Imran and Dong-Qing Wei", font: FONT, size: SZ })],
  alignment: AlignmentType.CENTER, spacing: { after: 60 },
}));
front.push(new Paragraph({
  children: [new TextRun({
    text: "State Key Laboratory of Microbial Metabolism, School of Life Sciences and Biotechnology, "
        + "Shanghai Jiao Tong University, Shanghai 200240, China",
    font: FONT, size: SZ_SMALL, italics: true })],
  alignment: AlignmentType.CENTER, spacing: { after: 300 },
}));

front.push(H1("Contents"));
[
  "Figure S1. Effect of adding the peptide sequence to TCR-IMAP.",
  "Figure S2. Per-complex per-residue AUPRC under each partitioning protocol.",
  "Figure S3. Interface characterization of the four simulated TCR–pMHC complexes.",
  "Figure S4. Sensitivity of the per-residue benchmark to the validation split.",
  "Figure S5. Sensitivity of the three primary tcr70 comparisons to cluster weighting.",
  "Figure S6. The comparator subset, per training seed.",
  "Figure S7. Legacy results on 28 MHC positions, documenting the candidate-set correction.",
  "Figure S8. Per-complex pairwise AUPRC grouped by MHC allele and by TRBV.",
  "Figure S9. Composition of the identity-threshold sweep.",
  "Table S1. Full benchmark: every method under every protocol.",
  "Table S2. Ablations of TCR-IMAP (tcr70, pairwise task).",
  "Table S3. Architectures, hyperparameters and parameter counts.",
  "Table S4. Partitioning protocols: fold composition.",
  "Table S5. Identity-threshold sweep: composition and results.",
  "Table S6. Molecular dynamics: system selection and per-replicate quality control.",
  "Table S7. Molecular dynamics: occupancy statistics.",
  "Table S8. Dataset descriptors.",
  "Table S9. Training-size-matched control: estimands, diagnostics and memberships.",
  "Table S10. Common-test-set and target-definition sensitivities.",
].forEach((t) => front.push(new Paragraph({
  children: rt(t, { size: SZ_SMALL }),
  spacing: { after: 40, line: 240, lineRule: LineRuleType.AUTO },
})));

front.push(new Paragraph({ children: [new PageBreak()] }));
front.push(H1("Supporting Figures"));

const siCaps = [
 ["Figure S1","Effect of adding the peptide sequence to TCR-IMAP: no detectable effect. Difference in pairwise macro-AUPRC under tcr70 between TCR-IMAP trained with and without the peptide sequence, with 95 % confidence intervals, in four scopes: all positions (+0.0100 [−0.0058, +0.0258], p = 0.196); all positions with the peptide sequence shuffled, a control (+0.0071 [−0.0072, +0.0222], p = 0.308); CDR1/2 rows only (+0.0107 [−0.0038, +0.0267], p = 0.133); and CDR3 rows only (+0.0082 [−0.0137, +0.0293], p = 0.416). Every 95 % CI includes zero (hollow markers). The macro-AUPRC without → with peptide is given beside each row. Statistics: paired cluster bootstrap over 69 TCR clusters, 2,000 resamples, recorded values. n = 174 complexes."],
 ["Figure S2","Per-complex per-residue AUPRC under each partitioning protocol. Distribution over test complexes of the per-residue AUPRC of TCR-IMAP, Struct-TCS, Seq-TCS, the within-fold position prior and nearest-neighbor (NN) transfer, under random (n = 174), double (n = 30), allele (n = 64), temporal (n = 60) and tcr70 (n = 174) partitioning. Boxes show the median and interquartile range (IQR), whiskers extend to 1.5 × IQR, and points are individual complexes. The means are those of Figure 4A. Per-complex AUPRC tracks a complex's number of contacts and is not a measure of difficulty; differences between methods are assessed with paired tests. Descriptive; no test."],
 ["Figure S3","Interface characterization of the four simulated TCR–pMHC complexes. Columns, systems (5yxn, 7n1f, 3vxs, 8eo8); colors, the three independent replicates. Every measure covers the analyzed window, 20–200 ns of each 200 ns production run, on trajectories made whole and free of periodic jumps. (A) Interface RMSD: RMSD of the Cα atoms of the canonical CDR loops after a least-squares fit on the MHC α1/α2 Cα atoms, relative to the crystal pose, every frame. (B) Buried solvent-accessible surface area of the TCR–pMHC interface, Shrake–Rupley with a 1.4 Å probe, sampled every 1 ns. (C) Per-loop Cα RMSF after the same fit. (D) Interface hydrogen bonds per frame, between the TCR and MHC α1/α2 and between the TCR and the peptide. 5yxn's disengagement appears jointly in these measures: interface RMSD exceeds 8 Å in more than a quarter of the frames of replicate 1 (maximum 10.7 Å), buried SASA falls in all three replicates, and CDR3β RMSF is elevated in replicates 1 and 2, where the β chain lost its contacts, but not in replicate 3, where the α chain disengaged. These measures are descriptive: with four systems they characterize these simulations only."],
 ["Figure S4","Sensitivity of the per-residue benchmark to the validation split used for model selection. (A) Per-residue macro-AUPRC of TCR-IMAP, Struct-TCS and Seq-TCS under each partitioning protocol, with early stopping on a random 15 % of the training complexes (random split validation; superseded) and on whole TCR clusters (cluster-disjoint validation; reported). The value above each protocol is the change in TCR-IMAP's macro-AUPRC from the random to the cluster-disjoint split: random −0.0541, double −0.0252, allele −0.0484, temporal −0.0122 and tcr70 −0.0165. The position prior and NN transfer use no validation subset and are identical under both splits. (B) Under tcr70, TCR-IMAP minus each comparator with 95 % confidence intervals under each split. The predeclared primary comparison, against Struct-TCS, is +0.0504 [+0.0188, +0.0828] (p = 0.001) under the random split and +0.0312 [−0.0049, +0.0650] (p = 0.088) under the cluster-disjoint split; against Seq-TCS, +0.0897 and +0.0799 (both p < 0.001); against the position prior, +0.0250 (p = 0.060) and +0.0085 (p = 0.546). (C) Number of TCR clusters with complexes on both sides of model selection, in each of the 18 folds under the random split: 5–15 (median 10.5). Under the cluster-disjoint split the number is 0 in every fold, by construction. Statistics: (B) paired cluster bootstrap over 69 TCR clusters, 20,000 resamples, seed 42."],
 ["Figure S5","Sensitivity of the three primary tcr70 comparisons to cluster weighting. Difference in macro-AUPRC between TCR-IMAP and the comparator with every test complex weighted equally (complex-macro, circles; the predeclared primary weighting) and with each of the 69 TCR clusters weighted equally (cluster-macro, squares). Per-residue vs Struct-TCS: +0.0312 [−0.0049, +0.0650], p = 0.088 and +0.0258 [−0.0016, +0.0518], p = 0.063. Per-residue vs the position prior: +0.0085 [−0.0162, +0.0371], p = 0.546 and +0.0164 [−0.0070, +0.0400], p = 0.174. Pairwise vs the pair-frequency prior: +0.0507 [+0.0215, +0.0841], p < 0.001 and +0.0360 [+0.0148, +0.0579], p = 0.001. The two weightings agree in every comparison. Filled markers, 95 % CI excludes zero. n = 174 complexes, 69 clusters."],
 ["Figure S6","The comparator subset, per training seed. The 88 structures shared with the comparator's dataset, at its 4.5 Å and 6 Å contact cutoffs, under a single random split of the kind the comparator used (n = 22 test complexes) and under tcr70 (n = 88). (A) TCR-IMAP − Struct-TCS per-residue macro-AUPRC. Circles, the reported three-seed ensemble with model selection on cluster-disjoint validation, with 95 % CI and p; crosses, the three training seeds individually; triangles, the same ensemble with model selection on a random validation subset (superseded). (B) Macro-AUPRC of each method. Under cluster-disjoint validation the 6 Å single-random-split cell reverses sign, from +0.0265 [+0.0008, +0.0506] to −0.0949 [−0.1335, −0.0517] (p < 0.001), and all three seeds agree on the sign (−0.060, −0.120 and −0.087). Under tcr70, TCR-IMAP leads at both cutoffs (+0.1057 and +0.0964, both p < 0.001). The single-random-split rows are reported for comparability with the published comparator setting, not as an evaluation we endorse: each rests on one fold of 22 test complexes, and only the structure subset and cutoff are shared."],
 ["Figure S7","Legacy results on 28 MHC positions, shown to document the correction. An earlier version of this analysis used 28 MHC positions selected from contact frequency pooled over all 174 complexes before the data were split into training and test folds. That is outcome-informed selection. It was replaced by the structurally defined 84-column α1/α2 window used throughout this paper, which does not depend on observed contacts. The legacy results (diamonds) are superseded and are shown only to document what the correction changed; the current values (circles) are those reported in the main text. Five headline comparisons, TCR-IMAP − comparator in macro-AUPRC: pairwise tcr70 vs the pair-frequency prior, +0.0398 (p = 0.010) → +0.0507 (p < 0.001); per-residue tcr70 vs Struct-TCS, +0.0273 (p = 0.001) → +0.0312 (p = 0.088); per-residue tcr70 vs the position prior, +0.0256 (p = 0.048) → +0.0085 (p = 0.546); identity sweep at < 0.70 vs the pair-frequency prior, +0.0127 (p = 0.429) → +0.0274 (p = 0.052); peptide input, CDR3 rows, +0.0380 (p < 0.001) → +0.0082 (p = 0.416). The per-residue legacy models were also selected on a random rather than a cluster-disjoint validation subset, so their change combines both corrections."],
 ["Figure S8","Per-complex pairwise AUPRC under tcr70, grouped by MHC allele and by TRBV. Points, individual complexes; lines, group means; TCR-IMAP (blue) and the within-fold pair-frequency prior (gray); n above each group. (A) By MHC allele. Allele annotation was parsed from the structure records and is missing for 108 of 174 complexes, shown as one separate group. (B) By TRBV subgroup. Groups with fewer than three complexes are omitted. This figure is descriptive: it is not a subgroup analysis and no test is performed. Groups are small and not independent of the TCR clusters used for evaluation, and per-complex AUPRC tracks a complex's number of contacts, so it is not a measure of difficulty."],
 ["Figure S9","Composition of the identity-threshold sweep. The sweep of Figure 3A,B removes, from each training fold, every complex whose TCR Levenshtein similarity ratio to any test complex of that fold reaches the threshold, on a fixed five-fold partition. (A) TCR-IMAP − pair-frequency prior, pairwise macro-AUPRC, with 95 % CI; filled markers, CI excludes zero; at 0.70, +0.0274, p = 0.052. (B) Complexes per fold. The test composition does not change: all 174 complexes are evaluated at every threshold, 34–35 per fold. The training composition does change, which is what the sweep manipulates: from 139–140 complexes per fold without restriction to 69–78 at 0.70. (C) The train/test similarity achieved: for each test complex, its maximum similarity to the training complexes its fold kept, summarized as the maximum and the mean over the 174 test complexes; dotted line, the threshold. At 0.70 the achieved maximum is 0.699 and the mean 0.649. The dash-dotted line is the mean_max_sim value recorded with the sweep, which is an upper bound on the achieved mean, not the achieved value."],
];
siCaps.forEach(([lab, cap], i) => {
  front.push(...figure(path.join(REPO, "figures/supplementary/Figure_S" + (i + 1) + ".png"), lab, cap, 6.4));
});

/* ================================ tables (landscape) ======================= */
const tablesBody = [H1("Supporting Tables")];
tablesBody.push(P("These tables are generated from the committed result tables by "
  + "`docs/manuscript/si_tables/build_si_tables.py`. Recomputed statistics use the pooled cluster bootstrap over the "
  + "69 TCR clusters (20,000 resamples, seed 42) on stored per-complex arrays and are asserted equal to the committed "
  + "values; the build stops on any mismatch. The corresponding CSV files accompany the deposited data.",
  { spacing: { after: 200, line: 260, lineRule: LineRuleType.AUTO } }));

["S1_benchmark", "S2_ablations", "S3_hyperparameters", "S4_protocol_folds",
 "S5_identity_sweep", "S6_md_selection_qc", "S7_md_occupancy", "S8_dataset_descriptors",
 "S9_size_matched_control", "S10_sensitivities"]
.forEach((name, i) => {
  if (i) tablesBody.push(new Paragraph({ children: [new PageBreak()] }));
  tablesBody.push(...convertTableDoc(path.join(REPO, "docs/manuscript/si_tables", name + ".md")));
});

const doc = new Document({
  creator: "Muhammad Imran",
  title: "Supporting Information - TCR-IMAP: Homology-Controlled Evaluation of TCR-MHC Contact Prediction",
  styles: { default: { document: { run: { font: FONT, size: SZ } } } },
  sections: [
    {
      properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 } } },
      children: front,
    },
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE },
          margin: { top: 1080, right: 1080, bottom: 1080, left: 1080 },
        },
      },
      children: tablesBody,
    },
  ],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(OUT, buf);
  console.log("written:", OUT, (buf.length / 1e6).toFixed(2), "MB");
});
