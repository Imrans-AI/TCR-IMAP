// Build v2 of the TCR-IMAP manuscript (main text). Supporting Information is a separate file,
// built by build_si_docx_v2.js. Numbers come from RESULTS_DRAFT.md, METHODS_DRAFT.md and
// FIGURE_CAPTIONS.md; none is recomputed.
// Numbers are taken from RESULTS_DRAFT.md, METHODS_DRAFT.md and FIGURE_CAPTIONS.md; none is recomputed.
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType, ImageRun,
  Bookmark, InternalHyperlink, PageBreak, TabStopType, TabStopPosition,
  Math: OMath, MathRun, MathFraction, MathSubScript, MathSuperScript, MathSubSuperScript,
  MathRadical, MathSum, MathRoundBrackets, MathSquareBrackets, MathAngledBrackets,
  convertInchesToTwip, LineRuleType,
} = require("docx");

const REPO = "/home/muhammad-imran/PhD_data/TCRMHC_Module_3_TCR-IMAP";
const OUT = path.join(REPO, "docs/manuscript/TCR-IMAP_manuscript_v3.docx");

/* ------------------------------------------------------------------ */
/* Figure 1 is the author-approved asset, not a script output.          */
/* It is installed at figures/main/Figure_1.png (6,480 x 4,230 px,      */
/* 7.2 x 4.7 in at 900 dpi) and documented in figures/main/             */
/* Figure_1.SOURCE.md. Current revision: Arial typography, 2 Nov 2026.  */
/* The editable source lives in                                         */
/* figures/source/fig1_arial_20261102/. scripts/70_fig1_workflow.py     */
/* is superseded and can no longer write here. The hash is verified on  */
/* every build so that a replaced or regenerated file cannot be shipped */
/* silently: the build aborts instead.                                  */
/* ------------------------------------------------------------------ */
const FIG1 = path.join(REPO, "figures/main/Figure_1.png");
const FIG1_SHA256 = "0684c8647a0b158c86af9229f2344f1eec4f8ae71b08e83c12cda0a26364d9a1";
{
  const got = require("crypto").createHash("sha256").update(fs.readFileSync(FIG1)).digest("hex");
  if (got !== FIG1_SHA256) {
    throw new Error(
      "figures/main/Figure_1.png is not the author-approved figure.\n" +
      "  expected sha256 " + FIG1_SHA256 + "\n" +
      "  found    sha256 " + got + "\n" +
      "Restore it from figures/source/fig1_arial_20261102/ (see figures/main/Figure_1.SOURCE.md) " +
      "before building. Do not update FIG1_SHA256 unless the author has approved a new figure."
    );
  }
  const b = fs.readFileSync(FIG1);
  const w = b.readUInt32BE(16), h = b.readUInt32BE(20);
  if (w !== 6480 || h !== 4230) throw new Error(`Figure 1 is ${w}x${h}, expected 6480x4230`);
}

const FONT = "Times New Roman";
const SZ = 22;        // 11 pt (half-points)
const SZ_SMALL = 20;  // 10 pt for captions
const LINE = { line: 360, lineRule: LineRuleType.AUTO }; // 1.5 lines

/* ------------------------------------------------------------------ */
/* rich-text mini-parser:  **bold**  *italic*  _{sub}  ^{sup}  [[1,2]] */
/* ------------------------------------------------------------------ */
function rt(text, base = {}) {
  const out = [];
  const re = /\*\*(.+?)\*\*|\*(.+?)\*|_\{(.+?)\}|\^\{(.+?)\}|\[\[([\d,\s–-]+)\]\]/g;
  let last = 0, m;
  const push = (t, extra = {}) => {
    if (t === "") return;
    out.push(new TextRun({ text: t, font: FONT, size: SZ, ...base, ...extra }));
  };
  while ((m = re.exec(text)) !== null) {
    push(text.slice(last, m.index));
    if (m[1] !== undefined) push(m[1], { bold: true });
    else if (m[2] !== undefined) push(m[2], { italics: true });
    else if (m[3] !== undefined) push(m[3], { subScript: true });
    else if (m[4] !== undefined) push(m[4], { superScript: true });
    else if (m[5] !== undefined) {
      const ids = m[5].split(",").map((s) => s.trim()).filter(Boolean);
      out.push(new TextRun({ text: "[", font: FONT, size: SZ, ...base }));
      ids.forEach((id, k) => {
        if (k) out.push(new TextRun({ text: ",", font: FONT, size: SZ, ...base }));
        out.push(new InternalHyperlink({
          anchor: "ref" + id.replace(/[^\d]/g, "") ,
          children: [new TextRun({ text: id, font: FONT, size: SZ, ...base, style: "Hyperlink" })],
        }));
      });
      out.push(new TextRun({ text: "]", font: FONT, size: SZ, ...base }));
    }
    last = re.lastIndex;
  }
  push(text.slice(last));
  return out;
}

const P = (text, opts = {}) => new Paragraph({
  children: rt(text),
  spacing: { after: 160, ...LINE },
  alignment: AlignmentType.JUSTIFIED,
  ...opts,
});

const H1 = (t, opts = {}) => new Paragraph({
  children: [new TextRun({ text: t, font: FONT, size: 26, bold: true })],
  spacing: { before: 320, after: 160 },
  heading: HeadingLevel.HEADING_1,
  ...opts,
});

const H2 = (t) => new Paragraph({
  children: [new TextRun({ text: t, font: FONT, size: 23, bold: true, italics: true })],
  spacing: { before: 240, after: 120 },
  heading: HeadingLevel.HEADING_2,
});

/* equation paragraph: centred maths + right-aligned number */
let eqN = 0;
const EQ = (children) => {
  eqN += 1;
  return new Paragraph({
    children: [
      new TextRun({ text: "\t", font: FONT, size: SZ }),
      new OMath({ children }),
      new TextRun({ text: "\t(" + eqN + ")", font: FONT, size: SZ }),
    ],
    tabStops: [
      { type: TabStopType.CENTER, position: TabStopPosition.MAX / 2 },
      { type: TabStopType.RIGHT, position: TabStopPosition.MAX },
    ],
    spacing: { before: 120, after: 160, ...LINE },
  });
};
const mr = (t) => new MathRun(t);
const msub = (b, s) => new MathSubScript({ children: [mr(b)], subScript: [mr(s)] });
const msup = (b, s) => new MathSuperScript({ children: [mr(b)], superScript: [mr(s)] });

/* figures */
function figure(file, label, caption, maxHeightIn = 7.1, exactAspect = false) {
  const img = fs.readFileSync(file);
  const dim = pngSize(img);
  const maxW = 6.5 * 96;
  const maxH = maxHeightIn * 96;
  let w = maxW, h = (dim.h / dim.w) * maxW;
  if (h > maxH) { h = maxH; w = (dim.w / dim.h) * maxH; }
  // Rounding both sides to whole 1/96 in changes the aspect ratio by up to ~0.1 %. Every figure
  // built from a script has been shipped that way since v1 and is left alone; the author-approved
  // Figure 1 is embedded with exactAspect so its proportions are preserved exactly.
  if (!exactAspect) { w = Math.round(w); h = Math.round(h); }
  return [
    new Paragraph({
      children: [new ImageRun({ data: img, type: "png", transformation: { width: w, height: h } })],
      alignment: AlignmentType.CENTER,
      spacing: { before: 200, after: 80 },
    }),
    new Paragraph({
      children: [
        new TextRun({ text: label + ". ", font: FONT, size: SZ_SMALL, bold: true }),
        ...rt(caption, { size: SZ_SMALL }),
      ],
      alignment: AlignmentType.JUSTIFIED,
      spacing: { after: 240, line: 240, lineRule: LineRuleType.AUTO },
    }),
  ];
}
function pngSize(buf) {
  return { w: buf.readUInt32BE(16), h: buf.readUInt32BE(20) };
}

/* references */
const REFS = [
  "Rossjohn, J.; Gras, S.; Miles, J. J.; Turner, S. J.; Godfrey, D. I.; McCluskey, J. T Cell Antigen Receptor Recognition of Antigen-Presenting Molecules. Annu. Rev. Immunol. 2015, 33, 169-200.",
  "La Gruta, N. L.; Gras, S.; Daley, S. R.; Thomas, P. G.; Rossjohn, J. Understanding the Drivers of MHC Restriction of T Cell Receptors. Nat. Rev. Immunol. 2018, 18, 467-478.",
  "Szeto, C.; Lobos, C. A.; Nguyen, A. T.; Gras, S. TCR Recognition of Peptide-MHC-I: Rule Makers and Breakers. Int. J. Mol. Sci. 2021, 22, 68.",
  "Krovi, S. H.; Kappler, J. W.; Marrack, P.; Gapin, L. Inherent Reactivity of Unselected TCR Repertoires to Peptide-MHC Molecules. Proc. Natl. Acad. Sci. U.S.A. 2019, 116, 22252-22261.",
  "Boughter, C. T.; Meier-Schellersheim, M. Conserved Biophysical Compatibility among the Highly Variable Germline-Encoded Regions Shapes TCR-MHC Interactions. eLife 2024, 12, RP90681.",
  "Dash, P.; Fiore-Gartland, A. J.; Hertz, T.; et al. Quantifiable Predictive Features Define Epitope-Specific T Cell Receptor Repertoires. Nature 2017, 547, 89-93.",
  "Hudson, D.; Fernandes, R. A.; Basham, M.; Ogg, G.; Koohy, H. Can We Predict T Cell Specificity with Digital Biology and Machine Learning? Nat. Rev. Immunol. 2023, 23, 511-521.",
  "Joglekar, A. V.; Li, G. T Cell Antigen Discovery. Nat. Methods 2021, 18, 873-880.",
  "Leem, J.; de Oliveira, S. H. P.; Krawczyk, K.; Deane, C. M. STCRDab: The Structural T-Cell Receptor Database. Nucleic Acids Res. 2018, 46, D406-D412.",
  "Gowthaman, R.; Pierce, B. G. TCR3d: The T Cell Receptor Structural Repertoire Database. Bioinformatics 2019, 35, 5323-5325.",
  "Jumper, J.; Evans, R.; Pritzel, A.; et al. Highly Accurate Protein Structure Prediction with AlphaFold. Nature 2021, 596, 583-589.",
  "Abramson, J.; Adler, J.; Dunger, J.; et al. Accurate Structure Prediction of Biomolecular Interactions with AlphaFold 3. Nature 2024, 630, 493-500.",
  "Wohlwend, J.; Corso, G.; Passaro, S.; Getz, N.; Reveiz, M.; Leidal, K.; Swiderski, W.; Atkinson, L.; Portnoi, T.; Chinn, I.; Silterra, J.; Jaakkola, T.; Barzilay, R. Boltz-1: Democratizing Biomolecular Interaction Modeling. bioRxiv 2024, 2024.11.19.624167 (preprint, v4 2025). DOI: 10.1101/2024.11.19.624167.",
  "Chai Discovery; Boitreaud, J.; Dent, J.; McPartlon, M.; Meier, J.; Reis, V.; Rogozhnikov, A.; Wu, K. Chai-1: Decoding the Molecular Interactions of Life. bioRxiv 2024, 2024.10.10.615955 (preprint). DOI: 10.1101/2024.10.10.615955.",
  "Bradley, P. Structure-Based Prediction of T Cell Receptor:Peptide-MHC Interactions. eLife 2023, 12, e82813.",
  "Yin, R.; Ribeiro-Filho, H. V.; Lin, V.; Gowthaman, R.; Cheung, M.; Pierce, B. G. TCRmodel2: High-Resolution Modeling of T Cell Receptor Recognition Using Deep Learning. Nucleic Acids Res. 2023, 51, W569-W576.",
  "Wu, F.; Zhao, Y.; Xiao, Y.; Qin, C.; Wang, F.; Wu, Z.; Huang, L.-K.; Liu, X.; Song, J.; He, B.; Rossjohn, J.; Yao, J. Fast and Accurate Modeling of TCR-Peptide-MHC Complexes Using tFold-TCR. bioRxiv 2025, 2025.01.12.632367 (preprint). DOI: 10.1101/2025.01.12.632367.",
  "Ascunce-Paris, A.; Farriol-Duran, R.; Romero-Durana, M.; Valencia, A.; Guallar, V. A Unified Framework for TCR-pMHC Structural Model Assessment. bioRxiv 2025, 2025.10.09.681411 (preprint). DOI: 10.1101/2025.10.09.681411.",
  "Lu, J.; Zhu, X.; Hu, X.; Zhang, C.; Feng, F. Benchmarking TCR-pMHC Structure Prediction: A Unified Evaluation and CDR3-Based Functional Insights. Brief. Bioinform. 2026, 27 (3), bbag289. DOI: 10.1093/bib/bbag289.",
  "Lin, Z.; Akin, H.; Rao, R.; et al. Evolutionary-Scale Prediction of Atomic-Level Protein Structure with a Language Model. Science 2023, 379, 1123-1130.",
  "Montemurro, A.; Schuster, V.; Povlsen, H. R.; et al. NetTCR-2.0 Enables Accurate Prediction of TCR-Peptide Binding by Using Paired TCR Alpha and Beta Sequence Data. Commun. Biol. 2021, 4, 1060.",
  "Lu, T.; Zhang, Z.; Zhu, J.; et al. Deep Learning-Based Prediction of the T Cell Receptor-Antigen Binding Specificity. Nat. Mach. Intell. 2021, 3, 864-875.",
  "Weber, A.; Born, J.; Rodriguez Martinez, M. TITAN: T-Cell Receptor Specificity Prediction with Bimodal Attention Networks. Bioinformatics 2021, 37, i237-i244.",
  "Sidhom, J.-W.; Larman, H. B.; Pardoll, D. M.; Baras, A. S. DeepTCR Is a Deep Learning Framework for Revealing Sequence Concepts within T-Cell Repertoires. Nat. Commun. 2021, 12, 1605.",
  "Peng, X.; Lei, Y.; Feng, P.; Jia, L.; Ma, J.; Zhao, D.; Zeng, J. Characterizing the Interaction Conformation between T-Cell Receptors and Epitopes with Deep Learning. Nat. Mach. Intell. 2023, 5, 395-407.",
  "Zhang, Y.; Wang, Z.; Jiang, Y.; Littler, D. R.; Gerstein, M.; Purcell, A. W.; Rossjohn, J.; Ou, H.-Y.; Song, J. Epitope-Anchored Contrastive Transfer Learning for Paired CD8+ T Cell Receptor-Antigen Recognition. Nat. Mach. Intell. 2024, 6, 1344-1358. DOI: 10.1038/s42256-024-00913-8.",
  "Gao, L.; Zhang, Y.; Ge, F.; Li, S.; Guo, Y.; Song, J.; Yu, D.-J. Structure-Directed Pan-Specific T-Cell Receptor-Peptide-Major Histocompatibility Complex Interaction Prediction. J. Chem. Inf. Model. 2025, 65, 4674-4686. DOI: 10.1021/acs.jcim.5c00055.",
  "Nielsen, M.; Eugster, A.; Jensen, M. F.; Goel, M.; Mayer, A.; Pelissier, A.; Valkiers, S.; Rodriguez Martinez, M.; Meynard-Piganeau, B.; Greiff, V.; Mora, T.; Walczak, A. M.; Croce, G.; Lopez Moreno, D.; Gfeller, D.; Meysman, P.; Barton, J. Lessons Learned from the IMMREP23 TCR-Epitope Prediction Challenge. ImmunoInformatics 2024, 16, 100045.",
  "Dens, C.; Laukens, K.; Bittremieux, W.; Meysman, P. The Pitfalls of Negative Data Bias for the T-Cell Epitope Specificity Challenge. Nat. Mach. Intell. 2023, 5, 1060-1062.",
  "Grazioli, F.; Moesch, A.; Machart, P.; et al. On TCR Binding Predictors Failing to Generalize to Unseen Peptides. Front. Immunol. 2022, 13, 1014256.",
  "Deng, L.; Ly, C.; Abdollahi, S.; Zhao, Y.; Prinz, I.; Bonn, S. Performance Comparison of TCR-pMHC Prediction Tools Reveals a Strong Data Dependency. Front. Immunol. 2023, 14, 1128326.",
  "Moris, P.; De Pauw, J.; Postovskaya, A.; et al. Current Challenges for Unseen-Epitope TCR Interaction Prediction and a New Perspective Derived from Image Classification. Brief. Bioinform. 2021, 22, bbaa318.",
  "Rajeshwar T., R.; et al. TCR-H: Explainable Machine Learning Prediction of T-Cell Receptor Epitope Binding on Unseen Datasets. Front. Immunol. 2024, 15, 1426173.",
  "Weber, A.; Pelissier, A.; Rodriguez Martinez, M. Reliable Interpretability of Biology-Inspired Deep Neural Networks and Challenges in TCR Binding Prediction. ImmunoInformatics 2024, 15, 100036.",
  "Whalen, S.; Schreiber, J.; Noble, W. S.; Pollard, K. S. Navigating the Pitfalls of Applying Machine Learning in Genomics. Nat. Rev. Genet. 2022, 23, 169-181.",
  "Kapoor, S.; Narayanan, A. Leakage and the Reproducibility Crisis in Machine-Learning-Based Science. Patterns 2023, 4, 100804.",
  "Walsh, I.; Fishman, D.; Garcia-Gasulla, D.; et al. DOME: Recommendations for Supervised Machine Learning Validation in Biology. Nat. Methods 2021, 18, 1122-1127.",
  "Teufel, F.; Gislason, M. H.; Almagro Armenteros, J. J.; Johansen, A. R.; Winther, O.; Nielsen, H. GraphPart: Homology Partitioning for Biological Sequence Analysis. NAR Genom. Bioinform. 2023, 5, lqad088.",
  "McDermott, M. B. A.; Hansen, L. H.; Zhang, H.; Angelotti, G.; Gallifant, J. A Closer Look at AUROC and AUPRC under Class Imbalance. Adv. Neural Inf. Process. Syst. 2024, 37.",
  "Lefranc, M.-P.; Pommie, C.; Ruiz, M.; et al. IMGT Unique Numbering for Immunoglobulin and T Cell Receptor Variable Domains and Ig Superfamily V-like Domains. Dev. Comp. Immunol. 2003, 27, 55-77.",
  "Cock, P. J. A.; Antao, T.; Chang, J. T.; et al. Biopython: Freely Available Python Tools for Computational Molecular Biology and Bioinformatics. Bioinformatics 2009, 25, 1422-1423.",
  "Joosten, R. P.; te Beek, T. A. H.; Krieger, E.; et al. A Series of PDB Related Databases for Everyday Needs. Nucleic Acids Res. 2011, 39, D411-D419.",
  "Mayer-Blackwell, K.; Schattgen, S.; Cohen-Lavi, L.; et al. TCR Meta-Clonotypes for Biomarker Discovery with tcrdist3 Enabled Identification of Public, HLA-Restricted Clusters of SARS-CoV-2 TCRs. eLife 2021, 10, e68605.",
  "Vaswani, A.; Shazeer, N.; Parmar, N.; et al. Attention Is All You Need. Adv. Neural Inf. Process. Syst. 2017, 30, 5998-6008.",
  "Loshchilov, I.; Hutter, F. Decoupled Weight Decay Regularization. International Conference on Learning Representations, 2019.",
  "Paszke, A.; Gross, S.; Massa, F.; et al. PyTorch: An Imperative Style, High-Performance Deep Learning Library. Adv. Neural Inf. Process. Syst. 2019, 32, 8024-8035.",
  "Pedregosa, F.; Varoquaux, G.; Gramfort, A.; et al. Scikit-Learn: Machine Learning in Python. J. Mach. Learn. Res. 2011, 12, 2825-2830.",
  "Holm, S. A Simple Sequentially Rejective Multiple Test Procedure. Scand. J. Stat. 1979, 6, 65-70.",
  "Field, C. A.; Welsh, A. H. Bootstrapping Clustered Data. J. R. Stat. Soc. B 2007, 69, 369-390.",
  "Huang, J.; Rauscher, S.; Nawrocki, G.; et al. CHARMM36m: An Improved Force Field for Folded and Intrinsically Disordered Proteins. Nat. Methods 2017, 14, 71-73.",
  "Jo, S.; Kim, T.; Iyer, V. G.; Im, W. CHARMM-GUI: A Web-Based Graphical User Interface for CHARMM. J. Comput. Chem. 2008, 29, 1859-1865.",
  "Jorgensen, W. L.; Chandrasekhar, J.; Madura, J. D.; Impey, R. W.; Klein, M. L. Comparison of Simple Potential Functions for Simulating Liquid Water. J. Chem. Phys. 1983, 79, 926-935.",
  "Abraham, M. J.; Murtola, T.; Schulz, R.; et al. GROMACS: High Performance Molecular Simulations through Multi-Level Parallelism from Laptops to Supercomputers. SoftwareX 2015, 1-2, 19-25.",
  "Hess, B. P-LINCS: A Parallel Linear Constraint Solver for Molecular Simulation. J. Chem. Theory Comput. 2008, 4, 116-122.",
  "Essmann, U.; Perera, L.; Berkowitz, M. L.; Darden, T.; Lee, H.; Pedersen, L. G. A Smooth Particle Mesh Ewald Method. J. Chem. Phys. 1995, 103, 8577-8593.",
  "Bussi, G.; Donadio, D.; Parrinello, M. Canonical Sampling through Velocity Rescaling. J. Chem. Phys. 2007, 126, 014101.",
  "Bernetti, M.; Bussi, G. Pressure Control Using Stochastic Cell Rescaling. J. Chem. Phys. 2020, 153, 114107.",
  "Hopkins, C. W.; Le Grand, S.; Walker, R. C.; Roitberg, A. E. Long-Time-Step Molecular Dynamics through Hydrogen Mass Repartitioning. J. Chem. Theory Comput. 2015, 11, 1864-1874.",
  "Knapp, B.; Ospina, L.; Deane, C. M. Avoiding False Positive Conclusions in Molecular Simulation: The Importance of Replicas. J. Chem. Theory Comput. 2018, 14, 6127-6138.",
  "Gowers, R. J.; Linke, M.; Barnoud, J.; et al. MDAnalysis: A Python Package for the Rapid Analysis of Molecular Dynamics Simulations. Proc. 15th Python Sci. Conf. 2016, 98-105.",
  "Shrake, A.; Rupley, J. A. Environment and Exposure to Solvent of Protein Atoms. Lysozyme and Insulin. J. Mol. Biol. 1973, 79, 351-371.",
  "Ma, J.; Singh, N. K.; Baker, B. M.; et al. Structural Basis of a Noncanonical T Cell Receptor Binding Geometry. Biochemistry 2020, 59, 4163-4175.",
 "Lee, J.; Cheng, X.; Swails, J. M.; et al.; Im, W. CHARMM-GUI Input Generator for NAMD, GROMACS, AMBER, OpenMM, and CHARMM/OpenMM Simulations Using the CHARMM36 Additive Force Field. J. Chem. Theory Comput. 2016, 12, 405-413.",
];

const refParagraphs = REFS.map((r, i) => new Paragraph({
  children: [
    new Bookmark({
      id: "ref" + (i + 1),
      children: [new TextRun({ text: "(" + (i + 1) + ")  ", font: FONT, size: SZ_SMALL })],
    }),
    new TextRun({ text: r, font: FONT, size: SZ_SMALL }),
  ],
  spacing: { after: 100, line: 240, lineRule: LineRuleType.AUTO },
  indent: { left: convertInchesToTwip(0.35), hanging: convertInchesToTwip(0.35) },
}));

/* ================================================================== */
/*  DOCUMENT BODY                                                      */
/* ================================================================== */
const body = [];

/* ---- title block ---- */
body.push(new Paragraph({
  children: [new TextRun({
    text: "TCR-IMAP: Homology-Controlled Evaluation of TCR\u2013MHC Contact Prediction",
    font: FONT, size: 34, bold: true,
  })],
  alignment: AlignmentType.CENTER,
  spacing: { after: 240, line: 300, lineRule: LineRuleType.AUTO },
}));
body.push(new Paragraph({
  children: [
    new TextRun({ text: "Muhammad Imran", font: FONT, size: SZ }),
    new TextRun({ text: "1", font: FONT, size: SZ, superScript: true }),
    new TextRun({ text: " and Dong-Qing Wei", font: FONT, size: SZ }),
    new TextRun({ text: "1,*", font: FONT, size: SZ, superScript: true }),
  ],
  alignment: AlignmentType.CENTER, spacing: { after: 100 },
}));
body.push(new Paragraph({
  children: [new TextRun({
    text: "1 State Key Laboratory of Microbial Metabolism, School of Life Sciences and Biotechnology, "
        + "Shanghai Jiao Tong University, Shanghai 200240, China",
    font: FONT, size: SZ_SMALL, italics: true,
  })],
  alignment: AlignmentType.CENTER, spacing: { after: 80 },
}));
body.push(new Paragraph({
  children: [new TextRun({ text: "* Corresponding author.", font: FONT, size: SZ_SMALL, italics: true })],
  alignment: AlignmentType.CENTER, spacing: { after: 320 },
}));

body.push(H1("Abstract"));
body.push(P(
"T-cell receptors (TCRs) recognize peptides presented by major histocompatibility complex (MHC) molecules, and the "
+ "residue pairs forming the TCR\u2013MHC interface are the level at which this recognition can be described "
+ "mechanistically. Because solved TCR\u2013peptide\u2013MHC structures remain few relative to the receptor "
+ "repertoire, sequence-based prediction of these contacts is an active area of research, but such methods are "
+ "evaluated on data sets assembled from families of near-identical receptors, and the contribution of that "
+ "redundancy to reported accuracy has not been quantified. In this study, we curated 174 "
+ "TCR\u2013peptide\u2013MHC class I structures, all with a human MHC class I heavy chain, which collapse into 69 TCR clusters at a sequence-similarity "
+ "ratio of 0.70, and defined contacts and the MHC candidate window from the crystal structures. We developed "
+ "TCR-IMAP, a dual-encoder cross-attention framework predicting the complete TCR-residue \u00d7 MHC-residue "
+ "contact map from sequence alone, together with a per-residue variant for comparison with published predictors. "
+ "Notably, we froze five partitioning protocols, declared the primary metric, protocol and comparison before "
+ "training, and measured each method against within-fold frequency priors. Our results show that apparent accuracy "
+ "depends on the evaluation design as strongly as on the model: restricting train/test similarity reduced the "
+ "advantage of TCR-IMAP over a frequency prior from +0.2750 to +0.0274 macro-AUPRC, and the ranking of two methods "
+ "reversed under a published single-split design. A training-size-matched control shows that this decline is not "
+ "explained by the smaller training sets that similarity filtering produces: drawing the same number of training "
+ "complexes at random instead left the advantage at +0.2296, and the paired contrast between the two matched-size "
+ "regimes was \u22120.2021 macro-AUPRC at the strictest threshold. Furthermore, under homology-controlled evaluation TCR-IMAP "
+ "exceeded the frequency prior and a reimplemented residue-level model on the pairwise task, whereas the "
+ "predeclared per-residue comparison remained unresolved at this sample size. In summary, we not only propose "
+ "TCR-IMAP for sequence-based TCR\u2013MHC contact map prediction but also provide the frozen protocols, the "
+ "characterized benchmark, the frequency-prior baselines and a model-free contact footprint as a reusable "
+ "evaluation framework for this task."));
body.push(P("**Keywords:** T-cell receptor; MHC class I; contact prediction; benchmark design; data leakage; "
+ "deep learning; molecular dynamics."));

body.push(H1("1. Introduction"));
body.push(P(
"The \u03b1\u03b2 T-cell receptor reads peptides presented by MHC molecules through an interface whose geometry is "
+ "conserved across otherwise unrelated receptors [[1,2]]. In complexes with MHC class I the receptor lies "
+ "diagonally across the peptide-binding groove, the germline-encoded CDR1 and CDR2 loops engage the \u03b11 and "
+ "\u03b12 helices of the heavy chain, and the somatically generated CDR3 loops engage the bound peptide [[3,4]]. "
+ "Which residue pairs actually form is the level at which recognition can be described mechanistically, and recent "
+ "structural analysis indicates that the repertoire-wide bias toward MHC is carried less by a fixed set of "
+ "conserved contacts than by biophysical compatibility between the germline loops and the MHC surface [[5]]. "
+ "Mapping TCRs to the antigens they recognize remains a central open problem in systems immunology [[6,7,8]], and "
+ "contact-level description is informative precisely where it is hardest to obtain."));
body.push(P(
"Those residue pairs are resolved experimentally by crystallography, and solved TCR\u2013pMHC structures number in "
+ "the hundreds against a repertoire of millions of distinct receptors [[9,10]]. Two routes have been taken to close "
+ "the gap. The first is structure prediction: AlphaFold2 and AlphaFold3, general co-folding models such as Boltz-1 "
+ "and Chai-1, and TCR-specific protocols including TCRdock, TCRmodel2 and tFold-TCR now generate complete complexes "
+ "[[11,12,13,14,15,16,17]]. Independent benchmarks nonetheless report that CDR3 conformation, docking orientation "
+ "and the TCR\u2013peptide interface remain the weakest components of those models [[18]]. On 70 complexes released "
+ "after the training cutoffs of the methods it evaluates, a benchmark of 13 models found MSA-based methods, and "
+ "AlphaFold3 in particular, the most accurate, and identified the pLDDT of the CDR3 region as an indicator of both "
+ "structural correctness and functional relevance [[19]]. Predicted structures are therefore increasingly reliable, "
+ "but the residue pairs they place in contact are inferred rather than observed, which is why the contacts analyzed "
+ "here are taken from crystal structures. The second route predicts the interaction directly from "
+ "sequence. Protein language models and sequence-based specificity predictors are now numerous [[20,21,22,23,24]], "
+ "and interface-level predictors have followed: TEIM predicts residue-level interaction patterns for "
+ "TCR\u2013epitope pairs [[25]], EPACT couples contrastive pretraining with structural fine-tuning and recovers "
+ "interchain distance matrices and contact residues [[26]], and Struct-TCS and Seq-TCS predict contact sites within "
+ "TCR\u2013pMHC complexes [[27]]."));
body.push(P(
"Establishing how well any of these performs is difficult, because the data on which they are evaluated are not "
+ "independent samples. Structures are deposited in families: the same receptor is re-solved with altered peptides, "
+ "with engineered mutations and in different crystal forms, so a randomly drawn test complex usually has a "
+ "near-identical counterpart in the training set. The consequences are documented in the adjacent task of "
+ "TCR\u2013epitope specificity prediction, where models that perform well on seen epitopes approach random "
+ "performance on unseen ones [[30,31,32,33]], where negative-set construction alone has been shown to carry much of "
+ "the reported signal [[29]], and where the IMMREP23 community challenge traced part of the measured performance to "
+ "the construction of the test data rather than to the methods [[28]]. The same failure mode is general to machine "
+ "learning in biology, and redundancy-aware partitioning together with explicit reporting of train/test similarity "
+ "are long-standing recommendations [[34,35,36,37,38]]."));
body.push(P(
"For TCR\u2013MHC contact prediction specifically, the question has not been asked. The principal comparator for "
+ "this task was evaluated on 94 structures under single train/test partitions, with contact-site splits as small as "
+ "41 training and 12 test complexes, and novelty was controlled by requiring that the concatenated CDR1\u2013CDR3 "
+ "sequence of each test chain be absent from training [[27]]. Because CDR1 and CDR2 are germline-encoded, that "
+ "criterion leaves individual loop sequences shared between training and test. Reported accuracy is therefore a "
+ "joint property of the method and the partitioning scheme, and the two have not been separated."));
body.push(P(
"This study separates them. We assembled 174 TCR\u2013pMHC class I structures with human MHC class I, defined contacts from the "
+ "crystal structures, and fixed a structurally defined MHC candidate window on the \u03b11/\u03b12 platform before "
+ "any model was trained. We then built TCR-IMAP, a dual-encoder cross-attention network that predicts the complete "
+ "TCR-residue \u00d7 MHC-residue contact map from amino-acid identities alone, together with a per-residue variant "
+ "for comparison with published residue-level predictors, and froze five partitioning protocols, four of them orderable by how much "
+ "train/test similarity they permit, with a declared primary metric, protocol and comparison. Four analyses then measure how much of the "
+ "apparent accuracy belongs to the evaluation rather than to the model: an identity-threshold sweep over the "
+ "training set, a correction of the validation split used for model selection, a docking-angle regression evaluated "
+ "under grouped partitions, and a comparison of method ordering under the comparator\u2019s own single-split "
+ "design. A molecular dynamics case study then asks whether high-scoring residue pairs that are non-contacts in the "
+ "crystal structure form transient contacts in solution. The outcome is an evaluation framework, and a set of "
+ "measurements made with it, rather than a claim of architectural novelty."));

body.push(...figure(FIG1, "Figure 1",
"Structural benchmark and sequence-based contact prediction with TCR-IMAP. (A) The benchmark and the evaluation "
+ "it feeds. STCRDab TCR–peptide–MHC class I complexes provide the structures, from which the crystal contact "
+ "labels and the aligned TCR and MHC sequence inputs are derived. The illustrated complex is 1AO7: the TCR "
+ "variable domains, the peptide and the MHC class I platform are shown, with TCRα in blue, TCRβ in purple, MHC "
+ "in green and the peptide in magenta; domains outside this interface view, including β2-microglobulin and the "
+ "α3 domain, are omitted. A contact is a minimum inter-residue heavy-atom distance of 4.5 Å. TCR-IMAP, the "
+ "reimplemented comparators and the within-fold frequency priors are then trained and evaluated under the five "
+ "frozen random and grouped partitioning protocols, compared as paired per-complex differences and assessed by "
+ "cluster bootstrap; the primary protocol is tcr70. Protocol definitions, candidate masks, the per-residue task "
+ "and the statistical procedure are given in Methods. (B) The pairwise TCR-IMAP architecture. The TCR side "
+ "embeds 75 aligned residue rows as residue, segment and position terms and the MHC side 84 aligned residue "
+ "columns as residue and position terms; separate Transformer encoders then produce residue representations. In "
+ "the cross-attention block the TCR representations act as queries (Q) and the MHC representations as keys (K) "
+ "and values (V), and a residual connection with layer normalization updates the TCR representations. The "
+ "bilinear head combines the updated TCR representations with the MHC representations, adds a learned "
+ "pair-position term, a scale and a bias, and a sigmoid (σ) maps the result to contact probabilities between 0 "
+ "and 1. Both maps are 1AO7: the labels in (A) are its crystal contacts and the map in (B) is its stored "
+ "out-of-fold prediction under the primary tcr70 protocol, not a drawn illustration. Grey cells are residue "
+ "pairs left unresolved by the structure, which are masked in training and in scoring, and both maps keep the "
+ "original row and column order. Structural coordinates define the labels and the mask and are not input "
+ "features of this sequence-based model. The architecture shown is the pairwise variant, which predicts the "
+ "complete TCR-residue × MHC-residue map; the predeclared primary comparison uses the separate per-residue "
+ "variant, which predicts for each TCR position whether it contacts any window residue and so matches the output "
+ "format of the published comparators (Methods). Descriptive throughout; no performance value appears in this "
+ "figure.", 7.1, true));

/* ---- methods ---- */
body.push(H1("2. Materials and Methods"));
body.push(P(
"Sections 2.1\u20132.4 describe the benchmark and the way it is partitioned, Sections 2.5\u20132.9 the model, its "
+ "variants and the baselines against which it is measured, Section 2.10 the metrics and the statistical analysis, "
+ "and Section 2.11 the molecular dynamics. Each partitioning protocol, the primary metric and the primary "
+ "comparison were fixed before the corresponding results were produced."));

body.push(H2("2.1. Structures, dataset construction and contact definition"));
body.push(P(
"TCR–peptide–MHC class I complexes were taken from the IMGT-numbered structure set of the structural T-cell "
+ "receptor database (STCRDab), retrieved on 15 August 2026 [[9]]. The download comprised 1,002 summary rows over 634 "
+ "unique PDB entries. Four criteria were applied simultaneously — an αβ TCR, a peptide antigen, MHC "
+ "class I and a human MHC heavy chain — leaving 275 rows, which were deduplicated by PDB entry to give 175 "
+ "complexes. One complex, 6uz1, was excluded on structural-scope grounds: it is a noncanonical single-chain A6 "
+ "variant whose contacts with MHC lie predominantly in the α3 domain (17 of its 31 contacts with resolved MHC "
+ "residues), outside the α1/α2 platform analyzed here [[62]]. The primary set therefore comprised 174 "
+ "complexes. MHC allele annotations were obtained by regular-expression parsing of the free-text compound field of "
+ "the STCRDab summary, which is why only 66 of the 175 indexed complexes carry one; this limits two of the "
+ "evaluation protocols (Section 2.3)."));
body.push(P(
"Contacts were defined from the crystal structures alone. Standard residues and heavy atoms were used, with the "
+ "default alternate location as parsed by Biopython [[41]]. For a TCR residue *i* and an MHC residue *j* the "
+ "interresidue distance was taken as the minimum over their heavy atoms,"));
body.push(EQ([
  msubsup("d", "ij", "min"), mr(" = "),
  new MathSubScript({ children: [mr("min")], subScript: [mr("u∈A(i), v∈A(j)")] }),
  mr("  d("), msub("x", "u"), mr(", "), msub("x", "v"), mr(")"),
]));
body.push(P(
"where *A*(·) denotes the heavy atoms of a residue, *x*_{u} the coordinates of atom *u*, and *d*(·,·) the "
+ "Euclidean distance. The pair was labeled a contact, *y*_{ij} = 1, when *d*^{min}_{ij} ≤ 4.5 Å and a "
+ "non-contact otherwise. A 6 Å "
+ "label set was built identically and used only for the comparator-subset analysis. A pair was treated as "
+ "*resolved* when both residues were present in the structure, and only resolved pairs entered training losses and "
+ "metrics. Across the 174 complexes the label array contained 893,928 resolved pairs, of which 5,281 (0.5908 %) were "
+ "contacts at 4.5 Å and 10,317 (1.1541 %) at 6 Å."));

body.push(H2("2.2. MHC candidate window and residue representation"));
body.push(P(
"MHC residues were restricted to a fixed, structurally defined window on the α1/α2 platform: IMGT positions "
+ "50–90 (α1) and 1050–1090 (α2) of the heavy chain [[40]]. The window was defined from IMGT domain "
+ "positions and not from observed contact frequencies. Every residue occupying a window position in any complex "
+ "received its own column; the insertion residues 1061A and 1072A, present in all 174 complexes, therefore form "
+ "separate columns, giving 84 columns. All 84 positions were resolved in all 174 complexes. Contacts were also "
+ "computed over the full 182-position α1/α2 G-domain to quantify what the window omits: the window captured "
+ "5,281 of the 5,294 G-domain contacts (99.75 %). The α3 domain and β2-microglobulin were outside the scope "
+ "of the contact task. For the descriptive footprint analysis, each G-domain position was assigned a "
+ "secondary-structure region from DSSP [[42]]: α1 or α2 helix when the majority code across structures was H, "
+ "G or I, β-sheet floor when it was E or B, and other G-domain otherwise, each call requiring at least 30 "
+ "structures resolving the position and a majority fraction above 0.5."));
body.push(P(
"Residues were keyed by IMGT number and insertion code throughout, in the contact labels, in the model inputs and in "
+ "every structural visualization, so that no residue replaced another sharing its IMGT number. On the TCR side, "
+ "candidate positions were IMGT 27–38 (CDR1), 56–65 (CDR2) and 105–124 (CDR3 together with the J-region "
+ "framework positions 118–124) of both chains. The TCR axis contained one row for every (segment, IMGT number, "
+ "insertion code) occupied in at least one complex, giving 75 rows: 68 occupied IMGT positions and 7 CDR3 insertion "
+ "positions. The contact labels are therefore an array of shape 174 × 75 × 84, with a matching array marking "
+ "resolved pairs. Model inputs were amino-acid identities only; no coordinates were supplied to TCR-IMAP at training "
+ "or prediction time. Tokens were placed by IMGT position and never packed, so the token in each row and the label "
+ "in that row refer to the same residue; this was verified against an independent per-residue manifest built from "
+ "the structure files, with no misalignment among 10,642 resolved TCR residues and 14,616 of 14,616 MHC tokens "
+ "matching."));

body.push(H2("2.3. Sequence similarity, clustering and partitioning protocols"));
body.push(P(
"Selection required a human MHC class I heavy chain but placed no condition on the receptor. The resulting receptors are "
+ "human in 165 of the 174 complexes; eight carry mouse \u03b1 and \u03b2 chains and one carries a mixed chimpanzee/human "
+ "\u03b1-chain annotation, in every case paired with human MHC class I."));
body.push(P(
"Sequence similarity, clustering and fold assignment used one string per complex, formed by concatenating the six "
+ "CDR segments in the order CDR1α, CDR2α, CDR3α, CDR1β, CDR2β, CDR3β over IMGT positions "
+ "27–38, 56–65 and 105–124, with one residue per IMGT position. Similarity between strings *a* and *b* "
+ "was the normalized indel similarity computed by Levenshtein.ratio, an edit-distance measure of the kind used to cluster TCR repertoires "
+ "[[43]],"));
body.push(EQ([
  mr("r(a, b) = "),
  new MathFraction({
    numerator: [msub("n", "a"), mr(" + "), msub("n", "b"), mr(" − "), msub("d", "I"), mr("(a, b)")],
    denominator: [msub("n", "a"), mr(" + "), msub("n", "b")],
  }),
]));
body.push(P(
"where *n*_{a} and *n*_{b} are the lengths of the two strings and *d*_{I} is the indel distance, the minimum number of insertions and deletions transforming one string into the other, so that a substitution costs two. This is an edit-distance similarity rather than gap-free alignment "
+ "identity. Complexes were clustered by single linkage at ratios of 0.70, 0.80 and 0.90, following the "
+ "homology-partitioning principle recommended for biological sequence data [[38]], giving 69, 82 and 89 "
+ "clusters. Five partitioning protocols were frozen before any model was trained under them and were never "
+ "regenerated: random (five folds of complexes; 174 complexes evaluated), allele (five folds of MHC alleles; 64), "
+ "an intersection-selected partition retained from earlier work and reported here as legacy (three folds; 30), temporal (a single "
+ "split with test complexes deposited on or after 1 January 2020; 60) and tcr70, the primary protocol (five folds of "
+ "clusters at a ratio of 0.70; 174). Folds with fewer than three test complexes were not evaluated. The allele and "
+ "intersection protocols cover only 64 and 30 of the 174 complexes because 108 complexes carry no parsed allele annotation "
+ "and the frozen protocol file assigns none of them to an allele or intersection test fold; one allele fold of two "
+ "complexes was additionally dropped by the three-complex minimum. The protocol file was not regenerated to repair "
+ "this. Under tcr70 the maximum train/test similarity ratio was 0.696 on the clustering strings, 0.707 over full "
+ "CDR-plus-framework sequences and 0.682 over canonical CDR sequences alone; no test complex exceeded 0.90 under any "
+ "definition, and duplicate-family leakage was zero in every fold. No cluster is shared between a test fold and the "
+ "complexes used to train or to validate under it: cluster disjointness was verified against the actual training and "
+ "early-stopping membership of each fold, not against the partition alone."));
body.push(P(
"Two of these protocols require qualification. The partition labeled double in earlier work is not a double holdout. "
+ "Its test folds were formed as the intersection of a held-out cluster group and a held-out allele group, but the "
+ "frozen protocol file records only the test complexes, and every consumer takes the training set to be all remaining "
+ "complexes. The held-out groups therefore remain in training: across its three folds, 13 of 30 test complexes share "
+ "a cluster with a training complex and 21 share an annotated allele, leaving 5 that are unseen in both senses. It is "
+ "consequently not a stricter protocol than the others, it is excluded from any ordering by strictness and from any "
+ "claim about jointly unseen clusters and alleles, and its results are reported only as a legacy partition. No "
+ "corrected replacement was constructed for this study. The allele protocol, in turn, holds out alleles by their "
+ "parsed annotation rather than by MHC sequence. The 84-residue \u03b11/\u03b12 window takes only 28 distinct values "
+ "across the 174 complexes, and for 32 of the 64 evaluated allele-test complexes an identical window is present in "
+ "training; in every one of those 32 cases the matching training complex is one of the 108 whose allele annotation "
+ "could not be parsed and which therefore never leaves the training set. An identical window does not establish an "
+ "identical allele, but it does mean the MHC input itself was not held out, so the protocol is described as holdout "
+ "by available allele annotation."));

body.push(H2("2.4. Identity-threshold sweep"));
body.push(P(
"To measure how pairwise per-complex macro-AUPRC depends on the maximum train/test similarity ratio, the pairwise "
+ "model and the pair-frequency prior were retrained on a fixed five-fold random partition of the 174 complexes "
+ "(seed 0, independent of the frozen protocols) after removing from each training fold every complex whose "
+ "similarity ratio to any test complex of that fold reached a threshold of 0.95, 0.90, 0.85, 0.80, 0.75 or 0.70, or "
+ "with no restriction. A fold was skipped if fewer than 30 training complexes remained. The test set was not "
+ "modified: all 174 complexes were evaluated at every threshold, so the sweep manipulates training-set composition "
+ "alone."));

body.push(H2("2.5. Model architecture"));
body.push(P(
"TCR-IMAP is a dual-encoder network with cross-attention that predicts a full TCR-row × MHC-column contact map "
+ "from sequence. TCR tokens were embedded (*d* = 128) and summed with a learned segment embedding over the six CDR "
+ "segments and a learned positional embedding per row, then encoded by a two-layer pre-norm Transformer encoder "
+ "(8 heads, feed-forward width 256, dropout 0.2) [[44]]. MHC tokens were embedded with their own token and "
+ "positional embeddings and encoded by a one-layer pre-norm Transformer of the same width. Each attention operation "
+ "is the scaled dot-product form"));
body.push(EQ([
  mr("Attention(Q, K, V) = softmax"),
  new MathRoundBrackets({ children: [
    new MathFraction({
      numerator: [mr("Q"), msup("K", "⊤")],
      denominator: [new MathRadical({ children: [msub("d", "k")] })],
    }),
  ]}),
  mr(" V"),
]));
body.push(P(
"with *Q*, *K* and *V* the query, key and value projections and *d*_{k} the per-head width. The TCR representations "
+ "attended to the MHC representations through multi-head cross-attention (8 heads, dropout 0.2), followed by a "
+ "residual connection and layer normalization,"));
body.push(EQ([
  msup("H", "′"), mr(" = LayerNorm"),
  new MathRoundBrackets({ children: [mr("H + MultiHead(H, G, G)")] }),
]));
body.push(P(
"where *H* ∈ ℝ^{75×d} and *G* ∈ ℝ^{84×d} are the encoded TCR and MHC representations. The "
+ "contact logit for TCR row *i* and MHC column *j* was"));
body.push(EQ([
  msub("L", "ij"), mr(" = s "),
  new MathFraction({
    numerator: [new MathAngledBrackets({ children: [msub("W", "t"), msub("h", "i"), mr(", "), msub("W", "m"), msub("g", "j")] })],
    denominator: [new MathRadical({ children: [mr("d")] })],
  }),
  mr(" + "), msub("P", "ij"), mr(" + b"),
]));
body.push(P(
"where *h*_{i} and *g*_{j} are the TCR and MHC representations of row *i* and column *j*, *W*_{t} and *W*_{m} are "
+ "learned bias-free projections, ⟨·,·⟩ is the inner product, *P* is a learned 75 × 84 positional "
+ "term initialized at zero, and *s* and *b* are learned scalars initialized at 1 and −4. Padding positions were "
+ "represented by their segment and positional embeddings only. TCR-IMAP has 529,310 trainable parameters at the "
+ "75 × 84 input size used for every result reported here."));

body.push(H2("2.6. Training"));
body.push(P(
"Models were trained full-batch for a fixed 400 epochs with AdamW (learning rate 10^{−3}, weight decay 0.01) and "
+ "gradient-norm clipping at 1.0 [[45,46]]. The loss was a class-weighted binary cross-entropy averaged over the set "
+ "*R* of resolved pairs,"));
body.push(EQ([
  mr("ℒ = − "),
  new MathFraction({ numerator: [mr("1")], denominator: [msub("N", "R")] }),
  new MathSum({ children: [msub("ℓ", "ij")], subScript: [mr("(i,j)∈R")] }),
]));
body.push(EQ([
  msub("ℓ", "ij"), mr(" = w "), msub("y", "ij"), mr(" log σ("), msub("L", "ij"),
  mr(") + (1 − "), msub("y", "ij"), mr(") log(1 − σ("), msub("L", "ij"), mr("))"),
]));
body.push(P(
"where *N*_{R} is the number of resolved pairs, σ is the logistic function and the positive-class weight "
+ "*w* is the ratio of negative to positive "
+ "resolved pairs in the training fold. The pairwise models use no validation subset: no early stopping and no "
+ "validation-based model selection were applied to them, so no test complex influenced training and the pairwise "
+ "results are unaffected by the validation-split correction described in Section 2.7. Three models were trained per "
+ "fold with seeds 42, 43 and 44 and their sigmoid outputs were averaged before scoring. All predictions reported are "
+ "out-of-fold."));

body.push(H2("2.7. Per-residue variant and cluster-disjoint model selection"));
body.push(P(
"For the per-residue task — whether a TCR row contacts any MHC window residue — a variant with the same "
+ "encoder–cross-attention design (*d* = 96, 6 heads, dropout 0.25) and a per-row two-layer output head was used "
+ "(291,073 parameters). A row was labeled positive if it contacted any window column and resolved if any of its "
+ "pairs was resolved. Per-residue models, including the comparators, were trained with Adam (learning rate "
+ "10^{−3}, weight decay 10^{−4}) for up to 300 epochs with the same weighted loss; validation AUPRC was "
+ "evaluated every five epochs, training stopped after ten evaluations without an increase, and the best-validation "
+ "weights were restored. Three seeds were averaged."));
body.push(P(
"Validation complexes were held out by whole TCR cluster. In each fold, clusters were drawn at random without "
+ "replacement from those present in the fold's training complexes, and all complexes of a drawn cluster were held "
+ "out, until at least 15 % of the fold's training complexes were in validation (achieved 15.1–20.2 %, median "
+ "16.3 %, across the 18 evaluated folds). No validation complex therefore shared a cluster with any complex "
+ "remaining in training. The same procedure, random seed and hyperparameters were used for all three per-residue "
+ "models, so model selection cannot favor one of them. The earlier scheme, early stopping on a random 15 % of the "
+ "non-test complexes, is reported as a sensitivity analysis because a validation subset drawn at random over "
+ "complexes is not independent of the complexes left in training when 57.5 % of complexes have an exact TCR twin."));

body.push(H2("2.8. Ablations and the peptide-aware variant"));
body.push(P(
"Under tcr70, eight ablations of the pairwise model were trained and evaluated exactly as the full model: masking "
+ "MHC residue identity, masking TCR residue identity, shuffling TCR residues among the occupied rows of each "
+ "segment, shuffling MHC residues across window columns, removing CDR3 and FR4 tokens, removing the germline CDR1 "
+ "and CDR2 tokens, removing cross-attention, and replacing the bilinear head with a plain dot product. Shuffles were "
+ "drawn once per complex with a fixed seed and kept across folds and training seeds. A secondary analysis added the "
+ "peptide sequence (up to 15 residues, packed, since peptides have no IMGT positions) through an embedding and a "
+ "second cross-attention step, with an otherwise identical model trained without the peptide as reference and a "
+ "shuffled-peptide control."));

body.push(H2("2.9. Baselines and comparators"));
body.push(P(
"Two within-fold frequency priors were fitted on the training labels of each fold only. The pair-frequency prior "
+ "predicts, for every (TCR row, MHC column) cell, the fraction of training complexes resolving that cell in which it "
+ "is a contact,"));
body.push(EQ([
  msub("f", "ij"), mr(" = "),
  new MathFraction({
    numerator: [new MathSum({ children: [msubsup("S", "ij", "c"), msubsup("y", "ij", "c")], subScript: [mr("c∈T")] })],
    denominator: [mr("max"), new MathRoundBrackets({ children: [
      new MathSum({ children: [msubsup("S", "ij", "c")], subScript: [mr("c∈T")] }), mr(", 1"),
    ]})],
  }),
]));
body.push(P(
"where *T* is the training fold and *S*^{c}_{ij} indicates that the pair is resolved in complex *c*. The "
+ "position-frequency prior is the same quantity for per-residue labels. Nearest-neighbor label transfer assigned "
+ "each test complex the per-residue labels of the training complex with the highest similarity ratio."));
body.push(P(
"The identity-threshold sweep confounds train/test similarity with training-set size, because removing similar "
+ "complexes also removes complexes. A training-size-matched control separates the two. For each of three thresholds "
+ "(0.95, 0.85 and 0.70) and each of the sweep\u2019s five folds, we counted how many training complexes the similarity "
+ "filter had retained and drew that same number uniformly without replacement from the fold\u2019s full non-test pool, "
+ "applying no similarity filter. The sweep\u2019s own fixed five-fold partition and its full set of 174 test complexes "
+ "were reused unchanged, as were the p84v2 representation, the pairwise architecture and the training configuration; "
+ "the pair-frequency prior was fitted on exactly the same subset as its model. Five draws were run with seeds fixed "
+ "in advance (1001\u20131005) and consumed in a declared order, each with the sweep\u2019s three training seeds, giving "
+ "225 model fits. The design, including the thresholds, the draw seeds, the metric, the comparison and the stopping "
+ "rules, was recorded and committed before any model was trained. The principal estimate is the mean over draws of "
+ "the per-complex difference between the model and the prior, bootstrapped over the 69 clusters and conditional on "
+ "the five realized draws; the alternative of averaging predictions across draws before scoring is reported "
+ "separately, since average precision is nonlinear and the two are not the same quantity. Draw-to-draw variation is "
+ "reported as a range over the five draws and is never pooled with the cluster-bootstrap interval. Subset "
+ "memberships, achieved train/test similarities, per-draw values and both estimands are given in Supporting "
+ "Table S9."));

body.push(P(
"Struct-TCS and Seq-TCS were reimplemented from the description of Gao et al. [[27]]. In our implementation, TCR and "
+ "MHC tokens are embedded and passed through two-layer LSTMs; the TCR branch then applies two multi-head "
+ "graph-attention blocks with per-head gating over a residue graph connecting TCR residues whose Cβ atoms (Cα for "
+ "glycine) lie within 8 Å in the crystal structure, and the MHC branch two self-attention blocks, followed by "
+ "cross-attention, a feed-forward layer and a per-residue sigmoid output. Struct-TCS therefore receives the complex's own TCR structure as input, whereas "
+ "TCR-IMAP does not. Seq-TCS is identical with the graph attention replaced by self-attention. The released "
+ "implementation could not be executed as published: its module layer imports a module absent from the repository, "
+ "its graph-attention encoder passes six positional arguments to a five-argument constructor, and its training loop "
+ "never steps the optimizer. The repository also contains no feature-construction code, graph builder, data set, "
+ "splits or trained weights, so both models were reimplemented. Model selection was additionally moved off the test "
+ "set, on which the released code early-stops. Both reimplemented variants are provided as source, together "
+ "with a provenance table recording every deviation from the published specification; no patch to the released "
+ "code was produced, because it was reimplemented rather than repaired. "
+ "Our reimplementation departs from the published specification in four respects that are stated here rather than "
+ "implied. The published models encode each residue with BLOSUM62 together with four physicochemical properties "
+ "— hydrophobicity, hydrophilicity, mass and isoelectric point — giving 24 features per residue; ours uses "
+ "learned embeddings trained from scratch. The published graph takes the Cα atom as each residue's "
+ "pseudoposition, ours the Cβ; the published edge radius is not stated numerically and no graph code is released, "
+ "so 8 Å is our choice. The published description applies graph attention to both the TCR and the MHC branch, "
+ "whereas the released code applies it to the TCR branch only, which is what we follow; the per-head gating in our "
+ "implementation is described in neither, and is taken from an unused class in the released repository; its effect "
+ "is quantified in Section 3.5. Finally, the "
+ "published contact threshold is 6 Å, whereas our primary analysis uses 4.5 Å and reports 6 Å as a secondary "
+ "cutoff. These are differences from an available specification, not consequences of the code defects. Both comparators were trained under the same protocols, folds, "
+ "candidate positions, labels, validation scheme and seeds as the TCR-IMAP per-residue variant. That study reports "
+ "94 structures in its STRUCT-CS data set, while its released identifier list contains 95 distinct entries; the "
+ "discrepancy is unexplained in the available material. Of those 95 listed identifiers, 88 are among our 174; on these the two methods were compared at the 4.5 Å and 6 "
+ "\u00c5 cutoffs under a single random split of the kind used in that study (22 test complexes) and under tcr70 "
+ "(88 test complexes). Only the structure subset and the cutoff are shared with the original evaluation, so this "
+ "analysis is not a reproduction of that study's own evaluation."));
body.push(P(
"TEIM-Res [[25]] was reimplemented for the pairwise task with one-dimensional convolutional encoders, an "
+ "outer-product interaction tensor projected to 32 channels, three dilated two-dimensional residual blocks and a "
+ "1 × 1 output convolution; the pretrained peptide autoencoder of the original model does not apply to "
+ "TCR–MHC contacts, so all encoders were trained from scratch. Training was not reproducible under default GPU "
+ "settings because the cuDNN kernel selected for the weight gradient of the dilated two-dimensional convolutions is "
+ "nondeterministic; this is a property of our reimplementation on this hardware and not of the published method. It "
+ "was removed by requesting deterministic algorithms, at a 12 % increase in runtime, under which two independent "
+ "runs reproduced every prediction bit for bit. These are architecture comparisons under common inputs, not reproductions of "
+ "published performance. No comparator is run on its own data set, features or splits. The original Struct-TCS study "
+ "already addresses TCR contacts with the MHC helices as well as with the peptide, across four single-split subsets "
+ "of 9 to 14 test complexes each, so the task is not new here; what differs is the residue representation, the "
+ "candidate set, the graph, the contact threshold, the data set and the partitioning. TEIM-Res predicts "
+ "CDR3\u03b2\u2013peptide contact maps in the original work and is given the full CDR\u2013MHC tensor here without "
+ "its pretrained sequence encoder, which has no counterpart for this input. None of the values reported here should "
+ "be compared with any value in the original publications."));
body.push(P(
"A separate sequence-based regression network, denoted AngleNet and distinct from TCR-IMAP, was trained to predict "
+ "the TCR docking angle recorded in STCRDab from the same token inputs (two-layer pre-norm Transformer TCR encoder, "
+ "*d* = 96, 6 heads, dropout 0.3, cross-attention to MHC embeddings, pooled regression head; AdamW, learning rate "
+ "10^{−3}, weight decay 0.05, smooth L1 loss, 400 epochs). It was evaluated under five-fold random, "
+ "allele-holdout and TRBV-holdout partitions formed within the analysis script, each with three seeds, on the 172 "
+ "complexes with a recorded docking angle."));

body.push(H2("2.10. Metrics and statistical analysis"));
body.push(P(
"For each test complex, the area under the precision–recall curve was computed over its resolved candidate pairs "
+ "(pairwise task) or resolved TCR rows (per-residue task) as the average precision [[47]],"));
body.push(EQ([
  msub("AP", "c"), mr(" = "),
  new MathSum({ children: [
    new MathRoundBrackets({ children: [msub("R", "n"), mr(" − "), msub("R", "n−1")] }),
    msub("P", "n"),
  ], subScript: [mr("n")] }),
  mr(",   macro-AUPRC = "),
  new MathFraction({ numerator: [mr("1")], denominator: [msub("N", "C")] }),
  new MathSum({ children: [msub("AP", "c")], subScript: [mr("c∈C")] }),
]));
body.push(P(
"where *P*_{n} and *R*_{n} are precision and recall at the *n*-th threshold, *C* is the set of evaluated complexes "
+ "and *N*_{C} its size. Complexes with no positive or no negative resolved candidate were undefined and omitted from that "
+ "comparison. The primary metric was the per-complex macro-AUPRC; per-complex ROC AUC was recorded as secondary, since the two summaries answer different questions about a ranking [[39]]. "
+ "Per-complex AUPRC correlates with a complex's number of contacts and is therefore not a measure of difficulty, so "
+ "statements about which complexes are easier or harder are based on paired differences between methods. The primary "
+ "metric, the primary protocol (tcr70) and the primary comparison (TCR-IMAP versus Struct-TCS on the per-residue "
+ "task) were declared before the analysis."));
body.push(P(
"Complexes are not independent, since complexes sharing a TCR cluster are correlated; no test was computed by "
+ "treating the 174 complexes as independent observations. Uncertainty for every comparison between two methods was "
+ "estimated from the paired per-complex differences by a pooled cluster bootstrap over the 69 clusters [[49]]: "
+ "clusters were resampled with replacement, all complexes of each drawn cluster were pooled, and the difference of "
+ "the two methods' complex-level means was recomputed,"));
body.push(EQ([
  msup("Δ", "(b)"), mr(" = "),
  new MathFraction({ numerator: [mr("1")], denominator: [msup("N", "(b)")] }),
  new MathSum({ children: [
    new MathRoundBrackets({ children: [msubsup("AP", "c", "A"), mr(" − "), msubsup("AP", "c", "B")] }),
  ], subScript: [mr("c∈"), msup("C", "(b)")] }),
]));
body.push(P(
"over *b* = 1, …, *B* resamples, where *C*^{(b)} is the *b*-th resampled set of complexes, *N*^{(b)} its size "
+ "and AP^{A}_{c} and AP^{B}_{c} the per-complex values of the two methods. Confidence intervals were 95 % "
+ "percentile intervals of the bootstrap distribution. Writing *p*_{lo} and *p*_{hi} for the proportions of resampled "
+ "differences at or below and at or above zero, the two-sided p value was"));
body.push(EQ([
  mr("p = 2 min( "), msub("p", "lo"), mr(" , "), msub("p", "hi"), mr(" ),   p ← max( p , 1/B )"),
]));
body.push(P(
"floored at 1/*B* and reported as p < 0.001 at the floor. We used *B* = 20,000 resamples with seed 42 wherever "
+ "per-complex arrays were available, including the identity-threshold sweep: that stage originally stored only "
+ "aggregates, and its per-complex values were regenerated under its own configuration and seeds, reproducing every "
+ "frozen aggregate to four decimal places before being used. The 2,000-resample values the sweep first recorded are "
+ "retained in the Supporting Information as the superseded record; the point estimates are identical under both "
+ "settings. For the peptide analysis no per-complex arrays were stored, and the values recorded at the time of the "
+ "run (2,000 resamples) are reported and labeled as such. As a sensitivity analysis, each comparison was repeated with the 69 clusters weighted equally. Point "
+ "estimates are independent of the bootstrap seed and the p values vary with it; across ten seeds, p ranged "
+ "0.0001–0.0002 for the pairwise comparison with the pair-frequency prior, 0.077–0.088 for the per-residue "
+ "comparison with Struct-TCS and 0.53–0.56 for the per-residue comparison with the position prior, so no "
+ "comparison crosses the two-sided 0.05 level across these seeds. The eight ablations were corrected by the Holm "
+ "step-down procedure [[48]]; Bonferroni-simultaneous 99.375 % intervals were taken from the same resampled "
+ "differences. Other comparisons are reported unadjusted, as secondary to the predeclared primary comparison. "
+ "Throughout, “significant” denotes the stated test at a two-sided level of 0.05."));
body.push(P(
"For the docking-angle regression, mean absolute error was compared with a median predictor and the coefficient of "
+ "determination was computed relative to the mean predictor,"));
body.push(EQ([
  msup("R", "2"), mr(" = 1 − "),
  new MathFraction({
    numerator: [new MathSum({ children: [
        new MathSuperScript({ children: [new MathRoundBrackets({ children: [msub("θ", "c"), mr(" − "), msubhat("θ", "c")] })], superScript: [mr("2")] })
      ], subScript: [mr("c")] })],
    denominator: [new MathSum({ children: [
        new MathSuperScript({ children: [new MathRoundBrackets({ children: [msub("θ", "c"), mr(" − "), mr("θ̄")] })], superScript: [mr("2")] })
      ], subScript: [mr("c")] })],
  }),
]));
body.push(P(
"where θ_{c} is the recorded docking angle of complex *c*, the circumflex denotes the prediction and the bar the "
+ "mean over complexes, so that *R*^{2} = 0 corresponds to predicting the mean angle. Mean absolute error is the "
+ "mean of the absolute differences between predicted and recorded angles."));

body.push(H2("2.11. Molecular dynamics"));
body.push(P(
"Simulations asked whether residue pairs to which TCR-IMAP assigns high pairwise scores, but that are non-contacts in "
+ "the crystal structure, form transient contacts in solution. The design, the system-selection criteria, the "
+ "simulation protocol and the primary statistic were committed to a dated document before any system was selected, "
+ "and the analysis definitions were appended, dated, before any occupancy was computed from a production trajectory. "
+ "It was predeclared that a low-occupancy outcome would be reported and that no system, seed, force field or cutoff "
+ "would be changed in search of a different result."));
body.push(P(
"Four systems were selected mechanically from the per-complex tcr70 differences in pairwise AUPRC between TCR-IMAP "
+ "and the pair-frequency prior, before any simulation: the two highest and the two lowest among complexes passing "
+ "predeclared filters (resolution ≤ 2.5 Å; distinct clusters; complete α1/α2 window; single peptide "
+ "of standard residues without covalent TCR–peptide linkage, engineered interchain disulfide or single-chain "
+ "construct). The selected systems were 5yxn and 3vxs (high difference) and 7n1f and 8eo8 (low difference), at "
+ "1.80–2.39 Å resolution. Each system contained the MHC heavy chain, β2-microglobulin, the peptide and "
+ "the TCR α and β chains as resolved; crystallographic waters, cryoprotectants, buffer ions and additives "
+ "were removed. Systems were built with CHARMM-GUI Solution Builder [[51,63]] using the CHARMM36m force field "
+ "[[50]] and the CHARMM-modified TIP3P water model [[52]] in a cubic box with a 10 Å buffer, neutralized at "
+ "150 mM KCl, and contained "
+ "359,000–416,000 atoms."));
body.push(P(
"Hydrogen masses were repartitioned, permitting a 4 fs timestep with LINCS constraints on bonds to hydrogen "
+ "[[54,58]]; electrostatics used particle-mesh Ewald [[55]]. After energy minimization, each replicate was "
+ "equilibrated for 125 ps under NVT at 300 K and for 1 ns under NPT at 300 K and 1 bar with position restraints, "
+ "followed by 200 ns of unrestrained production, using the velocity-rescaling thermostat and stochastic cell "
+ "rescaling [[56,57]]. Simulations were run with GROMACS 2023.2 [[53]]. Each system was simulated in three "
+ "replicates differing in initial velocity seed, giving 4 × 3 × 200 ns = 2,400 ns; all replicates were "
+ "pooled and no representative trajectory was selected [[59]]. Protein coordinates were written every 20 ps and the "
+ "first 20 ns of each production run was discarded, leaving 9,001 frames per replicate. Because CHARMM-GUI "
+ "topologies do not carry IMGT insertion codes, simulated residues were mapped to the label rows and columns by "
+ "sequence alignment and confirmed by coordinates after rigid superposition, with every TCR row and all 84 MHC "
+ "columns mapping to exactly one simulated residue."));
body.push(P(
"The occupancy of a residue pair was the fraction of retained frames in which the two residues' minimum heavy-atom "
+ "distance was at most 4.5 Å, the same definition as the crystal labels, computed under the minimum-image "
+ "convention with MDAnalysis [[60]],"));
body.push(EQ([
  msub("o", "ij"), mr(" = "),
  new MathFraction({ numerator: [msubsup("F", "ij", "contact")], denominator: [mr("F")] }),
]));
body.push(P(
"where *F* = 9,001 is the number of retained frames per replicate and *F*^{contact}_{ij} the number of those frames "
+ "in which *d*^{min}_{ij} ≤ 4.5 Å; system occupancy was the mean of the three replicate occupancies. The "
+ "predeclared primary statistic was the Spearman rank correlation between the TCR-IMAP score and occupancy among "
+ "crystal-negative pairs, computed per system and pooled, with uncertainty assessed by 10,000 within-system "
+ "permutations and labeled approximate because pairs are not independent. Every statistic was repeated with the "
+ "within-fold pair-frequency prior in place of TCR-IMAP scores; without this control, a positive correlation for "
+ "TCR-IMAP could not be distinguished from general positional contact propensity. A crystal-negative pair was "
+ "counted as meaningfully occupied at an occupancy of at least 0.25, and for each complex the *k* highest-scoring "
+ "crystal-negative pairs, with *k* equal to its number of crystal-positive pairs, were compared with the base rate "
+ "among all crystal-negative pairs. Interface characterization used backbone and interface RMSD, per-loop Cα "
+ "RMSF, buried solvent-accessible surface area computed by the Shrake–Rupley algorithm [[61]] and interface "
+ "hydrogen bonds. MD occupancy is force-field dependent and is not treated as ground truth; four systems constitute "
+ "a case-controlled analysis and not a benchmark."));

body.push(H2("2.12. Data and code availability"));
body.push(P(
"The STCRDab structures are publicly available. The frozen partitioning file, complex metadata, contact labels with row and column manifests, trained-model predictions, per-complex metric arrays and all analysis and figure scripts are prepared for public release. The code repository is https://github.com/Imrans-AI/TCR-IMAP (private during author review). No Zenodo deposit or archive DOI exists at the time of writing. A local release package has been assembled and tested in two parts — a code repository tree and a versioned data archive, both staged from the same source revision — containing the frozen partitioning file, the complex table, the contact labels with their row and column manifests, the stored per-complex and per-pair prediction arrays, every analysis script, the comparator reimplementations with a table recording every deviation from each published specification, an environment specification, SHA-256 checksums over every file and a worked example that reproduces the reported pairwise results from the stored arrays on a CPU in under a minute. Molecular dynamics input files, parameters, seeds, analysis scripts and the twelve output tables are included; the trajectories themselves are hundreds of gigabytes and are available from the corresponding author on reasonable request. Software: Python 3.10, PyTorch 2.10.0, NumPy 1.26.4, pandas 2.3.3, SciPy 1.15.2, Biopython 1.81, scikit-learn 1.7.2, Levenshtein 0.27.3 (over RapidFuzz 3.14.3), MDAnalysis 2.9.0, matplotlib 3.10.8, DSSP 4.2.2, GROMACS 2023.2 and CHARMM-GUI."));

/* ---- results ---- */
body.push(H1("3. Results"));

body.push(H2("3.1. Overview of TCR-IMAP and the homology-controlled evaluation framework"));
body.push(P(
"TCR-IMAP takes the amino-acid identities of the six CDR loops of both TCR chains and of a fixed structural window "
+ "on the MHC \u03b11/\u03b12 platform, and returns a score for every one of the 75 \u00d7 84 candidate residue "
+ "pairs (Figure 1B). Two encoders, a cross-attention block and a bilinear output head map the two sequences to the "
+ "contact map; no coordinates enter the model at training or prediction time. The inputs are nonetheless "
+ "position-indexed: each residue identity occupies the row or column that its IMGT number and insertion code "
+ "assign, so applying the model to a receptor outside this cohort requires that numbering step, which the "
+ "released package documents but does not itself implement. A per-residue variant of the same architecture "
+ "predicts, for each TCR "
+ "position, whether it contacts any window residue, which is the output format of the published comparators."));
body.push(P(
"The evaluation framework has three components, all fixed before the results below were produced. First, the 174 "
+ "complexes are clustered by TCR sequence and partitioned under five protocols, four of which can be ordered by how much "
+ "train/test similarity they permit, from "
+ "random folds of complexes to folds of whole TCR clusters (tcr70), with the fold assignments frozen. Second, the "
+ "primary metric (per-complex macro-AUPRC), the primary protocol (tcr70) and the primary comparison (TCR-IMAP "
+ "against Struct-TCS on the per-residue task) were declared in advance, and model selection is performed on "
+ "validation complexes drawn as whole clusters. Third, every method is measured against a within-fold frequency "
+ "prior fitted on training labels alone, so that any advantage is reported relative to positional contact "
+ "propensity rather than to chance. Uncertainty is estimated throughout by a bootstrap over the 69 TCR clusters, "
+ "because complexes sharing a cluster are not independent."));

body.push(H2("3.2. Sequence redundancy is pervasive and is not removed by deduplication"));
body.push(P(
"The benchmark comprised 174 TCR–peptide–MHC class I complexes, which formed 69 clusters at a TCR "
+ "similarity ratio of 0.70, and 82 and 89 clusters at ratios of 0.80 and 0.90 (Figure 2C). Of the "
+ "174 complexes, 100 (57.5 %) had an exact TCR twin elsewhere in the "
+ "set, and 123 of the 163 twin pairs (75.5 %) differed in peptide. V-gene usage was concentrated: TRAV12 occurred in "
+ "61 complexes (35.1 %) and TRBV6 in 49 (28.2 %) (Figure 2A). Redundancy was nearly invisible across pairs but "
+ "pervasive per complex: a similarity ratio of at least 0.90 held for 1.5 % of the 15,051 complex pairs but for the "
+ "nearest neighbor of 67.8 % of complexes (Figure 2B). Removing all 42 complexes belonging to duplicate families "
+ "left 132 complexes whose median nearest-neighbor similarity was still 0.967."));
body.push(P(
"The partitioning protocols differed accordingly in how close each test complex was to its training data (Figure 2D). "
+ "Under random partitioning the median test complex had a training complex at a similarity ratio of 0.983, and "
+ "60.3 % of test complexes, pooled over complexes, had one at 0.90 or above; under tcr70 the maximum over all 174 "
+ "test complexes was 0.696, and none reached 0.90. Over the three sequence definitions examined this bound spans "
+ "0.682–0.707, and no test complex exceeds 0.80 under any of them (Figure S9). Duplicate families were split "
+ "between training and test 7.4 times per fold under random partitioning and never under tcr70 (Figure 2E). The "
+ "benchmark is therefore dominated by near-duplicate TCRs, and of the five protocols only tcr70 keeps every test "
+ "complex below a similarity ratio of 0.90 to its training data."));
body.push(...figure(path.join(REPO, "figures/main/Figure_2.png"), "Figure 2",
"Sequence redundancy in the 174-complex benchmark and its leakage across partitioning protocols. "
+ "(A) V-gene usage: the eight most frequent TRAV and TRBV genes (most frequent, TRAV12 in 61 complexes and TRBV6 in "
+ "49). (B) Cumulative distributions of the TCR similarity ratio over all 15,051 complex pairs and for "
+ "each complex's nearest neighbor. A ratio ≥ 0.90 holds for 1.5 % of pairs but for the nearest neighbor of "
+ "67.8 % of complexes, and 100 of 174 complexes (57.5 %) have an exact TCR twin. (C) Complexes by cluster size at "
+ "70, 80 and 90 % similarity (69, 82 and 89 clusters). (D) For each test complex, the maximum similarity "
+ "ratio to any complex in its training folds, under each protocol. Under random partitioning the median test complex "
+ "has a training complex at 0.983 and 60.3 % have one at ≥ 0.90; under tcr70 the maximum is 0.696, and no test "
+ "complex reaches 0.90. (E) Duplicate families (identical CDR3α–CDR3β pair and peptide) split between "
+ "training and test, mean per fold: random 7.4, intersection (legacy) 0.67, allele 1.25, temporal 0 and tcr70 0. Complexes "
+ "evaluated under each protocol are given below each bar (random 174, intersection 30, allele 64, temporal 60 and tcr70 "
+ "174 of 174). The 108 complexes without an allele annotation fall outside every allele and intersection test fold, and "
+ "one allele fold of two complexes is below the minimum fold size. The similarity ratio is computed on the "
+ "concatenated CDR1, CDR2 and CDR3 slot sequences of both chains, one residue per IMGT position. Descriptive; no "
+ "test is performed. n = 174 complexes."));

body.push(H2("3.3. Apparent performance declines monotonically as train/test similarity is restricted"));
body.push(P(
"Apparent performance here denotes per-complex macro-AUPRC measured under a given partitioning scheme. In an "
+ "identity-threshold sweep on one fixed random five-fold partition, the pairwise advantage of TCR-IMAP over the "
+ "within-fold pair-frequency prior fell at every step, from +0.2750 (p < 0.001) without restriction to +0.0274 at a "
+ "All intervals and p values for this sweep use 20,000 bootstrap resamples at seed 42. The values first recorded "
+ "for this stage used the 2,000-resample default and gave [−0.0001, +0.0543], p = 0.052 at this threshold; the "
+ "point estimate is identical and the difference is Monte Carlo error in the bootstrap itself, which 2,000 "
+ "resamples do not resolve near an interval endpoint. The superseded values are retained in the Supporting "
+ "Information. The advantage at this threshold is small under either setting. "
+ "Over the "
+ "same range, TCR-IMAP's macro-AUPRC fell from 0.6313 to 0.3271 and the prior's from 0.3563 to 0.2997 (Figure 3A)."));
body.push(P(
"The test set was identical at every threshold: all 174 complexes were evaluated each time. What changed was the "
+ "training set, which shrank from 139–140 complexes per fold without restriction to 69–78 at 0.70; at that "
+ "threshold the train/test similarity actually achieved had a maximum of 0.699 and a mean of 0.649 (Figure S9). "
+ "TCR-IMAP's advantage over the pair-frequency prior fell roughly tenfold across the sweep while the prior's own "
+ "macro-AUPRC fell by 0.0566, so most of the decline reflects the model's loss. Statistical support weakened "
+ "progressively with homology control; at the strictest threshold the advantage was small and not significant."));
body.push(P(
"Restricting similarity also reduces the number of training complexes, from 139\u2013140 per fold without "
+ "restriction to 69\u201378 at the strictest threshold, so the sweep alone cannot separate homology from training "
+ "volume. A size-matched control separates them. For each of three thresholds and each fold, the same number of "
+ "training complexes the similarity filter had retained was instead drawn at random from the unrestricted pool, with "
+ "no similarity filter, on the same test complexes, over five draws with seeds fixed in advance and three training "
+ "seeds each. The principal control estimate is the mean over draws of the per-complex difference between the model "
+ "and the prior, with a cluster-bootstrap interval conditional on those five draws."));
body.push(P(
"Training on an average of 73.6 randomly chosen complexes per fold rather than 139\u2013140 left the advantage over "
+ "the pair-frequency prior at +0.2296 ([+0.1687, +0.2747]), against +0.2750 without restriction \u2014 a reduction "
+ "of about 0.0454. Reaching the same training size by removing the most similar complexes instead left it at "
+ "+0.0274 ([+0.0011, +0.0542]), a reduction of about 0.2476. Comparing the two regimes directly, per complex on the "
+ "same test complexes under the same cluster resamples, the similarity-filtered advantage lay below the size-matched "
+ "advantage by 0.1623 ([\u22120.2090, \u22120.1163]) at a threshold of 0.95, 0.1936 ([\u22120.2434, \u22120.1364]) "
+ "at 0.85 and 0.2021 ([\u22120.2455, \u22120.1471]) at 0.70, each p < 0.001, the gap widening as the threshold "
+ "tightened (Figure 3E). Smaller training sets alone therefore do not explain the decline. Across the five draws the "
+ "size-matched advantage ranged over 0.0170, 0.0129 and 0.0416 at the three thresholds; that spread is a range over "
+ "draws, not a confidence interval, because draws, training seeds and folds are not independent replicates."));
body.push(P(
"Two limits bound this result. Similarity filtering removes whole sequence families, so it changes the composition "
+ "and diversity of the training set as well as its redundancy with the test fold; the contrast measures the effect "
+ "of removing similar complexes and does not isolate redundancy from family composition. And the comparison against "
+ "the unrestricted value is a difference of point estimates, since the unrestricted per-complex values were not "
+ "regenerated; the quantity carrying an interval is the paired contrast between the two matched-size regimes. The "
+ "relative sizes of the two reductions are specific to this partition, model, filtering rule and set of draws, and "
+ "are not a general decomposition of the effect."));

body.push(...figure(path.join(REPO, "figures/main/Figure_3.png"), "Figure 3",
"Apparent performance depends on train/test redundancy. (A) Pairwise macro-AUPRC of TCR-IMAP and the "
+ "pair-frequency prior when every training complex whose TCR similarity ratio to any test complex of its "
+ "fold reaches the threshold is removed. One fixed random five-fold partition is used throughout, and all 174 "
+ "complexes are evaluated at every threshold. From no restriction to 0.70, TCR-IMAP falls by 0.3042 (0.6313 to "
+ "0.3271) and the prior by 0.0566 (0.3563 to 0.2997). (B) TCR-IMAP − prior with 95 % confidence intervals. The "
+ "advantage decreases at every step, from +0.2750 (p < 0.001) without restriction to +0.0274 at 0.70 "
+ "([+0.0011, +0.0542], p = 0.041). Filled markers, CI excludes zero. Intervals and p values use 20,000 bootstrap "
+ "resamples at seed 42; the 2,000-resample values first recorded for this stage gave [−0.0001, +0.0543], "
+ "p = 0.052 at this threshold from the same point estimate, a Monte Carlo difference in the bootstrap, and are "
+ "retained in the Supporting Information. (C) Mean absolute "
+ "error of the docking angles predicted by AngleNet, a separate sequence-based regression network trained from "
+ "scratch in each fold (not TCR-IMAP), under random, allele-holdout and TRBV-holdout partitions: three training "
+ "seeds each (points) and their mean (line). MAE 5.40–6.08° (random), 11.18–12.01° (allele "
+ "holdout) and 12.50–13.43° (TRBV holdout); dashed line, the median-predictor baseline (11.0°). Pearson "
+ "r across seeds is printed for each partition. (D) Docking-angle R² for the same runs: +0.42 to +0.56 under "
+ "random partitioning, −0.38 to −0.22 under allele holdout and −0.47 to −0.39 under TRBV holdout; "
+ "dashed line, R² = 0 (the mean predictor). (E) Training-size-matched control at the three thresholds for which "
+ "it was run. Restricting similarity also shrinks the training set, from 139–140 complexes per fold without "
+ "restriction to 91–96, 86–93 and 69–78 at thresholds of 0.95, 0.85 and 0.70, so the sweep alone cannot "
+ "separate homology from training volume. Each control draws that same number of training complexes at random from "
+ "the unrestricted pool, with no similarity filter, on the same test complexes; five draws with seeds fixed in "
+ "advance, three training seeds within each draw, 225 model fits in total. Squares, the similarity-filtered "
+ "advantage; small points, the five individual draws; large circles, the principal control estimate, which is the "
+ "mean over draws of the per-complex difference between the model and the prior and is conditional on those five "
+ "realized draws. The size-matched advantage is +0.2400 [+0.1827, +0.2832], +0.2364 [+0.1748, +0.2822] and "
+ "+0.2296 [+0.1687, +0.2747] against the filtered +0.0777, +0.0428 and +0.0274; the paired contrast between the "
+ "two regimes, computed per complex on the same test complexes under the same cluster resamples, is "
+ "−0.1623 [−0.2090, −0.1163], −0.1936 [−0.2434, −0.1364] and −0.2021 [−0.2455, −0.1471], all "
+ "p < 0.001 and printed beside each pair. Dotted line, the unrestricted advantage (+0.2750). Filled markers, 95 % "
+ "CI excludes zero. Error bars are cluster-bootstrap intervals; the five draws are plotted as points because their "
+ "spread is a range over draws and not a confidence interval, and the two are never pooled. Detailed subset "
+ "memberships, achieved similarities and per-draw values are in Supporting Table S9. Statistics: (B) and (E) "
+ "paired cluster bootstrap over 69 TCR clusters, 20,000 resamples, seed 42. n = 174 complexes (A, B, E); "
+ "172 of 174 (C, D)."));

body.push(H2("3.4. Homology between validation and training complexes inflated model selection"));
body.push(P(
"The per-residue models — TCR-IMAP's per-residue variant, Struct-TCS and Seq-TCS — were originally "
+ "early-stopped on a random 15 % of the non-test complexes, drawn without regard to TCR cluster. A validation "
+ "subset drawn at random over complexes is not independent of the complexes left in training when 57.5 % of "
+ "complexes have an exact TCR twin. Under that split, "
+ "5–15 clusters (median 10.5) had complexes on both sides of model selection in every one of the 18 evaluated "
+ "folds (Figure S4C). Drawing validation complexes as whole clusters instead reduced TCR-IMAP's per-residue "
+ "macro-AUPRC by 0.0541 under random partitioning but by only 0.0165 under tcr70 (Figure S4A). The inflation was "
+ "thus largest where train/test redundancy was largest."));
body.push(P(
"The change comes from the validation procedure rather than from nondeterminism. The position prior and "
+ "nearest-neighbor transfer, which use no validation subset, were identical to the last bit under both splits; the "
+ "deterministic-algorithm settings changed no value in any of the 15 protocol \u00d7 method cells; and both "
+ "configurations reproduced byte for byte. It cannot, however, be attributed to cluster disjointness alone: changing "
+ "the procedure also changes which complexes are held out, the validation fraction, the number of complexes left to "
+ "optimize on and which checkpoint early stopping selects, and these move together by construction. Drawing "
+ "validation as whole clusters until the 15 % target is met overshoots it, giving an achieved validation fraction of "
+ "15.1\u201320.2 % of each fold\u2019s non-test complexes against 7.2\u201317.1 % under the random draw, and leaving "
+ "91\u2013145 complexes for optimization against 96\u2013149. The correction removed statistical significance from "
+ "the predeclared primary comparison, per-residue TCR-IMAP versus Struct-TCS under tcr70, which moved from +0.0504 "
+ "(p = 0.001) under the random split to +0.0312 (p = 0.088) (Figure S4B). The primary comparison was not "
+ "re-designated after this result. A comparison that is significant under one procedure and not under the other is "
+ "not itself evidence that the procedures differ significantly; that is a separate question. As an additional "
+ "analysis, the stored per-complex results allow the two procedures to be compared directly, paired by complex under "
+ "the same cluster bootstrap. Correcting the procedure reduced TCR-IMAP\u2019s macro-AUPRC by 0.0541 "
+ "([\u22120.0734, \u22120.0316], p < 0.001) under random partitioning and by 0.0165 ([\u22120.0317, \u22120.0022], "
+ "p = 0.024) under tcr70, while Struct-TCS did not change significantly under tcr70 (+0.0027, [\u22120.0179, "
+ "+0.0251], p = 0.846). Under random partitioning all three learned models fell significantly, by 0.0444 to 0.0641, "
+ "so the inflation was not specific to one architecture. The same redundancy that inflates test performance "
+ "therefore also inflated model selection within our own pipeline."));

body.push(H2("3.5. Under homology control, the predeclared per-residue comparison is unresolved while the pairwise comparisons are supported"));
body.push(P(
"The predeclared primary comparison, per-residue macro-AUPRC under tcr70, gave TCR-IMAP 0.7281 and Struct-TCS "
+ "0.6969, a difference of +0.0312 (95 % CI [−0.0049, +0.0650], p = 0.088) (Figure 4B). This is insufficient "
+ "evidence of superiority: the interval spans values from a disadvantage of 0.0049 to an advantage of 0.0650. "
+ "TCR-IMAP exceeded Seq-TCS by +0.0799 ([+0.0455, +0.1161], p < 0.001) and the position prior by +0.0085 "
+ "([−0.0162, +0.0371], p = 0.546); a p value of 0.546 does not establish equivalence with the prior. The allele "
+ "and legacy intersection results below rest on 64 and 30 of the 174 complexes, for the reason given in Section 2.3. "
+ "The allele protocol withholds allele labels rather than MHC sequences, and the intersection partition is not a "
+ "double holdout; both qualifications are given in Section 2.3, and the intersection value is reported for "
+ "completeness only. Under the "
+ "other protocols the difference from Struct-TCS was −0.0030 (random, p = 0.798), +0.0696 (legacy intersection "
+ "partition, p = 0.015), "
+ "+0.0507 (allele, p = 0.001) and +0.0030 (temporal, p = 0.874) (Figure 4A, B)."));
body.push(P(
"The predeclared comparison is also sensitive to an implementation choice in the comparator itself. Our Struct-TCS "
+ "applies a per-head gate to the graph attention that the original description does not specify and the released "
+ "code does not use (Section 2.9). Retraining Struct-TCS under tcr70 without that gate, with everything else held, "
+ "gave a macro-AUPRC of 0.6927 and a difference from TCR-IMAP of +0.0354 ([+0.0021, +0.0668], p = 0.036), against "
+ "+0.0312 (p = 0.088) with it. The two comparator variants are not themselves distinguishable: paired over the same "
+ "complexes under the same resamples they differ by \u22120.0042 ([\u22120.0161, +0.0070], p = 0.441). A change to "
+ "the comparator too small to detect therefore moves the headline verdict across the conventional significance "
+ "threshold. We report the gated implementation, which is the one described in Section 2.9 and the one the "
+ "predeclaration refers to, and treat this as evidence that the comparison is fragile at this sample size rather "
+ "than as grounds for selecting whichever variant yields significance."));
body.push(P(
"On the secondary, pairwise task under tcr70, TCR-IMAP (macro-AUPRC 0.3598) exceeded the within-fold pair-frequency "
+ "prior (0.3090) by +0.0507 ([+0.0215, +0.0841], p < 0.001) and was higher on 112 of 174 complexes (Figure 4C). It "
+ "exceeded TEIM-Res (0.2725) by +0.0873 ([+0.0592, +0.1149], p < 0.001) and was higher on 143 of 174 complexes "
+ "(Figure 4D). Weighting the 69 clusters equally instead of the complexes gave the same verdicts (Figures 4E, S5): "
+ "the per-residue comparisons with Struct-TCS (+0.0258, p = 0.063) and with the position prior (+0.0164, p = 0.174) "
+ "remained nonsignificant, and the pairwise comparison with the prior remained significant (+0.0360, p = 0.001)."));
body.push(P(
"On the 88 structures shared with the Struct-TCS dataset, TCR-IMAP led under tcr70 at both of that study's contact "
+ "cutoffs: +0.1057 at 4.5 Å and +0.0964 at 6 Å (both p < 0.001; n = 88) (Figure 4F). Under a single random "
+ "split of the kind used in that study (22 test complexes), the difference was not significant at 4.5 Å "
+ "(+0.0269, p = 0.171) and reversed at 6 Å (−0.0949 [−0.1335, −0.0517], p < 0.001). The reversal "
+ "held in each of the three training seeds (−0.060, −0.120 and −0.087), but seed-to-seed variation within a "
+ "fixed split does not measure variation between splits, and its magnitude varied more than in any other cell "
+ "(Figure S6). Because these two cells differ in the partitioning protocol and also in which complexes are scored, "
+ "the tcr70 evaluation was restricted to exactly the 22 complexes the single random split had evaluated, which are a "
+ "subset of the 88. On those identical complexes the difference from Struct-TCS was +0.0912 ([+0.0298, +0.1514], "
+ "p = 0.001) at 4.5 Å and +0.0735 ([+0.0172, +0.1233], p = 0.010) at 6 Å, so the 6 Å reversal follows the "
+ "partitioning protocol and not the choice of test complexes. The 88-complex results above are unchanged. Because "
+ "the single-split cell rests on one fold of 22 complexes, it indicates that the ordering of the two methods is "
+ "unstable under the comparator's own single-split design, not that Struct-TCS is the better method. Only the "
+ "structure subset and cutoff are shared with that study, so this is not a reproduction of its evaluation."));
body.push(...figure(path.join(REPO, "figures/main/Figure_4.png"), "Figure 4",
"TCR-IMAP against published comparators. (A) Per-residue macro-AUPRC of TCR-IMAP, Struct-TCS, Seq-TCS, the "
+ "within-fold position prior and nearest-neighbor transfer under five partitioning protocols, the four orderable ones "
+ "arranged by how much train/test similarity they permit, "
+ "with the number of complexes evaluated. Each protocol is annotated with the TCR-IMAP and Struct-TCS means: random "
+ "0.806 and 0.809, intersection (legacy) 0.800 and 0.730, allele 0.776 and 0.726, temporal 0.694 and 0.691, tcr70 0.728 and 0.697. "
+ "(B) Paired difference TCR-IMAP − Struct-TCS per protocol with 95 % confidence intervals: random −0.0030 "
+ "(p = 0.798), intersection (legacy) +0.0696 (p = 0.015), allele +0.0507 (p = 0.001), temporal +0.0030 (p = 0.874) and tcr70, the "
+ "predeclared primary comparison, +0.0312 [−0.0049, +0.0650] (p = 0.088), which is not significant. (C, D) "
+ "Pairwise contact prediction under tcr70, one point per complex; dashed, identity line. (C) Against the within-fold "
+ "pair-frequency prior: TCR-IMAP higher on 112 of 174 complexes (64.4 %), pooled Δ = +0.0507 [+0.0215, +0.0841], "
+ "p < 0.001. (D) Against TEIM-Res: higher on 143 of 174 (82.2 %), pooled Δ = +0.0873 [+0.0592, +0.1149], "
+ "p < 0.001. (E) Cluster-macro sensitivity, tcr70: the three primary comparisons with complexes weighted equally "
+ "(the predeclared weighting) and with the 69 clusters weighted equally. The two weightings agree in every "
+ "comparison. (F) The 88 structures shared with the comparator's dataset, at its 4.5 Å and 6 Å cutoffs, "
+ "under a single random split of the comparator's design (n = 22 test complexes) and under tcr70 (n = 88). Only the "
+ "structure subset and cutoff are shared, so this is not a reproduction of the comparator's own evaluation. All "
+ "per-residue models are selected on cluster-disjoint validation. Statistics: paired cluster bootstrap over the 69 "
+ "TCR clusters, 20,000 resamples, seed 42, in every panel; filled markers, 95 % CI excludes zero."));

body.push(H2("3.6. Docking-angle prediction falls below the mean predictor under grouped evaluation"));
body.push(P(
"A separate sequence-based regression network (AngleNet), trained from scratch in each fold and distinct from "
+ "TCR-IMAP, predicted the TCR docking angle for the 172 complexes with a recorded angle (Figure 3C, D). Under random "
+ "partitioning it explained part of the variance (R² +0.420 to +0.556 across three seeds; mean absolute error "
+ "5.40–6.08°). Under allele holdout R² fell to −0.382 to −0.216 "
+ "(11.18–12.01°), and under TRBV holdout to −0.468 to −0.394 (12.50–13.43°). A negative "
+ "R² means that the predictions were worse than predicting the mean angle, and the grouped mean absolute errors "
+ "also exceeded that of the median predictor (11.0°). A docking-angle regressor that performs well under random "
+ "partitioning therefore performs worse than a constant under grouped evaluation."));

body.push(H2("3.7. The observed contact footprint concentrates on the α1/α2 helices and avoids the β-sheet floor"));
body.push(P(
"This analysis is descriptive and involves no model and no data split. The 84-column α1/α2 window captured "
+ "5,281 of the 5,294 contacts observed across the 182-position G-domain (99.75 %) (Figure 5D). Of these, 88.3 % lay "
+ "on helix positions, and the β-sheet floor carried no contact in 595,952 resolved pairs. CDR2β contacted "
+ "no α2-helix position in 40,560 resolved pairs. Three CDR3β insertion rows were resolved in fewer than "
+ "five complexes and were excluded from frequency displays and from the identification of frequently contacted "
+ "positions, because a frequency over one to three complexes is not a population estimate. With them excluded, the "
+ "most frequent pair was CDR1β 36 × α1 76 (frequency 0.667), and 14 MHC positions exceeded 0.30 of the "
+ "maximum summed frequency (Figure 5A, B). Each loop preferred one helix over the 69 clusters (Figure 5E): "
+ "CDR1α, CDR2α and CDR3β preferred α2, and CDR3α, CDR1β and CDR2β preferred "
+ "α1, with Holm-adjusted p = 0.003 for every loop. The footprint therefore concentrates on the α1/α2 "
+ "helices, places the germline loops of the two chains over opposite helices, and does not reach the β-sheet "
+ "floor."));
body.push(...figure(path.join(REPO, "figures/main/Figure_5.png"), "Figure 5",
"Observed structural contact footprint on the MHC α1/α2 platform. (A) Contact frequency (contacts / "
+ "complexes in which both residues are resolved; minimum heavy-atom distance ≤ 4.5 Å) for 75 TCR rows "
+ "× the 84-column α1/α2 window. The strip marks DSSP helix (dark) and non-helical boundary (light) "
+ "columns. Boxes mark the helix each loop prefers (see E); solid outline, germline loops; dashed, CDR3; color, TCR "
+ "chain. Three CDR3β insertion rows resolved in fewer than five complexes are hatched and kept off the color "
+ "scale. With these rows excluded, the most frequent pair is CDR1β 36 × α1 76 (0.667). Square-root "
+ "color scale. (B) The MHC platform of 1ao7 colored by each position's contact frequency summed over the TCR rows "
+ "resolved in at least five complexes, relative to its own maximum. Sticks mark the 14 positions above 0.30 of the "
+ "maximum. (C) The 1ao7 TCR–pMHC complex with its CDR loops. (D) Contact rate (%) by loop and MHC region over "
+ "the 182-position G-domain. The β-sheet floor carries 0 contacts in 595,952 resolved pairs; the 84-column "
+ "window captures 5,281 of the 5,294 G-domain contacts (99.75 %). (E) Per-cluster α1 − α2 contact-rate "
+ "difference with 95 % bootstrap intervals over the 69 TCR clusters. α-chain germline loops prefer α2 and "
+ "β-chain germline loops α1; CDR3α prefers α1 and CDR3β α2. CDR2β contacts no "
+ "α2-helix position (0 of 40,560 resolved pairs). Statistics: (E) cluster bootstrap, 2,000 resamples, seed 42; "
+ "after Holm correction across the six loops every loop has p = 0.003. Descriptive: no model and no data split. "
+ "n = 174 complexes."));

body.push(H2("3.8. Residue shuffling and germline-loop removal reduce performance after multiplicity correction"));
body.push(P(
"Of eight ablations of the pairwise model under tcr70 (full model 0.3598), three reduced macro-AUPRC after Holm "
+ "correction, with differences given as full − ablated: shuffling TCR residues within loops (+0.0503, Holm "
+ "p < 0.001), removing the germline CDR1 and CDR2 loops (+0.0461, Holm p < 0.001) and shuffling MHC residues across "
+ "positions (+0.0283, Holm p = 0.014) (Figure 6A). Masking TCR residue identity (+0.0337, Holm p = 0.104) or MHC "
+ "residue identity (+0.0180, Holm p = 0.160) did not survive correction, nor did removing CDR3 (+0.0121, Holm "
+ "p = 0.483). Removing cross-attention (−0.0008) or replacing the bilinear head with a dot product "
+ "(−0.0039) showed no detectable effect after correction (Holm p = 1.000 for both). None of the five ablations "
+ "that failed to survive correction is evidence that the component or the information it removes is unnecessary: "
+ "with 174 complexes in 69 clusters and a correction across eight tests, an effect of the size of the surviving "
+ "three could go undetected. What the three surviving ablations establish is that shuffling residues across "
+ "positions on either side of the interface, and removing the germline CDR1 and CDR2 loops, reduce pairwise "
+ "macro-AUPRC. This does not show that the model learns physical residue–residue interactions."));
body.push(...figure(path.join(REPO, "figures/main/Figure_6.png"), "Figure 6",
"Sequence dependence of TCR-IMAP, and the complexes on which it gains and loses most against the pair-frequency "
+ "prior (tcr70, pairwise task). (A) Change in pairwise macro-AUPRC when one input or architectural component is "
+ "removed (full − ablated; full TCR-IMAP 0.3598). Bars are Bonferroni-simultaneous 99.375 % intervals across "
+ "the eight ablations; raw and Holm-adjusted p values are given beside each row. Three ablations lowered performance "
+ "after correction: TCR residues shuffled within loops (+0.0503, Holm p < 0.001), CDR1/2 removed (+0.0461, Holm "
+ "p < 0.001) and MHC residues shuffled across positions (+0.0283, Holm p = 0.014). Masking TCR residue identity "
+ "(+0.0337, Holm p = 0.104) or MHC residue identity (+0.0180, Holm p = 0.160) did not survive correction. Removing "
+ "cross-attention (−0.0008) or replacing the bilinear head with a dot product (−0.0039) had no detectable "
+ "effect. (B) Per-complex difference between TCR-IMAP and the within-fold pair-frequency prior, ranked: TCR-IMAP is "
+ "higher on 112 of 174 complexes (64.4 %); pooled Δ = +0.0507 [+0.0215, +0.0841], p < 0.001. (C, D) Observed "
+ "contacts (left) and TCR-IMAP's k highest-scoring residue pairs (right), where k is the complex's number of "
+ "observed contacts, from out-of-fold tcr70 predictions. View onto the MHC α1/α2 platform from the TCR. "
+ "(C) 3qfj, the largest gain over the prior (pairwise AUPRC 0.566 vs 0.287). (D) 3qeq, the largest loss to it "
+ "(0.293 vs 0.452). Overlap counts equal precision at k for these two complexes and are illustrative, not a "
+ "performance estimate. Statistics: paired cluster bootstrap over 69 TCR clusters, 20,000 resamples, seed 42. "
+ "n = 174 complexes."));

body.push(H2("3.9. In molecular dynamics the frequency prior ranks solution-phase contacts better than the model"));
body.push(P(
"Four complexes selected before any simulation for a high (5yxn, 3vxs) or low (7n1f, 8eo8) per-complex gain of "
+ "TCR-IMAP over the prior were simulated in three 200 ns replicates each, and occupancy was measured over "
+ "20–200 ns (Figure 7). High-scoring crystal-negative pairs were enriched for solution-phase contact: among the "
+ "k highest-scoring crystal-negative pairs (k = 22–36, each complex's crystal-positive count), 1–5 per "
+ "system reached occupancy ≥ 0.25 for TCR-IMAP and 1–7 for the prior, a fraction of 3.3–16.1 % and "
+ "3.3–22.6 % against base rates of 0.12–0.26 %; this is a 28- to 88-fold enrichment, but it rests on "
+ "1–7 pairs per system (Figure 7C)."));
body.push(P(
"Over all 19,957 crystal-negative pairs, of which 381 had any contact, occupancy correlated more strongly with the "
+ "prior than with TCR-IMAP (Spearman ρ 0.441 versus 0.216), in every system and in every replicate (Figure 7A). "
+ "The difference was driven by the > 8 Å stratum of crystal-frame distance (0.208 versus 0.103), which held "
+ "19,414 of the pairs, only 94 of which had any contact and none of which reached 0.25. The two scores were similar "
+ "at 6–8 Å (0.309 versus 0.303; 423 pairs). At 4.5–6 Å (120 pairs) TCR-IMAP was numerically "
+ "higher (0.241 versus 0.179), but the difference between the scores was not tested. The within-system permutation p "
+ "values are not quoted: every system-level and pooled value is at the test's floor, and because the "
+ "4,598–5,429 crystal-negative pairs of each system share residues and are not independent, a value at the "
+ "floor carries no information beyond the sign of ρ. The two scores were only moderately correlated with each "
+ "other (Spearman ρ 0.33–0.39 per system, 0.36 pooled), so they are genuinely different rankings."));
body.push(P(
"Crystal-positive pairs retained contact moderately: system mean occupancy was 0.33–0.60 (Figure 7D). In 5yxn "
+ "the TCR disengaged in all three replicates, retaining 21.0 %, 29.6 % and 34.4 % of its starting TCR–pMHC "
+ "contacts at 200 ns against 57.6–85.2 % in the other systems, but the chain that lost contact differed: "
+ "TCRβ in replicates 1 and 2, TCRα in replicate 3. Backbone RMSD stayed within 1.5–7.7 Å in every "
+ "replicate, so this was a change of interface rather than a loss of fold (Figure 7E, Figure S3). These simulations "
+ "do not show that TCR-IMAP's crystal-negative predictions are transient contacts; across all pairs the frequency "
+ "prior ranks solution-phase occupancy better."));
body.push(...figure(path.join(REPO, "figures/main/Figure_7.png"), "Figure 7",
"Solution-phase contact of crystal-negative residue pairs, against TCR-IMAP scores and the pair-frequency prior. "
+ "Four TCR–pMHC complexes, selected before any simulation for a high (5yxn, 3vxs) or low (7n1f, 8eo8) "
+ "per-complex gain of TCR-IMAP over the prior, were simulated in three independent 200 ns replicates each. The "
+ "occupancy of a residue pair is the fraction of the 9,001 frames between 20 and 200 ns in which its minimum "
+ "heavy-atom distance is ≤ 4.5 Å, averaged over the three replicates. (A) Spearman ρ between score and "
+ "occupancy over crystal-negative resolved pairs, grouped by the pair's minimum heavy-atom distance in the crystal "
+ "structure, with the four systems pooled. Over all pairs (n = 19,957, 381 non-zero), occupancy correlates more "
+ "strongly with the prior than with TCR-IMAP (ρ = 0.44 vs 0.22), in every system and every replicate. The "
+ "difference arises in the > 8 Å stratum (0.21 vs 0.10), which holds 19,414 pairs, of which 94 have any contact "
+ "and none reaches 0.25. The two scores are similar in the 6–8 Å stratum (0.31 vs 0.30; 423 pairs). In the "
+ "4.5–6 Å stratum (120 pairs) TCR-IMAP is numerically higher (0.24 vs 0.18), but no test of the difference "
+ "between scores was made. (B) Score against occupancy for every crystal-negative pair, colored by stratum. "
+ "(C) Fraction of the k highest-scoring crystal-negative pairs with occupancy ≥ 0.25 against the base rate, per "
+ "system; k is the system's crystal-positive count (22–36). Both scores enrich high-occupancy pairs 28- to "
+ "88-fold over the base rate, on 1–7 pairs per system. (D) Occupancy of the crystal-positive pairs per "
+ "replicate and per system. Above each replicate, the percentage of frames that keep at least 50 % of the starting "
+ "TCR–pMHC contacts, read two ways: by contact count and by contact identity. (E) Backbone RMSD over "
+ "20–200 ns for each replicate; dashed line, "
+ "system mean (5yxn 4.25, 7n1f 3.98, 3vxs 3.12 and 8eo8 3.80 Å). Backbone RMSD stays between 1.5 and 7.7 Å "
+ "in every replicate: that is a change of interface, not a loss of fold."));

body.push(H2("3.10. Adding the peptide sequence has no detectable effect"));
body.push(P(
"Adding the peptide sequence to TCR-IMAP changed pairwise macro-AUPRC under tcr70 by +0.0100 over all positions "
+ "([−0.0058, +0.0258], p = 0.196), +0.0071 with the peptide shuffled ([−0.0072, +0.0222], p = 0.308), "
+ "+0.0107 over CDR1/2 rows (p = 0.133) and +0.0082 over CDR3 rows ([−0.0137, +0.0293], p = 0.416) (Figure S1). "
+ "Every 95 % confidence interval included zero, so no effect of the peptide input was detected in any scope."));

/* ---- discussion ---- */
body.push(H1("4. Discussion"));
body.push(P(
"Interface-level prediction is attractive for T-cell recognition because the contacts carry the mechanism, and "
+ "because structures exist for only a small fraction of the repertoire [[7,9]]. Both routes to closing that gap are "
+ "advancing: co-folding models now build complete TCR\u2013pMHC complexes [[12,16,17]], and sequence models recover "
+ "residue-level interaction patterns and contact maps [[25,26,27]]. The measurement has not kept pace. Structural "
+ "datasets for this task are families of near-identical receptors, and the comparator study for contact prediction "
+ "used a single partition and could not exclude that the CDR1 and CDR2 sequences of its test complexes had been "
+ "seen during training [[27]]. Where this has been examined in the adjacent specificity task, evaluation design "
+ "accounted for a substantial part of the reported performance [[28,29,31]]."));
body.push(P(
"We therefore built TCR-IMAP as an instrument for measuring that dependence rather than as an end in itself: a "
+ "dual-encoder cross-attention network predicting the complete contact map from amino-acid identities alone, "
+ "evaluated under five protocols frozen in advance, with a declared primary metric, protocol and comparison, and "
+ "uncertainty estimated by a bootstrap over TCR clusters rather than over complexes."));
body.push(P(
"Four analyses converge on the same conclusion. Restricting the maximum train/test similarity ratio while holding "
+ "the test set fixed reduced the advantage over a within-fold frequency prior at every step, from +0.2750 to "
+ "+0.0274. A size-matched control shows that smaller training sets alone do not explain that decline: removing the "
+ "same number of training complexes at random left the advantage at +0.2296, and the paired contrast between the "
+ "two matched-size regimes was −0.2021 ([−0.2455, −0.1471], p < 0.001). Similarity filtering also removes "
+ "whole sequence families, so it changes the composition of the training set as well as its redundancy with the "
+ "test fold, and the contrast should be read as the effect of removing similar complexes rather than as a "
+ "decomposition into separable causes. "
+ "Drawing model-selection complexes as whole clusters rather than at random cost 0.0541 of per-residue macro-AUPRC "
+ "under random partitioning against 0.0165 under cluster-disjoint partitioning, and moved the predeclared "
+ "comparison from p = 0.001 to p = 0.088. A docking-angle regressor explained part of the variance under random "
+ "partitioning and fell below the mean predictor under grouped partitions, placing the same dependence on a scalar "
+ "target rather than a contact label. On the 88 structures shared with the comparator, the ordering of the two "
+ "methods reversed at the 6 \u00c5 cutoff when the comparator\u2019s single-split design was used, which indicates "
+ "that this ordering is unstable under that design rather than that the comparator is the better method. Three of "
+ "the four analyses were measured on our own pipeline. Together they establish that apparent accuracy is inflated "
+ "under the evaluated protocols on this benchmark; they do not measure how far the same inflation extends beyond "
+ "it."));
body.push(P(
"Under homology control the pairwise comparisons resolve and the predeclared per-residue comparison does not. "
+ "TCR-IMAP exceeded the pair-frequency prior by +0.0507 and the reimplemented TEIM-Res by +0.0873, each at "
+ "p < 0.001, with the same verdicts when clusters rather than complexes were weighted equally. The predeclared "
+ "comparison with Struct-TCS gave +0.0312 (95 % CI [\u22120.0049, +0.0650], p = 0.088): the interval runs from a "
+ "disadvantage to an advantage more than twice the estimated difference, so it is unresolved in either direction at "
+ "this sample size. The primary comparison was not re-designated afterwards and the pairwise task was not promoted "
+ "in its place, since selecting the primary analysis after the outcome is the error examined here. Two "
+ "architectural components, cross-attention and the bilinear output head, were not shown to be necessary, and "
+ "adding the peptide sequence produced no detected effect in any scope. The comparator that recurs throughout is a "
+ "frequency prior rather than a network: the per-residue difference from the position prior was +0.0085 "
+ "(p = 0.546), and in molecular dynamics the pair-frequency prior ranked solution-phase occupancy of "
+ "crystal-negative pairs above TCR-IMAP overall (Spearman \u03c1 0.441 against 0.216), the difference arising in "
+ "pairs more than 8 \u00c5 apart in the crystal frame. Those simulations neither validate nor refute the model; "
+ "they indicate that a positional frequency prior already captures much of what separates contacting from "
+ "non-contacting positions on this platform, and that it is the baseline methods of this class must clear."));
body.push(P(
"Several limitations bound these conclusions. The benchmark is 174 complexes in 69 clusters and every interval "
+ "resamples those clusters. The predeclared per-residue difference is +0.0312 with a 95 % confidence interval of "
+ "[−0.0049, +0.0650], an interval that spans values from a disadvantage of 0.0049 to an advantage of 0.0650; at "
+ "this sample size the comparison is unresolved, and we do not quantify a threshold below which a difference would "
+ "be undetectable, since no power analysis supporting such a bound was performed. The tcr70 protocol clusters on TCR sequence, so it constrains TCR similarity by construction "
+ "while MHC and peptide similarity are not quantities it partitions on. MHC redundancy is addressed instead by the "
+ "allele protocol, which evaluates 64 of the 174 complexes and holds out alleles by their parsed annotation rather "
+ "than by MHC sequence: the 84-residue window takes only 28 distinct values across the cohort, and for 32 of those "
+ "64 complexes an identical window is present in training, in every case contributed by a complex whose allele "
+ "annotation could not be parsed. Allele holdout therefore does not guarantee that the MHC window itself is "
+ "unseen. Per-complex AUPRC correlates with a complex\u2019s "
+ "number of contacts and is not a measure of difficulty, which is why every comparison is a paired difference "
+ "between methods on the same complexes [[39]]. The molecular dynamics comprises four complexes in three replicates "
+ "each and is a case study rather than a benchmark. The scope is human MHC class I and the \u03b11/\u03b12 "
+ "platform, and sequence identity throughout is an edit-distance similarity rather than gap-free alignment "
+ "identity."));
body.push(P(
"Two things the work does not support should be stated plainly: it does not establish that TCR-IMAP predicts "
+ "contacts accurately enough to guide the engineering of TCR sequences, no design experiment having been performed "
+ "here, and on the per-residue task it does not establish an advantage over a positional frequency prior. What it "
+ "does provide is the measurement apparatus: five partitioning protocols with frozen fold assignments, so that a "
+ "new method is evaluated on the same partitions reported here; a benchmark of 174 complexes whose redundancy is "
+ "characterized by cluster count, exact twins and per-complex nearest-neighbor similarity; within-fold "
+ "frequency-prior baselines; and a model-free contact footprint of the \u03b11/\u03b12 platform that uses no model "
+ "and no data split. These serve three purposes. They give sequence-based contact predictors, and the interface "
+ "predictions of co-folding systems [[12,16,26]], a declared and redundancy-controlled protocol under which to be "
+ "compared. They provide a template for evaluation design in adjacent tasks, where similarity between partitions "
+ "has the same effect and where community challenges have reached the same conclusion [[28,37]]. And the footprint "
+ "itself, being independent of any model, is directly usable in structural analysis of TCR\u2013pMHC interfaces."));

body.push(H1("5. Conclusions"));
body.push(P(
"Apparent accuracy in sequence-based TCR\u2013MHC contact prediction is a joint property of the method and the "
+ "evaluation design. On this benchmark, controlling train/test sequence similarity, making the model-selection "
+ "split cluster-disjoint, grouping the partitions of a docking-angle regressor and re-running a published "
+ "single-split design each change the measured result, in one case reversing the ordering of two methods. Under "
+ "homology control TCR-IMAP exceeds a within-fold frequency prior and a reimplemented residue-level comparator on "
+ "the pairwise task, while the predeclared per-residue comparison remains unresolved at this sample size. The "
+ "frozen protocols, the characterized benchmark, the frequency-prior baselines and the model-free footprint are "
+ "released so that the conditionality of method rankings in this task is visible and testable."));

body.push(H1("Data and Code Availability"));
body.push(P(
"A release package containing all code, the frozen partitioning file, the complex table, the contact labels with their row and column manifests, the stored per-complex and per-pair prediction arrays, the comparator reimplementations, an environment specification, SHA-256 checksums over every file and a worked example has been assembled and tested from a clean copy outside the working repository. The code and data packages are staged from the same source revision. The code repository is https://github.com/Imrans-AI/TCR-IMAP (private during author review); a public release and a versioned Zenodo archive are planned. No archive DOI has been assigned. The worked example reproduces the reported pairwise results from the stored arrays on a CPU in under a minute; regenerating those predictions by retraining requires a GPU and is documented separately. The structures analyzed are publicly available from STCRDab [[9]] and the Protein Data Bank. Four cleaned coordinate files (3VXS, 5YXN, 7N1F and 8EO8) are included as MD inputs; the remaining source coordinates are referenced through retrieval instructions. The comparator reimplementations and a provenance table documenting adaptations, execution defects and unavailable components are included; the original authors’ repositories are not redistributed. Molecular dynamics input files, parameters, seeds, analysis scripts and the twelve output tables are included, so the occupancy analysis can be repeated independently; the trajectories are hundreds of gigabytes and are available from the corresponding author on reasonable request. Original code is licensed under MIT; original derived data and documentation are licensed under CC BY 4.0. Third-party material retains its own terms."));

body.push(H1("Funding and Acknowledgments"));
body.push(P("D.-Q.W. is supported by the National Key R&D Program of China (2024YFA1306400 and 2023YFE0199200), the National Natural Science Foundation of China (32570803 and 32030063), and the Intergovernmental International Scientific and Technological Innovation Cooperation Program. The computations were performed in part at Pengcheng Laboratory and the Center for High-Performance Computing, Shanghai Jiao Tong University."));

body.push(H1("Supporting Information"));
body.push(P(
"Supporting Information is provided as a separate document and contains Supporting Figures S1\u2013S9 (peptide "
+ "input; per-complex AUPRC distributions; interface characterization of the simulated systems; validation-split "
+ "sensitivity; cluster-weighting sensitivity; the comparator subset per training seed; legacy results documenting "
+ "the candidate-set correction; per-allele and per-TRBV breakdown; and the composition of the identity-threshold "
+ "sweep) together with Supporting Tables S1\u2013S10, which include the training-size-matched control's estimands, "
+ "diagnostics and subset memberships, the comparator subset evaluated on a common set of test complexes, and the "
+ "canonical-CDR-only scoring sensitivity."));

/* ---- references ---- */
body.push(H1("References", { pageBreakBefore: true }));
body.push(...refParagraphs);

/* helper used above (declared late is fine for function decls) */
function msubsup(b, sub, sup) {
  return new MathSubSuperScript({ children: [mr(b)], subScript: [mr(sub)], superScript: [mr(sup)] });
}
function msubhat(b, sub) {
  return new MathSubScript({ children: [mr(b + "̂")], subScript: [mr(sub)] });
}

const doc = new Document({
  creator: "Muhammad Imran",
  title: "TCR-IMAP: Homology-Controlled Evaluation of TCR-MHC Contact Prediction",
  styles: {
    default: { document: { run: { font: FONT, size: SZ } } },
    paragraphStyles: [
      { id: "Hyperlink", name: "Hyperlink", basedOn: "Normal",
        run: { color: "0563C1", underline: {} } },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children: body,
  }],
});

Packer.toBuffer(doc).then((buf) => {
  fs.writeFileSync(OUT, buf);
  console.log("written:", OUT, (buf.length / 1e6).toFixed(2), "MB", "| equations:", eqN);
});
