"""Second M3 design: no dataset counts; straight architecture flow.

Run with dl_env Python. Protein artwork is rendered separately by render_structure.py.
Output dimensions are fixed, without bbox cropping: 7.2 x 4.7 inches, 600 dpi.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Ellipse, Arc
from PIL import Image

OUT=Path(__file__).resolve().parents[1]
ROOT=Path('/home/muhammad-imran/PhD_data/TCRMHC_Module_3_TCR-IMAP')
z=np.load(ROOT/'results/tables/pairmaps_tcr70_p84v2.npz',allow_pickle=True)
idx=list(z['pdb']).index('1ao7')
md=pd.read_csv(ROOT/'data/processed/complex_metadata.tsv',sep='\t')
assert z['Y'].shape==(174,75,84)
assert len(md)==174 and md.clu70.nunique()==69
assert bool(z['covered'][idx]) and np.isfinite(z['pred_model'][idx]).all()
mask=z['S'][idx]>0
observed=np.ma.array(z['Y'][idx],mask=~mask)
prediction=np.ma.array(z['pred_model'][idx],mask=~mask)

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7.4,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
 'axes.linewidth':0.6,'text.color':'#000000'})
fig=plt.figure(figsize=(7.2,4.7),facecolor='white')
ax=fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,70); ax.axis('off')
C={'blue':'#538EB4','bluefill':'#E7F0F7','purple':'#8E80B7','green':'#86AA74',
   'greenfill':'#EAF2E6','orange':'#E9BE97','orangefill':'#FAEBDD',
   'ink':'#26313A','line':'#5B626A','grey':'#F1F3F5','muted':'#68717B'}
def text(x,y,s,size=7.4,weight='normal',color=None,ha='center',va='center',**kw):
 # Consistent black typography, as requested; colours remain in structures and outlines.
 return ax.text(x,y,s,fontsize=size,fontweight=weight,color='#000000',ha=ha,va=va,**kw)
def box(x,y,w,h,s,fc='white',ec='#A9B3BC',size=7.1,weight='normal',radius=.9):
 p=FancyBboxPatch((x,y),w,h,boxstyle=f'round,pad=0.02,rounding_size={radius}',
  facecolor=fc,edgecolor=ec,linewidth=.8,zorder=3)
 ax.add_patch(p); text(x+w/2,y+h/2,s,size,weight,zorder=4,linespacing=1.35)
 return p
def arrow(points,color=None,lw=.9,style='-|>',dashed=False):
 for a,b in zip(points[:-2],points[1:-1]):
  ax.plot([a[0],b[0]],[a[1],b[1]],color=color or C['line'],lw=lw,zorder=2,
          linestyle='--' if dashed else '-')
 ax.add_patch(FancyArrowPatch(points[-2],points[-1],arrowstyle=style,mutation_scale=8,
   linewidth=lw,color=color or C['line'],linestyle='--' if dashed else '-',zorder=2,
   shrinkA=0,shrinkB=1))
def subaxis(x,y,w,h):
 return fig.add_axes([x/100,y/70,w/100,h/70])
def heatmap(x,y,w,h,arr,pred=False):
 a=subaxis(x,y,w,h)
 cm=plt.get_cmap('Blues').copy(); cm.set_bad('#DDE1E5')
 a.imshow(arr,cmap=cm,vmin=0,vmax=1,aspect='auto',interpolation='nearest',rasterized=False)
 a.set_xticks([]);a.set_yticks([])
 for s in a.spines.values():s.set_edgecolor('#98A3AD');s.set_linewidth(.6)
 return a
def tokens(x,y,w,h,cols):
 n=len(cols)
 for i,col in enumerate(cols):
  ax.add_patch(Rectangle((x+i*w/n,y),w/n,h,facecolor=col,edgecolor='white',lw=.65,zorder=3))
def stack(x,y,w,h,col):
 for i in [2,1,0]:
  ax.add_patch(Rectangle((x+i*.65,y+i*.55),w,h,facecolor=col,edgecolor='#82919E',lw=.65,zorder=3))

# A: an overview, not a record of every processing or audit step.
text(1.5,67.2,'A',12,'bold',ha='left')


# Compact white database cylinder with two internal curved dividers.
ax.add_patch(Rectangle((3.7,53.4),8.0,5.6,fc='white',ec='none',zorder=3))
for xx in (3.7,11.7):
 ax.plot([xx,xx],[53.4,59.0],color='#7891A5',lw=.9,zorder=4)
for yy in (53.4,55.3,57.2):
 ax.add_patch(Arc((7.7,yy),8.0,1.6,theta1=180,theta2=360,ec='#7891A5',lw=.9,zorder=4))
ax.add_patch(Ellipse((7.7,59.0),8.0,1.6,fc='white',ec='#7891A5',lw=.9,zorder=4))
text(7.7,62.0,'STCRDab',7.3,zorder=5)
text(7.7,47.8,'TCR–pMHC\nstructures',7.0,linespacing=1.5)
arrow([(12.5,56.1),(19,56.1)])

# Original molecular image; only variable domains and MHC platform displayed.
mol=subaxis(19.4,43.0,25.5,22.0)
mol.imshow(Image.open(OUT/'assets/1ao7_interface.png'));mol.axis('off')
text(23.0,63.3,'TCRα',7.0,color=C['blue'])
text(39.8,63.3,'TCRβ',7.0,color=C['purple'])
text(31.9,43.5,'MHC-I platform',7.0,color='#61864F')
text(39.5,49.9,'peptide',6.7,color='#B46A9F',ha='left')
ax.plot([33.5,38.8],[51.7,50.5],color='#B46A9F',lw=.65)

# Branch: structures supply labels and aligned sequence inputs.
arrow([(44.5,56),(47,56),(47,59.7),(50,59.7)])
arrow([(47,56),(47,48),(50,48)])
heatmap(51.0,55.6,12.8,7.3,observed)
text(57.4,64.0,'Contact labels',7.3)
text(57.4,53.8,'4.5 Å cutoff',6.5,color=C['muted'])
tokens(51,49.5,12.8,1.75,[C['blue']]*3+[C['purple']]*3)
tokens(51,46.5,12.8,1.75,[C['green']]*6)
text(57.4,44.4,'Aligned sequences',7.0)

# Evaluation is visible without turning the overview into a checklist.
arrow([(64.5,59.6),(69.8,59.6)])
arrow([(64.5,48.1),(67.2,48.1),(67.2,52.0),(69.8,52.0)])
box(70,46.3,28.0,16.2,'',fc='white',ec='#B0BBC5',radius=1.7)
text(84,60.4,'Train and evaluate',8.0,'normal')
text(84,56.7,'Random and grouped splits',7.0)
text(84,53.4,'TCR-IMAP, comparators and priors',6.7)
text(84,49.8,'Paired cluster bootstrap',7.0)
ax.plot([73,95],[58.7,58.7],color='#D5DCE2',lw=.6,zorder=5)

# B: faithful pairwise architecture. Cross-attention updates TCR only.
ax.plot([1.5,98.5],[40.8,40.8],color='#DCE2E7',lw=.7)
text(1.5,38.2,'B',12,'bold',ha='left')


text(8.0,33.3,'TCR residues',7.2,'normal')
tokens(2.6,28.7,11.0,2.1,[C['blue']]*4+[C['purple']]*4)
text(8.1,26.7,'Aligned residues',6.5,color=C['muted'])
text(8.0,21.0,'MHC residues',7.2,'normal')
tokens(2.6,16.4,11.0,2.1,[C['green']]*8)
text(8.1,14.4,'Aligned residues',6.5,color=C['muted'])

box(18,26.0,15.0,7.5,'Residue + segment\n+ position\nembeddings',fc='white',ec=C['blue'],size=6.5)
box(18,13.7,15.0,7.5,'Residue + position\nembeddings',fc='white',ec=C['green'],size=6.5)
arrow([(13.8,29.75),(17.9,29.75)])
arrow([(13.8,17.45),(17.9,17.45)])

box(37,26.0,14.0,7.5,'Transformer\nencoder',fc='white',ec=C['blue'],size=7.0)
box(37,13.7,14.0,7.5,'Transformer\nencoder',fc='white',ec=C['green'],size=7.0)
arrow([(33.1,29.75),(36.9,29.75)])
arrow([(33.1,17.45),(36.9,17.45)])

box(56,24.4,14.0,9.1,'Cross-\nattention',fc='white',ec=C['purple'],size=7.4)
text(63.0,22.4,'Residual + norm',6.3,color=C['muted'])
arrow([(51.1,29.75),(55.9,29.75)])
text(53.4,31.4,'Q',6.0,color=C['muted'])
arrow([(51.1,17.45),(53.5,17.45),(53.5,26.3),(55.9,26.3)])
text(52.6,23.0,'K, V',5.7,color=C['muted'],ha='right')

box(74,24.4,10.0,9.1,'Bilinear\nhead',fc='white',ec='#CE9B6C',size=7.0)
arrow([(70.1,28.9),(73.9,28.9)])
arrow([(51.1,17.45),(79.0,17.45),(79.0,24.3)])
text(66.5,15.2,'MHC representations',6.1,color=C['muted'])

# Direct, left-to-right output rather than an arrow around the page edge.
arrow([(84.1,28.9),(88.1,28.9)])
text(86.0,30.7,'σ',7.2,color=C['muted'])
heatmap(88.3,24.0,10.0,9.4,prediction,True)
text(93.3,20.4,'Contact\nprobabilities',7.0,linespacing=1.3)
# Input and output dimensions are documented in the caption.
cb=subaxis(89.3,13.5,8.0,.65)
cb.imshow(np.linspace(0,1,128)[None,:],cmap='Blues',vmin=0,vmax=1,aspect='auto')
cb.axis('off')
text(88.0,13.8,'0',5.6,color=C['muted'])
text(98.6,13.8,'1',5.6,color=C['muted'])

# Concise visual key; masks shown consistently in both maps.
ax.plot([2.6,98.5],[9.7,9.7],color='#E4E8EC',lw=.65)
text(2.6,6.4,'Model inputs: amino-acid identities',6.7,ha='left',color=C['muted'])
ax.add_patch(Rectangle((51.0,5.6),1.6,1.6,fc='#DDE1E5',ec='#B1B8BE',lw=.5))
text(53.5,6.4,'Unresolved pairs masked',6.4,ha='left',color=C['muted'])
text(98.3,6.4,'Maps: 1AO7',6.4,ha='right',color=C['muted'])
text(50.0,2.7,'Q: TCR queries   ·   K: MHC keys   ·   V: MHC values   ·   σ: sigmoid (0–1 probabilities)',
     6.5,color=C['ink'])

stem=OUT/'M3_Figure1_final'
fig.savefig(str(stem)+'.pdf',dpi=600)
fig.savefig(str(stem)+'.svg',dpi=600)
fig.savefig(str(stem)+'_900dpi.png',dpi=900)
fig.savefig(OUT/'M3_Figure1_final_preview.png',dpi=180)
plt.close(fig)
im=Image.open(str(stem)+'_900dpi.png')
assert im.size==(6480,4230),im.size
print('Exported',stem,'6480 × 4230 pixels; 7.2 × 4.7 inches at 900 dpi')

caption='''Figure 1. Structural benchmark and sequence-based contact prediction with TCR-IMAP.
(A) STCRDab complexes provide the structural benchmark, aligned sequence inputs and crystal-contact
labels. The 1AO7 illustration shows the TCR variable domains, peptide and MHC-I platform; domains
outside this interface view are omitted. TCRα is blue, TCRβ purple, MHC green and peptide magenta.
Contacts are defined by a minimum inter-residue heavy-atom distance of 4.5 Å. Methods and within-fold
priors are compared under frozen random and grouped partitioning protocols using paired per-complex
metrics and cluster-bootstrap inference; the primary protocol is tcr70. Protocol definitions,
candidate masks and per-residue comparisons are specified in Methods.
(B) The pairwise TCR-IMAP variant embeds 75 aligned TCR rows and 84 MHC columns. Separate
Transformer encoders produce residue representations; TCR queries attend to MHC keys and values.
A residual connection and layer normalization update the TCR representations. The bilinear head
combines these with the MHC representations, together with a learned positional term, scale and
bias, and a sigmoid yields contact probabilities. The displayed prediction is the stored tcr70
out-of-fold prediction for 1AO7, not a fabricated illustration. Grey map cells indicate unresolved
pairs; both maps retain the original row and column order. Structural coordinates define labels
and masks and are not input features for this sequence-based model. The separate per-residue
variant is described in Methods.'''
(OUT/'Figure1_caption.txt').write_text(caption)
(OUT/'DESIGN_NOTES.txt').write_text('''Review draft, not installed in the M3 manuscript.
The supplied comparator figures informed the separation of study overview and model architecture,
restrained colours, short labels and molecular context. No comparator artwork was copied.
The original molecular illustration was ray-traced from the M3 STCRDab 1AO7 structure with PyMOL.
All other diagram elements are editable vectors in SVG/PDF; the original molecular rendering is
embedded raster artwork. The PyMOL session and rendering script are supplied for structural edits.
This alternative omits all dataset counts and input/output dimensions and simplifies the architecture to a direct horizontal flow.
The 900-dpi PNG contains 6480 × 4230 pixels. Panel letters A and B are bold; descriptive panel headings are omitted. Train and evaluate, TCR residues and MHC residues use regular weight. PDF/SVG preserve vector text and connectors.
Scientific source: scripts/19_train_pairmap.py PairNet; frozen pairmaps_tcr70_p84v2.npz;
complex_metadata.tsv; data/structures/stcrdab/imgt/1ao7.pdb.
No manuscript files or M3 repository files were modified. No training or MD analysis was performed.
''')
