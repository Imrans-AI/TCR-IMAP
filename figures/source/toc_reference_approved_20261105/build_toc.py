from pathlib import Path
import os,json,hashlib
os.environ.setdefault('MPLCONFIGDIR','/tmp/m3-ga-matplotlib')
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap,Normalize
from matplotlib.patches import Polygon,Rectangle
from PIL import Image
O=Path(__file__).resolve().parent
R=O.parents[2]
M=O/'complex.png'
P=R/'results/tables/pairmaps_tcr70_p84v2.npz';S=R/'results/tables/similarity_curve_p84v2_authoritative.csv'
z=np.load(P,allow_pickle=True);i=list(z['pdb']).index('1ao7');v=z['S'][i]>0;pred=z['pred_model'][i]
sweep=pd.read_csv(S).sort_values('threshold',ascending=False)
assert v.sum()==5208 and z['covered'][i] and pred.shape==(75,84)
assert np.isclose(sweep['diff'].iloc[0],.2750) and np.isclose(sweep['diff'].iloc[-1],.0274)
plt.rcParams.update({'font.family':'Arial','font.size':6,'pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','axes.linewidth':.65})
W,H=82.5,44.5
fig=plt.figure(figsize=(W/25.4,H/25.4),facecolor='white')
a=fig.add_axes([0,0,1,1],xlim=(0,W),ylim=(0,H));a.axis('off')
def text(x,y,s,bold=False,color='#202020',ha='center',family='Arial'):
 return a.text(x,y,s,ha=ha,va='center',fontsize=6,fontfamily=family,fontweight='bold' if bold else 'normal',color=color)
def inset(x,y,w,h):return fig.add_axes([x/W,y/H,w/W,h/H])
def arrow(x1,y1,x2,y2,half=.25,head=.9,tip=1.4):
 d=np.array([x2-x1,y2-y1],float);d/=np.linalg.norm(d);n=np.array([-d[1],d[0]]);p=np.array([x1,y1]);q=np.array([x2,y2]);b=q-tip*d
 a.add_patch(Polygon([p-half*n,b-half*n,b-head*n,q,b+head*n,b+half*n,p+half*n],fc='#555555',ec='none'))
# Molecule-to-sequence visual motif follows the supplied reference.
m=inset(1.0,19.7,22.5,22.5);m.axis('off');m.imshow(Image.open(M),interpolation='none')
arrow(12.25,19.3,12.25,16.9,half=.2,head=.65,tip=.85)
# Representative real sequence fragments, drawn as flat residue tiles.
fragments=[('DRGSQS','#56B4E9'),('MNHEY','#D55E00'),('PWIEQE','#BFBFBF')]
for row,(seq,col) in enumerate(fragments):
 width=len(seq)*2.65;left=12.25-width/2;y=14.0-row*2.85
 for k,aa in enumerate(seq):
  x=left+k*2.65
  a.add_patch(Rectangle((x,y),2.35,2.35,facecolor='white',edgecolor=col,linewidth=.65))
  text(x+1.175,y+1.175,aa,color=col if row!=2 else '#666666',family='Courier New')
text(12.25,4.0,'Sequence input',bold=True)
arrow(25.0,25.0,28.8,25.0)
# Actual stored out-of-fold prediction, now given a full-sized scientific panel.
h=inset(32.0,15.0,21.0,23.0)
cmap=LinearSegmentedColormap.from_list('contacts',['#F7F7F7','#08306B'],N=256);cmap.set_bad('white')
h.imshow(np.ma.array(pred,mask=~v),cmap=cmap,norm=Normalize(0,1),interpolation='none',aspect='auto')
h.set_xticks([]);h.set_yticks([])
h.set_xlabel('MHC residues',fontsize=6,labelpad=2)
h.set_ylabel('TCR residues',fontsize=6,labelpad=2)
for sp in h.spines.values():sp.set_color('#202020');sp.set_linewidth(.65)
text(42.5,40.1,'Predicted contacts')
text(42.5,4.0,'TCR-IMAP',bold=True)
arrow(54.5,25.0,58.2,25.0)
# Evaluation uses the real threshold sweep, not the unrelated reference's R-squared.
c=inset(60.3,15.0,19.0,23.0)
x=np.arange(len(sweep));y=sweep['diff'].to_numpy()
c.plot(x,y,color='#0072B2',lw=1.15,marker='o',ms=2.6,mec='white',mew=.45)
c.axhline(0,color='#777777',lw=.5,ls=(0,(2,2)))
c.set_xlim(-.2,6.2);c.set_ylim(-.01,.32)
c.set_xticks([0,6]);c.set_xticklabels(['',''])
c.set_yticks([0,.3]);c.set_yticklabels(['0','+0.3'],fontsize=6)
c.yaxis.tick_left();c.tick_params(axis='y',length=1.2,pad=1,width=.5);c.tick_params(axis='x',length=1.2,width=.5)
for sp in c.spines.values():sp.set_color('#202020');sp.set_linewidth(.65)
text(69.8,40.1,'Δ macro-AUPRC')
c.text(.5,.91,'vs frequency prior',transform=c.transAxes,ha='center',va='center',fontsize=6)
text(60.3,12.9,'unrestricted',ha='left');text(80.3,12.9,'0.70',ha='right')
text(70.0,9.0,'maximum train/test\nsimilarity')
text(69.8,4.0,'Evaluation',bold=True)
fig.canvas.draw()
assert sum(t.get_fontweight()=='bold' for t in a.texts)==3
renderer=fig.canvas.get_renderer()
for axis in [a,h,c]:
    for t in (list(axis.texts)+([*axis.get_xticklabels(),*axis.get_yticklabels(),axis.xaxis.label,axis.yaxis.label] if axis.axison else [])):
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        assert b.x0>=0 and b.x1<=fig.bbox.width and b.y0>=0 and b.y1<=fig.bbox.height,(t.get_text(),b)

for ext in ['pdf','svg']:fig.savefig(O/f'TCR_IMAP_TOC_reference_style.{ext}',dpi=900)
fig.savefig(O/'TCR_IMAP_TOC_reference_style_900dpi.png',dpi=900)
fig.savefig(O/'TCR_IMAP_TOC_reference_style_preview.png',dpi=400)
plt.close(fig)
Image.open(O/'TCR_IMAP_TOC_reference_style_900dpi.png').convert('RGB').save(O/'TCR_IMAP_TOC_reference_style_900dpi.tif',dpi=(900,900),compression='tiff_lzw')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
(O/'PROVENANCE.json').write_text(json.dumps({'size_cm':[8.25,4.45],'raster_dpi':900,'reference':'user-supplied descriptor / machine-learning / prediction composition; no reference data copied','fonts':'Arial 6 pt; Courier New 6 pt residue tokens; only bottom stage labels bold','structure_source':str(M),'structure_sha256':sha(M),'structure':'actual complex, intact geometry; molecular context, not a structural model input','sequence_fragments':fragments,'sequence_scope':'representative resolved alpha CDR1, beta CDR1 and MHC window fragments from same example; not all model inputs shown','map_source':str(P),'map_sha256':sha(P),'map_transform':'unmodified stored predictions, unresolved cells white','curve_source':str(S),'curve_sha256':sha(S),'curve_values':y.tolist(),'thresholds':sweep.threshold.tolist(),'curve_scope':'full-cohort evaluation; threshold settings equally spaced; no synthetic regression data','legends':False,'older_versions_preserved':True},indent=2)+'\n')
print('Created reference-style TOC using the existing molecular rendering, stored predictions and measured evaluation.')
