from pathlib import Path
import numpy as np
from pymol import cmd

OUT=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[3]   # repository root, so the script works from any clone
cmd.reinitialize()
cmd.load(str(ROOT/'data/structures/stcrdab/imgt/1ao7.pdb'),'complex')
cmd.remove('not polymer')
cmd.remove('not ((chain D+E and resi 1-128) or (chain A and (resi 1-90 or resi 1001-1090)) or chain C)')
def center(s):return cmd.get_coords(s).mean(axis=0)
mhc=center('chain A and name CA');tcr=center('chain D+E and name CA')
x=center('chain E and name CA')-center('chain D and name CA');x/=np.linalg.norm(x)
y=tcr-mhc;y-=x*np.dot(x,y);y/=np.linalg.norm(y)
z=np.cross(x,y);z/=np.linalg.norm(z)
cmd.load_coords((cmd.get_coords('all')-mhc)@np.stack([x,y,z],axis=1),'all')
cmd.hide('everything');cmd.show('cartoon')
for name,rgb,sel in [('tcra',[86/255,180/255,233/255],'chain D'),('tcrb',[213/255,94/255,0],'chain E'),
                     ('mhc',[191/255]*3,'chain A'),('peptide',[255/255,217/255,47/255],'chain C')]:
    cmd.set_color(name,rgb);cmd.color(name,sel)
cmd.cartoon('tube','chain C')
cmd.set('cartoon_tube_radius',.6,'chain C')
cmd.set('cartoon_fancy_helices',1);cmd.set('cartoon_flat_sheets',1)
cmd.set('cartoon_loop_radius',.23)
cmd.set('orthoscopic',1);cmd.set('ray_shadows',0)
cmd.set('specular',0);cmd.set('ambient',1);cmd.set('direct',0)
cmd.set('ray_trace_mode',1);cmd.set('ray_trace_color',[.32,.32,.32])
cmd.set('ray_trace_gain',.08)
cmd.set('ray_opaque_background',0);cmd.bg_color('white')
cmd.set_view((1,0,0,0,1,0,0,0,1,0,0,-180,0,12,0,100,260,-20))
cmd.zoom('all',buffer=1.5)

# Retain all six CDR regions without moving atoms.
for chain in ['D','E']:
    for number,span in enumerate(['27-38','56-65','105-117'],1):
        sel=f'chain {chain} and resi {span}'
        assert cmd.count_atoms(f'({sel}) and name CA')>0
cmd.set('cartoon_loop_radius',.34)
cmd.set('depth_cue',0)
cmd.set('ray_trace_fog',0)
# Separate the two helical walls in projection; de-emphasize the underlying sheet.
cmd.create('floor','chain A and not (resi 50-90 or resi 1050-1090)')
cmd.remove('complex and chain A and not (resi 50-90 or resi 1050-1090)')
cmd.set('cartoon_transparency',.75,'floor')

# Rotate the intact complex about the alpha/beta separation axis: domain centres stay level.
cmd.turn('x',20)
cmd.zoom('all',buffer=1.5)
cmd.png(str(OUT/'1ao7_flat_cartoon.png'),width=2100,height=2500,dpi=900,ray=1)
print('Rendered intact complex, level TCR domain centres, exposed helical groove.')
