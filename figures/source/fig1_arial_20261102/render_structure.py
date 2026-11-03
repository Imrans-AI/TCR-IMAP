"""Render original M3 structural input; no comparator artwork is reused."""
from pathlib import Path
import numpy as np
from pymol import cmd

OUT = Path(__file__).resolve().parent / 'assets'
# Repository root, derived from this file's location: figures/source/<dir>/render_structure.py
ROOT = Path(__file__).resolve().parents[3]
cmd.reinitialize()
cmd.load(str(ROOT/'data/structures/stcrdab/imgt/1ao7.pdb'), 'complex')
cmd.remove('not polymer')
cmd.remove('not ((chain D+E and resi 1-128) or (chain A and (resi 1-90 or resi 1001-1090)) or chain C)')
# View basis puts the TCR above the platform and separates the two TCR chains.
def center(s): return cmd.get_coords(s).mean(axis=0)
mhc = center('chain A and name CA')
tcr = center('chain D+E and name CA')
y = tcr-mhc; y /= np.linalg.norm(y)
x = center('chain E and name CA') - center('chain D and name CA')
x -= y*np.dot(x,y); x /= np.linalg.norm(x)
z = np.cross(x,y); z /= np.linalg.norm(z)
rot = np.stack([x,y,z],axis=1)
coords = cmd.get_coords('all'); cmd.load_coords((coords-mhc)@rot,'all')
cmd.hide('everything')
cmd.show('cartoon')
cmd.set_color('tcra', [0.36,0.64,0.82])
cmd.set_color('tcrb', [0.58,0.54,0.77])
cmd.set_color('mhc', [0.56,0.72,0.50])
cmd.set_color('pep', [0.83,0.42,0.68])
cmd.color('tcra','chain D'); cmd.color('tcrb','chain E')
cmd.color('mhc','chain A'); cmd.color('pep','chain C')
cmd.show('sticks','chain C'); cmd.set('stick_radius',0.18)
cmd.set('cartoon_fancy_helices',1); cmd.set('cartoon_flat_sheets',1)
cmd.set('cartoon_loop_radius',0.20)
cmd.set('orthoscopic',1); cmd.set('ray_shadows',0)
cmd.set('specular',0.15); cmd.set('ambient',0.65); cmd.set('direct',0.35)
cmd.set('ray_opaque_background',0); cmd.bg_color('white')
cmd.set_view((1,0,0,0,1,0,0,0,1,0,0,-180,0,12,0,100,260,-20))
cmd.zoom('all', buffer=1.5)
cmd.png(str(OUT/'1ao7_interface.png'),width=1500,height=1650,dpi=600,ray=1)
cmd.save(str(OUT/'1ao7_interface.pse'))
print('Rendered', OUT/'1ao7_interface.png')
