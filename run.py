#!/usr/bin/env python3
import vasplot as vpl

# Input parameters
target  = [ [ 'Mn', 'd' ], [ 'Ta', 'd' ] ]
# [ [ atom ], [ orbital ], [ color ] ]
# atom    : Species or number (count from 0)
# orbital : s, p, d, t2g, eg, s:0, py:1, pz:2, px:3, dxy:4, dyz:5, dz2:6, dxz:7, x2-y2:8 (default : all)
# color   : (default : all) (cf. ['C', 'all', 'red'] )
# If doesn't declare target, or target = '', [], None, or doesn't include atomic information,
# only plot band or total dos
# Example : None | 'Mn'  |  ['Mn','d']  |  ['Mn', 'all', 'red'] |
#           [ ['C'], ['C', 'p'], ['C', 'all', 'red'], ['C', 'px'], [0], [ [1,2], 'd']
#           , [ 'C', [1,2,3] ], ['C', ['px, 'py'] ], [ ['C', 'Si'], ['s', 'p'] ]
e_range = [-4., 2.]

# Input files
Fdos  = '../vasprun.xml'
Fband = './vasprun.xml'
Fkpts = './KPOINTS'
Fsave = './ygband.png' # None, 'Show', [ 'Show', './ygband.pdf', 'ygband.png' ]
spin_proj = None       # 1, 'x', 2, 'y', 3, 'z'
x_range   = None       # [0,150] | [0, [50, 100], 150] normal : choose range to see, if outer->KPOINTS
e_ticks   = None       # [[-4., 0.], [-2., 2.], [0., 4.], [2., 6.]] (vs fermi level, label)
fermi     = None       # Set external fermi level

vpl.vasplot(target=target, Fdos=Fdos, Fband=Fband, Fkpts=Fkpts, e_range=e_range, Fsave=Fsave, spin_proj=spin_proj, x_range=x_range, e_ticks=e_ticks, fermi = fermi)
#vpl.vasplot(target=target, Fdos=Fdos, Fband=Fband, Fkpts=Fkpts, e_range=e_range, Fsave=Fsave, spin_proj=spin_proj, x_range=x_range, e_ticks=e_ticks, fermi = fermi, plot_mode = None)
