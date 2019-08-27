#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap
import matplotlib.colors as mcolors
import set_constant as cst

transpose   = lambda a: list(map(list,zip(*a)))
Reigenval   = lambda tag, f, x: transpose([ [ float(b.text.split()[0]) - f for b in tag[k] ]
                                           for k in range(x[0],x[1])])
#Reigenval   = lambda tag, f, x: transpose([ [ float(k[b].text.split()[0]) - f for b in range(x[0],x[1]) ] for k in tag])
Rprocar     = lambda tag, w, x: [ [ [ [ float(o) for o in i.text.split() ] for i in tag[k][b] ]
                                   for b in range(w[0],w[1]) ] for k in range(x[0],x[1])]
combine_fat = lambda tag, t: transpose([ [ sum(sum( b[i][o] for o in t.orb)
                                       for i in t.ion) for b in k ] for k in tag ])
normalize   = lambda fat, norm: [ cst.max_width*f/norm for f in fat ]
maxfat      = lambda fat: max( max(b) for b in fat )

class Band:
  def __init__(self, proj_tag, fermi, ispin, kptlist, e_range, x_range):
    Band.proj_tag = proj_tag
    Band.ispin    = ispin
    Band.x_range  = x_range
    self.Reigen(fermi)
    if Band.ispin == 2: print("  Spin up")
    Band.eigen, Band.window = self.select_band(Band.eigen, kptlist, e_range)
    if Band.ispin == 2:
      print("  Spin down")
      Band.eigendn, Band.windowdn = self.select_band(Band.eigendn, kptlist, e_range)

# Read Eigenval
  def Reigen(self, fermi):
    print("======= Read EIGENVAL =======")
    eigenval_tag = Band.proj_tag.find("eigenvalues").find("array").find("set")
    Band.eigen = Reigenval(eigenval_tag[0], fermi, Band.x_range)
    if Band.ispin == 2: Band.eigendn = Reigenval(eigenval_tag[1], fermi, Band.x_range)

    return

  def select_band(self, eigen, kptlist, e_range):
    bandmin = []
    bandmax = []
    for b in eigen:
      bandmin.append(min(b))
      bandmax.append(max(b))

    state = 'under'
    window = []
    fermi_range  = []
    for i in range(len(bandmin)):
      if state == 'under' and bandmax[i] > e_range[0]:
        state = 'vb'
        window.append(i)
      if state == 'vb' and bandmax[i] > 0.:
        state = 'fermi'
        fermi_range.append(i)
      if state == 'fermi' and bandmin[i] > 0.:
        state = 'cb'
        fermi_range.append(i)
      if state == 'cb' and bandmin[i] > e_range[1]:
        state = 'over'
        window.append(i)
        break

    if state != 'over': window.append(len(bandmin))

    print("Band index in energy window : {0:d} ~ {1:d}".format(window[0],window[1]))
    if fermi_range[0] == fermi_range[1]:
      i = fermi_range[0]-1
      k = eigen[i].index(bandmax[i])
      print('VBM : {0:6.3f}:{1:s}'.format(bandmax[i],kptlist[k].text))
      i = fermi_range[1]
      k = eigen[i].index(bandmin[i])
      print('CBM : {0:6.3f}:{1:s}'.format(bandmin[i],kptlist[k].text))
    else:
      for i in range(fermi_range[0], fermi_range[1]):
        print("fermi {0:6.3f} ~ {0:6.3f}".format(bandmin[i],bandmax[i]))

    return eigen[window[0]:window[1]], window

  def plotband(self, ax, Kdist):
# Plot spin band structure
    if Band.ispin == 2:
      for b in Band.eigen:
        ax.plot(Kdist, b, color='tab:red',linewidth=1.)
      for b in Band.eigendn:
        ax.plot(Kdist, b, color='tab:blue',linewidth=1.)

# Plot band structure
    else:
      for b in Band.eigen:
        ax.plot(Kdist, b, color='grey',linewidth=1.)

class FatBand(Band):
  procar = None
  def __init__(self, target):
    if FatBand.procar == None:
      print("======== Read PROCAR ========")
      procar_tag = Band.proj_tag.find("array").find("set")
# procar[k][b][i][o]
      FatBand.procar = Rprocar(procar_tag[0], Band.window, Band.x_range)
      if Band.ispin == 2: FatBand.procardn = Rprocar(procar_tag[1], Band.window, Band.x_range)
    self.target = target
    self.fat = combine_fat(FatBand.procar, target)
    self.fmax = maxfat(self.fat)
    if Band.ispin == 2:
      self.fatdn = combine_fat(FatBand.procardn, target)
      self.fmax = max(self.fmax, maxfat(self.fatdn))

  def plotfatband(self, ax, Kdist, norm):

# Plot fat band
    nbands = len(Band.eigen)
    lines = []
    for e, f in list(zip(Band.eigen, self.fat)):
      line = self.plotfillband(ax, Kdist, e, normalize(f, norm))
    if Band.ispin == 2:
      for e, f in list(zip(Band.eigendn, self.fatdn)):
        self.plotfillband(ax, Kdist, e, normalize(f, norm))
    return line

  def plotfillband(self, ax, Kdist, eigen, fat):
    datup = []
    datdn = []
    for bk, fk in list(zip(eigen, fat)):
      datup.append(bk+fk)
      datdn.append(bk-fk)
    line=ax.fill_between(Kdist, datup, datdn, label=self.target.name
                   ,edgecolor=self.target.color, facecolor=self.target.color, alpha=0.33)
    return line

class SpinProjBand(Band):
  combine_fat = lambda self, tag: transpose([ [ sum(sum( o for o in i) for i in b) for b in k ] for k in tag ])
  minfat      = lambda self, fat: max( -min(b) for b in fat )
  def __init__(self, proj_tag, fermi, kptlist, e_range, spin_proj, x_range):
    Band.proj_tag = proj_tag
    Band.ispin    = 4
    self.Reigen(fermi)
    Band.eigen, Band.window = self.select_band(Band.eigen, kptlist, e_range)
    procar_tag = Band.proj_tag.find("array").find("set")
# procar[k][b][i][o]
    procar = Rprocar(procar_tag[spin_proj], Band.window, Band.x_range)
    self.fat = self.combine_fat(procar)
    self.fmax = max(maxfat(self.fat),self.minfat(self.fat))

  def plotband(self, ax, Kdist):
    import numpy as np
    red  = mcolors.to_rgb('tab:red')
    grey = mcolors.to_rgb('tab:grey')
    blue = mcolors.to_rgb('tab:blue')
    N = 128
    vals = []
    for i in range(N):
      temp = []
      for j in range(3):
        temp.append((blue[j]*(N-i)+grey[j]*i)/N)
      temp.append(1.)
      vals.append(temp)
    for i in range(N):
      temp = []
      for j in range(3):
        temp.append((grey[j]*(N-i)+red[j]*i)/N)
      temp.append(1.)
      vals.append(temp)
    cmap = ListedColormap(vals)
    for b, f in list(zip(Band.eigen, self.fat)):
      segments = []
      old = [Kdist[0], b[0]]
      for k, e in list(zip(Kdist[1:], b[1:])):
        segments.append([old, [k,e]])
        old = [k,e]
#     points = np.array([Kdist, b]).T.reshape(-1, 1, 2)
#     segments = np.concatenate([points[:-1], points[1:]], axis=1)

      lc = LineCollection(segments, cmap=cmap, norm=plt.Normalize(-self.fmax, self.fmax))
      lc.set_array(np.array(f))
      lc.set_linewidth(1.)
#     ax = plt.gca()
      ax.add_collection(lc)
