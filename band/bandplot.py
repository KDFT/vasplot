#!/usr/bin/env python3
import xml.etree.cElementTree as ET
import re
try: from band.band import Band, FatBand, SpinProjBand, transpose
except: from band import Band, FatBand, SpinProjBand, transpose

class bandplot:
  def __init__(self, ax, Fband, Fkpts, target, e_range, fermi, ispin, lgd, spin_proj, x_range, e_ticks, plot_mode):
    print("\n==== Band plotting start ====")
# Set plot environments
    root = ET.parse(Fband).getroot()
    Kdist, kptlist = self.Rkdist(root)
    kticks, klabel = self.Rlabel(Fkpts)
    kticks, klabel, Kdist, x_range = self.fixlabel(root, x_range, Kdist, kticks, klabel)

    self.set_env(ax, Kdist, kticks, klabel, e_range, e_ticks)

    proj_tag = root.find("calculation")
#   if 'projected' in [i.tag for i in proj_tag.getchildren()]:
    if any(child.tag == 'projected' for child in proj_tag):
      proj_tag = root.find("calculation").find("projected")
    band, fatband, norm = self.bandrepair(proj_tag, target, kptlist, e_range, fermi, ispin, spin_proj, x_range)

    if plot_mode == 'line': lgd = False

    legends = []
    for fb in fatband:
      line = fb.plotfatband(ax, Kdist, norm, plot_mode)
      legends.append(line)
# Set legend
    if lgd: self.set_legend(ax, target, legends)

    if plot_mode != 'line': band.plotband(ax, Kdist)

    print("===== Band plotting end =====")

# Read KPOINTS (Kdist)
  def Rkdist(self, root):
    print("======= Get Kdistance =======")
    rec_tag = root.find('structure[@name="finalpos"]').find("crystal").find('varray[@name="rec_basis"]')
# reciprocal lattice vectors
    rec_basis = []
    for recV in rec_tag:
      rec_vec = []
      for i in recV.text.split():
        rec_vec.append(float(i))
      rec_basis.append(rec_vec)

    kptlist = root.find("kpoints").find("varray")

# Kdist[kpt]
    Kdist =[0.]
    kpt2 = None
    distold = None
#   for k in kptlist:
    for x,k in enumerate(kptlist):
      kpt1 = [0., 0., 0.]
      for i,ktemp in enumerate(k.text.split()):
        for j in range(3):
          kpt1[j] += rec_basis[i][j] * float(ktemp)
      if kpt2 != None:
        if kpt1 != kpt2:
          dsum = 0.
          for i in range(3):
            dsum += (kpt2[i]-kpt1[i])**2
          dist = dsum**0.5
          if distold == None:
            Kdist.append(Kdist[-1] + dist)
          elif abs((distold - dist)/distold) < 2. :
            Kdist.append(Kdist[-1] + dist)
          else:
            Kdist.append(Kdist[-1])
# Add k-distance along band line
          distold = dist
# For discontinued band
        else:
          Kdist.append(Kdist[-1])
      kpt2 = kpt1

    return Kdist, kptlist

# Read KPOINTS (tag)
  def Rlabel(self, Fkpts):
    if Fkpts == None: return [], []
    print("======= Read KPOINTS ========")
    f     = open(Fkpts)
    lines = f.readlines()
    f.close()

    for i, l in enumerate(lines):
      if l.split()[0].isdigit():
        num_grd = int(l.split()[0])
        kmode   = lines[i+1][0].lower()
        lines   = lines[i+2:]
        break

    if kmode == 'l':
      for i, l in enumerate(lines):
        if 'reciprocal' in l.lower():
          lines = lines[i+1:]
          break
      labels = []
      for l in lines:
        line = l.split()
        if len(line) > 3: labels.append(line[-1])
      for i in range(len(labels)):        # G to Gamma
        if 'G' in labels[i]: labels[i] = '$\Gamma$'

# kticks[special point]
# klabel[special point]
      klabel = [ labels[0] ]
      kticks = [0,num_grd-1]
      for i in range(int(len(labels)/2)-1):
        if labels[2*i+1] == labels[2*i+2]:
          klabel.append(labels[2*i+1])
# For discontinued band
        else:
          klabel.append(labels[2*i+1] + '|' + labels[2*i+2])
        kticks.append((i+2)*num_grd-1)
      klabel.append(labels[-1])
      return kticks, klabel

    else:
      print("There are no line mode KPOINTS")

  def fixlabel(self, root, x_range, Kdist, kticks, klabel):
# Fake
    temp = []
    startk = 0
    for i, w in enumerate(root.find("kpoints").find('varray[@name="weights"]')):
      if float(w.text) == 0.: temp.append(i)
    if temp != []:
      if isinstance(x_range,list):
        outer = [i for i in x_range if isinstance(i,int)]
        inner = [i for i in x_range if isinstance(i,list)]
        if len(inner) == 0 and len(outer) == 0:
          startk = min(temp)
          Kdist = Kdist[startk:max(temp)]
          x_range = [startk,startk +len(Kdist)]
        elif (len(inner) == 0 and len(outer) >= 2):
          startk = min(outer)
          x_range = [startk, max(outer)]
          Kdist   = Kdist[startk:max(outer)]
        elif (len(inner) >= 2 and len(outer) == 0):
          startk = min(inner)
          x_range = [startk, max(inner)]
          Kdist   = Kdist[startk:max(inner)]
        elif (len(inner) >= 2 and len(outer) >= 2):
          startk = min(outer)
          x_range = [min(inner), max(inner)]
          Kdist   = Kdist[startk:max(outer)]
        elif len(outer) == 1:
          startk = min(outer)
          if kticks != []: Kdist = Kdist[startk:startk+kticks[-1]+1]
          else           : Kdist = Kdist[startk:]
          if len(inner) < 2:
            x_range = [startk, startk+len(Kdist)]
          else:
            x_range = [min(inner),max(inner)]
      else:
        startk = min(temp)
        Kdist = Kdist[startk:max(temp)]
        x_range = [startk, startk+len(Kdist)]
    if x_range == None: x_range = [0, len(Kdist)]
    newticks = []
    newlabel = []
    if startk-1 in kticks:
      newticks.append(Kdist[0])
      newlabel.append(klabel[kticks.index(startk-1)])
    for i in range(len(Kdist)):
      if i+startk in kticks:
        newticks.append(Kdist[i])
        newlabel.append(klabel[kticks.index(i+startk)])
    return newticks, newlabel, Kdist, x_range

# Set plot environments
  def set_env(self, ax, Kdist, kticks, klabel, e_range, e_ticks):
    ax.set_xticks(kticks)
    ax.set_xticklabels(klabel)
    ax.set_xlim(0,Kdist[-1])
    ax.set_ylim(e_range[0], e_range[1])
    ax.set_ylabel(r"Energy vs E$_{\rm{F}}$ (eV)")
    if e_ticks != None:
      e_ticks, e_label = transpose(e_ticks)
      ax.set_yticks(e_ticks)
      ax.set_yticklabels(e_label)
    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1)
    ax.grid(True)
    ax.set(aspect = (Kdist[-1]-Kdist[0])/(e_range[1]-e_range[0]))

    return

# Band Data Process
  def bandrepair(self, proj_tag, target, kptlist, e_range, fermi, ispin, spin_proj, x_range):
    print("======== Repair band ========")

    band = Band(proj_tag, fermi, ispin, kptlist, e_range, x_range)
    if spin_proj != None: band = SpinProjBand(proj_tag, fermi, kptlist, e_range, spin_proj, x_range)

    fatband = []
    for t in target:
      fatband.append(FatBand(t))

# Normalize fat bands width
    norm = 0.
    for fb in fatband:
      if fb.fmax > norm: norm = fb.fmax

    print(' norm  : {0:6.3f}'.format(norm))

    return band, fatband, norm

  def set_legend(self, ax, target, legends):
    maxlgd = 0
    for t in target:
      if '$' in t.name:
        lgd_num = t.name.find('$') + int((len(t.name) - t.name.find('$')-3)*2/3)
      else            : lgd_num = len(t.name)
      if lgd_num > maxlgd: maxlgd = lgd_num
    if   maxlgd < 6 : ncol = 4
    elif maxlgd < 10: ncol = 3
    elif maxlgd < 17: ncol = 2
    else            : ncol = 1
    ax.legend(handles=legends,bbox_to_anchor=(0., -0.2, 1., .102), loc=2,
                         ncol=ncol, mode="expand", borderaxespad=0.)

    return
