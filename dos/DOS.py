#!/usr/bin/env python3
import set_constant as cst

pi = 3.141592
e  = 2.718281 #8284590452353602874

# pdos[ion][energy][orbital] = partial DOS
Rpdos     = lambda tag, s: [ [ [ float(o) for o in e.text.split()[1:] ] for e in i[s] ] for i in tag ]
Gaussian  = lambda x: 1/(cst.s*pow(2*pi, 1/2)) * pow(e, -1/2 * (x/cst.s)**2)
transpose = lambda a: list(map(list,zip(*a)))
combine   = lambda pdos, t: [ sum(sum( e[i][o] for o in t.orb) for i in t.ion) for e in transpose(pdos) ]

class DOS:
  weight = None
  def __init__(self, dos_tag, fermi, ispin, e_range):
    DOS.dos_tag = dos_tag
    DOS.ispin   = ispin
    self.Rtotdos(fermi, ispin)
    self.Choose_Erange(e_range)

# Read total dos
  def Rtotdos(self, fermi, ispin):
    dostot_tag = DOS.dos_tag.find('total').find('array').find('set')

# energy[nedos]
    energy   = []

# dos[nedos]
    dos   = []
    for i in dostot_tag[0]:                           # energy
      dosE = i.text.split()
      energy.append(float(dosE[0])-fermi)
      dos.append(float(dosE[1]))
    DOS.energy = energy
    self.dos = dos
    if ispin == 2:
      dosdn = []
      for i in dostot_tag[1]:                           # energy
        dosdn.append(-float(i.text.split()[1]))
      self.dosdn = dosdn
    return

# Choose energy range in e_range
  def Choose_Erange(self, e_range):
    DOS.de = (DOS.energy[-1]-DOS.energy[0])/(len(DOS.energy)-1)

    e_range_index = []
    for i, e in enumerate(DOS.energy):
      if len(e_range_index) == 0 and e > e_range[0]:
        e_range_index.append(i-1)
      if len(e_range_index) == 1 and e > e_range[1]:
        e_range_index.append(i+1)
        break

    if e_range_index[0] < 1                : e_range_index[0] = 1
    if len(e_range_index) == 1: e_range_index.append(len(DOS.energy) - 2)
    if e_range_index[1] > len(DOS.energy)-2: e_range_index[1] = len(DOS.energy) - 2
    DOS.energy = DOS.energy[e_range_index[0]:e_range_index[1]]
    DOS.e_range_index = e_range_index

# Gaussian broadening
  def BroaDOS(self):
    if DOS.weight == None: DOS.weight, DOS.width = self.norm_Gaussian()

    self.dos = self.broadening(self.dos)
    if DOS.ispin == 2: self.dosdn = self.broadening(self.dosdn)
    return

  def norm_Gaussian(self):
    tol   = 0.0001
    w = []
    for i in range(len(DOS.energy)):
      temp = Gaussian(i*DOS.de) * DOS.de
      if temp > tol:
        w.append(temp)
      else:
        width = i-1
        break
    weight = w[:]
    w.reverse()

# weight[-width,width]
    weight.extend(w[:-1])

# normalized Gaussian function
    norm = 0.
    for i in weight:
      norm += i
    for i in range(len(weight)):
      weight[i] = weight[i]/norm

    return weight, width

  def broadening(self, dos):
    max_range = len(dos)
    newdos = []
    for i in range(DOS.e_range_index[0], DOS.e_range_index[1]):
      temp = 0.
      for j in range(-DOS.width, DOS.width+1):
# Except end value
        if i+j < 1          : continue
        if i+j > max_range-2: continue
        temp += dos[i+j] * DOS.weight[j]
      newdos.append(temp)
    return newdos

  def get_dosmax(self):
    dosmax = max(self.dos)
    if DOS.ispin == 2: dosmax = max(dosmax, -min(self.dosdn))
    return dosmax

  def plotDOS(self, ax):
# Plot spin DOS
    if DOS.ispin == 2:
      ax.plot(self.dos,DOS.energy,color='tab:red', linewidth=0.5)
      ax.plot(self.dosdn,DOS.energy,color='tab:blue', linewidth=0.5)

# Plot DOS
    else:
      ax.plot(self.dos,DOS.energy,color='grey', linewidth=0.5)

class PDOS(DOS):
  pdos = None
  occ  = None
  def __init__(self, target):
    self.target = target
    if PDOS.pdos == None:
      partial_tag = DOS.dos_tag.find('partial').find('array').find('set')
# pdos[ion][energy][orbital] = partial DOS
      PDOS.pdos = Rpdos(partial_tag, 0)
      if DOS.ispin == 2:
        PDOS.pdosdn = Rpdos(partial_tag, 1)
    self.combine_PDOS()

  def combine_PDOS(self):
    self.dos = combine(PDOS.pdos, self.target)
    if DOS.ispin == 2:
      self.dosdn = [ -i for i in combine(PDOS.pdosdn, self.target) ]
    return

# Plot PDOS
  def plotPDOS(self, ax):
    if PDOS.occ == None:
      PDOS.occ = []
      for i in DOS.energy:
        PDOS.occ.append(i<0.)

    ax.plot(self.dos, DOS.energy, linewidth=0.5, color=self.target.color)
    if DOS.ispin == 2:
      ax.plot(self.dosdn, DOS.energy, linewidth=0.5, color=self.target.color)
      line=ax.fill_betweenx(DOS.energy, self.dos, self.dosdn, where=PDOS.occ
                      ,facecolor=self.target.color, alpha=0.5, label=self.target.name)
    else:
      line=ax.fill_betweenx(DOS.energy, self.dos, where=PDOS.occ
                      ,facecolor=self.target.color, alpha=0.5, label=self.target.name)

    return line
