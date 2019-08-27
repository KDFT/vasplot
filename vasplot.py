#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import random
from target import set_target
from dos.DOSplot import DOSplot
from band.bandplot import bandplot

class vasplot:
  def __init__(self, target=None, Fdos=None, Fband=None, Fkpts=None, e_range=[-3,3], Fsave=None, spin_proj=None, x_range=None, e_ticks=None):
# Read basic information
    ispin, fermi, species = self.read_parm(Fdos, Fband)
    if spin_proj != None:
      if spin_proj.lower() == 'x': spin_proj = 1
      if spin_proj.lower() == 'y': spin_proj = 2
      if spin_proj.lower() == 'z': spin_proj = 3

# Set target
    target, title, lgd = self.pretreatment(target, species, ispin, spin_proj)

# Get figure
    self.get_plot(target, Fdos, Fband, Fkpts, e_range, fermi, ispin, title, lgd, spin_proj, x_range, e_ticks)

# Save figure
    self.save_plot(Fsave)

# Read basic information
  def read_parm(self, Fdos=None, Fband=None):
    if Fdos != None   : root = ET.parse(Fdos).getroot()
    elif Fband != None: root = ET.parse(Fband).getroot()
    else:
      print("Please give me a vasprun.xml")
      exit()

    parmElec_tag = root.find("parameters").find('separator[@name="electronic"]')
    spin_tag = parmElec_tag.find('separator[@name="electronic spin"]')
    if spin_tag.find('i[@name="LNONCOLLINEAR"]').text == ' F  ':
      ispin=int(spin_tag.find('i[@name="ISPIN"]').text)
    else: ispin = 4

    fermi = float(root.find("calculation").find("dos").find("i").text)   # fermi level

# Species[ion] = atomic symbol
    species = []
    for i in root.find('atominfo').find('array').find('set'):
      species.append(i[0].text.strip())

    return ispin, fermi, species

  def pretreatment(self, target, species, ispin, spin_proj):
    title, spec_index = self.get_title(species, ispin)  # Make title

    target = self.target2list(target, spec_index)       # Pretreat target str
    if target == []: return [], title, False

    color = self.get_color(target, ispin, spin_proj)    # Make color list

# Get target information and make name for legend
    newtarget = []
    for t,c in list(zip(target, color)):
      newtarget.append(set_target(t, c, ispin, species, spec_index))
    self.print_target(newtarget)

    return newtarget, title, True

  def get_title(self, species, ispin):
# spec_index[atomic symbol] = [ions]
    spec_index = {}
    for i,s in enumerate(species):
      if s in spec_index.keys(): spec_index[s].append(i)
      else                     : spec_index[s] = [i]

    title = ''
    for s in spec_index.keys():
      title += s + '$_{' + str(len(spec_index[s]))+  '}$'

    return title, spec_index

  def target2list(self, target, spec_index):
# For single input
    if not isinstance(target,list): target = [ target ]

# For single set input
    for t in target:
      if not isinstance(t,list) and not isinstance(t,int) and t not in spec_index.keys():
        target = [ target ]

# For mixed dimension input
    for i in range(len(target)):
      if not isinstance(target[i],list):
        target[i] = [ target[i] ]

# delete taget if the first term is not atom
    for i,t in enumerate(target):
      if not isinstance(t[0],list): temp = [ t[0] ]
      else                        : temp = t[0]
      for j in temp:
        if not isinstance(j,int) and j not in spec_index.keys():
          del target[i]
    return target

  def get_color(self, target, ispin, spin_proj):
# Choose project color
    color=list(mcolors.TABLEAU_COLORS.keys())
    if ispin==2 or ispin==4 and spin_proj != None:
      color.remove('tab:blue')
      color.remove('tab:red')

    if len(target) > len(color):
      newcolor = list(mcolors.CSS4_COLORS.keys())
      for i in range(len(color), len(target)):
        num = random.randint(0, len(newcolor)-1)
        color.append(newcolor.pop(num))

    for t in target:
      if len(t) == 3:
        for i,c in enumerate(color):
          if mcolors.same_color(t[2],c):
            del color[i]

    for i, t in enumerate(target):
      if len(t) == 3:
        color.insert(i,t[2])

    return color

# Show legend name
  def print_target(self, target):
# For remove $_{ }$
    name = []
    for t in target:
      if '$' in t.name:
        temp = t.name.find('$')
        name.append(t.name[:temp] + '_' + t.name[temp+3:-2])
      else:
        name.append(t.name)

    lentn = max([ len(str(n)) for n in name ])
    lenti = max([ len(str(t.ion)) for t in target ])
    lento = max([ len(str(t.orb)) for t in target ])
    lentc = max([ len(str(t.color)) for t in target ])

    if lentn < 7: lentn = 7
    if lenti < 6: lenti = 6
    if lento < 8: lento = 8
    if lentc < 6: lentc = 6

    print(' ' * (lentn- 6) + 'Target' + ' ' * (lenti - 4) + 'atoms'
          + ' ' * (lento - 6) + 'orbital' + ' ' * (lentc - 4) + 'color')
    for n, t in list(zip(name, target)):
      data =  ' ' * (lentn - len(str(n)) ) + str(n) + ' '
      data += ' ' * (lenti - len(str(t.ion)) ) + str(t.ion) + ' '
      data += ' ' * (lento - len(str(t.orb)) ) + str(t.orb) + ' '
      data += ' ' * (lentc - len(str(t.color)) ) + str(t.color)
      print(data)

    return

  def get_plot(self, target, Fdos, Fband, Fkpts, e_range, fermi, ispin, title, lgd, spin_proj, x_range, e_ticks):
    if Fdos == None:
      fig, ax = plt.subplots(1,1)
      bandplot(ax, Fband, Fkpts, target, e_range, fermi, ispin, lgd, spin_proj, x_range, e_ticks)

    elif Fband == None:
      fig, ax = plt.subplots(1,1)
      DOSplot(ax, Fdos, target, e_range, fermi, ispin, lgd, e_ticks)

    else:
      fig, ax = plt.subplots(1,2, gridspec_kw = {'width_ratios':[3, 1], 'wspace':0.1})
      bandplot(ax[0], Fband, Fkpts, target, e_range, fermi, ispin, lgd, spin_proj, x_range, e_ticks)
      DOSplot(ax[1], Fdos, target, e_range, fermi, ispin, False, e_ticks)

    fig.suptitle(title, y=0.92)

    return

  def save_plot(self, Fsave):
    if Fsave == None:
      print("Show")
      plt.show()
      return

    if not isinstance(Fsave,list): Fsave = [ Fsave ]

    for i in Fsave:
      if i.lower() == 'show':
        print("Show")
        plt.show()
      else:
        print("Save : " + i)
        plt.savefig(i,dpi=300, bbox_inches='tight')

    return
