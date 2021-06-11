#!/usr/bin/env python3
import xml.etree.cElementTree as ET
from dos.DOS import DOS, PDOS, transpose

class DOSplot:
  def __init__(self, ax, Fdos, target, e_range, fermi, ispin, lgd=False, e_ticks=None):
    print("\n==== PDOS plotting start ====")
    dos_tag = ET.parse(Fdos).getroot().find("calculation").find("dos")
    dos, pdos, dosmax = self.dosrepair(dos_tag, target, e_range, fermi, ispin)

# Set plot environments
    self.set_env(ax, e_range, dosmax, ispin, e_ticks)

    dos.plotDOS(ax)
    legends = []
    for p in pdos:
      line = p.plotPDOS(ax)
      legends.append(line)

# If PDOS only, set legend
    if lgd == True:
      ax.legend(handles=legends,bbox_to_anchor=(0., -0.23, 1., .102)
                , loc=2, ncol=1, mode="expand", borderaxespad=0.)

    print("===== PDOS plotting end =====")

# Read DOSCAR
  def dosrepair(self, dos_tag, target, e_range, fermi, ispin):
    print("======== Read DOSCAR ========")
    dos = DOS(dos_tag, fermi, ispin, e_range)
# Sum for target PDOS
    print("======== Choose DOS =========")
    pdos = []
    for t in target:
      pdos.append(PDOS(t))

# total DOS broadening and cutting
    print("==== Gaussian Broadening ====")
    dos.BroaDOS()
# PDOS broadening and cutting
    for p in pdos:
      p.BroaDOS()

    dosmax = []
    if target == []: dosmax = dos.get_dosmax()
    else:
      for p in pdos: dosmax.append(p.get_dosmax())
      dosmax = max(dosmax)

    print('dosmax : {0:6.3f}'.format(dosmax))

    return dos, pdos, dosmax

  def set_env(self, ax, e_range, dosmax, ispin, e_ticks):
    ax.set_xlabel('DOS[states/eV]')
    if ispin == 2: ax.set_xlim([-dosmax,dosmax])
    else         : ax.set_xlim([0,dosmax])
    ax.set_ylim(e_range)
    if e_ticks != None:
      e_ticks, e_label = transpose(e_ticks)
      ax.set_yticks(e_ticks)
      ax.set_yticklabels(e_label)
    ax.yaxis.tick_right()

    if ispin == 2: ax.set(aspect = 6. * dosmax/(e_range[1]-e_range[0]))
    else         : ax.set(aspect = 3. * dosmax/(e_range[1]-e_range[0]))

    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1)
    ax.axvline(x=0, color='black', linewidth=1)
    return
