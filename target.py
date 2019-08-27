#!/usr/bin/env python3

# Translate target information
orbit_name = {'all':[0,1,2,3,4,5,6,7,8], 's':[0], 'p':[1,2,3], 'd':[4,5,6,7,8]
              , 't2g':[4,5,7], 'eg':[6,8], 'py':[1], 'pz':[2], 'px':[3]
              , 'dxy':[4], 'dyz':[5], 'dz2':[6], 'dxz':[7], 'x2-y2':[8] }

class set_target:
  def __init__(self, target=None, color=None, ispin=None, species=None, spec_index=None):
    self.ion   = self.get_target_ion(target, species)
    self.orb   = self.get_target_orb(target)
    self.name  = self.get_target_name(target, species, spec_index)
    self.color = color

  def get_target_ion(self, target, species):
    ion = set()
    if not isinstance(target[0],list): target[0] = [ target[0] ]
    for i in target[0]:
      if isinstance(i,int): ion.add(i)
      else:
        for j, s in enumerate(species):
          if s == i: ion.add(j)
    return list(ion)

  def get_target_orb(self, target):
    if len(target) == 1: target.append('all')
    orb = set()
    if not isinstance(target[1],list): target[1] = [ target[1] ]
    for i in target[1]:
      if isinstance(i,int): orb.add(i)
      else:
        if i.lower() in orbit_name.keys(): orb.update(orbit_name[i])
        else:
          print("No orbital tag")
          exit()

    return list(orb)

# Make legend name from low data
  def get_target_name(self, target, species, spec_index):
# atom's name for legend
# ion index[atomic symbol] = [ions]
    ion_index = {}
    for i in self.ion:
      if species[i] in ion_index.keys(): ion_index[species[i]].append(i)
      else                             : ion_index[species[i]] = [i]

    name = ''
    for i in ion_index.keys():
      name += i
# if ion contains all ions of a species, target name -> species
# else, target name -> species(ion)
      if ion_index[i] != spec_index[i]:
        name += '('
        for j in ion_index[i]:
          name += str(j) + ','
        name = name[:-1] + ')'

# orbital's name for legend
    tempo = self.orb[:]
# if targeto contains all orbitals, target name -> (only) species (ion)

# else. target name -> _orbitals
    if tempo != orbit_name['all']:
      name += '$_{'
# From large classification, if targeto included by orbit_name,
# orb - orbit_name's orbitals, and name += orbit_name.
      for i,n in orbit_name.items():
        if n == orbit_name['all']: continue
        if len(tempo) == 0       : break
        flag = True
        for l in n:
          if l not in tempo: flag = False
        if flag == True:
          name += i + ','
          for l in n:
            tempo.remove(l)
      name = name[:-1]
      name += '}$'
    return name
