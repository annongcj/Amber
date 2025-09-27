#!/usr/bin/env python

import parmed as pmd

print(pmd.__version__)

# change modifed names to conventional name so gbneck2nu can correctly assign parameters

parm = pmd.load_file('../prmtop')
print(set(res.name for res in parm.residues))

amap = {
        'CFZ' : 'C',
        'GF2' : 'G',
        'AF2' : 'A',
        'UFT' : 'U',
}

for res in parm.residues:
    if res.name in amap:
        res.name = amap[res.name]

print(set(res.name for res in parm.residues))
parm.save('fake_prmtop.parm7')
