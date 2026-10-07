Remove water residues from the protein structure

from vmd menu, open plugins -> tkconsole to remove water residues, typing
`set a [atomselect top protein]`

and create a new file with `$a writepdb filename.pdb`

then exit the tkconsole and you can open the new file

script: ubq.pgn
```
package require psfgen

; provide the topology list of the top 27 protein and lipid
topology top_all27_prot_lipid.inp

; istidine in different protonation states, so alias for easy use
; indistinguishable in cristallography. the most common is HSE
pdbalias residue HIS HSE

; isoleucine
pdbalias atom ILE CD1 CD

; generate chain of atom based on the pdb file
segment U {pdb ubqp.pdb}

; assign coordinates to atoms against the pdb file
coordpdb ubqp.pdb U

; guess coordinates for missing atoms
guesscoord

;create new files because we added new atoms (in particular coordinates for missing)
writepdb ubq.pdb
writepsf ubq.psf
```

---

#### top_all27_prot_lipid.inp

describe topology of amminoacids' building atoms, with definitions of atoms (matching the ones in the pdb file), their types, partial charges.
Can be used in any pdb, so I don't have to describe the same stuff each time in different pdb files.

e.g. alanine (ALA, the first in the file), the hydrogens hb1/2/3 have the same type HA despite having different names.

bonds between atoms are defined. + before the atom name.

impr for improper torsions

coordinates part denoted by IC:
> distance first-second, angle between first-second-third, angle between second-third-fourth, ..., distance n-1 n.

---

Import the script file sourcing in vkconsole `source ubq.pgn`

psf file: atom number, chain name, residue name, atom name, atom type, charge

- nbond section: defines the bonds between atoms by their number

can open pdb with its psf in the same time with `vmd psffile.psf pdbfile.pdb`

add water: in tkconsole with `package require solvate` and then `solvate psffile.psf pdbfile.pdb -t 5 -o ubq_wb`
- -t 5: add a layer of at least 5 Å of water around the protein
- -o: output file basename (generates psf and pdb files)

> apply filter to visualize better (e.g. newcartoon)

> in periodic boundary conditions, so we have distance between atoms of 5 + 5 = 10 Å, and tipically this is the cutoff distance.

create new files to compensate protein charge via menu: modelling -> add ions
