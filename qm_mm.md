# Preparation to Make a QM/MM Simulation with GROMACS and DFTB+

First, run a free MM simulation (see `gmx_md.md` for instructions).
Copy the .gro-structure file and the `topol.top` file.

## Preparation of the Topology File:

Create atom type LA (link atom) in the forcefield before `[moleculetype]`:

```
[ atomtypes ]
; name	at.num	mass	charge	ptype	sigma		epsilon
 LA	1	0.0000	0.00000	A	0.000000000000	0.00000	
```

Set charges of the QM atoms, CA and HA in `[atoms]` section to zero and distribute
the excess charge on the backbone N- and C-atom. The net charge should be zero.
An example for a cysteine residue could look like this:

```
; residue  22 CYM rtp CYX  q  0.0
   344          N     22    CYM      N    344    -0.3585      14.01
   345          H     22    CYM      H    345     0.2719      1.008
   346         CT     22    CYM     CA    346     0.0000      12.01
   347         H1     22    CYM     HA    347     0.0000      1.008
   348         CT     22    CYM     CB    348     0.0000      12.01 ; QM
   349         H1     22    CYM    HB1    349     0.0000      1.008 ; QM
   350         H1     22    CYM    HB2    350     0.0000      1.008 ; QM
   351          S     22    CYM     SG    351     0.0000      32.06 ; QM
   352          C     22    CYM      C    352     0.6545      12.01
   353          O     22    CYM      O    353    -0.5679         16   ; qtot 3
```
	
Add link atoms at the end of the `[atoms]` section. Continue the numbering.
Atom type: LA, mass and charge: zero

```
; Link atoms
  3250         LA     22    XXX     LA   3250     0.0000       0.00
  3251         LA      2    XXX     LA   3251     0.0000       0.00
  3252         LA     27    XXX     LA   3252     0.0000       0.00
```

Set the bond type in the `[bonds]` section from 1 to 5 for QM/QM bonds:

```
  346   348     1 ; QM/MM - CA/CB
  346   352     1 
  348   349     5 ; QM - CB/HB1 
  348   350     5 ; QM - CB/HB2 
  348   351     5 ; QM - CB/S 
  351  1687     5 ; QM -  S/S
```

Create `[virtual_sites2]` section after the `[bonds]` section to place the link atoms along
the CA-CB bond in C-H distance to CB (value a):

```
[virtual_sites2]
; Link-atom    QM-atom    MM-atom    funct.    a
  3250          348        346       1         0.72  ; LA-CB-CA (CYX22)
  3251         1684       1682       1         0.72  ; LA-CB-CA (CYX2)
  3252         2140       2138       1         0.72  ; LA-CB-CA (CYX27)
```

Create `[constraints]` section after the `[virtual_sites2]` section to constrain the
CA-CB bond length to 0.153 nm:

```
[constraints]
; QM-atom    MM-atom    func.    b0
  348        346        2        0.153
  1684       1682       2        0.153
  2140       2138       2        0.153
```

## Preparation of the Structure File (.gro):

Add the number of link atoms to the total number of atoms right at the beginning
of the .gro file.

Add the link atoms in the same order as in the topol.top file, so after the last
residue of the protein. This could look something like this:

```
96TRP    HE3 3245  12.000   1.488   5.646
96TRP    CD2 3246  12.072   1.326   5.531
96TRP      C 3247  12.287   1.378   5.819
96TRP    OC1 3248  12.233   1.481   5.868
96TRP    OC2 3249  12.389   1.379   5.745
22XXX     LA 3250   0.000   0.000   0.000
 2XXX     LA 3251   0.000   0.000   0.000
27XXX     LA 3252   0.000   0.000   0.000
106SOL     OW 3250  10.312   5.497   2.429
106SOL    HW1 3251  10.350   5.412   2.406
106SOL    HW2 3252  10.316   5.500   2.524
```

Create an index file. In the menu, select `0` for the whole system and `q` to
save and quit:

```
gmx_d make_ndx -f *.gro -o index.ndx
```

Create the section `[QM_region]` at the end of the index file and add the
numbers of the QM atoms and link atoms, for example:

```
[QM_region]
348 349 350 351
1684 1685 1686 1687
2140 2141 2142 2143
3250 3251 3252
```

Copy the files `dftb_in.hsd` and `single_submit.sh` into the working directory.
The total charge of the QM region needs to be set at the beginning of the
`dftb_in.hsd` file at `Charge = `

Edit the `single_submit.sh` file:

- replace NAME with a name for the simulation

- comment / uncomment the export options

- comment / uncomment the options for plumed

- if plumed file is used: insert `plumed.dat` to `cp dftb_in.hsd md.tpr
  plumed.dat $workdir` in line 23

Choose a fitting .mdp file:

- p01: 50 ps in 0.5 fs steps, no continuation, generates velocities, berendsen
  pressure coupling; -> for initial short simulation / equilibration of the
  system

- p02: 100 ps in 0.5 fs steps, no continuation, generates velocities,
  Parinello-Rahman pressure coupling

- r01: 1 ns in 1 fs steps, continuation, doesn't generate velocities,
  Parinello-Rahman pressure coupling

- r02: 1 ns in 1 fs steps, no continuation, generates velocities,
  Parinello-Rahman pressure coupling

**Note** that in the .mdp file you also need to set the charge of the QM region,
just like in the `dftb_in.hsd` file.

Create the .tpr-file for the QM/MM simulation WITHOUT continuation:

```
gmx_d grompp -f *.mdp -c *.gro -p topol.top -o md.tpr -n index.ndx -maxwarn 1
```

Create the .tpr-file for the QM/MM simulation WITH continuation:
```
gmx_d grompp -f *.mdp -c *.gro -p topol.top -o md.tpr -n index.ndx -t *.cpt -maxwarn 1
```

Submit the simulation:

```
qsub single_submit.sh
```
