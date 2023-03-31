# Usage

Usage of Hornussen is described.  
This document is being updated.

Last update: March 28, 2023

<br/>

## Inputs & Outputs

### Inputs:
  - lime.conf  
  - energies.csv  
  - ver_atomic_forces.csv  
  - init.xyz : it should be same with *.xyz

**template/**  
  - apdft.conf  
  - imp_mod_cli1.sh  
  - imp_mod_cli2.sh  

**inputs/**  
  - *.xyz : specified by lime.conf  
  - target_molecules.inp : specified by lime.conf  

<br/>

### Outputs:
  - design_opt.dat  
  - elapsed_time.dat  
  - geom_opt_hist/  

<br/>

## Restart option

The following outputs of the terminated design should be retained:  
  - design_opt.dat  
  - restart_design.dat  
  - geom_opt_hist/  

<br/>

## Preprocessing

1. Prepare *.xyz for initial molecular geometry of a reference molecule.

2. Specify candidate molecules in chemical space.  
  To get unique candidate molecules, the following parameters should be
  given in apdft.conf:
    - apdft_maxdz,  
    - conf.apdft_specifytargets,  
    - conf.apdft_targetpositions,  
    - apdft_readtargetpath  

  apdft.conf  
  - apdft_includeonly      : it specifies mutated atoms.  
  - apdft_maxdz            : for accurate calculations, it should be
                             the number of mutated atoms.  
  - apdft_specifytargets   : whether to specify target molecules or not.  
  - apdft_targetpositions  : target atom positions to be mutated.  
  - apdft_readtargetpath   : path including the file name.  

    e.g., when a reference molecule is benzene  
    benzene.xyz:  
    C         -2.09726        2.41992        0.00000  
    C         -0.69947        2.47902       -0.00000  
    C          0.05061        1.29805       -0.00000  
    C         -0.59710        0.05797       -0.00000  
    C         -1.99490       -0.00113        0.00000  
    C         -2.74498        1.17984        0.00000  
    H         -0.19838        3.43838       -0.00000  
    H          1.13198        1.34377       -0.00000  
    H         -0.01682       -0.85566       -0.00000  
    H         -2.49598       -0.96049        0.00000  
    H         -3.82635        1.13412        0.00000  
    H         -2.67755        3.33356        0.00000  

    apdft_maxdz = 6  
    apdft_specifytargets = True  
    apdft_targetpositions = 0,1,2,3,4,5  
    apdft_readtargetpath = None  

3. Make an output target_molecules.inp of a list of target molecules which do not have
   equivalent molecules.

4. Perform geometry optimization of a reference molecule by using an electronic structure
   method which will be combined with APDFT and make mol.xyz with optimized geometry.

<br/>

## Design

1. Read inputs  
  All inputs of Lime is read by using a configuration file, lime.conf.  
    design_init_mol_geom_path  : path including the file name for an initial
    molecular geometry,  
    e.g., /home/test/benzene.xyz  
    design_target_mol_path     : path including the file name for a list of
    target molecules,  
    e.g., /home/test/target_molecules.inp.  
    It can be generated using modified APDFT.

2. Generate participation coefficients and normalized participation coefficients

3. ...

<br/>

## How to run efficient calculations?

  Current version of Hornussen may need some code modulations

  1. If reading guess is needed, change pyscf2.py in APDFT/src/apdft/calculator/templates.
  2. Move required CSV files.
  3. Intra-node parallerization of APDFT: flag_ap_smp = True for APDFT
     calculations in physics.py in APDFT/src/apdft.
  4. Inter-node parallerization of QM: set num_mpi_proc in ase_apdft.py in APDFT/src/apdft/ase.
  5. Remove profiler options of cli.py in APDFT/src
  6. Prepare lime.conf
  7. Prepare apdft.conf

  3 can not be used for general property design with optimality criteria.

<br/>

## Applicability

  In lime.conf
  - design_method: standard or optimality_criteria
  - design_target_property: atomization_energy, total_energy, or ele_dipole  
    - Atomization energy (atomization_energy) is maximized  
    - Total energy (total_energy) is minimized  
    - Electrictric dipole strength (ele_dipole) is maximized

  Only single equilibrium geometry can be predicted and accounted in the target
  properties in the current version.