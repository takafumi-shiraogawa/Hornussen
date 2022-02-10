Usage
===

Usage of Lime is described.


Flowchart
-------------------

1. Preprocessing
  1.1. Get a list of target molecules (target_molecules.inp) by modified APDFT.
  1.2. Geometry optimization of the reference molecule by electronic structure theory
       to be used in modified APDFT.


Preprocessing
-------------------

1. Prepare *.xyz for initial molecular geometry of a reference molecule.

2. Specify target molecules in chemical space.

  In apdft.conf,
    apdft_maxdz,
    conf.apdft_specifytargets,
    conf.apdft_targetatom,
    conf.apdft_targetpositions,
    apdft_readtargetpath
  should be specified.
  We can get target molecules which do not have equivalent molecules.
  This strategy does not affect the number of QM calculations.
  To reduce the number of QM calculations, some developments are needed.

  apdft.conf
    apdft_includeonly      : it specifies mutated atoms.
    apdft_maxdz            : for accurate calculations, it should be
                             the number of mutated atoms.

    apdft_specifytargets   : whether to specify target molecules or not.
    apdft_targetatom       : a target atom type.
    apdft_targetpositions  : target atom positions to be mutated.
    apdft_readtargetpath   : path including the file name.

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
  apdft_targetatom = 6
  apdft_targetpositions = 0,1,2,3,4,5
  apdft_readtargetpath = None

  Note for further developments:
    If apdft_targetpositions does not cover all the atoms of a reference molecule,
    the number of QM calculations is larger than the required ones.
    It may be possible to use apdft_includeonly to specify atoms to be mutated,
    which affects target molecules and QM calculations, to reduce the cost.
    However, the "energies_geometries" mode does not correspond to "apdft_includeonly".

3. Make an output target_molecules.inp of a list of target molecules which do not have
   equivalent molecules.

4. Perform geometry optimization of a reference molecule by using an electronic structure
   method which will be combined with APDFT and make mol.xyz with optimized geometry.


Design
-------------------

1. Read inputs
  All inputs of Lime is read by using a configuration file, lime.conf.
    design_init_mol_geom_path  : path including the file name for an initial molecular geometry,
                                 e.g., /home/test/benzene.xyz
    design_target_mol_path     : path including the file name for a list of target molecules,
                                 e.g., /home/test/target_molecules.inp. It can be generated using
                                 modified APDFT.

2. Generate participation coefficients and normalized participation coefficients