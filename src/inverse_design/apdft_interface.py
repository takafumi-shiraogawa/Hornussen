import os
import csv
import numpy as np
import shutil
from ase import Atoms
# from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
import apdft as APDFTtool
import apdft.ase.ase_apdft as APDFT
from apdft.ase.ase_opt import ASE_OPT
import inverse_design.design as Inverse_Design


# Conversion factor from Angstrom to Bohr
ang_to_bohr = 1 / 0.52917721067
# Conversion factor from hartree to eV
har_to_ev = 27.21162

# Hartree / Bohr to eV / Angstrom
hb_to_ea = har_to_ev * ang_to_bohr


def get_property_values(property_name, dict, num_mol, apdft_order = 1):
  """ Read property values """
  property_values = np.zeros(num_mol)
  for i, row in enumerate(dict):
    property_values[i] = np.array(
        row["%s%s%s" % (str(property_name), "_order", str(apdft_order))], dtype=np.float64)

  return property_values


class APDFT_Proc():
  """ APDFT processors for inverse design """

  def __init__(self, num_mol, num_atom):
    self._num_mol = num_mol
    self._num_atom = num_atom

  def read_potential_energies(self, path_potential_energies):
    """ Read potential energy of target molecules
    Args:
      path_potential_energies  : A string of path of APDFT potential energies, e.g., /home/test/energies.csv
    Returns:
      potential_energies       : A (the number of molecules) array of potential energies of target molecules. [Hartree]
    """
    file_total_energies = open(path_potential_energies, "r")
    dict_total_energies = csv.DictReader(file_total_energies)
    potential_energies = get_property_values("total_energy", dict_total_energies, self._num_mol)
    file_total_energies.close()

    return potential_energies

  def read_atomic_forces(self, path_atomic_forces):
    """ Read atomic forces of target molecules
    Args:
      path_atomic_forces  : A string of path of APDFT atomic forces, e.g., /home/test/ver_atomic_forces.csv
    Returns:
      atomic_forces       : A (the number of molecules, the number of atoms, 3) array of atomic forces of target molecules. [Hartree / Bohr]
    """
    atomic_forces = np.zeros((self._num_mol, self._num_atom, 3))
    for i in range(self._num_atom):
      for didx, dim in enumerate('xyz'):
        file_atomic_forces = open(path_atomic_forces, "r")
        dict_atomic_forces = csv.DictReader(file_atomic_forces)
        try:
          atomic_forces[:, i, didx] = get_property_values(
              "ver_atomic_force_%s_%s" % (str(i), str(dim)), dict_atomic_forces, self._num_mol)
        # For one-dimensional calculation (when only z-component is given)
        except KeyError:
          # For z-Cartesian component
          if didx == 2:
            file_atomic_forces.close()
            file_atomic_forces = open(path_atomic_forces, "r")
            dict_atomic_forces = csv.DictReader(file_atomic_forces)
            atomic_forces[:, i, 2] = get_property_values(
              "ver_atomic_force_%s" % str(i), dict_atomic_forces, self._num_mol)
          # For x- and y-Cartesian components
          else:
            atomic_forces[:, i, didx] = 0.0
        file_atomic_forces.close()

    return atomic_forces


class ASE_APDFT_Interface(APDFT.mod_APDFT):
  """ APDFT-ASE calculators interface of APDFT for Lime's inverse design. """

  def __init__(self, norm_part_coeff):
    # Conduct __init__ of APDFT.mod_APDFT
    super().__init__()

    # Read normalized participation coefficients
    self._norm_part_coeff = norm_part_coeff

    # Set the number of target molecules
    self._num_target_mol = len(self._norm_part_coeff)

    # Set the number of atoms in a molecule
    self._num_atom = len(self.nuclear_numbers)


  def calc_weight_energy_and_atom_forces(self, path):
    """ Calculate weight energy and atomic forces by reading
        It is a partly modified version of calc_weight_energy_and_atom_forces
        of Lime's inverse_design.design.Inverse_Design

    Args:
      path            : A string for a path of energies.csv and ver_atomic_forces.csv
                        generated by modified APDFT, not including a file name.
      norm_part_coeff : A (the number of molecules) array of normalized participation
                        coefficients.
    Rerturns:
      weight_energy      : A scalar
      weight_atom_forces : A (the number of molecules) array
    """
    apdft_proc = APDFT_Proc(self._num_target_mol, self._num_atom)

    # Read potential energies of target molecules
    energies = apdft_proc.read_potential_energies("%s/energies.csv" % path)

    # # Check
    # print("energies")
    # print(energies)
    # print("")

    # Read atomic forces of target molecules
    atomic_forces = apdft_proc.read_atomic_forces(
        "%s/ver_atomic_forces.csv" % path)

    # # Check
    # print("atomic_forces")
    # print(atomic_forces[:2])
    # print("")

    # Calculate weighted energy by normalized participation participati
    weight_energy = Inverse_Design.Design_Tools.get_weight_property(
        energies, self._norm_part_coeff)

    # # Check
    # print("weight_energy")
    # print(weight_energy)
    # print("")

    # Calculate weighted atomic forces by normalized participation participati
    weight_atomic_forces = Inverse_Design.Design_Tools.get_weight_atomic_forces(
        atomic_forces, self._norm_part_coeff)

    # # Check
    # print("weight_atomic_forces")
    # print(weight_atomic_forces)
    # print("")

    return weight_energy, weight_atomic_forces


  def read_results(self):
    """ Read calculated energy and atomic forces. """
    path = 'work/temp'

    # TODO: now can handle the full-cartesian atomic force calculation only.
    #       need to be generalized to handle the one-dimensional calculation.
    pot_energy, atom_forces = ASE_APDFT_Interface.calc_weight_energy_and_atom_forces(
        self, path)

    # Save results of geometry optimization
    with open('./work/geom_opt.dat', 'a') as f:
      print("APDFT results:", self.num_opt_step, pot_energy, file = f)
      for i in range(self._num_atom):
        print("APDFT geometry:", self.num_opt_step, *self.atoms.positions[i, :], file = f)
      print("APDFT geometry: -----", file = f)
      print("", file = f)

    self.results = {'energy': pot_energy * har_to_ev,
                    'forces': atom_forces * hb_to_ea,
                    'stress': np.zeros(6),
                    'dipole': np.zeros(3),
                    'charges': np.zeros(self._num_atom),
                    'magmom': 0.0,
                    'magmoms': np.zeros(self._num_atom)}
    APDFT.handle_APDFT.save_results(self.num_opt_step)


class ASE_OPT_Interface(ASE_OPT):
  """ APDFT-ASE geometry optimization interface of APDFT for Lime's inverse design. """

  def imp_ase_opt(nuclear_numbers, coordinates, norm_part_coeff, fmax_au = 0.005):
    """ Implement ASE geometry optimization.

    How to perform geometry optimization?
      import inverse_design.apdft_interface as APDFT_Interface
      APDFT_Interface.ASE_OPT_Interface.imp_ase_opt()

    Args:
      norm_part_coeff  :  A (the number of target molecules) arrray

    Returns:
      MOL._get_positions() : A (the number of atoms, 3) arrray of optimized geometry
                             of a molecule
    """
    # # coordinates of init.xyz should be in Angstrom.
    # nuclear_numbers, coordinates = APDFTtool.read_xyz("init.xyz")

    molstring = ASE_OPT.get_molstring(nuclear_numbers)

    MOL = Atoms(molstring,
              positions=coordinates,
              calculator=ASE_APDFT_Interface(norm_part_coeff))

    # Remove an old results of geometry optimization
    if os.path.isfile('BFGSLineSearch.dat'):
      os.remove('BFGSLineSearch.dat')

    # dyn = BFGS(MOL)
    dyn = BFGSLineSearch(MOL, logfile="BFGSLineSearch.dat")
    dyn.run(fmax=fmax_au * hb_to_ea)

    # Move the output of geometry optimization into the working directory.
    shutil.move("BFGSLineSearch.dat", "./work/")

    return MOL._get_positions()