import csv
import numpy as np
import time
from ase import Atoms
# from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
import apdft as APDFTtool
import apdft.ase.ase_apdft as APDFT
from apdft.ase.ase_opt import ASE_OPT


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
        except:
          # TODO: exciption handling is required
          pass
        file_atomic_forces.close()

    return atomic_forces


class ASE_APDFT_Interface(APDFT.mod_APDFT):
  """ APDFT-ASE calculators interface of APDFT for Lime's inverse design. """


class ASE_OPT_Interface(ASE_OPT):
  """ APDFT-ASE geometry optimization interface of APDFT for Lime's inverse design. """

  def imp_ase_opt():
    """ Implement ASE geometry optimization.

    How to perform geometry optimization?
      import inverse_design.apdft_interface as APDFT_Interface
      APDFT_Interface.ASE_OPT_Interface.imp_ase_opt()

    Args:
    """

    # Conversion factor from Angstrom to Bohr
    ang_to_bohr = 1 / 0.52917721067
    # Conversion factor from hartree to eV
    har_to_ev = 27.21162

    # Hartree / Bohr to eV / Angstrom
    hb_to_ea = har_to_ev * ang_to_bohr


    start = time.time()

    # coordinates of init.xyz should be in Angstrom.
    nuclear_numbers, coordinates = APDFTtool.read_xyz("init.xyz")

    molstring = ASE_OPT.get_molstring(nuclear_numbers)

    MOL = Atoms(molstring,
              positions=coordinates,
              calculator=ASE_APDFT_Interface())

    # dyn = BFGS(MOL)
    dyn = BFGSLineSearch(MOL)
    dyn.run(fmax=0.005 * hb_to_ea)

    elapsed_time = time.time() - start
    print("elapsed_time:{0}".format(elapsed_time) + "[sec]")

    with open('elapsed_time.dat', 'w') as tfile:
      tfile.write("elapsed_time:{0}".format(elapsed_time) + "[sec]")