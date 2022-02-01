import csv
import numpy as np


def get_property_values(property_name, dict, num_mol, apdft_order = 1):
  """ Read property values """
  property_values = np.zeros(num_mol)
  for i, row in enumerate(dict):
    property_values[i] = np.array(
        row["%s%s%s" % (str(property_name), "_order", str(apdft_order))], dtype=np.float64)

  return property_values


class APDFT_Proc():
  """ APDFT processors for inverse design """

  def __init__(self, num_mol):
    self._num_mol = num_mol

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