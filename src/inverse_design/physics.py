import numpy as np


class Calc_Prop():
  """ Calculate properties of molecules in chemical space """

  def calc_sum_free_atom_energies(mol_target_list, free_atom_energies):
    """ Calculate sums of free atom energies of each molecule in chemical space
    Args:
      mol_target_list    : A (the number of molecules, the number of atoms) array
                           of a list of target molecules.
      free_atom_energies : A ('atom', 'atom_energy') dictionary including H, B, C, N, O
                           free atom energies. [nondimensional, Hartree]
    Returns:
      sum_free_atom_energies : A (the number of molecules) array of sums of
                               free atom energies of each molecule in chemical space. [Hartree]
    """
    # Get the number of molecules
    num_mol = len(mol_target_list)

    sum_free_atom_energies = np.zeros(num_mol)

    for i, mol in enumerate(mol_target_list):
      for j, atom in enumerate(mol):
        sum_free_atom_energies[i] += free_atom_energies['atom_energy'][free_atom_energies['atom'].index(atom)]

    return sum_free_atom_energies


  def calc_atomization_energies(energies, sum_free_atom_energies):
    """ Calculate atomization energies of each molecule in chemical space
    Args:
      energies               : A (the number of molecules) array of
                               potential energies of each molecule in chemical space. [Hartree]
      sum_free_atom_energies : A (the number of molecules) array of sums of
                               free atom energies of each molecule in chemical space. [Hartree]
    Returns:
      atomization_energies : A (the number of molecules) array of atomization energies
    """
    # Calculate atomization energies
    # Definition of atomization energy is (sum of atom energies - potential energy of a molecule)
    atomization_energies = sum_free_atom_energies - energies

    return atomization_energies