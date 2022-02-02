#!/usr/bin/env python
import configparser
import os
import csv
import apdft
import apdft.math as apm

class Option:

  def get_input_paths():
    """ Get paths including input file names from the configuration file.

    Returns:
      init_mol_geom_path : A string of the path for the molecular geometry.
      target_mol_path    : A string of the path for the list of target molecules.
    """
    is_file = os.path.isfile('lime.conf')
    if not is_file:
      raise FileExistsError("lime.conf does not excist.")

    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    init_mol_geom_path = lime_conf['design']['design_init_mol_geom_path']
    target_mol_path = lime_conf['design']['design_target_mol_path']

    if init_mol_geom_path == "":
      raise ValueError("design_init_mol_geom_path must be given in lime.conf.")

    if target_mol_path == "":
      raise ValueError("design_target_mol_path must be given in lime.conf.")

    return init_mol_geom_path, target_mol_path


  def get_input_design():
    """ Get an input for design from the configuration file.

    Returns:
      design_target_property : A string of the target property to be designed.
    """
    is_file = os.path.isfile('lime.conf')
    if not is_file:
      raise FileExistsError("lime.conf does not excist.")

    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    design_target_property = lime_conf['design']['design_target_property']

    if design_target_property not in ['atomization_energy']:
      raise ValueError("design_target_property must be atomization_energy in lime.conf.")

    return design_target_property


  def get_free_atom_energies():
    path_data = os.path.dirname(__file__).replace('src/inverse_design', 'data')
    path_free_atom_energies = "%s%s" % (str(path_data), "/free_atom_energies.csv")
    file_free_atom_energies = open(path_free_atom_energies, "r")

    dict_free_atom_energies = csv.DictReader(file_free_atom_energies)
    dict_dict_free_atom_energies = {}
    dict_dict_free_atom_energies['atom'] = []
    dict_dict_free_atom_energies['atom_energy'] = []

    # TODO: 'G16_CCSD_highspin' need to correspond to the calculation method
    type_atom_energy = 'G16_CCSD_highspin'

    for i, row in enumerate(dict_free_atom_energies):
      dict_dict_free_atom_energies['atom'].append(int(row['atom']))
      dict_dict_free_atom_energies['atom_energy'].append(float(row[type_atom_energy]))

    return dict_dict_free_atom_energies


  def get_inputs():
    """ Get inputs for inverse design

    Returns:
      geom_coordinate : A (the number of molecules, 3) array of the molecular geometry. [Angstrom]
      mol_target_list : A (the number of molecules, the number of atoms) array of the list of target molecules.
                        It can be generated from modified APDFT.
      free_atom_energies : A ('atom', 'atom_energy') dictionary. [nondimension, Hartree]
    """
    init_mol_geom_path, target_mol_path = Option.get_input_paths()

    # Read an initial geometry of a molecule in chemical space
    try:
      nuclear_numbers, geom_coordinate = apdft.read_xyz(init_mol_geom_path)
    except FileNotFoundError:
      print('Unable to open input file "%s".' % init_mol_geom_path, level="error")

    # Read target molecules
    mol_target_list = apm.IntegerPartitions.read_target_molecules(target_mol_path)

    if len(nuclear_numbers) != len(mol_target_list[0]):
      raise ValueError("A list of target molecules needs to have a reference molecule in the first line.")
    for i in range(len(nuclear_numbers)):
      if nuclear_numbers[i] != mol_target_list[0][i]:
        raise ValueError("A list of target molecules needs to have a reference molecule in the first line.")

    # Read free atom energies
    free_atom_energies = Option.get_free_atom_energies()

    return geom_coordinate, mol_target_list, free_atom_energies