#!/usr/bin/env python
import configparser
import apdft
import apdft.math as apm

class Option:

  def get_input_paths():
    """ Get paths including input file names from the configuration file.

    Returns:
      init_mol_geom_path : A string of the path for the molecular geometry.
      target_mol_path    : A string of the path for the list of target molecules.
    """
    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    init_mol_geom_path = lime_conf['design']['design_init_mol_geom_path']
    target_mol_path = lime_conf['design']['design_target_mol_path']

    if init_mol_geom_path == "":
      raise ValueError("design_init_mol_geom_path must be given in lime.conf.")

    if target_mol_path == "":
      raise ValueError("design_target_mol_path must be given in lime.conf.")

    return init_mol_geom_path, target_mol_path


  def get_inputs():
    """ Get inputs for inverse design

    Returns:
      geom_coordinate : A (the number of molecules, 3) array of the molecular geometry. [Angstrom]
      mol_target_list : A (the number of molecules, the number of atoms) array of the list of target molecules.
                        It can be generated from modified APDFT.
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

    return geom_coordinate, mol_target_list