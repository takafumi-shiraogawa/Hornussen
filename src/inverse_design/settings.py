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
      raise FileExistsError("lime.conf does not exist.")

    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    init_mol_geom_path = lime_conf['design']['design_init_mol_geom_path']
    target_mol_path = lime_conf['design']['design_target_mol_path']

    if init_mol_geom_path == "":
      raise ValueError("design_init_mol_geom_path must be given in lime.conf.")

    if target_mol_path == "":
      raise ValueError("design_target_mol_path must be given in lime.conf.")

    return init_mol_geom_path, target_mol_path


  def get_input_params():
    """ Get input parameters. """

    is_file = os.path.isfile('lime.conf')
    if not is_file:
      raise FileExistsError("lime.conf does not exist.")

    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    try:
      perturb_ampli = lime_conf['design']['design_perturb_ampli']
    except:
      perturb_ampli = 0.05

    try:
      max_design_opt_iter = lime_conf['design']['design_max_opt_iter']
    except:
      max_design_opt_iter = 1000

    try:
      design_opt_criter = lime_conf['design']['design_opt_criter']
    except:
      design_opt_criter = 0.0001

    return float(perturb_ampli), int(max_design_opt_iter), abs(float(design_opt_criter))


  def get_input_design():
    """ Get an input for design from the configuration file.

    Returns:
      design_target_property : A string of the target property to be designed.
      design_calc_level      : A string that specifies a calculation level of free atom energies.
    """
    is_file = os.path.isfile('lime.conf')
    if not is_file:
      raise FileExistsError("lime.conf does not exist.")

    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    design_target_property = lime_conf['design']['design_target_property']

    try:
      design_restart = lime_conf['design']['design_restart']
      if design_restart == "True":
        design_restart = bool(design_restart)
      elif design_restart == "False":
        design_restart = bool("")
    except:
      design_restart = bool("")

    try:
      design_flag_scale_gradient = lime_conf['design']['design_scale_gradient']
      if design_flag_scale_gradient == "True":
        design_flag_scale_gradient = bool(design_flag_scale_gradient)
      elif design_flag_scale_gradient == "False":
        design_flag_scale_gradient = bool("")
    except:
      design_flag_scale_gradient = bool("")

    try:
      design_geom_optimizer = lime_conf['design']['design_geom_optimizer']
      if design_geom_optimizer not in ['BFGS', 'BFGSLineSearch', 'STEEPEST_DESCENT']:
        raise ValueError(
          "design_geom_optimizer must be BFGS or BFGSLineSearch or STEEPEST_DESCENT in lime.conf."
        )
    except:
      design_geom_optimizer = 'BFGSLineSearch'

    if design_target_property not in ['atomization_energy', 'total_energy']:
      raise ValueError(
          "design_target_property must be atomization_energy or total_energy in lime.conf.")

    if design_target_property == 'atomization_energy':
      try:
        design_calc_level = lime_conf['design']['design_calc_level']
      except:
        raise ValueError(
            "design_calc_level for atomization energy must be given in lime.conf.")
    else:
      design_calc_level = None

    return design_target_property, design_restart, design_calc_level, design_flag_scale_gradient, design_geom_optimizer


  def get_free_atom_energies(design_calc_level):
    """ Get free atom energies.

    Args:
      design_calc_level  : A string that specifies a calculation level of free atom energies.
    Returns:
      free_atom_energies : A ('atom', 'atom_energy') dictionary. [nondimension, Hartree]
    """
    if design_calc_level.lower() == 'ccsd/def2-tzvp' or design_calc_level.lower() == 'ccsd/def2tzvp':
      type_atom_energy = 'pyscf_CCSD_highspin'
    elif design_calc_level.lower() == 'hf/def2-tzvp' or design_calc_level.lower() == 'hf/def2tzvp':
      type_atom_energy = 'pyscf_HF_highspin'
    elif design_calc_level.lower() == 'pbe0/def2-tzvp' or design_calc_level.lower() == 'pbe0/def2tzvp':
      type_atom_energy = 'pyscf_PBE0_highspin'
    elif design_calc_level.lower() == 'b3lyp/def2-tzvp' or design_calc_level.lower() == 'b3lyp/def2tzvp':
      type_atom_energy = 'pyscf_B3LYP_highspin'
    else:
      raise ValueError("Free atom energies in %s have not been calculated yet." % design_calc_level.lower())

    path_data = os.path.dirname(__file__).replace('src/inverse_design', 'data')
    path_free_atom_energies = "%s%s" % (str(path_data), "/free_atom_energies.csv")
    file_free_atom_energies = open(path_free_atom_energies, "r")

    dict_free_atom_energies = csv.DictReader(file_free_atom_energies)
    dict_dict_free_atom_energies = {}
    dict_dict_free_atom_energies['atom'] = []
    dict_dict_free_atom_energies['atom_energy'] = []

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

    # if len(nuclear_numbers) != len(mol_target_list[0]):
    #   raise ValueError("A list of target molecules needs to have a reference molecule in the first line.")
    # for i in range(len(nuclear_numbers)):
    #   if nuclear_numbers[i] != mol_target_list[0][i]:
    #     raise ValueError("A list of target molecules needs to have a reference molecule in the first line.")

    return geom_coordinate, mol_target_list


  def get_debug_params():
    """ Get debug parameters. """

    is_file = os.path.isfile('lime.conf')
    if not is_file:
      raise FileExistsError("lime.conf does not exist.")

    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    try:
      debug_design = lime_conf['debug']['debug_design']
      if debug_design == "True":
        debug_design = bool(debug_design)
      elif debug_design == "False":
        debug_design = bool("")
      else:
        raise ValueError("debug_design is invalid.")
    except:
      debug_design = bool("")

    return debug_design
