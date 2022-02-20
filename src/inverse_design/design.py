#!/usr/bin/env python
import os
import shutil
import numpy as np
import random
from inverse_design.physics import Calc_Prop
from inverse_design.opt import line_searcher, optimality_criteria
from inverse_design.apdft_interface import APDFT_Proc, ASE_OPT_Interface

class Design_Tools():
  """ Tools of inverse design """

  def norm_part_coeff(part_coeff):
    """ Normalize participation coefficients """

    # Calculate a sum of double of each participation coefficients
    sum_double_part_coeff = np.sum(np.square(part_coeff))

    # Calculate normalized participation coefficients
    return np.square(part_coeff) / sum_double_part_coeff


  def sin_norm_part_coeff(part_coeff):
    """ Normalize participation coefficients """

    # Calculate a sum of double of sin of each participation coefficients
    sum_double_sin_part_coeff = np.sum(np.square(np.sin(part_coeff)))

    # Calculate normalized participation coefficients
    norm_part_coeff = np.square(np.sin(part_coeff)) / sum_double_sin_part_coeff

    return norm_part_coeff


  def gener_local_part_coeff(num_mol, target_mol):
    """ Generate localized participation coefficients """

    # Set localized participation coefficients
    local_part_coeff = np.zeros(num_mol)
    local_part_coeff[target_mol] = 10

    return local_part_coeff


  def gener_random_part_coeff(random_seed, num_mol):
    """ Generate randomized participation coeffcients """

    # Set initial participation coefficients
    init_part_coeff = np.zeros(num_mol)

    # Set random seed
    random.seed(random_seed)

    # Generate randomized participation coefficients
    for i in range(num_mol):
      init_part_coeff[i] = random.random()

    return init_part_coeff


  def perturb_part_coeff(part_coeff, perturb_param = None):
    """ Perturb participation coefficients """

    # If the perturbation parameter is not given
    if perturb_param == None:
      # Get the perturbation parameter
      # 1% of the participation coefficients is distributed to all candidates
      perturb_param = np.sum(part_coeff) * 0.01 / len(part_coeff)

    # Perturb participation coefficients
    new_part_coeff = part_coeff[:] + perturb_param

    # Normalize new participation coefficients
    new_norm_part_coeff = Design_Tools.norm_part_coeff(new_part_coeff)

    return new_part_coeff, new_norm_part_coeff


  def redistr_part_coeff(part_coeff, perturb_ampli = None):
    """ Re-distribute participation coefficients """

    # Get indexes of nonzero components of paticipation coefficients
    idx_nonzero = np.nonzero(part_coeff)

    # Calculate a sum of participation coefficients
    sum_part_coeff = np.sum(part_coeff)

    perturb_param = 0.0
    for i, idx in enumerate(idx_nonzero):
      if perturb_ampli == None:
        perturb_param += part_coeff[idx] * 0.05
      else:
        perturb_param += part_coeff[idx] * perturb_ampli

    new_part_coeff = part_coeff
    for i in range(len(part_coeff)):
      if i not in idx_nonzero:
        new_part_coeff[i] += perturb_param
      else:
        new_part_coeff[i] -= perturb_param

    # Normalize new participation coefficients
    # new_norm_part_coeff = Design_Tools.norm_part_coeff(new_part_coeff)
    if len(idx_nonzero) != 1:
      raise ValueError("Initial participation coefficients are not valid.")
    new_norm_part_coeff = part_coeff / sum_part_coeff

    return new_part_coeff, new_norm_part_coeff


  def get_change_norm_part_coeff(old_norm_part_coeff, new_norm_part_coeff):
    """ Estimate a change of normalized participation coefficients """

    # Calculate a sum of changes of normalized participation coefficients
    return np.sum(np.abs(old_norm_part_coeff - new_norm_part_coeff))


  def update_part_coeff(part_coeff, gradient, scale_gradient = 1.0):
    """ Update participation coefficients by a line search """

    # Perform the steepest descent line search
    new_part_coeff = line_searcher.steepest_descent(
        part_coeff, gradient, scale_gradient)

    # Normalize new participation coefficients
    norm_new_part_coeff = Design_Tools.norm_part_coeff(new_part_coeff)

    return new_part_coeff, norm_new_part_coeff


  def calc_change_ratio(part_coeff, gradient, scale_gradient):
    """ Calculate a sum of changes of normalized participation coefficients """

    # Normalize participation coeffients
    norm_part_coeff = Design_Tools.norm_part_coeff(part_coeff)

    # Perform the steepest descent line search
    new_part_coeff, norm_new_part_coeff = Design_Tools.update_part_coeff(
        part_coeff, gradient, scale_gradient)

    # Get the number of molecules in chemical space
    num_mol = len(part_coeff)

    # Calculate a sum of changes of normalized participation coefficients
    change_ratio = 0.0
    for i in range(num_mol):
      change_ratio += abs(norm_new_part_coeff[i] - norm_part_coeff[i])

    return change_ratio


  def get_scale_gradient(part_coeff, gradient):
    """ Get a scale factor which leads to small change (<0.1) of the molecule in the line search """

    # Calculate a scale factor of gradients in the line search
    for i in range(10000):
      # Set a scale factor of gradients in the line search
      # At first iteration, scale_gradient is 1.0.
      scale_gradient = 1.0 * (0.99 ** i)

      # Calculate a change ratio of normalized participation coefficients
      # by the linear search update
      change_ratio = Design_Tools.calc_change_ratio(part_coeff, gradient, scale_gradient)

      # If the change ratio of normalized participation coefficients is smaller than 0.1
      if change_ratio <= 0.1:
        break

    if i == 9999:
      raise ValueError("Gradients maybe too large and lead to large molecular change in the update.")

    # Write scale
    with open("scale_gradient.dat", mode='w') as fh:
      print("iteration,", i, ", scale factor,", scale_gradient, file=fh)

    return scale_gradient


  def get_weight_property(properties, norm_part_coeff, penalty_factor = 1.0):
    """ Calculate weighted property """
    weight_property = np.sum(np.multiply(properties, np.power(norm_part_coeff, penalty_factor)))

    return weight_property


  def get_weight_atomic_forces(atomic_forces, norm_part_coeff):
    """ Calculate weighted atomic forces """
    num_atom = len(atomic_forces[0])
    weight_atomic_force = np.zeros((num_atom, 3))
    for i in range(num_atom):
      for j in range(3):
        weight_atomic_force[i, j] = np.multiply(atomic_forces[:, i, j], norm_part_coeff).sum()

    return weight_atomic_force


  def get_weight_property_gradient(properties, part_coeff):
    """ Calculate weighted property gradient
    Args:
      properties : A (the number of molecules) array of properties of target molecules.
      part_coeff : A (the number of molecules) array of participation coefficients, not normalized one.
    Returns:
      weight_property_gradient : A (the number of molecules) array of gradients of properties
                                 with respect to participation coefficients, not normalized one.
    """
    num_mol = len(properties)
    weight_property_gradient = np.zeros(num_mol)

    double_sum_double_part_coeff = (np.sum(np.square(part_coeff))) ** 2.0

    for i in range(num_mol):

      double_weight_property_diff = np.sum(np.multiply(properties[i] - properties, np.square(part_coeff)))

      weight_property_gradient[i] = -2.0 * part_coeff[i] * \
          (double_weight_property_diff / double_sum_double_part_coeff)

    return weight_property_gradient


class Geom_OPT_Tools():
  """ Tools of inverse design """

  def save_geom_opt_hist(opt_step):
    """ Save a geometry optimization history. """

    path = "geom_opt_hist/geom_opt-%s" % opt_step

    shutil.copytree('./work', path)


class Inverse_Design():
  """ Inverse design based on chemical space of geometrically relaxed molecules """

  def __init__(
        self,
        geom_coordinate = None,
        mol_target_list = None,
        design_target_property = None,
        free_atom_energies = None
    ):

    self._geom_coordinate = geom_coordinate
    self._mol_target_list = mol_target_list
    self._design_target_property = design_target_property

    # Get the number of target molecules
    self._num_target_mol = len(self._mol_target_list)

    # Get the number of atoms of a molecule
    self._num_atom = len(self._geom_coordinate)

    # Calculate sums of free atom energies of each molecule in chemical space
    self._sum_free_atom_energies = Calc_Prop.calc_sum_free_atom_energies(
        self._mol_target_list, free_atom_energies)


  def calc_weight_energy_and_atom_forces(self, path, norm_part_coeff):
    """ Calculate weight energy and atomic forces by reading
        the results of

    Args:
      path            : A string for a path of energies.csv and ver_atomic_forces.csv
                        generated by modified APDFT, not including a file name.
      norm_part_coeff : A (the number of molecules) array of normalized participation
                        coefficients.
    Rerturns:
      energies           : A (the number of molecules) array
      atomic_forces      : A (the number of molecules, 3) array
      weight_energy      : A scalar
      weight_atom_forces : A (the number of molecules) array
    """
    apdft_proc = APDFT_Proc(self._num_target_mol, self._num_atom)

    # Read potential energies of target molecules
    energies = apdft_proc.read_potential_energies("%s/energies.csv" % path)

    # Read atomic forces of target molecules
    atomic_forces = apdft_proc.read_atomic_forces("%s/ver_atomic_forces.csv" % path)

    # Calculate weighted energy by normalized participation participati
    weight_energy = Design_Tools.get_weight_property(energies, norm_part_coeff)

    # Calculate weighted atomic forces by normalized participation participati
    weight_atomic_forces = Design_Tools.get_weight_atomic_forces(atomic_forces, norm_part_coeff)

    return energies, atomic_forces, weight_energy, weight_atomic_forces


  def calc_atomization_energies_and_gradients(energies, sum_free_atom_energies, norm_part_coeff, part_coeff):
    """" Calculate atomization energies and its sign-inverted gradients """
    # Calculate atomization energies
    atomization_energies = Calc_Prop.calc_atomization_energies(
        energies, sum_free_atom_energies)

    # Calculate weighted atomization energy
    weight_atomization_energy = Design_Tools.get_weight_property(
        atomization_energies, norm_part_coeff)

    # Calculate sign-inverted gradients of atomization energies with respect to participation coefficients
    weight_atomization_energy_gradient = Design_Tools.get_weight_property_gradient(
        atomization_energies, part_coeff)

    return atomization_energies, weight_atomization_energy, weight_atomization_energy_gradient


  def update_output(self, w_opt_step, norm_part_coeff, weight_atomization_energy, \
    weight_atomization_energy_gradient, file_name = 'design_opt.dat'):
    """ Make and update an output of the design. """
    with open(str(file_name), 'a') as f:
      print('Step', w_opt_step, file=f)
      # Molecule
      for i in range(self._num_target_mol):
        print("Lime molecule:", 'step%i' % (w_opt_step),
              'molecule%i' % i, norm_part_coeff[i], file=f)

      # Atomization energy
      print("Lime atomization energy:", 'step%i' % (w_opt_step),
            weight_atomization_energy, file=f)

      # Atomization energy gradients
      for i in range(self._num_target_mol):
        print("Lime atomization energy gradients:", 'step%i' % (w_opt_step),
              'molecule%i' % i, weight_atomization_energy_gradient[i], file=f)

      # Molecular geometry
      for i in range(self._num_atom):
        print("Lime geometry:", 'step%i' % (w_opt_step),
              'atom%i' % i, *self._geom_coordinate[i, :], file=f)
      print("Lime geometry: -----", file = f)

      print("", file = f)


  def interpolation(self, idx_two_mols, type_interp, geom_opt, num_div = 10):
    """ Perform interpolation between two real molecules.

    Args:
      idx_two_mols : A (2) arrays of indexes that specify two molecules to be interpolated.
      type_interp  : A string of a type of the interpolation.
      geom_opt     : A boolean for geometry optimization.
      num_div      : A scalar of the number of divided componenets in the interpolation.
    """

    if type_interp not in ['w', 'a']:
      raise ValueError("type_interp should be w or a in interpolation")

    if type(geom_opt) is not bool:
      raise ValueError("geom_opt should be boolean in interpolation")

    # Remove an old results of the design
    if os.path.isfile('interpolation.dat'):
      os.remove('interpolation.dat')

    # Remove an old directory for saving geometry optimization histories.
    if os.path.isdir("geom_opt_hist/"):
      shutil.rmtree("geom_opt_hist/")

    # Exception handling
    if len(idx_two_mols) != 2:
      raise ValueError("idx_two_mols should be two molecules in interpolation.")

    if type_interp == 'w':
      part_coeff_mol1 = Design_Tools.gener_local_part_coeff(self._num_target_mol, idx_two_mols[0])
      part_coeff_mol2 = Design_Tools.gener_local_part_coeff(self._num_target_mol, idx_two_mols[1])
    else:
      part_coeff = Design_Tools.gener_local_part_coeff(self._num_target_mol, idx_two_mols[0])
      norm_part_coeff_mol1 = Design_Tools.norm_part_coeff(part_coeff)
      part_coeff = Design_Tools.gener_local_part_coeff(self._num_target_mol, idx_two_mols[1])
      norm_part_coeff_mol2 = Design_Tools.norm_part_coeff(part_coeff)

    change = 1.0 / num_div

    for idx_num_div in range(num_div + 1):
      # Perform interpolation of two molecules
      if type_interp == 'w':
        interp_part_coeff = (1.0 - change * idx_num_div) * \
            part_coeff_mol1 + change * idx_num_div * part_coeff_mol2
        interp_norm_part_coeff = Design_Tools.norm_part_coeff(interp_part_coeff)
        # This does not make sense
        part_coeff = interp_part_coeff
      else:
        interp_norm_part_coeff = (1.0 - change * idx_num_div) * \
            norm_part_coeff_mol1 + change * idx_num_div * norm_part_coeff_mol2

      if geom_opt:
        # Perform geometry optimization
        self._geom_coordinate = ASE_OPT_Interface.imp_ase_opt(
            self._mol_target_list[0], self._geom_coordinate, interp_norm_part_coeff)

      ### Calculate properties
      # Read energies
      # Read atomic forces
      # Calculate weighted properties
      # Calculate weighted potential energy
      # Calculate weighted atomic forces
      if geom_opt:
        path_inputs = './work/temp'
      else:
        path_inputs = '.'
      energies, atomic_forces, weight_energy, weight_atomic_forces = Inverse_Design.calc_weight_energy_and_atom_forces(
            self, path_inputs, interp_norm_part_coeff)

      # If the target property to be designed is atomization energy
      if self._design_target_property == 'atomization_energy':

        # Calculate atomization energies
        # Calculate weighted atomization energy
        # Calculate gradients of atomization energies with respect to participation coefficients
        atomization_energies, weight_atomization_energy, weight_atomization_energy_gradient = Inverse_Design.calc_atomization_energies_and_gradients(
            energies, self._sum_free_atom_energies, interp_norm_part_coeff, part_coeff)

      ### Save results
      if geom_opt:
        # Save work/ for geometry optimization
        Geom_OPT_Tools.save_geom_opt_hist(idx_num_div + 1)

      # Save results of the design
      Inverse_Design.update_output(self, idx_num_div + 1, interp_norm_part_coeff,
                                   weight_atomization_energy, weight_atomization_energy_gradient, 'interpolation.dat')

    if geom_opt:
      shutil.rmtree("work/")


  def design(self, perturb_ampli, flag_debug_design, flag_design_restart):
    """ Perform inverse design

    Args:
      perturb_ampli : an amplitude of perturbation for participation coefficients.
      flag_debug_design : a boolean of debug
      flag_design_restart : a boolean of restart
    """

    # Restart design
    if flag_design_restart:
      if os.path.isdir("work/"):
        shutil.rmtree("work/")

      if not os.path.isdir("geom_opt_hist/"):
        raise Exception("To restart design, geom_opt_hist/ is required.")

      if os.path.isfile('design_opt.dat'):
        raise Exception("To restart design, design_opt.dat is required.")

      if os.path.isfile('elapsed_time.dat'):
        raise Exception("To restart design, elapsed_time.dat is required.")

    ### 1. Generate participation coefficients

    # Get localized participation coefficients to a reference molecule
    part_coeff = Design_Tools.gener_local_part_coeff(self._num_target_mol, 0)

    # Normalize participation coefficients
    norm_part_coeff = Design_Tools.norm_part_coeff(part_coeff)


    ### 2. Calculate properties
    ### 2.1. Read energies
    ### 2.2. Read atomic forces
    ### 3. Calculate weighted properties
    ### 3.1. Calculate weighted potential energy
    ### 3.2. Calculate weighted atomic forces
    energies, atomic_forces, weight_energy, weight_atomic_forces = Inverse_Design.calc_weight_energy_and_atom_forces(
        self, '.', norm_part_coeff)

    # If the target property to be designed is atomization energy
    if self._design_target_property == 'atomization_energy':

      ### 3.3. Calculate atomization energies
      ### 3.4. Calculate weighted atomization energy
      ### 3.5. Calculate gradients of atomization energies with respect to participation coefficients
      atomization_energies, weight_atomization_energy, weight_atomization_energy_gradient = Inverse_Design.calc_atomization_energies_and_gradients(
        energies, self._sum_free_atom_energies, norm_part_coeff, part_coeff)

      # Remove an old results of the design
      if os.path.isfile('design_opt.dat'):
        os.remove('design_opt.dat')

      # Save results of the design
      Inverse_Design.update_output(
          self, 0, norm_part_coeff, weight_atomization_energy, weight_atomization_energy_gradient)


    ### 4. Perturb the molecule and calculate weighted properties
    ### 4.1. Calculate weighted potential energy
    # Save for check change of normalized participation coefficients
    temp_norm_part_coeff = np.copy(norm_part_coeff)

    # part_coeff, norm_part_coeff = Design_Tools.perturb_part_coeff(part_coeff)
    part_coeff, norm_part_coeff = Design_Tools.redistr_part_coeff(part_coeff, perturb_ampli)

    # Perform geometry optimization
    print("Perform geometry optimization")

    # Remove an old directory for saving geometry optimization histories.
    if os.path.isdir("geom_opt_hist/"):
      shutil.rmtree("geom_opt_hist/")

    # Set a maximum number of the molecular species
    # TODO: change it from 1. 1 is given for checking performance.
    if flag_debug_design:
      max_w_opt_step = 1
    else:
      max_w_opt_step = 1000

    # Loop for design
    for w_opt_step in range(max_w_opt_step):
      # Note that self._mol_target_list[0] is not used in geometry optimization.
      self._geom_coordinate = ASE_OPT_Interface.imp_ase_opt(
          self._mol_target_list[0], self._geom_coordinate, norm_part_coeff)

      ### Calculate properties
      # Read energies
      # Read atomic forces
      # Calculate weighted properties
      # Calculate weighted potential energy
      # Calculate weighted atomic forces
      energies, atomic_forces, weight_energy, weight_atomic_forces = Inverse_Design.calc_weight_energy_and_atom_forces(
          self, './work/temp', norm_part_coeff)

      # If the target property to be designed is atomization energy
      if self._design_target_property == 'atomization_energy':

        # Calculate atomization energies
        # Calculate weighted atomization energy
        # Calculate gradients of atomization energies with respect to participation coefficients
        atomization_energies, weight_atomization_energy, weight_atomization_energy_gradient = Inverse_Design.calc_atomization_energies_and_gradients(
            energies, self._sum_free_atom_energies, norm_part_coeff, part_coeff)


      ### Save results
      # Save work/ for geometry optimization
      Geom_OPT_Tools.save_geom_opt_hist(w_opt_step + 1)

      # Save results of the design
      Inverse_Design.update_output(
          self, w_opt_step + 1, norm_part_coeff, weight_atomization_energy, weight_atomization_energy_gradient)

      # Save a data to restart design
      # design step, participation coefficients, molecular geometry
      with open("restart_design.dat", mode='w') as fh:
        print("Design step:", file=fh)
        print(w_opt_step, file=fh)

        print("Participation coefficients:", file=fh)
        for idx in range(self._num_target_mol):
          print(part_coeff[idx], file=fh)

        print("Molecular geometry:", file=fh)
        for idx in range(self._num_atom):
          print(*self._geom_coordinate[idx, :], file=fh)

      ### Update molecular species
      part_coeff, norm_part_coeff = Design_Tools.update_part_coeff(
          part_coeff, weight_atomization_energy_gradient)

    shutil.rmtree("work/")