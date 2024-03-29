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
    new_norm_part_coeff = Design_Tools.norm_part_coeff(new_part_coeff)

    if len(idx_nonzero) != 1:
      raise ValueError("Initial participation coefficients are not valid.")

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


  def calc_change_ratio(part_coeff, norm_part_coeff, gradient, scale_gradient):
    """ Calculate a sum of changes of normalized participation coefficients """

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


  def get_scale_gradient(part_coeff, norm_part_coeff, gradient, design_step):
    """ Get a scale factor which leads to small change (<0.1) of the molecule in the line search """

    # Calculate a scale factor of gradients in the line search
    for i in range(100000):
      # Set a scale factor of gradients in the line search
      # At first iteration, scale_gradient is 1.0.
      scale_gradient = 1.0 * (0.99 ** i)

      # Calculate a change ratio of normalized participation coefficients
      # by the linear search update
      change_ratio = Design_Tools.calc_change_ratio(
          part_coeff, norm_part_coeff, gradient, scale_gradient)

      # If the change ratio of normalized participation coefficients is smaller than 0.05
      if change_ratio <= 0.01:
        break

    if i == 99999:
      raise ValueError("Gradients maybe too large and lead to large molecular change in the update.")

    # Write scale
    with open("scale_gradient.dat", mode='a') as fh:
      print("design step: ", design_step, "iteration: ", i, "change ratio: ", change_ratio, file=fh)
      print("design step: ", design_step, "iteration: ", i, "scale factor: ", scale_gradient, file=fh)
      print("", file=fh)

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
        weight_atomic_force[i, j] = np.sum(np.multiply(atomic_forces[:, i, j], norm_part_coeff))

    return weight_atomic_force


  def get_weight_property_gradient(properties, part_coeff, sign_inv=False):
    """ Calculate weighted property gradient
    Args:
      properties : A (the number of molecules) array of properties of target molecules.
      part_coeff : A (the number of molecules) array of participation coefficients, not normalized one.
      sign_inv   : A boolean for sign inversion of gradients
    Returns:
      weight_property_gradient : A (the number of molecules) array of gradients of properties
                                 with respect to participation coefficients, not normalized one.
    """
    num_mol = len(properties)
    weight_property_gradient = np.zeros(num_mol)

    double_sum_double_part_coeff = (np.sum(np.square(part_coeff))) ** 2.0

    if sign_inv:
      sign_param = -1.0
    else:
      sign_param = 1.0

    for i in range(num_mol):

      double_weight_property_diff = np.sum(np.multiply(properties[i] - properties, np.square(part_coeff)))

      weight_property_gradient[i] = sign_param * -2.0 * part_coeff[i] * \
          (double_weight_property_diff / double_sum_double_part_coeff)

    return weight_property_gradient


class Geom_OPT_Tools():
  """ Tools of inverse design """

  def save_geom_opt_hist(opt_step):
    """ Save a geometry optimization history. """

    path = "geom_opt_hist/geom_opt-%s" % str(opt_step)

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

    if self._design_target_property == 'atomization_energy':
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


  def get_energy_from_file(self, path):
    # Read potential energies of target molecules
    apdft_proc = APDFT_Proc(self._num_target_mol, self._num_atom)
    energies = apdft_proc.read_potential_energies("%s/energies.csv" % path)

    return energies


  def get_ele_dipole_from_file(self, path):
    # Read electric dipoles of target molecules
    apdft_proc = APDFT_Proc(self._num_target_mol, self._num_atom)
    ele_dipoles = apdft_proc.read_ele_dipoles("%s/dipoles.csv" % path)

    return ele_dipoles


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

    return weight_atomization_energy, weight_atomization_energy_gradient, atomization_energies


  def calc_total_energies_and_gradients(energies, norm_part_coeff, part_coeff):
    """" Calculate total energies and its sign-inverted gradients """
    # Calculate weighted total energy
    weight_total_energy = Design_Tools.get_weight_property(
        energies, norm_part_coeff)

    # Calculate sign-inverted gradients of atomization energies with respect to participation coefficients
    weight_total_energy_gradient = Design_Tools.get_weight_property_gradient(
        energies, part_coeff, True)

    return weight_total_energy, weight_total_energy_gradient


  def update_output(self, w_opt_step, norm_part_coeff, weight_design_property, \
    file_name = 'design_opt.dat'):
    """ Make and update an output of the design. """
    with open(str(file_name), 'a') as f:
      print('Step', w_opt_step, file=f)
      # Molecule
      for i in range(self._num_target_mol):
        print("Lime molecule:", 'step%s' % str(w_opt_step),
              'molecule%i' % i, norm_part_coeff[i], file=f)

      # designed property
      print("Lime design property:", 'step%s' % str(w_opt_step),
            weight_design_property, file=f)

      # # Designed property gradients
      # for i in range(self._num_target_mol):
      #   print("Lime design property gradients:", 'step%s' % str(w_opt_step),
      #         'molecule%i' % i, weight_design_property_gradient[i], file=f)

      # Molecular geometry
      for i in range(self._num_atom):
        print("Lime geometry:", 'step%s' % str(w_opt_step),
              'atom%i' % i, *self._geom_coordinate[i, :], file=f)
      print("Lime geometry: -----", file = f)

      print("", file = f)


  def gener_restart_file(self, design_step, part_coeff):
    """ Generate a restart file. """
    with open("restart_design.dat", mode='w') as fh:
      print("Design step:", file=fh)
      print(design_step, file=fh)

      print("Participation coefficients:", file=fh)
      for idx in range(self._num_target_mol):
        print(part_coeff[idx], file=fh)

      print("Molecular geometry:", file=fh)
      for idx in range(self._num_atom):
        print(*self._geom_coordinate[idx, :], file=fh)


  def read_restart_file(self):
    """ Generate a restart file. """
    with open("restart_design.dat", mode='r') as fh:
      lines = fh.readlines()

    part_coeff = []
    coordinates = []

    flag_part_coeff = False
    flag_mol_geom = False

    for idx_line, line in enumerate(lines):
      line = line.strip()

      if ":" in line:
        if line == "Participation coefficients:":
          idx_part_coeff = idx_line
          flag_part_coeff = True
        elif line == "Molecular geometry:":
          idx_mol_geom = idx_line
          flag_mol_geom = True
        continue

      # Read a design step
      if idx_line == 1:
        design_step = int(line)

      if flag_part_coeff:
        # Read participation coefficients
        if idx_line > idx_part_coeff and idx_line <= idx_part_coeff + self._num_target_mol:
          part_coeff.append(float(line))

      if flag_mol_geom:
        # Read nuclear coordinates
        if idx_line > idx_mol_geom:
          parts = line.split()
          coordinates.append([float(_) for _ in parts[0:3]])

    return design_step, np.array(part_coeff), np.array(coordinates)


  def interpolation(self, idx_two_mols, type_interp, geom_opt, design_geom_optimizer, num_div = 10):
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
    elif type_interp == 'a':
      part_coeff = Design_Tools.gener_local_part_coeff(self._num_target_mol, idx_two_mols[0])
      norm_part_coeff_mol1 = Design_Tools.norm_part_coeff(part_coeff)
      part_coeff = Design_Tools.gener_local_part_coeff(self._num_target_mol, idx_two_mols[1])
      norm_part_coeff_mol2 = Design_Tools.norm_part_coeff(part_coeff)
    else:
      raise NotImplementedError("Interpolation type should be 'w' or 'a'.")

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
            self._mol_target_list[0], self._geom_coordinate, interp_norm_part_coeff, design_geom_optimizer)

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
        weight_design_property, weight_design_property_gradient, atomization_energies = Inverse_Design.calc_atomization_energies_and_gradients(
            energies, self._sum_free_atom_energies, interp_norm_part_coeff, part_coeff)

      elif self._design_target_property == 'total_energy':

        # Calculate total energies
        # Calculate weighted total energy
        # Calculate gradients of total energies with respect to participation coefficients
        weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
          energies, interp_norm_part_coeff, part_coeff)

      elif self._design_target_property == 'ele_dipole':
        ele_dipoles = self.get_ele_dipole_from_file(path_inputs)
        str_ele_dipoles = np.linalg.norm(ele_dipoles, axis=1)
        weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
          str_ele_dipoles, interp_norm_part_coeff, part_coeff)

      ### Save results
      if geom_opt:
        # Save work/ for geometry optimization
        Geom_OPT_Tools.save_geom_opt_hist(idx_num_div + 1)

      # Save results of the design
      Inverse_Design.update_output(self, idx_num_div + 1, interp_norm_part_coeff,
                                   weight_design_property, 'interpolation.dat')

    if geom_opt:
      shutil.rmtree("work/")


  def design(self, perturb_ampli, max_design_opt_iter, design_opt_criter, flag_debug_design, \
    flag_design_restart, flag_scale_gradient, design_geom_optimizer, design_method):
    """ Perform inverse design

    Args:
      perturb_ampli : an amplitude of perturbation for participation coefficients.
      max_design_opt_iter : an scalar of the maximum iteration number of design optimization.
      design_opt_criter : an scaler of the design optimization criterion.
      flag_debug_design : a boolean of debug
      flag_design_restart : a boolean of restart
    """

    # Restart design
    if flag_design_restart:
      if os.path.isdir("work/"):
        shutil.rmtree("work/")

      if not os.path.isdir("geom_opt_hist/"):
        raise Exception("To restart design, geom_opt_hist/ is required.")

      if not os.path.isfile('design_opt.dat'):
        raise Exception("To restart design, design_opt.dat is required.")

      if not os.path.isfile('restart_design.dat'):
        raise Exception("To restart design, restart_design.dat is required.")

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
    # energies, atomic_forces, weight_energy, weight_atomic_forces = Inverse_Design.calc_weight_energy_and_atom_forces(
    #     self, '.', norm_part_coeff)
    energies = self.get_energy_from_file('.')

    # If the target property to be designed is atomization energy
    if self._design_target_property == 'atomization_energy':

      ### 3.3. Calculate atomization energies
      ### 3.4. Calculate weighted atomization energy
      ### 3.5. Calculate gradients of atomization energies with respect to participation coefficients
      weight_design_property, weight_design_property_gradient, atomization_energies = Inverse_Design.calc_atomization_energies_and_gradients(
        energies, self._sum_free_atom_energies, norm_part_coeff, part_coeff)

    elif self._design_target_property == 'total_energy':

      ### 3.3. Calculate total energies
      ### 3.4. Calculate weighted total energy
      ### 3.5. Calculate gradients of total energies with respect to participation coefficients
      weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
        energies, norm_part_coeff, part_coeff)

    elif self._design_target_property == 'ele_dipole':
      ele_dipoles = self.get_ele_dipole_from_file('.')
      str_ele_dipoles = np.linalg.norm(ele_dipoles, axis=1)
      weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
        str_ele_dipoles, norm_part_coeff, part_coeff)

    if not flag_design_restart:
      # Remove an old results of the design
      if os.path.isfile('design_opt.dat'):
        os.remove('design_opt.dat')
      if os.path.isfile('scale_gradient.dat'):
        os.remove('scale_gradient.dat')

      # Save results of the design
      Inverse_Design.update_output(
          self, 0, norm_part_coeff, weight_design_property)


    ### 4. Perturb the molecule and calculate weighted properties
    ### 4.1. Calculate weighted potential energy
    # # Save for check change of normalized participation coefficients
    # temp_norm_part_coeff = np.copy(norm_part_coeff)

    # part_coeff, norm_part_coeff = Design_Tools.perturb_part_coeff(part_coeff)
    part_coeff, norm_part_coeff = Design_Tools.redistr_part_coeff(part_coeff, perturb_ampli)

    # Perform geometry optimization
    print("Perform geometry optimization")

    # Set a maximum number of the molecular species
    # TODO: change it from 1. 1 is given for checking performance.
    if flag_debug_design:
      max_w_opt_step = 1
    else:
      max_w_opt_step = max_design_opt_iter

    if flag_design_restart:
      init_w_opt_step, part_coeff, self._geom_coordinate = Inverse_Design.read_restart_file(self)
      norm_part_coeff = Design_Tools.norm_part_coeff(part_coeff)
      if os.path.isdir("geom_opt_hist/geom_opt-%s/" % str(init_w_opt_step + 1)):
        shutil.rmtree("geom_opt_hist/geom_opt-%s/" % str(init_w_opt_step + 1))
    else:
      init_w_opt_step = 0
      # Remove an old directory for saving geometry optimization histories.
      if os.path.isdir("geom_opt_hist/"):
        shutil.rmtree("geom_opt_hist/")

    # Loop for design
    for w_opt_step in range(init_w_opt_step, max_w_opt_step):
      # Note that self._mol_target_list[0] is not used in geometry optimization.
      self._geom_coordinate = ASE_OPT_Interface.imp_ase_opt(
          self._mol_target_list[0], self._geom_coordinate, norm_part_coeff, design_geom_optimizer)

      ### Calculate properties
      # Read energies
      # Read atomic forces
      # Calculate weighted properties
      # Calculate weighted potential energy
      # Calculate weighted atomic forces
      # energies, atomic_forces, weight_energy, weight_atomic_forces = Inverse_Design.calc_weight_energy_and_atom_forces(
      #     self, './work/temp', norm_part_coeff)
      energies = self.get_energy_from_file('./work/temp')

      # If the target property to be designed is atomization energy
      if self._design_target_property == 'atomization_energy':

        # Calculate atomization energies
        # Calculate weighted atomization energy
        # Calculate gradients of atomization energies with respect to participation coefficients
        weight_design_property, weight_design_property_gradient, atomization_energies = Inverse_Design.calc_atomization_energies_and_gradients(
            energies, self._sum_free_atom_energies, norm_part_coeff, part_coeff)

      elif self._design_target_property == 'total_energy':

        # Calculate total energies
        # Calculate weighted total energy
        # Calculate gradients of total energies with respect to participation coefficients
        weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
          energies, norm_part_coeff, part_coeff)

      elif self._design_target_property == 'ele_dipole':
        ele_dipoles = self.get_ele_dipole_from_file('./work/temp')
        str_ele_dipoles = np.linalg.norm(ele_dipoles, axis=1)
        weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
          str_ele_dipoles, norm_part_coeff, part_coeff)

      ### Save results
      # Save work/ for geometry optimization
      Geom_OPT_Tools.save_geom_opt_hist(w_opt_step + 1)

      # Save results of the design
      Inverse_Design.update_output(
          self, w_opt_step + 1, norm_part_coeff, weight_design_property)


      ### Update molecular species
      if design_method == 'standard':
        # Get a scale factor for the gradient
        if not flag_scale_gradient:
          scale_factor_gradient = 1.0
        else:
          scale_factor_gradient = Design_Tools.get_scale_gradient(
              part_coeff, norm_part_coeff, weight_design_property_gradient, w_opt_step + 1)

        # Update participation coefficients and normalized ones
        part_coeff, norm_part_coeff = Design_Tools.update_part_coeff(
            part_coeff, weight_design_property_gradient, scale_factor_gradient)

      elif design_method == 'optimality_criteria':
        old_norm_part_coeff = np.copy(norm_part_coeff)
        if self._design_target_property == 'atomization_energy':
          norm_part_coeff = optimality_criteria.do_optimality_criteria(
              norm_part_coeff, atomization_energies)

        elif self._design_target_property == 'total_energy':
          norm_part_coeff = optimality_criteria.do_optimality_criteria(
              norm_part_coeff, -energies)

        elif self._design_target_property == 'ele_dipole':
          # Here True is for minimization of the strength of the electric dipole
          # moment. If we specify False, it is a design of a molecule with
          # the dipole moment with the large strength.
          norm_part_coeff = optimality_criteria.do_optimality_criteria(
              norm_part_coeff, str_ele_dipoles, False)

      # Save data for restarting design
      # design step, participation coefficients, molecular geometry
      Inverse_Design.gener_restart_file(self, w_opt_step + 1, part_coeff)


      ### Convergence
      if design_method == 'standard':
        if np.max(np.abs(weight_design_property_gradient)) < design_opt_criter:
          flag_design_opt_conv = True
          break
      elif design_method == 'optimality_criteria':
        # In the initial 11 steps, the design does not terminate.
        if w_opt_step == 0:
          criter_file = open(file='criter.dat', mode='a')
        print(w_opt_step, np.sum(np.abs(norm_part_coeff - old_norm_part_coeff)), file=criter_file)
        if w_opt_step > 10:
          if np.sum(np.abs(norm_part_coeff - old_norm_part_coeff)) < 0.01:
            flag_design_opt_conv = True
            break

        # if np.sum(np.abs(norm_part_coeff - old_norm_part_coeff)) < 0.01:
        #   flag_design_opt_conv = True
        #   break

      # Check
      if w_opt_step == max_w_opt_step - 1:
        flag_design_opt_conv = True
        break

    if os.path.isdir("work/"):
      shutil.rmtree("work/")


    ### Calculate a designed molecule
    # If design optimization converges
    if flag_design_opt_conv:
      # Round off normalized participation coefficients
      norm_part_coeff = np.round(norm_part_coeff)

      # Perform geometry optimization
      self._geom_coordinate = ASE_OPT_Interface.imp_ase_opt(
          self._mol_target_list[0], self._geom_coordinate, norm_part_coeff, design_geom_optimizer)

      # Read energies
      # Read atomic forces
      # Calculate weighted properties
      # Calculate weighted potential energy
      # Calculate weighted atomic forces
      # energies, atomic_forces, weight_energy, weight_atomic_forces = Inverse_Design.calc_weight_energy_and_atom_forces(
      #     self, './work/temp', norm_part_coeff)
      energies = self.get_energy_from_file('./work/temp')

      # If the target property to be designed is atomization energy
      if self._design_target_property == 'atomization_energy':

        # Calculate atomization energies
        # Calculate weighted atomization energy
        # Calculate gradients of atomization energies with respect to participation coefficients
        weight_design_property, weight_design_property_gradient, atomization_energies = Inverse_Design.calc_atomization_energies_and_gradients(
            energies, self._sum_free_atom_energies, norm_part_coeff, norm_part_coeff)

      elif self._design_target_property == 'total_energy':

        # Calculate total energies
        # Calculate weighted total energy
        # Calculate gradients of total energies with respect to participation coefficients
        weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
          energies, norm_part_coeff, part_coeff)

      elif self._design_target_property == 'ele_dipole':
        ele_dipoles = self.get_ele_dipole_from_file('./work/temp')
        str_ele_dipoles = np.linalg.norm(ele_dipoles, axis=1)
        weight_design_property, weight_design_property_gradient = Inverse_Design.calc_total_energies_and_gradients(
          str_ele_dipoles, norm_part_coeff, part_coeff)

      # Save work/ for geometry optimization
      Geom_OPT_Tools.save_geom_opt_hist('real')

      # Save results of the design
      Inverse_Design.update_output(
          self, 'real', norm_part_coeff, weight_design_property)

    # If design optimization does not converge
    else:
      with open('design_opt.dat', 'a') as f:
        print('Design optimization did not converge.', file=f)