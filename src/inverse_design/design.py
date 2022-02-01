#!/usr/bin/env python
import numpy as np
import random
from inverse_design.opt import line_searcher

class Design_Tools():
  """ Tools of inverse design """

  def norm_part_coeff(part_coeff):
    """ Normalize participation coefficients """

    # Get the number of molecules in chemical space
    num_mol = len(part_coeff)

    # Calculate a sum of double of each participation coefficients
    sum_double_part_coeff = 0
    for i in range(num_mol):
      sum_double_part_coeff += part_coeff[i] ** 2.0

    # Set normalized participation coefficients
    norm_part_coeff = np.zeros(num_mol)

    # Calculate normalized participation coefficients
    for i in range(num_mol):
      norm_part_coeff[i] = (part_coeff[i] ** 2.0) / sum_double_part_coeff

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


  def update_part_coeff(part_coeff, gradient, scale_gradient):
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


class Inverse_Design():
  """ Inverse design based on chemical space of geometrically relaxed molecules """

  def __init__(
        self,
        geom_coordinate = None,
        mol_target_list = None
    ):

    self._geom_coordinate = geom_coordinate
    self._mol_target_list = mol_target_list


  def design(self):
    """ Perform inverse design """

    # Get the number of target molecules
    num_target_mol = len(self._mol_target_list)

    # Check
    print("")
    print("the number of target molecules")
    print(num_target_mol)
    print("")

    ### 1. Generate participation coefficients

    # Get localized participation coefficients to a reference molecule
    part_coeff = Design_Tools.gener_local_part_coeff(num_target_mol, 0)

    # Check
    print("initial part_coeff")
    print(part_coeff)
    print("")

    # Normalize participation coefficients
    norm_part_coeff = Design_Tools.norm_part_coeff(part_coeff)

    # Check
    print("norm_part_coeff")
    print(norm_part_coeff)
    print("")