#!/usr/bin/env python
import numpy as np
import random

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


class Inverse_Design():
  """ Inverse design based on chemical space of geometrically relaxed molecules """