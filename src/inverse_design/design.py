#!/usr/bin/env python
import numpy as np
import random

class Design_Tools():
  """ Tools of inverse design """

  def gener_init_part_coeff(random_seed, num_mol):
    """ Generate inital participation coeffcients """

    # Set initial participation coefficients
    init_part_coeff = np.zeros(num_mol)

    # Set random seed
    random.seed(random_seed)

    # Generate randomized initial participation coefficients
    for i in range(num_mol):
      init_part_coeff[i] = random.random()

    return init_part_coeff


class Inverse_Design():
  """ Inverse design based on chemical space of geometrically relaxed molecules """