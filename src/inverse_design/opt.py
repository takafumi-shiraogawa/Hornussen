#!/usr/bin/env python
import numpy as np

class line_searcher():
  """ Perform line searches for sequential optimization methods."""

  # Perform a line search by a steepest descent minimization manner
  def steepest_descent(variable, gradient, grad_scale_factor = 1.0):
    """ Calculates updated variables.

		Args:
			variable:		A N array of variables
      gradient:      A N array of derivatives of objective function with respect to given variables
			grad_scale_factor:			A scale factor of gradients
		Returns:
			next_variable: A N array of updated variables
		"""
    next_variable = np.zeros(len(variable))
    for i in range(len(variable)):
      next_variable[i] = variable[i] - (grad_scale_factor * gradient[i])

    return next_variable


class optimality_criteria():
  """ Perform optimality criteria method. """

  def calc_scale_factor(variables, object_functions, multiplier, penalty_factor):
    """ Calculate scale factor D
    Args:
      variables        : A N array of variables.
      object_functions : A N array of objective functions.
      multiplier       : A scalar of a Lagrangian multiplier.
      penalty_factor   : A scalar of a penalty factor.
    Returns:
      A N array of the scale factor D.
    """
    return np.multiply(variables ** (penalty_factor - 1.0), object_functions) * penalty_factor / multiplier