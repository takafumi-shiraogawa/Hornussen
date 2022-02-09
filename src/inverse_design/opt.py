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
    return variable - grad_scale_factor * gradient


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


  def calc_scaled_variables(variables, scale_factors_D, damp_factor = 0.5):
    """ Calculate scaled variables by the scale factors D with the damping coefficient. """
    if damp_factor == 0.5:
      return np.multiply(np.sqrt(scale_factors_D), variables)
    else:
      return np.multiply(scale_factors_D ** damp_factor, variables)


  def update_variables(variables, scaled_variables):
    """ Update variables. """
    # Set allowable change in each variable
    # This value is obtained from O. Sigmund, A 99 line topology optimization code
    # written in Matlab, 2001.
    allow_change = 0.2

    for i in range(len(variables)):
      if scaled_variables[i] <= max(0.0, variables[i] - allow_change):
        variables[i] = max(0.0, variables[i] - allow_change)
      elif scaled_variables[i] >= min(1.0, variables[i] + allow_change):
        variables[i] = min(1.0, variables[i] + allow_change)
      elif scaled_variables[i] > max(0.0, variables[i] - allow_change) and scaled_variables[i] < min(1.0, variables[i] + allow_change):
        variables[i] = scaled_variables[i]
      else:
        raise Exception("Strange behavior!")

    return variables