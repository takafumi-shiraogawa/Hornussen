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