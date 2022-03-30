#!/usr/bin/env python
from functools import wraps
import time
import inverse_design.design as ds
import inverse_design.settings as iconf
# from apdft.settings import Configuration


def stop_watch(func):
  """ Measure time """
  @wraps(func)
  def wrapper(*args, **kargs):
    start = time.time()

    result = func(*args, **kargs)

    elapsed_time = time.time() - start

    print("")
    print("elapsed_time:{0}".format(elapsed_time) + "[sec]")

    with open('elapsed_time.dat', 'w') as tfile:
      tfile.write("elapsed_time:{0}".format(elapsed_time) + "[sec]")

    return result
  return wrapper


@stop_watch
def ignition_design():
  geom_coordinate, mol_target_list = iconf.Option.get_inputs()
  perturb_ampli, max_design_opt_iter, design_opt_criter = iconf.Option.get_input_params()
  design_target_property, flag_design_restart, design_calc_level, flag_scale_gradient, \
    design_geom_optimizer, design_method = iconf.Option.get_input_design()
  if design_target_property == 'atomization_energy':
    free_atom_energies = iconf.Option.get_free_atom_energies(design_calc_level)

  flag_debug_design = iconf.Option.get_debug_params()

  if design_target_property == 'atomization_energy':
    derivatives = ds.Inverse_Design(
        geom_coordinate, mol_target_list, design_target_property, free_atom_energies)
  elif design_target_property == 'total_energy' or design_target_property == 'ele_dipole':
    derivatives = ds.Inverse_Design(
        geom_coordinate, mol_target_list, design_target_property)

  derivatives.design(perturb_ampli, max_design_opt_iter,
                     design_opt_criter, flag_debug_design, flag_design_restart,
                     flag_scale_gradient, design_geom_optimizer, design_method)


@stop_watch
def ignition_interpolation(idx_two_mol, type_interp, geom_opt, num_div):
  geom_coordinate, mol_target_list, free_atom_energies = iconf.Option.get_inputs()
  design_target_property, design_restart = iconf.Option.get_input_design()

  derivatives = ds.Inverse_Design(
      geom_coordinate, mol_target_list, design_target_property, free_atom_energies)

  derivatives.interpolation(idx_two_mol, type_interp, geom_opt, num_div)
