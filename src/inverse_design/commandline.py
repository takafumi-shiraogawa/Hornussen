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
  geom_coordinate, mol_target_list, free_atom_energies = iconf.Option.get_inputs()
  design_target_property = iconf.Option.get_input_design()

  derivatives = ds.Inverse_Design(
      geom_coordinate, mol_target_list, design_target_property, free_atom_energies)

  derivatives.design()