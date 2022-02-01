#!/usr/bin/env python
import inverse_design.design as ds
import inverse_design.settings as iconf
# from apdft.settings import Configuration

def ignition():
  geom_coordinate, mol_target_list = iconf.Option.get_inputs()

  derivatives = ds.Inverse_Design(geom_coordinate, mol_target_list)

  derivatives.design()