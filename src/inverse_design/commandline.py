#!/usr/bin/env python
import inverse_design.design as ds
import inverse_design.settings as iconf
# from apdft.settings import Configuration

def ignition():
  geom_coordinate, mol_target_list, free_atom_energies = iconf.Option.get_inputs()

  derivatives = ds.Inverse_Design(geom_coordinate, mol_target_list, free_atom_energies)

  derivatives.design()