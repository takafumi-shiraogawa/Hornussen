#!/usr/bin/env python
import inverse_design.design as ds
import inverse_design.settings as iconf
# from apdft.settings import Configuration

def ignition():
  conf1, conf2 = iconf.Option.get_options()

  derivatives = ds.Inverse_Design(conf1, conf2)

  derivatives.design()