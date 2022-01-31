#!/usr/bin/env python
import configparser

class Option:

  def get_options():
    lime_conf = configparser.ConfigParser()
    lime_conf.read('lime.conf')

    return lime_conf['design']['design_init_mol_geom_path'], lime_conf['design']['design_target_mol_path']