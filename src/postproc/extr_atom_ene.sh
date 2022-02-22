# Atomization energy
grep "Lime atomization energy:" design_opt.dat | \
  awk '{print $NF}' > extr_atom_energy.dat
