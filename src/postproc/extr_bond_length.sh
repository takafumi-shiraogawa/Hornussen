# Normalized bond lengths
grep "Lime geometry: step" design_opt.dat | \
  awk '{print $NF}' > extr_bond_length.dat
