# Atomization energy
grep "Lime design property:" design_opt.dat | \
  awk '{print $NF}' > extr_design_prop.dat
