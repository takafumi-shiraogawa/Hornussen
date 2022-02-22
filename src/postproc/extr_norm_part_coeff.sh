# Normalized participation coefficients
grep "Lime molecule:" design_opt.dat | \
  awk '{print $NF}' | sed 's/molecule//' > extr_norm_part_coeff.dat
