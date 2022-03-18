
# How to use?
# . extr_geom.sh 3

echo "Get step$1 geometry!"

# Write the number of atoms
grep "Lime geometry:" design_opt.dat | grep -v "\----" | grep "step$1" \
  | awk '{print $(NF-2), $(NF-1),$NF}' | wc -l | tr -d ' ' > geom_step$1.xyz

echo "step$1 geometry" >> geom_step$1.xyz

# Write the geometry
grep "Lime geometry:" design_opt.dat | grep -v "\----" | grep "step$1" \
  | awk '{print $(NF-2), $(NF-1),$NF}' >> geom_step$1.xyz
