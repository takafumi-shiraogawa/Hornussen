import numpy as np
import postproc.postprocessing as pp

# How to use?:
# Perform
# extr_atom_ene.sh
# extr_norm_part_coeff.sh
# extr_bond_length.sh
# and
# import postproc.anal_design as pa
# pa.analyze_design.make_figure_opt_hist_atom_ene()
# pa.analyze_design.make_figure_opt_hist_norm_part_coeff(2)
# pa.analyze_design.make_figure_opt_hist_bond_length(2)

class analyze_design():

  def make_figure_opt_hist_atom_ene():

    # After extr_atom_ene.sh

    atom_ene = []
    step = []

    data = open('extr_atom_energy.dat', mode='r')

    count = 0

    for row in data:
      count += 1
      step.append(count)

      row = row.rstrip('\n').split()
      atom_ene.append(float(row[0]))

    data.close()

    step = np.array(step)
    atom_ene = np.array(atom_ene)

    label_data = ['Atomization energy']
    label_x = 'Design step'
    label_y = 'Atomization energy / Hartree'

    range_x = [0, 400, 800, 1200]
    range_y = [0.15, 0.25, 0.35, 0.45]
    lim_x = [-100, 1300]
    lim_y = [0.15, 0.45]

    pp.figure_proc.make_figure(step, atom_ene, label_data, label_x,
                            label_y, range_x, range_y, lim_x, lim_y, pic_name='opt_hist_atom_ene')


  def make_figure_opt_hist_norm_part_coeff(num_target_mol):

    # After extr_norm_part_coeff.sh

    norm_part_coeff = []
    for i in range(num_target_mol):
      norm_part_coeff.append([])
    step = []

    data = open('extr_norm_part_coeff.dat', mode='r')

    count = 0
    real_count = 0

    for row in data:
      count += 1
      if count == num_target_mol:
        real_count += 1
        step.append(real_count - 1)

      row = row.rstrip('\n').split()
      norm_part_coeff[count - 1].append(float(row[0]))

      if count == num_target_mol:
        count = 0

    data.close()

    step = np.array(step)
    norm_part_coeff = np.array(norm_part_coeff)

    label_data = ['$a_{BF}$', '$a_{CO}$']
    label_x = 'Design step'
    label_y = '$a_{i}$'

    range_x = [0, 400, 800, 1200]
    range_y = [-0.25, 0, 0.25, 0.5, 0.75, 1, 1.25]
    lim_x = [-100, 1300]
    lim_y = [-0.25, 1.25]

    pp.figure_proc.make_figure(step, norm_part_coeff, label_data, label_x,
                            label_y, range_x, range_y, lim_x, lim_y, pic_name='opt_hist_norm_part_coeff')


  def make_figure_opt_hist_bond_length(num_atom):

    # After extr_bond_length.sh

    atom_pos = []
    for i in range(num_atom):
      atom_pos.append([])
    step = []

    data = open('extr_bond_length.dat', mode='r')

    count = 0
    real_count = 0

    for row in data:
      count += 1
      if count == num_atom:
        real_count += 1
        step.append(real_count - 1)

      row = row.rstrip('\n').split()
      atom_pos[count - 1].append(float(row[0]))

      if count == num_atom:
        count = 0

    data.close()

    step = np.array(step)
    atom_pos = np.array(atom_pos)

    bond_length = np.zeros(np.shape(atom_pos)[1])
    for i in range(np.shape(atom_pos)[1]):
      bond_length[i] = abs(atom_pos[0, i] - atom_pos[1, i])

    label_data = 'Bond length / Å'
    label_x = 'Design step'
    label_y = 'Bond length / Å'

    range_x = [0, 400, 800, 1200]
    range_y = [1.0, 1.1, 1.2, 1.3, 1.4]
    lim_x = [-100, 1300]
    lim_y = [1.0, 1.4]

    pp.figure_proc.make_figure(step, bond_length, label_data, label_x,
                            label_y, range_x, range_y, lim_x, lim_y, pic_name='opt_hist_bond_length')