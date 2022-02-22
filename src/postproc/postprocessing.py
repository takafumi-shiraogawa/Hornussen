import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np

class figure_proc():

  def make_figure(data_x, data_y, label_data, label_x, label_y, \
    range_x=None, range_y=None, lim_x=None, lim_y=None, pic_name='output'):
    """ Make a picture. """

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    plt.rcParams["xtick.major.size"] = 8
    plt.rcParams["ytick.major.size"] = 8

    fig, ax = plt.subplots(figsize=(6, 4), tight_layout=True)

    if data_y.ndim == 1:
      ax.plot(data_x, data_y, marker="o")
    else:
      for i in range(len(data_y)):
        if i < 2:
          spec_maker = "o"
        if i > 1:
          spec_maker = "v"
        if i > 3:
          spec_maker = "^"
        ax.plot(data_x, data_y[i], label=label_data[i], marker=spec_maker)

    x_range = range_x
    y_range = range_y

    xticklabels = ax.get_xticklabels()
    yticklabels = ax.get_yticklabels()
    xlabel = "a"
    ylabel = "b"

    ax.set_xlabel(label_x, fontsize=18, fontname='Arial')
    ax.set_ylabel(label_y, fontsize=18, fontname='Arial')

    if x_range is not None and y_range is not None:
      ax.set_xticks(x_range)
      ax.set_yticks(y_range)
      ax.set_xticklabels(x_range, fontsize=16, fontname='Arial', position=(0.0, -0.02))
      ax.set_yticklabels(y_range, fontsize=16, fontname='Arial', position=(-0.02, 0.0))

    fmt = '%.2f'
    yticks = mtick.FormatStrFormatter(fmt)
    ax.yaxis.set_major_formatter(yticks)
    if not isinstance(data_x[0], np.int64):
      fmt = '%.1f'
      xticks = mtick.FormatStrFormatter(fmt)
      ax.xaxis.set_major_formatter(xticks)

    fig.show()

    if data_y.ndim != 1:
      ax.legend(loc = "upper right", fontsize=12)

    if lim_x != [] and lim_y != []:
      plt.xlim(lim_x)
      plt.ylim(lim_y)

    plt.savefig("%s.png" % pic_name)