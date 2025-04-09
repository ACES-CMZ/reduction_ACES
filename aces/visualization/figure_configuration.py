import numpy as np
from astropy import units as u
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab as pl

from aces.visualization.compare_to_other import format_ax

pl.rcParams['figure.facecolor'] = 'w'
pl.rcParams['figure.figsize'] = (10, 8)
pl.rcParams['font.size'] = 16
pl.rcParams['xtick.direction'] = 'in'
pl.rcParams['ytick.direction'] = 'in'

colors1 = pl.cm.gray_r(np.linspace(0., 1, 128))
#colors2 = pl.cm.inferno(np.linspace(0, 1, 128))
colors2 = pl.cm.hot(np.linspace(0, 1, 128))
colors2b = pl.cm.afmhot(np.linspace(0, 1, 128))

colors = np.vstack((colors1, colors2))
mymap = mcolors.LinearSegmentedColormap.from_list('my_colormap', colors)

colors_stack2 = np.vstack((colors1, colors2b))
mymap2 = mcolors.LinearSegmentedColormap.from_list('my_afmcolormap', colors_stack2)

# from Gravity paper?
# https://ui.adsabs.harvard.edu/abs/2019A%26A...625L..10G/abstract says 8.178
# https://ui.adsabs.harvard.edu/abs/2021A%26A...647A..59G/abstract 8.275
distance = 8.275*u.kpc

# mymap_gray is just the lower half of my_colormap
colors_gray = np.vstack((colors1, [colors1[-1]]*128))
mymap_gray = mcolors.LinearSegmentedColormap.from_list('my_colormap_gry', colors_gray)