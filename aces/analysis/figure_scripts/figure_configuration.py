import numpy as np
from astropy import units as u
import matplotlib.colors as mcolors
import pylab as pl
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

distance = 8.1*u.kpc

# mymap_gray is just the lower half of my_colormap
colors_gray = np.vstack((colors1, [colors1[-1]]*128))
mymap_gray = mcolors.LinearSegmentedColormap.from_list('my_colormap_gry', colors_gray)