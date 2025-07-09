import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import seaborn as sns

# Function to load FITS files
def load_fits(fits_file):
    # Open the FITS file and return the data
    with fits.open(fits_file) as hdul:
        dat = hdul[0].data
        datmax = np.nanmax(dat)
        dat_norm = dat / datmax
        dat_norm[np.isnan(dat_norm)] = 0.0
        return dat_norm

cubepath = '..'


######## Read the continum
contfilename = '../../continuum/12m_continuum_commonbeam_circular_reimaged_mosaic_MUSTANGfeathered.fits'

######## Read the mom0
molnamelist = ['HCOP_mopra','HNCO_7m12mTP','SiO21', 'SO21', 'H13COp', 'H13CN', 'HN13C', 'HC15N','CS21','H40a','SO32','HC3N','CH3CHO','NSplus']

mom0filename = [f'{cubepath}/{molname}_CubeMosaic_masked_hlsig_dilated_mom0.fits' for molname in molnamelist]

#fits_files = np.append(mom0filename, contfilename)
fits_files = [contfilename] + mom0filename

maps = [load_fits(file) for file in fits_files]

# Function to flatten 2D arrays (maps) to 1D for correlation calculation
def flatten_map(map_data):
    return map_data.flatten()

# Compute the 2D correlation coefficients between pairs of maps
n_maps = len(maps)
correlation_matrix = np.zeros((n_maps, n_maps))

# Calculate correlation coefficients between maps
for i in range(n_maps):
    for j in range(i, n_maps):  # Use symmetry: Correlation matrix is symmetric
        # Flatten the 2D arrays and compute the correlation coefficient
        flattened_i = flatten_map(maps[i])
        flattened_j = flatten_map(maps[j])
        corr = np.corrcoef(flattened_i, flattened_j)[0, 1]
        correlation_matrix[i, j] = corr
        correlation_matrix[j, i] = corr  # Symmetric matrix

#molstringlist = [r'SiO$^{\phantom{x}}$', r'SO$^{\phantom{x}}$', r'H$^{13}$CO$^+$', r'H$^{13}$CN', r'HN$^{13}$C', r'HC$^{15}$N']
molstringlist = [r'HCO$^+$',r'HNCO',r'SiO', r'SO (2$_2$-1$_1$)', r'H$^{13}$CO$^+$', r'H$^{13}$CN', r'HN$^{13}$C', r'HC$^{15}$N',r'CS',r'H40$\alpha$',r'SO (2$_3$-1$_2$)', r'HC$_3$N', r'CH$_3$CHO', r'NS$^+$']
tracer_list = ['Continuum'] + molstringlist

# Plot the correlation matrix
plt.figure(figsize=(10, 9))

#mask = np.triu(np.ones_like(correlation_matrix, dtype=bool))

half_corr = correlation_matrix[1:,:-1]
mask = np.triu(np.ones_like(half_corr, dtype=bool), k=1)

#sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", xticklabels=tracer_list, yticklabels=tracer_list,vmax=1, vmin=-1, mask=mask)
sns.heatmap(half_corr, annot=True, cmap="coolwarm", xticklabels=tracer_list[:-1], yticklabels=tracer_list[1:],vmax=1, vmin=-1, mask=mask)

plt.savefig('2Dcorr_matrix_more.pdf', dpi=100, bbox_inches='tight')
