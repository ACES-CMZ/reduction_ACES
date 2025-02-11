import glob
import warnings
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import spectral_cube
import astropy.units as u
from astropy.stats import mad_std

warnings.filterwarnings("ignore", module="spectral_cube")


def set_plot_params():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.weight": "bold",
        "axes.labelweight": "bold",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 1.5,
        "ytick.major.width": 1.5,
        "axes.linewidth": 1.5,
        "font.size": 14,
        "axes.labelsize": 16,
    })


def load_spectral_data(file_path):
    """
    Load spectral data from a FITS file and extract valid spectral measurements.

    Args:
        file_path: Path to the FITS file containing spectral cube data.

    Returns:
        dict: Contains spectrum array, spectral axis, validity mask, and index mappings
             between masked and original arrays.
    """
    cube = spectral_cube.SpectralCube.read(file_path).to(u.K)
    spectrum = cube.mean(axis=(1, 2))
    spectral_axis = cube.spectral_axis.to("MHz")
    spectrum_arr = np.array(spectrum)
    valid_mask = ~np.isnan(spectrum_arr)

    return {
        "spectrum": spectrum_arr[valid_mask],
        "spectral_axis": spectral_axis[valid_mask],
        "mask": valid_mask,
        "original_to_masked": np.cumsum(valid_mask) - 1,
        "masked_to_original": np.arange(len(spectrum_arr))[valid_mask],
    }


def get_masked_ranges(line_free_ranges, mask, original_to_masked):
    """
    Convert spectral ranges to their masked array equivalents.

    Args:
        line_free_ranges: List of tuples containing start and end indices of line-free regions.
        mask: Boolean array indicating valid data points.
        original_to_masked: Array mapping original indices to masked array indices.

    Returns:
        list: Tuples of start and end indices in the masked array space.
    """
    masked_ranges = []
    for start, end in line_free_ranges:
        if start < len(mask) and mask[start]:
            masked_start = original_to_masked[start]
            end_channel = min(end, len(mask) - 1)
            while end_channel > start and not mask[end_channel]:
                end_channel -= 1
            if end_channel > start:
                masked_ranges.append((masked_start, original_to_masked[end_channel]))
    return masked_ranges


def auto_select_line_free_ranges_sigma_clip(spectrum_arr, min_range_length=10, sigma_threshold=2.0, max_iter=10):
    """
    Automatically identify line-free regions using iterative sigma clipping.

    Args:
        spectrum_arr: Array of spectral data.
        min_range_length: Minimum length of a valid line-free range.
        sigma_threshold: Number of standard deviations for clipping.
        max_iter: Maximum number of sigma-clipping iterations.

    Returns:
        list: Tuples of start and end indices for line-free regions.
    """
    mask = np.ones_like(spectrum_arr, dtype=bool)
    previous_mask = np.zeros_like(spectrum_arr, dtype=bool)
    iteration = 0

    while iteration < max_iter and not np.array_equal(mask, previous_mask):
        previous_mask = mask.copy()
        median_val = np.median(spectrum_arr[mask])
        sigma = mad_std(spectrum_arr[mask])
        mask = spectrum_arr < (median_val + sigma_threshold * sigma)
        iteration += 1

    ranges = []
    i = 0
    while i < len(mask):
        if mask[i]:
            start = i
            while i < len(mask) and mask[i]:
                i += 1
            if (i - start) >= min_range_length:
                ranges.append((start, i - 1))
        else:
            i += 1

    return ranges


def prepare_line_free_data(spectrum_arr, spectral_axis, masked_ranges):
    """
    Prepare line-free spectral data for analysis.

    Args:
        spectrum_arr: Array of spectral data.
        spectral_axis: Array of frequency values.
        masked_ranges: List of tuples containing start and end indices of line-free regions.

    Returns:
        dict: Contains concatenated spectrum segments, adjusted spectral axis, and offset value.
    """
    spectrum_segments = [spectrum_arr[start:end + 1] for start, end in masked_ranges]
    spectral_segments = [spectral_axis[start:end + 1] for start, end in masked_ranges]
    offset_value = spectral_axis[masked_ranges[0][0]]

    return {
        "spectrum": np.concatenate(spectrum_segments),
        "spectral_axis": np.concatenate(spectral_segments) - offset_value,
        "offset_value": offset_value,
    }


def fit_single_sin(t, y, f_guess):
    """
    Fit a single sinusoid to data.

    Args:
        t: Time/frequency values.
        y: Amplitude values.
        f_guess: Initial frequency guess for the fit.

    Returns:
        dict: Contains fit parameters, fit function, and covariance matrix.
    """
    A_guess = (np.max(y) - np.min(y)) / 2.0

    def sinfunc(t_val, A, f, phi):
        return A * np.sin(2 * np.pi * f * t_val + phi)

    popt, pcov = scipy.optimize.curve_fit(sinfunc, t, y, p0=[A_guess, f_guess, 0.0], maxfev=10000)
    A, f, phi = popt

    def fitfunc(t_val):
        return A * np.sin(2 * np.pi * f * t_val + phi)

    return {
        "A": A,
        "f": f,
        "phi": phi,
        "fitfunc": fitfunc,
        "popt": popt,
        "pcov": pcov
    }


def iterative_fit_sinusoids(t, y, threshold_frac=0.05, max_iter=20):
    """
    Iteratively fit multiple sinusoids to data using Fourier analysis.

    Args:
        t: Time/frequency values.
        y: Amplitude values.
        threshold_frac: Fractional threshold for peak identification.
        max_iter: Maximum number of fitting iterations.

    Returns:
        dict: Contains offset, individual fits, complete model, residual, and model function.
    """
    t, y = np.array(t), np.array(y)
    dt = t[1] - t[0]
    offset = np.mean(y)
    residual = y - offset
    model_sum = np.zeros_like(y)
    fits = []

    freqs = np.fft.rfftfreq(len(y), d=dt)
    threshold_power = threshold_frac * np.max(np.abs(np.fft.rfft(residual)[1:]) ** 2)

    for _ in range(max_iter):
        FT = np.fft.rfft(residual)
        power = np.abs(FT) ** 2
        if len(power) <= 1:
            break

        idx_peak = np.argmax(power[1:]) + 1
        if power[idx_peak] < threshold_power:
            break

        fit_result = fit_single_sin(t, residual, freqs[idx_peak])
        sine_fit = fit_result["fitfunc"](t)
        model_sum += sine_fit
        residual -= sine_fit
        fits.append(fit_result)

    def model_func(t_input):
        return offset + sum(
            fit["A"] * np.sin(2 * np.pi * fit["f"] * t_input + fit["phi"]) for fit in fits
        )

    return {
        "offset": offset,
        "fits": fits,
        "model": offset + model_sum,
        "residual": residual,
        "model_func": model_func
    }


def create_full_spectrum_plots(full_spectrum, full_spectral_axis, iter_fit, offset_value, output_prefix, masked_ranges):
    """
    Create plots showing original, corrected, and residual spectra.

    Args:
        full_spectrum: Complete spectral data array.
        full_spectral_axis: Complete frequency axis array.
        iter_fit: Dictionary containing fit results.
        offset_value: Frequency offset value.
        output_prefix: Prefix for output filenames.
        masked_ranges: List of tuples containing line-free regions.
    """
    full_axis_corr = full_spectral_axis - offset_value
    baseline_full = iter_fit["model_func"](full_axis_corr.value)
    corrected_full_spectrum = full_spectrum - baseline_full

    fig, (ax_orig, ax_corr, ax_resid) = plt.subplots(1, 3, figsize=(18, 6))
    ymin = min(np.min(full_spectrum), np.min(corrected_full_spectrum))
    ymax = max(np.max(full_spectrum), np.max(corrected_full_spectrum))
    yrange = ymax - ymin
    ylims = (ymin - 0.05 * yrange, ymax + 0.05 * yrange)

    for ax, data, title in [
        (ax_orig, full_spectrum, "Original Spectrum"),
        (ax_corr, corrected_full_spectrum, "Corrected Full Spectrum"),
        (ax_resid, full_spectrum - corrected_full_spectrum, "Residual (Original - Corrected)")
    ]:
        ax.plot(full_axis_corr, data, color="cornflowerblue", linewidth=1.5)
        ax.set_ylim(ylims)
        ax.grid(True, alpha=0.3)
        ax.set_title(title)
        ax.set_xlabel("Frequency Offset [MHz]")
        ax.set_ylabel("Intensity [K]")

    for start, end in masked_ranges:
        ax_orig.axvspan(full_axis_corr[start].value, full_axis_corr[end].value, color="gray", alpha=0.2)

    plt.tight_layout()
    fig.savefig(f"{output_prefix}_full_spectrum_correction_iterative_sigma.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def main():
    """
    Main function to process FITS files and perform baseline corrections.
    Skips specific spectral windows and previously corrected files.
    SPWs 21 and 23 are skipped as this doesn't work very well for them,
    likely due to the combo of broad lines + narrow bandwidth.
    """
    fits_files = glob.glob("*.fits")
    if not fits_files:
        print("No .fits files found in the current directory")
        return

    set_plot_params()
    for file_path in fits_files:
        if any(x in file_path for x in ["spw21", "spw23", "baseline_corrected"]):
            continue

        try:
            data = load_spectral_data(file_path)
            masked_ranges = auto_select_line_free_ranges_sigma_clip(data["spectrum"], min_range_length=20, sigma_threshold=3.0, max_iter=100)
            line_free_data = prepare_line_free_data(data["spectrum"], data["spectral_axis"], masked_ranges)

            iter_fit = iterative_fit_sinusoids(line_free_data["spectral_axis"].value,
                                               line_free_data["spectrum"],
                                               threshold_frac=0.05,
                                               max_iter=20)

            create_full_spectrum_plots(data["spectrum"].copy(),
                                       data["spectral_axis"].copy(),
                                       iter_fit,
                                       line_free_data["offset_value"],
                                       file_path.replace(".fits", ""),
                                       masked_ranges)

            print(f"Successfully processed {file_path}")
        except Exception as err:
            print(f"Error processing {file_path}: {err}")


if __name__ == "__main__":
    main()