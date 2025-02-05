import numpy as np
from astropy import units as u
from astropy.io import ascii
import string


def exp_to_tex(st):
    if st == 'nan':
        return '-'
    elif 'e' in st:
        pt1, pt2 = st.split('e')
        return "{0}\\ee{{{1:d}}}".format(pt1, int(pt2))
    return st


def format_float(st):
    return exp_to_tex("{0:0.2g}".format(st))


latexdict = ascii.latex.latexdicts['AA']
latexdict['tabletype'] = 'table*'
latexdict['tablealign'] = 'htp'


def rounded(value, error, extra=1):
    """
    Return the value and error both rounded to the error's first digit
    """

    if error == 0:
        return (0, 0)

    if hasattr(value, 'unit'):
        value = value.value

    digit = int(np.ceil(-np.log10(error))) + extra
    assert np.round(error, digit) != 0
    return np.round(value, digit), np.round(error, digit)#, digit


def rounded_arr(value, error, extra=1):
    return np.array(list(map(lambda x, y: rounded(x, y, extra=extra)[0],
                             value, error)))


def round_to_n(x, n):
    if np.isnan(x):
        return np.nan
    elif x == 0:
        return 0
    else:
        return round(x, -int(np.floor(np.log10(np.abs(x)))) + (n - 1))


def strip_trailing_zeros(x, nsf=2):
    if '.' in x:
        y = x.rstrip("0")
        if len(y.rstrip('.')) < nsf or (sum(y.count(x) for x in string.digits[1:]) < nsf):
            # if y = 1, for example
            return y + "0"
        return y.rstrip(".")
    else:
        return x
