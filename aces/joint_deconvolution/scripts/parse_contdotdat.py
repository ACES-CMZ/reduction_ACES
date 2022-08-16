import numpy as np
import string

def parse_contdotdat(filepath):

    selections = []

    with open(filepath, 'r') as fh:
        for line in fh:
            if "LSRK" in line:
                selections.append(line.split()[0])


    return ";".join(selections)


from casatools import quanta 

qq = quanta()

def contchannels_to_linechannels(contsel, freqslist, return_fractions=False):
    """
    Parameters
    ----------
    contsel : str
        A CASA selection string with assumed units of frequency and no assumed
        spectral windows.
    freqslist : dict
        A dictionary of frequency arrays, where the key is the spectral window
        number and the value is a numpy array of frequencies
    Returns
    -------
    channel_selection : str
        A comma-separated string listing the *channels* corresponding to lines.
        Each section will be labeled by the appropriate SPW.  For example, you
        might get: "0:1~15;30~40,1:5~10,15~20"
    """

    new_sel = []

    line_fraction = {}

    for spw,freq in freqslist.items():
        fmin, fmax = np.min(freq), np.max(freq)
        if fmin > fmax:
            raise ValueError("this is literally impossible")
        selected = np.zeros_like(freq, dtype='bool')
        for selstr in contsel.split(";"):
            lo, hi = selstr.strip(string.ascii_letters).split("~")
            unit = selstr.lstrip(string.punctuation + string.digits)
            flo = qq.convert({'value':float(lo), 'unit':unit}, 'Hz')['value']
            fhi = qq.convert({'value':float(hi), 'unit':unit}, 'Hz')['value']
            if flo > fhi:
                flo,fhi = fhi,flo

            # only include selections that are at least partly in range
            if fmin < fhi < fmax or fmax > flo > fmin:
                selected |= (freq > flo) & (freq < fhi)
            # but also allow for the case where EVERYTHING is included
            elif fhi > fmax and flo < fmin:
                selected[:] = True

        # invert from continuum to line
        invselected = ~selected

        line_fraction[spw] = invselected.sum() / float(invselected.size)
        if line_fraction[spw] == 0:
            continue

        # get the indices where we swap between selected and not
        chans = np.where(invselected[1:] != invselected[:-1])[0].tolist()

        if invselected[0]:
            # if the first index is 'True', then we start with selected
            chans = [0] + chans
        if invselected[-1]:
            chans = chans + [len(freq)-1]

        if len(chans) % 2 > 0:
            raise ValueError("Found an odd number of channel endpoints in "
                             "line inclusion for spw {0}. ".format(spw))

        selchan = ("{0}:".format(spw) +
                   ";".join(["{0}~{1}".format(lo,hi)
                             for lo, hi in zip(chans[::2], chans[1::2])]))

        new_sel.append(selchan)

    if return_fractions:
        return ",".join(new_sel), line_fraction
    else:
        return ",".join(new_sel)

