import numpy as np
import string

def parse_contdotdat(filepath):

    selections = []

    with open(filepath, 'r') as fh:
        for line in fh:
            if "LSRK" in line:
                selections.append(line.split()[0])


    return ";".join(selections)

try:
    try:
        from __casac__.quanta import quanta
        from taskinit import msmdtool
        from taskinit import mstool
    except (ImportError, ModuleNotFoundError):
        from casatools import quanta
        from casatools import msmetadata as msmdtool
        from casatools import ms as mstool
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
                # the code below doesn't know how to handle the case where
                # no lines are selected; simplest is to simply *not flag anything*
                # in those channels (there are no line channels in that window)
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

    def freq_selection_overlap(ms, freqsel, spw=0):
        """
        For a given frequency selection string (e.g., '215~216GHz;900~950GHz'),
        find the subset of the `;`-separated entries that overlap with the
        measurement set and return those.

        Parameters
        ----------
        ms : str
            The name of a measurement set to compare against
        freqsel : str
            A frequency-based selection string
        spw : int
            The spectral window number, default to 0
        """

        # use the default (topo or lsrk or whatever) frequency, don't cvel it,
        # because the inputs are assumed to be in the same frame.
        # it's important to be exact here, since we're using the fmin/fmax
        # to fill in when one of the ends of a selection is out of range.
        msmd = msmdtool()
        msmd.open(ms)
        frequencies_in_ms = msmd.chanfreqs(spw)
        msmd.close()

        fmin, fmax = frequencies_in_ms.min(), frequencies_in_ms.max()

        new_sel = []

        for selstr in freqsel.split(";"):
            lo, hi = map(float, selstr.strip(string.ascii_letters).split("~"))
            unit = selstr.lstrip(string.punctuation + string.digits)
            if lo > hi:
                lo, hi = hi, lo

            flo = qq.convert({'value':float(lo), 'unit':unit}, 'Hz')['value']
            fhi = qq.convert({'value':float(hi), 'unit':unit}, 'Hz')['value']

            if ((fhi < fmax) and (fhi > fmin)) and ((flo > fmin) and (flo < fmax)):
                # if the whole thing is in range...
                new_sel.append(selstr)
            elif ((fhi < fmax) and (fhi > fmin)) and not ((flo > fmin) and (flo < fmax)):
                # if the lower end is out of range, but the upper is in
                new_selstr = "{0}~{1}GHz".format(fmin / 1e9, fhi / 1e9)
                new_sel.append(new_selstr)
            elif (not ((fhi < fmax) and (fhi > fmin))) and ((flo > fmin) and (flo < fmax)):
                # if the upper end is out of range, but the lower is in
                new_selstr = "{0}~{1}GHz".format(flo / 1e9, fmax / 1e9)
                new_sel.append(new_selstr)

        return "{0}:".format(spw) + ";".join(new_sel)

    def cont_channel_selection_to_contdotdat(cont_channel_selection, msname,
                                             spw_mapping=None):
        """
        Create a cont.dat file from continuum channel selections like 0:10~20, etc.


        Convert the result to a selection string:
            fselstr = ",".join(str(x)+":"+ ";".join(freqsel[x]) for x in freqsel)
        """

        ms.open(msname)

        freqsels = {}

        for spwsel in cont_channel_selection.split(","):
            spwn = int(spwsel.split(":")[0])
            if spw_mapping is not None and spwn in spw_mapping:
                spw = spw_mapping[spwn]
            elif spw_mapping is not None:
                continue
            print("spectral window = {spw}".format(spw=spw))
            freqs = ms.cvelfreqs(spw)

            chansel = spwsel.split(":")[1]
            freqsels[spw] = []
            for chs in chansel.split(";"):
                lo,hi = map(int, chs.split("~"))
                freqsels[spw].append("{0}~{1}GHz".format(freqs[lo]/1e9,
                                                         freqs[hi]/1e9))

        ms.close()

        return freqsels

except ModuleNotFoundError:
    # if using this not in a casa6 environment
    pass


    """
    flagchannels='0:0~60;180~300;2790~2880;3280~3360;3460~3490;3830~3839,1:60~130;200~250;320~420;580~650;1000~1040;1200~1360;1420~1460;1720~1790;1860~1919,2:40~300;630~700;800~1000;1440~1640;1780~1919,3:100~150;470~540;640~820;920~980;1220~1260;1370~1420;1710~1780,4:0~60;180~300;2790~2880;3280~3360;3460~3490;3830~3839,5:60~130;200~250;320~420;580~650;1000~1040;1200~1360;1420~1460;1720~1790;1860~1919,6:40~300;630~700;800~1000;1440~1640;1780~1919,7:100~150;470~540;640~820;920~980;1220~1260;1370~1420;1710~1780,8:0~60;180~300;2790~2880;3280~3360;3460~3490;3830~3839,9:60~130;200~250;320~420;580~650;1000~1040;1200~1360;1420~1460;1720~1790;1860~1919,10:40~300;630~700;800~1000;1440~1640;1780~1919,11:100~150;470~540;640~820;920~980;1220~1260;1370~1420;1710~1780,12:0~60;180~300;900~1050;1860~1950;2100~2140;2230~2280;2790~2880;3050~3100;3280~3360;3460~3490;3590~3650;3830~3839,13:60~130;200~250;265~285;320~420;435~460;580~650;670~700;760~810;1000~1040;1200~1360;1420~1460;1720~1790;1800~1840;1860~1919,14:40~300;630~700;800~1000;1440~1640;1780~1919,15:100~150;470~540;640~820;920~980;1170~1190;1220~1260;1370~1420;1710~1780'
    freqsel = cont_channel_selection_to_contdotdat(flagchannels, ('./science_goal.uid___A001_X1290_X44/group.uid___A001_X1290_X45/member.uid___A001_X1290_X46/calibrated/calibrated_final.ms/'), spw_mapping={0:25,1:27,2:29,3:31})
    fselstr = ";".join( ";".join(freqsel[x]) for x in freqsel)
    ms.open('./science_goal.uid___A001_X1290_X44/group.uid___A001_X1290_X45/member.uid___A001_X1290_X46/calibrated/calibrated_final.ms/')
    freqs = {25: ms.cvelfreqs(25), 27: ms.cvelfreqs(27), 29:ms.cvelfreqs(29), 31:ms.cvelfreqs(31)}
    ms.close()
    contsel = contchannels_to_linechannels(fselstr, freqs)
    contfreqsel = cont_channel_selection_to_contdotdat(contsel, ('./science_goal.uid___A001_X1290_X44/group.uid___A001_X1290_X45/member.uid___A001_X1290_X46/calibrated/calibrated_final.ms/'), spw_mapping={0:25,1:27,2:29,3:31})
    print("Field: Sgr_B2_DS")
    for entry in contfreqsel:
        print()
        print(f"SpectralWindow: {entry}")
        for row in contfreqsel[entry]:
            print(f"{row} LSRK")
    """
