import os
import glob
import pandas as pd
from pathlib import Path
from astropy.table import Table


def process_string(input_string):
    return ''.join(char.lower() for char in input_string if char not in ' -()')


def read_table(filename):
    table = Table.read(filename, format='csv')
    return {
        process_string(str(line)): {
            "mol_12m_spw": str(int(spw_12m)),
            "mol_7m_spw": str(int(spw_7m)),
            "mol_TP_spw": str(int(spw_tp)),
            "restfreq": restfreq
        }
        for line, spw_12m, spw_7m, spw_tp, restfreq in zip(
            table['Line'], table['12m SPW'], table['7m SPW'], table['TP SPW'], table['Rest (GHz)']
        )
    }


def perform_mstransform(obs_dir, seven_m_cube, line_spws, line, phasecenter, restfreq, vel_0):
    output_vis = str(obs_dir / (Path(seven_m_cube).stem + '_1chan.ms'))
    common_params = {
        'vis': seven_m_cube,
        'outputvis': output_vis,
        'regridms': True,
        'phasecenter': phasecenter,
        'mode': 'frequency',
        'spw': line_spws[line]['mol_7m_spw'],
        'width': 1,
        'start': vel_0,
        'nchan': 1,
        'restfreq': str(restfreq),
        'outframe': 'LSRK',
        'veltype': 'radio'
    }

    for datacolumn in ['corrected', 'data']:
        try:
            mstransform(datacolumn=datacolumn, **common_params)
            print(f"[INFO] Successfully transformed {seven_m_cube} using datacolumn='{datacolumn}'")
            return True
        except Exception as e:
            print(f"[ERROR] Failed to transform {seven_m_cube} using datacolumn='{datacolumn}': {e}")

    print(f"[ERROR] Aborting transformation for {seven_m_cube} after all attempts failed.")
    return False


def mstrans(aces_workdir, aces_data, aces_rootdir, line_spws, line, phasecenter, vel_0):
    sb_names = pd.read_csv(aces_rootdir / 'aces/data/tables/aces_SB_uids.csv')

    visfiles = []
    for seven_m_mous_id in sb_names['7m MOUS ID']:
        obs_dir = aces_workdir / 'seven_m_single_chans' / f'member.uid___A001_{seven_m_mous_id}'
        obs_dir.mkdir(parents=True, exist_ok=True)

        ms_files = glob.glob(str(aces_data / f'member.uid___A001_{seven_m_mous_id}/calibrated/working/*.ms'))
        for ms_file in ms_files:
            visfiles.append((ms_file, obs_dir))

    successful_transformations = 0
    total_transformations = len(visfiles)

    for seven_m_cube, obs_dir in visfiles:
        output_ms = obs_dir / (Path(seven_m_cube).stem + '_1chan.ms')
        print(f'Transforming: {Path(seven_m_cube).name} to {output_ms}')
        if perform_mstransform(obs_dir, seven_m_cube, line_spws, line, phasecenter, line_spws[line]['restfreq'], vel_0):
            successful_transformations += 1

    print(f"[INFO] Completed {successful_transformations} out of {total_transformations} transformations successfully.")


def main():
    aces_rootdir = Path(os.getenv('ACES_ROOTDIR'))
    aces_workdir = Path(os.getenv('ACES_WORKDIR'))
    aces_data = Path(os.getenv('ACES_DATA'))

    phasecenter = 'ICRS 17h46m03.5s -28d50m02.6s'
    vel_0 = '87.9106 GHz'

    line_table = aces_rootdir / 'aces/data/tables/linelist.csv'
    lines = ['hnco43']

    line_spws = read_table(line_table)
    for line in lines:
        mstrans(aces_workdir, aces_data, aces_rootdir, line_spws, line, phasecenter, vel_0)


if __name__ == "__main__":
    main()