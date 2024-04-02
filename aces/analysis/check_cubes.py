from aces.analysis.statcont_cubes import check_cube, check_fits_file
import os
import glob


def main():
    dirnames = sorted(glob.glob("/orange/adamginsburg/ACES/data/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_*/calibrated/working/"))

    if os.getenv('SLURM_ARRAY_TASK_ID') is not None:
        slurm_array_task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))

        for ii, dirname in enumerate(dirnames):
            if ii == slurm_array_task_id:
                print(dirname)
                files = glob.glob(f"{dirname}/*.statcont.contsub.fits")

                for fn in files:
                    print(fn)
                    check_fits_file(fn, remove=True)
                    if os.path.exists(fn):
                        check_cube(fn, remove=True)


if __name__ == "__main__":
    main()
