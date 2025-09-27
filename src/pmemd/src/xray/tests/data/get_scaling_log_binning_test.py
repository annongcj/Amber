from cctbx import crystal, miller
from cctbx.array_family import flex


def main():
    # Initialize unit cell dimensions
    uc_dimensions = [10., 13, 17, 90, 90, 90]
    # uc_dimensions = [10., 20, 30, 90, 90, 90]

    # Initialize Miller indices set
    N = 20
    indices = [(i + 1, j + 1, k + 1) for i in range(N) for j in range(N) for k in range(N)]
    # Initialize Miller set object
    ms = miller.set(crystal_symmetry=crystal.symmetry(unit_cell=uc_dimensions, space_group='P 1'),
                    anomalous_flag=False,
                    indices=flex.miller_index(indices))
    # Initialize f_obs array, values do not matter
    f_obs = ms.array(flex.double([i + 1 for i in range(len(indices))]))
    f_obs.set_observation_type_xray_amplitude()

    # Divide resolution range into bins
    selections = f_obs.log_binning()

    # Dump resolution ranges of bins
    with open('generated/test_scaling_log_binning.txt', 'w') as f:
        f.write(str(len(selections)) + '\n')
        for s in selections[::-1]:
            f.write('%20.15e %20.15e\n' % (f_obs.select(s).resolution_range()[::-1]))


if __name__ == "__main__":
    main()
