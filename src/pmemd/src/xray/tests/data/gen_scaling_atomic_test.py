import numpy as np
from cctbx import crystal, miller
from cctbx.array_family import flex

from scaling_common import SFmanager, ScalingManager, \
    dump_sfs_and_k_arrays, f_obs_from_f_calc_f_mask_1, f_obs_from_f_calc_f_mask_2


def main():
    scatterers_coordinates = [(0.0, 2, 4), (5, 3, 4), (6, 10, 5), (1, 1, 2)]
    scatterers_types = [2, 3, 2, 5]
    # Initialize unit cell dimensions
    uc_dimensions = [10., 13, 17, 90, 90, 90]
    # uc_dimensions = [10., 20, 30, 90, 90, 90]

    # Initialize Miller indices set
    # IMPORTANT: it is not a complete set!
    # This affects the expansion of the values in zones onto the individual reflections!
    n_h = 7
    n_k = 9
    n_l = 11
    indices = [(i + 1, j + 1, k + 1) for i in range(n_h) for j in range(n_k) for k in range(n_l)]
    # Initialize Miller set object
    # Beware that the reflections might have non-unique resolution values due to unit cell dimensions
    # This affects resolution binning
    ms = miller.set(crystal_symmetry=crystal.symmetry(unit_cell=uc_dimensions, space_group='P 1'),
                    anomalous_flag=False,
                    indices=flex.miller_index(indices))
    # indices = ms.indices()

    # Initialize R-free values array with every 10th 'free' value
    flags = flex.bool([(i % 10 == 0) for i in range(len(indices))])

    # Prepare test cases for extensive coverage
    sf_set_cases = []

    # CASE 1: Initialize fake f_obs, f_calc and f_mask arrays
    f_obs = ms.array(flex.double([i + 1 for i in range(len(indices))]))
    f_obs.set_observation_type_xray_amplitude()
    f_calc = miller.array(
        miller_set=f_obs,
        data=flex.complex_double([(sf[1] + 2) * (np.exp(1j * np.sin(i + 1))) for i, sf in enumerate(f_obs)])
    ).set_info(f_obs.info()).set_observation_type(f_obs)
    f_mask = [miller.array(
        miller_set=f_obs,
        data=flex.complex_double([0.8 * (sf[1] + 2) * (np.exp(1j * np.cos(i + 2))) for i, sf in enumerate(f_obs)])
    ).set_info(f_obs.info()).set_observation_type(f_obs)]  # 0.8 is needed to avoid cancelling out with f_calc
    sf_set_cases.append((f_obs, f_calc, f_mask, flags))

    # CASES 2 & 3: A bit more realistic f_obs structure factors
    # sf_set_cases.append((f_obs_from_f_calc_f_mask_1(f_calc, f_mask[0]), f_calc, f_mask, flags))
    # sf_set_cases.append((f_obs_from_f_calc_f_mask_2(f_calc, f_mask[0]), f_calc, f_mask, flags))

    # CASES 4, 5 & 6: Realistic structure factors f_calc, f_mask and fake + "realistic" f_obs
    structure_factors = SFmanager(ms, scatterers_types, scatterers_coordinates)
    f_calc, f_mask = structure_factors.compute_SFs()
    sf_set_cases.append((f_obs, f_calc, f_mask, flags))
    sf_set_cases.append((f_obs_from_f_calc_f_mask_1(f_calc, f_mask[0]), f_calc, f_mask, flags))
    sf_set_cases.append((f_obs_from_f_calc_f_mask_2(f_calc, f_mask[0]), f_calc, f_mask, flags))

    # Initialize random k-arrays
    n_bins = len(f_obs.log_binning())
    n_of_k_in_tests = 5
    k1_in = [np.ones(len(indices))]
    k2_in_bins_in = [np.ones(n_bins)]
    k3_in = [np.zeros(len(indices))]
    k4_in = [np.ones(len(indices))]
    for i in range(n_of_k_in_tests - 1):
        np.random.seed(i * 4)
        k1_in.append(np.random.random_sample(len(indices)))
        np.random.seed(i * 4 + 1)
        k2_in_bins_in.append(np.random.random_sample(n_bins))
        np.random.seed(i * 4 + 2)
        k3_in.append(np.random.random_sample(len(indices)))
        np.random.seed(i * 4 + 3)
        k4_in.append(np.random.random_sample(len(indices)))

    # scaling_manager = ScalingManager(f_obs, f_calc, f_mask, ms.array(flags),
    #                                  k1_in[1],
    #                                  k2_in[1],
    #                                  k3_in[1],
    #                                  k4_in[1],
    #                                  )
    # d_f = scaling_manager.f_obs.d_spacings().select(flags).sort().data()
    # d_w = scaling_manager.f_obs.d_spacings().select(~flags).sort().data()
    # for i, s in enumerate(-1.0/flex.pow2(d_f)/4.0):
    #     if i < 5 or i > len(d_f) - 6:
    #         print(i+1, s)
    # for i, s in enumerate(-1.0/flex.pow2(d_w)/4.0):
    #     if i < 5 or i > len(d_w) - 6:
    #         print(i+1, s)

    # Run scaling tests
    for i, case in enumerate(sf_set_cases):
        f_obs, f_calc, f_mask, flags = case

        for func_name in [
            'update_k_iso_exp',
            'update_k_bulk_k_iso',
            'update_k_bulk_k_iso_via_cubic_eq',
            'update_k_aniso',
            'full'
        ]:
            ks_in = []
            ks_out = []

            filename = 'generated/test_scaling_%s_%02d.txt' % (func_name, i + 1)
            for test_number in range(n_of_k_in_tests):
                print("Test #", test_number + 1, "(file:", filename, ")")

                scaling_manager = ScalingManager(f_obs, f_calc, f_mask, ms.array(flags),
                                                 k1_in[test_number],
                                                 k2_in_bins_in[test_number],
                                                 k3_in[test_number],
                                                 k4_in[test_number],
                                                 )
                ks_in.append(scaling_manager.extract_scaling_coefficients())
                ks_out.append(getattr(scaling_manager, func_name)())
                print()


            dump_sfs_and_k_arrays(filename,
                                  f_obs, f_calc, f_mask, flags,
                                  [(ki, kj) for ki, kj in zip(ks_in, ks_out)])


if __name__ == "__main__":
    main()
