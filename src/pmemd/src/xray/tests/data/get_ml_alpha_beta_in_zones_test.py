import numpy as np
from cctbx import crystal, miller
from cctbx.array_family import flex
from mmtbx.max_lik.maxlik import alpha_beta_est_manager


def main():
    # Initialize unit cell dimensions
    uc_dimensions = [10., 13, 17, 90, 90, 90]
    # uc_dimensions = [10., 20, 30, 90, 90, 90]

    # Initialize Miller indices set
    # IMPORTANT: it is not a complete set!
    # This affects the expansion of the values in zones onto the individual reflections!
    N = 20
    indices = [(i + 1, j + 1, k + 1) for i in range(N) for j in range(N) for k in range(N)]
    # Initialize Miller set object
    # Beware that the reflections might have non-unique resolution values due to unit cell dimensions
    # This affects resolution binning
    ms = miller.set(crystal_symmetry=crystal.symmetry(unit_cell=uc_dimensions, space_group='P 1'),
                    anomalous_flag=False,
                    indices=flex.miller_index(indices))
    # indices = ms.indices()

    # Initialize f_obs and f_calc arrays
    f_obs = ms.array(flex.double([i + 1 for i in range(len(indices))]))
    f_obs.set_observation_type_xray_amplitude()

    f_calc = miller.array(
        miller_set=f_obs,
        data=flex.double([abs(sf[1] * (1 + 0.1 * np.sin(i + 1))) for i, sf in enumerate(f_obs)])
    ).set_info(f_obs.info()).set_observation_type(f_obs)

    # Initialize R-free values array with every 10th 'free' value
    flags = flex.bool([(i % 10 == 0) for i in range(len(indices))])
    # flags = f_obs.generate_r_free_flags().data()

    # Estimate alpha and beta of ML-target
    estimator = alpha_beta_est_manager(
        f_obs=f_obs,
        f_calc=f_calc,
        free_reflections_per_bin=140,
        flags=flags,
        interpolation=True,  # default Phenix option
        epsilons=f_obs.epsilons().data().as_double()  # epsilons==1 in P1 group
    )

    # Print out alpha and beta in zones into a file and into stdout
    alpha_in_zones, beta_in_zones = estimator.alpha_in_zones, estimator.beta_in_zones
    with open('generated/test_ml_alpha_beta_in_zones.txt', 'w') as f:
        f.write(str(len(alpha_in_zones)) + '\n')
        f.write(' '.join(['%20.15e' % a for a in alpha_in_zones]) + '\n')
        f.write(' '.join(['%20.15e' % a for a in beta_in_zones]) + '\n')


if __name__ == "__main__":
    main()
