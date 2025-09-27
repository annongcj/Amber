import numpy as np
from cctbx import crystal, miller
from cctbx.array_family import flex
from mmtbx.max_lik.maxlik import alpha_beta_est_manager

import boost_adaptbx.boost.python as bp

ext = bp.import_ext("cctbx_miller_ext")
from cctbx_miller_ext import binning


class my_alpha_beta_est_manager(alpha_beta_est_manager):
    # modified version which uses test set binner limits
    def my_alpha_beta_for_each_reflection(self, f_obs=None):
        if f_obs is None: f_obs = self.f_obs
        alpha = flex.double(f_obs.size())
        beta = flex.double(f_obs.size())
        # f_obs.setup_binner(n_bins= len(self.alpha_in_zones))

        # only the following 5 lines of code are added
        # retrieve test set binner limits
        limits = self.f_obs_test.binner().limits()
        # cover the whole resolution range
        d_max, d_min = f_obs.resolution_range()
        limits[0] = 1 / (d_max + 1) ** 2  # f_obs.binner().limits()[0]
        limits[-1] = 1 / (max(d_min - 1, 0.01)) ** 2  # f_obs.binner().limits()[-1]
        # apply f_obs_test limits
        f_obs.use_binning(binning=binning(f_obs.unit_cell(), limits))

        binner = f_obs.binner()
        if (self.interpolation == True):
            # no need to do smoothing again
            az = flex.double(self.alpha_in_zones)
            bz = flex.double(self.beta_in_zones)
            alpha = binner.interpolate(az, 0)
            beta = binner.interpolate(bz, 0)
        elif (self.interpolation == False):
            for i_bin, az, bz in zip(binner.range_used(), self.alpha_in_zones,
                                     self.beta_in_zones):
                sel = binner.selection(i_bin)
                alpha.set_selected(sel, az)
                beta.set_selected(sel, bz)
        alpha = miller.array(miller_set=f_obs, data=alpha)
        beta = miller.array(miller_set=f_obs, data=beta)
        return alpha, beta


def main():
    # Initialize unit cell dimensions
    uc_dimensions = [10., 13, 13, 90, 90, 90]
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
    estimator = my_alpha_beta_est_manager(
        f_obs=f_obs,
        f_calc=f_calc,
        free_reflections_per_bin=140,
        flags=flags,
        interpolation=True,  # default Phenix option
        epsilons=f_obs.epsilons().data().as_double()  # epsilons==1 in P1 group
    )

    # Print out individual alpha and beta based on "our" binning into a file
    alpha, beta = estimator.my_alpha_beta_for_each_reflection()

    with open('generated/test_ml_alpha_beta_individual.txt', 'w') as f:
        f.write(str(alpha.size()) + '\n')
        f.write(' '.join(['%20.15e' % a[1] for a in alpha]) + '\n')
        f.write(' '.join(['%20.15e' % a[1] for a in beta]) + '\n')


if __name__ == "__main__":
    main()
