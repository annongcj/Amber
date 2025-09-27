from cctbx import miller
from cctbx import xray
from cctbx.array_family import flex
from cctbx.xray.structure_factors import from_scatterers
from mmtbx import masks

import numpy as np

import cctbx_scaler as scaler


def try_and_return_scale(self,
              k_isotropic_exp=None,
              k_isotropic=None,
              k_mask=None,
              k_anisotropic=None,
              selection=None):
    import mmtbx.arrays
    from mmtbx import bulk_solvent

    if(k_isotropic_exp is None): k_isotropic_exp = self.core.k_isotropic_exp
    if(k_isotropic is None):     k_isotropic     = self.core.k_isotropic
    if(k_mask is None):          k_mask          = self.core.k_mask()
    if(k_anisotropic is None):   k_anisotropic   = self.core.k_anisotropic
    c = mmtbx.arrays.init(
        f_calc          = self.core.f_calc,
        f_masks         = self.core.f_mask(),
        k_isotropic_exp = k_isotropic_exp,
        k_isotropic     = k_isotropic,
        k_anisotropic   = k_anisotropic,
        k_masks         = k_mask)
    sel = self.selection_work.data()
    if(selection is not None): sel = selection & sel
    scale = bulk_solvent.scale(self.f_obs.data(), c.f_model.data(), sel)
    return bulk_solvent.r_factor(self.f_obs.data(), c.f_model.data(), sel, scale), scale


def bulk_solvent_scaling(self, r_start):
    from mmtbx import bulk_solvent

    k_mask = flex.double(self.f_obs.size(), -1)
    k_mask_bin = flex.double()
    k_total = self.core.k_anisotropic * self.core.k_isotropic_exp

    def get_k_mask_trial_range(x, shift=0.05, grid_size=11):
        assert grid_size % 2 == 1, "Grid size must be odd in order to include `shift` exactly"
        x = min(x, 1.0)
        return flex.double([z for z in np.linspace(x-shift, x+shift, grid_size) if 0 <= z <= 1])

    for i_cas, cas in enumerate(self.cores_and_selections):
        selection, core, selection_use, sel_work = cas
        f_obs = self.f_obs.select(selection).data()
        f_calc = core.f_calc.data() * k_total.select(selection)
        f_mask = core.f_mask().data() * k_total.select(selection)

        r = flex.double()
        k = flex.double()
        #
        if (self.bss_result.k_mask_bin_orig is not None):
            x0 = self.bss_result.k_mask_bin_orig[i_cas]
            k_mask.set_selected(selection, x0)
            r0 = self.try_scale(k_mask=k_mask, selection=selection)
            r.append(r0)
            k.append(x0)
        #
        fmv = flex.min(flex.abs(f_mask).select(selection_use))
        if (abs(fmv) > 1.e-9):
            obj1 = bulk_solvent.overall_and_bulk_solvent_scale_coefficients_analytical(
                f_obs=f_obs,
                f_calc=f_calc,
                f_mask=f_mask,
                selection=selection_use)
            k_mask.set_selected(selection, obj1.x_best)
            k.append(obj1.x_best)
        else:
            k_mask.set_selected(selection, 0)
            k.append(0)
        r.append(self.try_scale(k_mask=k_mask, selection=selection))
        #
        s = flex.sort_permutation(r)
        x = k.select(s)[0]
        # fine-sample k_mask around minimum of LS to fall into minimum of R
        k_mask_bin_, k_isotropic_bin_ = \
            bulk_solvent.k_mask_and_k_overall_grid_search(
                f_obs,
                f_calc,
                f_mask,
                get_k_mask_trial_range(x=x),
                selection_use)
        k_mask_bin.append(k_mask_bin_)
        k_mask.set_selected(selection, k_mask_bin_)

        # k_mask_bin.append(x)
        # k_mask.set_selected(selection, x)\
    #
    k_mask_bin_smooth = self.smooth(k_mask_bin)
    k_mask = populate_bin_to_individual_k_mask_linear_interpolation(self,
                                                                    k_mask_bin=k_mask_bin_smooth)
    k_isotropic = self._k_isotropic_as_scale_k1(r_start=r_start, k_mask=k_mask)

    r_try, scale_try = try_and_return_scale(self, k_mask=k_mask, k_isotropic=k_isotropic)
    if (r_try < r_start):
        self.core = self.core.update(k_isotropic=k_isotropic*scale_try, k_masks=k_mask)
        # self.bss_result.k_mask_bin_orig = k_mask_bin
        # self.bss_result.k_mask_bin_smooth = k_mask_bin_smooth
        # self.bss_result.k_mask = k_mask
        self.bss_result.k_isotropic = k_isotropic*scale_try
    r = self.r_factor()
    return r


def populate_bin_to_individual_k_mask_linear_interpolation(self, k_mask_bin):
    from mmtbx import bulk_solvent

    assert len(k_mask_bin) == len(self.cores_and_selections)

    def linear_interpolation(x1, x2, y1, y2):
        k = 0
        if (x1 != x2): k = (y2 - y1) / (x2 - x1)
        b = y1 - k * x1
        return k, b

    result1 = flex.double(self.f_obs.size(), -1)
    result2 = flex.double(self.f_obs.size(), -1)
    result = flex.double(self.f_obs.size(), -1)
    for i, cas in enumerate(self.cores_and_selections):
        selection, zzz, zzz, selection_work = cas  # fix for WORK only reflexes
        ss = self.ss.select(selection_work)        # fix for WORK only reflexes
        x1, x2 = flex.min(ss), flex.max(ss)        # fix for WORK only reflexes
        y1 = k_mask_bin[i]
        if (i == len(k_mask_bin) - 1):
            y2 = k_mask_bin[i]  # fix last bin interpolation
        else:
            y2 = k_mask_bin[i + 1]
        k, b = linear_interpolation(x1=x1, x2=x2, y1=y1, y2=y2)
        bulk_solvent.set_to_linear_interpolated(self.ss, k, b, selection, result1)
        result2.set_selected(selection, y1)
        r1 = self.try_scale(k_mask = result1, selection=selection) # XXX inefficient
        r2 = self.try_scale(k_mask = result2, selection=selection) # XXX inefficient
        if (r1 < r2):
            bulk_solvent.set_to_linear_interpolated(self.ss, k, b, selection, result)
        else:
            result.set_selected(selection, y1)
    assert (result < 0).count(True) == 0
    return result


def k_mask_grid_search(self, r_start):
    from mmtbx import bulk_solvent

    # k_mask_trial_range = flex.double([i/1000. for i in range(0,650,50)])
    # k_mask_trial_range = flex.double([i / 1000. for i in range(0, 1010, 10)])
    k_mask_trial_range = flex.double(np.linspace(0, 1, 101))
    k_mask = flex.double(self.f_obs.size(), 0)
    k_mask_bin = flex.double()
    # k_iso_bin = flex.double()
    k_isotropic = flex.double(self.f_obs.size(), 0)
    k_total = self.core.k_anisotropic * self.core.k_isotropic_exp
    for i_cas, cas in enumerate(self.cores_and_selections):
        selection, core, selection_use, sel_work = cas
        f_obs = self.f_obs.select(selection)
        k_total_ = k_total.select(selection)

        f_calc = core.f_calc.data() * k_total_
        f_mask = core.f_mask().data() * k_total_

        k_mask_bin_, k_isotropic_bin_ = \
            bulk_solvent.k_mask_and_k_overall_grid_search(
                f_obs.data(),
                f_calc,
                f_mask,
                k_mask_trial_range,
                selection_use)

        k_mask_bin.append(k_mask_bin_)
        # k_iso_bin.append(k_isotropic_bin_)
        k_mask.set_selected(selection, k_mask_bin_)
        k_isotropic.set_selected(selection, k_isotropic_bin_)

    k_mask_bin_smooth = self.smooth(k_mask_bin)
    k_mask = populate_bin_to_individual_k_mask_linear_interpolation(self,
        k_mask_bin=k_mask_bin_smooth)
    # k_isotropic = self._k_isotropic_as_scale_k1(r_start=r_start, k_mask=k_mask)

    # update unconditionally to init cycles 2+ of the 'full' loop
    self.bss_result.k_mask_bin_orig = k_mask_bin

    r_try, scale_try = try_and_return_scale(self, k_mask=k_mask, k_isotropic=k_isotropic)
    if (r_try < r_start - 1e-10):  # the tiny shift is for tests numerical stability
        self.core = self.core.update(k_masks=k_mask, k_isotropic=k_isotropic*scale_try)
        # ???? in fact, only k_mask_bin_orig is used
        self.bss_result.k_mask_bin_orig = k_mask_bin
        self.bss_result.k_mask_bin_smooth = k_mask_bin_smooth
        self.bss_result.k_mask = k_mask
        self.bss_result.k_isotropic = k_isotropic*scale_try
    r = self.r_factor()
    if (self.verbose):
        print("      r_final: %6.4f (r_low: %6.4f)" % (r, self._r_low()))
    return r


def set_k_isotropic_exp(self, r_start, b_lower_limit=-100):
    import mmtbx.arrays
    import scitbx.math
    from mmtbx import bulk_solvent

    k_iso   = flex.double(self.core.k_isotropic.size(), 1) # Done at start only!
    k_aniso = flex.double(self.core.k_isotropic.size(), 1) # Done at start only!
    arrays = mmtbx.arrays.init(
        f_calc=self.core.f_calc,
        f_masks=self.core.f_mask(),
        k_isotropic=k_iso,
        k_anisotropic=k_aniso,
        k_masks=self.core.k_mask())
    sel = self.selection_work.data()
    #
    # At least in one example this gives more accurate answer but higher R than start!
    #
    rf = scitbx.math.gaussian_fit_1d_analytical(
        x=flex.sqrt(self.ss).select(sel),
        y=self.f_obs.data().select(sel),
        z=abs(arrays.f_model).data().select(sel))
    if (rf.b < b_lower_limit): return r_start

    k1 = rf.a * flex.exp(-self.ss * rf.b)
    r1 = self.try_scale(k_isotropic_exp=k1, k_anisotropic=k_aniso)

    #
    # At least in one example this gives less accurate answer but lower R than start!
    #
    o = bulk_solvent.f_kb_scaled(
        f1 = self.f_obs.data().select(sel),
        f2 = flex.abs(arrays.f_model.data()).select(sel),
        b_range = flex.double(range(-100,100,1)),
        ss = self.ss.select(sel))
    k2 = o.k() * flex.exp(-self.ss * o.b())
    r2 = self.try_scale(k_isotropic_exp = k2, k_anisotropic=k_aniso)

    if r2 + 1e-5 < r1:  # the tiny shift is for tests numerical stability
        r1 = r2
        k1 = k2

    if r1 < r_start + 1e-5:  # the tiny shift is for tests numerical stability
        self.core = self.core.update(k_isotropic_exp = k1)
    # Note: since k_iso_exp is fitted with k_aniso=1, generally one needs to update the r-factor, accordingly.
    # In the case of normal scaling, this assumption is always true.

    return self.r_factor()


# This is the modified anisotropic_scaling function without forced "force_to_use_expmin" fit,
# if polynomial fit ends up with better r-factor, but negative coefficients
def anisotropic_scaling(self, r_start, use_highres):
    from mmtbx import bulk_solvent
    from cctbx import adptbx
    import boost_adaptbx.boost.python as bp
    ext = bp.import_ext("mmtbx_f_model_ext")

    r_expanal, r_poly, r_expmin = None, None, None
    k_anisotropic_expanal, k_anisotropic_poly, k_anisotropic_expmin = None, None, None
    scale_matrix_expanal, scale_matrix_poly, scale_matrix_expmin = None, None, None
    sel         = self.selection_work.data()

    if(use_highres):
        sel_ = self.f_obs.d_spacings().data() < self.d_hilo
        sel = sel & sel_

    f_model_abs = flex.abs(self.core.f_model_no_aniso_scale.data().select(sel))
    f_obs       = self.f_obs.data().select(sel)
    mi          = self.f_obs.indices().select(sel)
    uc          = self.f_obs.unit_cell()
    mi_all      = self.f_obs.indices()
    # try exp_anal
    if(self.try_expanal):
        obj = bulk_solvent.aniso_u_scaler(
            f_model_abs    = f_model_abs,
            f_obs          = f_obs,
            miller_indices = mi,
            adp_constraint_matrix = self.adp_constraints.gradient_sum_matrix())
        u_star = self.adp_constraints.all_params(tuple(obj.u_star_independent))
        scale_matrix_expanal = adptbx.u_as_b(adptbx.u_star_as_u_cart(uc, u_star))
        k_anisotropic_expanal = ext.k_anisotropic(mi_all, u_star)
        r_expanal = self.try_scale(k_anisotropic = k_anisotropic_expanal)
        if(self.verbose):
            print("      r_expanal: %8.6f"%r_expanal)
    # try poly
    if(self.try_poly):
        obj = bulk_solvent.aniso_u_scaler(
            f_model_abs    = f_model_abs,
            f_obs          = f_obs,
            miller_indices = mi,
            unit_cell      = uc)
        scale_matrix_poly = obj.a
        k_anisotropic_poly = ext.k_anisotropic(mi_all, obj.a, uc)
        r_poly = self.try_scale(k_anisotropic = k_anisotropic_poly) # + 1e-6 ##  for tests only
        if(self.verbose):
            print("      r_poly   : %8.6f"%r_poly)

    # EDIT: instead of force_to_use_expmin optimization, use k_anisotropic_poly only if all coefficients are positive
    if(k_anisotropic_poly is not None and (k_anisotropic_poly<=0).count(True)>0):
        k_anisotropic_poly = None
        r_poly = None
        if(self.verbose):
            print("      r_poly rejected since there is a negative coefficient")
    # END of EDIT

    # select best
    r = [(r_expanal, k_anisotropic_expanal, scale_matrix_expanal),
         (r_poly,    k_anisotropic_poly,    scale_matrix_poly),
         (r_expmin,  k_anisotropic_expmin,  scale_matrix_expmin)]
    if self.verbose:
        print(r_expanal, r_poly, r_start)
    r_best = r_start
    k_anisotropic_best = None
    scale_matrix_best = None
    for result in r:
        r_factor, k_anisotropic, scale_matrix = result
        if(r_factor is not None and r_factor < r_best):
            r_best = r_factor
            k_anisotropic_best = k_anisotropic.deep_copy()
            scale_matrix_best = scale_matrix[:]
    if(scale_matrix_best is None):
        if(self.verbose):
            print("      result rejected due to r-factor increase")
    else:
        self.scale_matrices = scale_matrix_best
        self.core = self.core.update(k_anisotropic = k_anisotropic_best)
        r_aniso = self.r_factor()
        if(self.verbose):
            self.format_scale_matrix()
            print("      r_final  : %8.6f"%r_aniso)
    return r_best


class SFmanager:
    def __init__(self, ms_set, scatterers_types, scatterers_cart_sites):
        assert len(scatterers_types) == len(scatterers_cart_sites)
        self.ms = ms_set
        self.dummy_array = self.ms.array(flex.double(ms_set.size(), 1))

        self.x = xray.structure(crystal_symmetry=self.ms.crystal_symmetry())
        for s in scatterers_types:
            if s == 2:
                s_c = 'C'
            elif s == 3:
                s_c = 'N'
            elif s == 4:
                s_c = 'O'
            elif s == 5:
                s_c = 'S'
            else:
                s_c = 'H'
            self.x.add_scatterer(xray.scatterer(label=s_c))
        self.x.set_sites_cart(flex.vec3_double(scatterers_cart_sites))

        self.manager = masks.manager(self.dummy_array, self.x)

    def compute_SFs(self, new_sites=None):
        if new_sites:
            self.x.set_sites_cart(flex.vec3_double(new_sites))
        return from_scatterers(miller_set=self.ms)(self.x, self.ms, 'direct').f_calc(), \
               self.manager.shell_f_masks(self.x, force_update=True)


class ScalingManager:
    def __init__(self, f_obs, f_calc, f_mask, r_free_flags, k_iso_exp, k_iso_in_bin, k_masks, k_aniso):
        self.f_obs = f_obs
        self.result = scaler.run(
            f_obs=self.f_obs,
            f_calc=f_calc,
            f_mask=f_mask,
            r_free_flags=r_free_flags,
            ss=1. / flex.pow2(self.f_obs.d_spacings().data()) / 4.,  # division by 4.0 is by design of cctbx
            # number_of_cycles=100,  # as in scaler.run class
            ### We can try varying 'try_*' options
            verbose=False)

        expanded_k_iso = flex.double(self.f_obs.size(), -1)
        for i, cas in enumerate(self.result.cores_and_selections):
            selection, _, _, _ = cas
            expanded_k_iso.set_selected(selection, k_iso_in_bin[i])

        self.result.core.update(k_isotropic_exp=flex.double(k_iso_exp),
                                k_isotropic=flex.double(expanded_k_iso),
                                k_masks=flex.double(k_masks),
                                k_anisotropic=flex.double(k_aniso))

    def extract_scaling_coefficients(self):
        return [self.result.core.k_isotropic_exp,
                self.result.core.k_isotropic,
                self.result.core.k_anisotropic,
                self.result.core.k_masks[0]]

    def update_k_aniso(self):
        r0 = self.result.r_factor()
        r1 = anisotropic_scaling(self.result, r0, use_highres=False)
        print("%8.6f" % r0)
        print("%8.6f" % r1)
        return self.extract_scaling_coefficients()

    def update_k_bulk_k_iso_via_cubic_eq(self):
        r0 = self.result.r_factor()
        self.result.bss_result.k_mask_bin_orig = None
        r1 = bulk_solvent_scaling(self.result, r0)
        print("%8.6f" % r0)
        print("%8.6f" % r1)
        return self.extract_scaling_coefficients()

    def update_k_iso_exp(self):
        r0 = self.result.r_factor()
        r1 = set_k_isotropic_exp(self.result, r0)
        print("%8.6f" % r0)
        print("%8.6f" % r1)
        return self.extract_scaling_coefficients()

    def update_k_bulk_k_iso(self):
        r0 = self.result.r_factor()
        r1 = k_mask_grid_search(self.result, r0)
        print("%8.6f" % r0)
        print("%8.6f" % r1)
        return self.extract_scaling_coefficients()

    def full(self):
        self.result.bss_result.k_mask_bin_orig = None
        number_of_cycles = 20
        verbose = self.result.verbose
        for cycle in range(number_of_cycles):
            r_start = self.result.r_factor()
            r_start0 = r_start
            if verbose:
                print("  cycle %d:" % cycle)
                print("    r(start): %6.4f" % r_start)
                # bulk-solvent and overall isotropic scale
            if self.result.bulk_solvent:
                if cycle == 0:
                    print("%8.6f" % r_start0)
                    for mic in [1, 2]:
                        r_start = set_k_isotropic_exp(self.result, r_start=r_start)
                        r_start = k_mask_grid_search(self.result, r_start=r_start)
                        r_start = set_k_isotropic_exp(self.result, r_start=r_start)
                else:
                    r_start = bulk_solvent_scaling(self.result, r_start=r_start)
                    if verbose:
                        print("    r(bulk_solvent_scaling): %6.4f" % r_start)
                # anisotropic scale
            r_start = anisotropic_scaling(self.result, r_start=r_start, use_highres=False)
            if r_start > r_start0 - self.result.auto_convergence_tolerance:
                break
        print("%8.6f" % r_start)
        return self.extract_scaling_coefficients()


def f_obs_from_f_calc_f_mask_1(f_calc, f_mask):
    return miller.array(
        miller_set=f_calc,
        data=flex.double([abs(sf_p + sf_m * 0.35 * np.exp(-ss / 4 * 46.0))
                          for sf_p, sf_m, ss in zip(f_calc.data(),
                                                    f_mask.data(),
                                                    1. / flex.pow2(f_calc.d_spacings().data()))])
    ).set_info(f_calc.info()).set_observation_type(f_calc)


def f_obs_from_f_calc_f_mask_2(f_calc, f_mask):
    return miller.array(
        miller_set=f_calc,
        data=flex.double([abs(sf_p) + abs(sf_m) * 0.35 * np.exp(-ss / 4 * 46.0)
                          for sf_p, sf_m, ss in zip(f_calc.data(),
                                                    f_mask.data(),
                                                    1. / flex.pow2(f_calc.d_spacings().data()))])
    ).set_info(f_calc.info()).set_observation_type(f_calc)


def dump_sfs_and_k_arrays(filename,
                          f_obs, f_calc, f_mask, flags,
                          ks_in_ks_out
                          ):
    order = ['k_iso_exp', 'k_iso', 'k_anisotropic', 'k_mask']

    with open(filename, 'w') as f:
        f.write('# number of reflections\n')
        f.write('%d\n' % f_obs.data().size())

        f.write('# h  array\n')
        for h, _, _ in f_obs.indices():
            f.write('%d ' % h)
        f.write('\n')
        f.write('# k  array\n')
        for _, k, _ in f_obs.indices():
            f.write('%d ' % k)
        f.write('\n')
        f.write('# l  array\n')
        for _, _, l in f_obs.indices():
            f.write('%d ' % l)
        f.write('\n')

        f.write('# f_obs array\n')
        for sf in f_obs.data():
            f.write('%20.15e ' % sf)
        f.write('\n')

        f.write('# f_calc array\n')
        for sf in f_calc.data():
            f.write('(%20.15e,%20.15e) ' % (sf.real, sf.imag))
        f.write('\n')

        f.write('# f_mask array\n')
        for sf in f_mask[0].data():
            f.write('(%20.15e,%20.15e) ' % (sf.real, sf.imag))
        f.write('\n')

        f.write('# r_flags array\n')
        for flag in flags:
            f.write('T ' if flag else 'F ')
        f.write('\n')

        f.write('# number of scaling tests\n')
        f.write('%d\n' % len(ks_in_ks_out))

        for t, (k_in, k_out) in enumerate(ks_in_ks_out):
            for i in range(4):
                ki_in = k_in[i]
                ki_out = k_out[i]
                f.write('# test %d, %s in\n' % (t + 1, order[i]))
                for e in ki_in:
                    f.write('%20.15e ' % e)
                f.write('\n')
                f.write('# test %d, %s out\n' % (t + 1, order[i]))
                for e in ki_out:
                    f.write('%20.15e ' % e)
                f.write('\n')
