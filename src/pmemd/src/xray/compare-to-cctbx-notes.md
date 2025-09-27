### Correspondence between cctbx and `xray_*_module`s

[comment]: <> (TODO)

### Discrepancy between cctbx and `xray_*_module`s
 

### In `mmtbx.bulk_solvent.scaler.run.k_mask_grid_search`:

 - Avoid division `Fobs` on `k_total` [gh#663](https://github.com/cctbx/cctbx_project/issues/663)
 - [*] Exclude `k_iso` from `k_total`
 - We build grid of `[0, 1]` tiny bit differently than cctbx. It manifests itself as numerical discrepancies in subsequent optimization cycles 

### In `mmtbx.bulk_solvent.scaler.run.bulk_solvent_scaling`:

- [*] Exclude `k_iso` from `k_total`
- [*] We use regular grid intersected with `[0, 1]` unlike cctbx

### Issues to file

1. Fitting done with `k_iso=1` and `k_aniso`, but tested with current values
   [🔗](https://github.com/cctbx/cctbx_project/blob/f476ac7af391357632757bf698f269333d1756e3/mmtbx/bulk_solvent/scaler.py#L247)
   [gh#655](https://github.com/cctbx/cctbx_project/issues/655)

2. `r_best` should be calculated over `selection`, not whole set
   [🔗](https://github.com/cctbx/cctbx_project/blob/f476ac7af391357632757bf698f269333d1756e3/mmtbx/bulk_solvent/bulk_solvent.h#L1151)
   [gh#654](https://github.com/cctbx/cctbx_project/issues/654)

3. `k_overall_best` should be initialized in accordance with `r_best`
[🔗](https://github.com/cctbx/cctbx_project/blob/f476ac7af391357632757bf698f269333d1756e3/mmtbx/bulk_solvent/bulk_solvent.h#L1150-L1151)
   [gh#654](https://github.com/cctbx/cctbx_project/issues/654)

4. Interpolation values for the last (best) resolution bin makes no sense
   [🔗](https://github.com/cctbx/cctbx_project/blob/f476ac7af391357632757bf698f269333d1756e3/mmtbx/bulk_solvent/scaler.py#L337)

5. `k_total` should use only work reflexes in
   [[1]](https://github.com/cctbx/cctbx_project/blob/a15fcd9b500288fe471a47b6e32eabf89f78701e/mmtbx/bulk_solvent/scaler.py#L372) and
   [[2]](https://github.com/cctbx/cctbx_project/blob/a15fcd9b500288fe471a47b6e32eabf89f78701e/mmtbx/bulk_solvent/scaler.py#L490)