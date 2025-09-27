#ifndef AMBER_XRAY_BULK_MASK_GPU_H
#define AMBER_XRAY_BULK_MASK_GPU_H

#include "xray/BulkMask.h"
#include <thrust/device_vector.h>
#ifdef AMBER_PLATFORM_AMD
#include <hipfft.h>

#define CUFFT_SUCCESS HIPFFT_SUCCESS
#define CUFFT_D2Z HIPFFT_D2Z
#define CUFFT_Z2D HIPFFT_Z2D
#define CUFFT_R2C HIPFFT_R2C
#define CUFFT_C2R HIPFFT_C2R
#define cufft hipfft
#define cufftComplex hipfftComplex
#define cufftDestroy hipfftDestroy
#define cufftDoubleComplex hipfftDoubleComplex
#define cufftExecC2R hipfftExecC2R
#define cufftExecD2Z hipfftExecD2Z
#define cufftExecR2C hipfftExecR2C
#define cufftExecZ2D hipfftExecZ2D
#define cufftHandle hipfftHandle
#define cufftPlan3d hipfftPlan3d

#else
#include <cufft.h>
#endif

namespace xray {

class BulkMaskGPU : public BulkMask {
public:
  BulkMaskGPU() = default;
  ~BulkMaskGPU() override;
  void init(xray::UnitCell unit_cell, int n_atom, int *mask_grid_size,
            double *reciprocal_norms, double *mask_cutoffs,
            int *mask_bs_grid, double shrink_r,
            int n_hkl, complex_double* f_mask, int* hkl) override;
  void update_grid(int n_atom, const double *frac) override;

  void calc_f_bulk() override;

private:
  void shrink();
  thrust::device_vector<int> m_dev_mask_grid;
  thrust::device_vector<int> m_dev_mask_grid_not_shrank;
  thrust::device_vector<GridPoint> m_dev_shrink_neighbours;
  thrust::device_vector<double> m_dev_mask_cutoffs;
  thrust::device_vector<double> m_dev_frac;
  thrust::device_ptr<Sym33> m_dev_metric_tensor;

  thrust::device_vector<double> m_dev_fft_in;
  thrust::device_vector<cufftDoubleComplex> m_dev_fft_out;

  std::vector<cufftDoubleComplex> m_fft_out;

  cufftHandle m_plan;
};
} // namespace xray

#endif // AMBER_XRAY_BULK_MASK_GPU_H
