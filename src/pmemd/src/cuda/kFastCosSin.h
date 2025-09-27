#include "copyright.i"

//---------------------------------------------------------------------------------------------
// AMBER NVIDIA CUDA GPU IMPLEMENTATION: PMEMD VERSION
//
// July 2017, by Scott Le Grand, David S. Cerutti, Daniel J. Mermelstein, Charles Lin, and
//               Ross C. Walker
//---------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------------------------
// This special version of sincos is designed for |a| < 6*PI. On a GTX 285 it is about 25%
// faster than sincos from the CUDA math library.  Also uses 8 fewer registers than the CUDA
// math library's sincos.  Maximum observed error is 2 ulps across range stated above.
// Infinities and negative zero are not handled according to C99 specifications.  NaNs are
// handled fine.
//---------------------------------------------------------------------------------------------
{
  PMEDouble t, u, s, c, j, a2;
  int i;

  i = __double2int_rn(a * (PMEDouble)6.3661977236758138e-1);
  j = (PMEDouble)i;
  a = __fma_rn(-j, (PMEDouble)1.57079632679489660e+000, a); // PIO2_HI
  a = __fma_rn(-j, (PMEDouble)6.12323399573676600e-017, a); // PIO2_LO
  a2 = a * a;
  u =                 (PMEDouble)-1.136788825395985E-011;
  u = __fma_rn(u, a2, (PMEDouble)2.087588480545065E-009);
  u = __fma_rn(u, a2, (PMEDouble)-2.755731555403950E-007);
  u = __fma_rn(u, a2, (PMEDouble)2.480158729365970E-005);
  u = __fma_rn(u, a2, (PMEDouble)-1.388888888888074E-003);
  u = __fma_rn(u, a2, (PMEDouble)4.166666666666664E-002);
  u = __fma_rn(u, a2, (PMEDouble)-5.000000000000000E-001);
  u = __fma_rn(u, a2, (PMEDouble)1.000000000000000E+000);
  t =                 (PMEDouble)1.5896230157221844E-010;
  t = __fma_rn(t, a2, (PMEDouble)-2.5050747762850355E-008);
  t = __fma_rn(t, a2, (PMEDouble)2.7557313621385676E-006);
  t = __fma_rn(t, a2, (PMEDouble)-1.9841269829589539E-004);
  t = __fma_rn(t, a2, (PMEDouble)8.3333333333221182E-003);
  t = __fma_rn(t, a2, (PMEDouble)-1.6666666666666630E-001);
  t = t * a2;
  t = __fma_rn(t, a, a);
  if (i & 1) {
    s = u;
    c = t;
  }
  else {
    s = t;
    c = u;
  }
  if (i & 2) {
    s = -s;
  }
  i++;
  if (i & 2) {
    c = -c;
  }
  *sptr = s;
  *cptr = c;
}
