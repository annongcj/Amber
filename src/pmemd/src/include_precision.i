! If no precision is defined then use DPDP
#ifdef pme_float
#undef pme_float
#endif

#ifdef pme_double
#undef pme_double
#endif

#ifdef real_if_SPDP
#undef real_if_SPDP
#endif

#define pme_float double precision
#define pme_double double precision
#define real_if_SPDP(x) (x)

#ifdef pmemd_SPDP
#undef pme_float
#undef pme_double
#undef real_if_SPDP
#define pme_float real
#define pme_double double precision
#define real_if_SPDP(x) (real((x)))
#endif

#ifdef pmemd_SPSP
#undef pme_float
#undef pme_double
#undef real_if_SPDP
#define pme_float real
#define pme_double real
#define real_if_SPDP(x) (real((x)))
#endif
