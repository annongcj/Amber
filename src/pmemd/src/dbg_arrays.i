#ifdef DBG_ARRAYS
#define DBG_ARRAYS_TIME_STEP(nstep)  call timestep_dbg_arrays(nstep)
#define DBG_ARRAYS_DUMP_3DBLE(tag, id, dble3, N) call DBG_ARRAYS_DUMP_3DBLE_FUN(tag, id, dble3, N)
#define DBG_ARRAYS_DUMP_CRD(tag, id, dble3, N, wrap) call DBG_ARRAYS_DUMP_CRD_FUN(tag, id, dble3, N, wrap)

#else
#define DBG_ARRAYS_TIME_STEP(nstep)
#define DBG_ARRAYS_DUMP_3DBLE(tag, id, dble3, N)
#define DBG_ARRAYS_DUMP_CRD(tag, id, dble3, N, wrap)

#endif


