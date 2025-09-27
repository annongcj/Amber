#ifdef GTI
#  ifndef GTI_KERNEL_ATTRIBUTE
#    define GTI_KERNEL_ATTRIBUTE

// HIP-TODO: Remove after investigation why some kernels don't work correctly with large blocks
#if !defined(AMBER_PLATFORM_AMD)
static const int MAX_THREADS_PER_BLOCK = 1024;
#define _kPlainHead_ __global__ void
#else
static const int MAX_THREADS_PER_BLOCK = 256;
#define _kPlainHead_ __global__ void __LAUNCH_BOUNDS__(MAX_THREADS_PER_BLOCK, 1)
#endif

#define _kReduceFrcHead_  __global__ void __launch_bounds__(REDUCEFORCES_THREADS_PER_BLOCK, 1)

#  endif /* GTI_KERNEL_ATTRIBUTE */
#endif /* GTI */
