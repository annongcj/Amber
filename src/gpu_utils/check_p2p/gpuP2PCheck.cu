/*
 * This code prints out information about all available GPUs,
 * and tests each possible pair for P2P compatibility. It also
 * checks that the compute mode is correctly set.
 *
 * By Robin M. Betz and Ross C. Walker
 * May 13, 2014
 */

#include <stdio.h>
#include <assert.h>
#include <cuda_runtime_api.h>
#include <helper_cuda.h>

#define MAX_DEVICES          8
#define PROCESSES_PER_DEVICE 1
#define DATA_BUF_SIZE        4096

#ifdef __linux
#include <unistd.h>
#include <sched.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <linux/version.h>

void getDeviceCount()
{
    // We can't initialize CUDA before fork() so we need to spawn a new process

    pid_t pid = fork();

    if (0 == pid)
    {
        int i,j;
        int count, uvaCount = 0;
        int uvaOrdinals[MAX_DEVICES];
        checkCudaErrors(cudaGetDeviceCount(&count));

        printf("CUDA-capable device count: %i\n", count);
        for (i=0; i<count; ++i) {
            cudaDeviceProp prop;
            checkCudaErrors(cudaGetDeviceProperties(&prop, i));

            if (prop.computeMode != cudaComputeModeDefault) {
                printf(" GPUs must be in Compute Mode Default to run\n");
                printf(" Please run the following as root to change the Compute Mode: \n");
                printf("   nvidia-smi -c 0\n"); 
                exit(EXIT_FAILURE);
            }
 
            if (prop.unifiedAddressing) {
                uvaOrdinals[uvaCount] = i;
                printf("   GPU%d \"%15s\"\n", i, prop.name);
                uvaCount += 1;
            } else
                printf("   GPU%d \"%15s\"     NOT UVA capable\n", i, prop.name);
       }

        if (uvaCount < 2) {
          printf("Fewer than 2 UVA-capable GPUs found.\n");
          printf("Skipping peer it peer access test.\n");
          exit(EXIT_FAILURE);
        }

        // Check possibility for peer accesses, relevant to our tests
        int canAccessPeer_ij, canAccessPeer_ji;

        printf("\nTwo way peer access between:\n");
        for (i=0; i<uvaCount; ++i) {
          for (j=i+1; j<uvaCount; ++j) {
            checkCudaErrors(cudaDeviceCanAccessPeer(&canAccessPeer_ij, uvaOrdinals[i], uvaOrdinals[j]));
            checkCudaErrors(cudaDeviceCanAccessPeer(&canAccessPeer_ji, uvaOrdinals[j], uvaOrdinals[i]));

            if (canAccessPeer_ij*canAccessPeer_ji) {
                printf("   GPU%d and GPU%d: YES\n", uvaOrdinals[i], uvaOrdinals[j]);
            } else {
                printf("   GPU%d and GPU%d: NO\n", uvaOrdinals[i], uvaOrdinals[j]);
            }
          }
        }

        exit(EXIT_SUCCESS);
    } else
    {
        int status;
        waitpid(pid, &status, 0);
  
        if (status) {
          exit(EXIT_SUCCESS);
        }
    }
}

#endif
int main(int argc, char **argv)
{
#if CUDART_VERSION >= 4010 && defined(__linux)
    // Check this is built 64 bit
    #if !defined(__x86_64) && !defined(AMD64) && !defined(_M_AMD64)
        printf("%s is only supported on 64-bit Linux OS and the application must be built as a 64-bit target. Test is being waived.\n", argv[0]);
        exit(EXIT_WAIVED);
    #endif

    // Check if CUDA_VISIBLE_DEVICES is set
    char *cudaVis = getenv("CUDA_VISIBLE_DEVICES");
    if (cudaVis!=NULL) {
      printf("CUDA_VISIBLE_DEVICES=\"%s\"\n", cudaVis);
    } else {
      printf("CUDA_VISIBLE_DEVICES is unset.\n");
    }

    // Check device info, status, and P2P communication
    getDeviceCount();

    // Shut down
    cudaDeviceReset();

#else // Using CUDA 4.0 and older or non Linux OS
    printf("This test requires CUDA 4.1 and Linux to build and run, waiving testing\n\n");
    exit(EXIT_WAIVED);
#endif
}
