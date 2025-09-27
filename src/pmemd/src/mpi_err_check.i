! Place after the MPI calls that return error status. error status message variable name should be err_code_mpi. Uses variables from parallel_dat module 
#ifdef CHECK_MPI_ERRORS
#define MPI_ERR_CHK(message) \
  IF (err_code_mpi .ne. MPI_SUCCESS) then; \
    write(0,*) "### MPI COMMUNICATION ERROR in ", message; \
    call MPI_Error_string(err_code_mpi, mpi_err_msg_str, mpi_error_msg_len, mpi_err_tmp); \
    write(0,*) mpi_err_msg_str(1:mpi_error_msg_len); \
  End if
#else
#define MPI_ERR_CHK(message)
#endif
