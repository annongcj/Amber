!*******************************************************************************
!
! Module: mdin_emil_dat_mod
!
! Description: read in (mostly optional) extra config information for EMIL calcs
!              
!*******************************************************************************

#ifdef EMIL

module mdin_emil_dat_mod


  use file_io_dat_mod, only : mdin, mdout
  use file_io_mod,     only : nmlsrc  
  use gbl_constants_mod

  implicit none
  private

    
  ! Public variables
  character(256), public   :: emil_paramfile, emil_logfile, emil_model_infile, &
                              emil_model_outfile

  ! Subroutines
  public   init_emil_dat
  private  validate_and_log_filename 

  contains

  
!*******************************************************************************
!
! Subroutine:  init_emil_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_emil_dat()
   implicit none

   ! Local variables:
   integer :: ifind
   namelist /emil_cntrl/       emil_paramfile, emil_logfile, emil_model_infile, &
                               emil_model_outfile

   !let the emil module set the defaults if we have no preferences ourselves
   emil_paramfile     = ''
   emil_logfile       = ''
   emil_model_infile  = ''
   emil_model_outfile = ''

   ! Read input in namelist format:
   rewind(mdin)                  ! Search the mdin file from the beginning
   call nmlsrc('emil_cntrl', mdin, ifind)

   if (ifind .ne. 0) then        

      ! Namelist found. Read it:
      read(mdin, nml = emil_cntrl)
      write(mdout, '(a, a)') info_hdr,       'Found an "emil_cntrl" namelist'
      call validate_and_log_filename(emil_paramfile, "emil_paramfile")
      call validate_and_log_filename(emil_logfile,   "emil_logfile")
      call validate_and_log_filename(emil_model_infile, "emil_model_infile")
      call validate_and_log_filename(emil_model_outfile, "emil_model_outfile")
      write(mdout, '(a)') ''
   else
      write(mdout, '(a, a)') info_hdr,  &
        'Did not find an "emil_cntrl" namelist: assuming default emil i/o files'
      
   end if

end subroutine init_emil_dat

subroutine validate_and_log_filename(filename, id_string)

   implicit none

   !arguments
   character(*), intent(in) :: id_string
   character(*), intent(in) :: filename
 
   !local variables
   integer                  :: string_length

   string_length = len_trim(filename)
   if( string_length .ge. 255 ) then
           write(mdout,'(a,a)') warn_hdr, &
                        id_string // ' is >= max length! Probably cropped.'
   end if
   if( string_length .gt. 0 ) then
           write(mdout,'(a,a)') extra_line_hdr, &
                        id_string // ' set to: ' // trim(filename)
   else
           write(mdout,'(a,a)') extra_line_hdr, &
                        id_string // ' not set in namelist, using default.'
   end if
       
end subroutine validate_and_log_filename

end module mdin_emil_dat_mod
#endif


