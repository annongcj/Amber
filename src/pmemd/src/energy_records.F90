
! Module that holds the energy records for PME and GB to avoid cyclic
! dependencies.

module energy_records_mod

  implicit none

  !!! IMPORTANT !!!
  ! IF YOU CHANGE THIS DERIVED TYPE, YOU ****MUST**** UPDATE THE STRUCT IN
  ! gputypes.h OR TRANSFERRING BETWEEN THE C AND F90 CODES WILL BE SUBTLY WRONG
  !!! IMPORTANT !!!
  type gb_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot
    double precision    :: elec_tot
    double precision    :: gb
    double precision    :: surf
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
    double precision    :: angle_ub !CHARMM Urey-Bradley Energy
    double precision    :: imp      !CHARMM Improper Energy
    double precision    :: cmap     !CHARMM CMAP
    double precision    :: dvdl
    double precision    :: amd_boost ! AMD boost energy
    double precision    :: gamd_boost ! GAMD boost energy
    double precision    :: emap ! EMAP restraint energy
    double precision    :: nfe  ! restraint energy in NFE module
    double precision    :: phmd ! virtual particles in continuous CpHMD
  end type gb_pot_ene_rec

  integer, parameter    :: gb_pot_ene_rec_size = 20

  type(gb_pot_ene_rec), parameter      :: null_gb_pot_ene_rec = &
    gb_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                   0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

  ! Potential energies, with breakdown, from pme.  This is intended to be the
  ! external interface to potential energy data produced by this module, in
  ! particular by subroutine pme_force().

  !!! IMPORTANT !!!
  ! IF YOU CHANGE THIS DERIVED TYPE, YOU ****MUST**** UPDATE THE STRUCT IN
  ! gputypes.h OR TRANSFERRING BETWEEN THE C AND F90 CODES WILL BE SUBTLY WRONG
  !!! IMPORTANT !!!
  type pme_pot_ene_rec
    sequence
    double precision    :: total
    double precision    :: vdw_tot      ! total of dir, recip
    double precision    :: vdw_dir
    double precision    :: vdw_recip
    double precision    :: elec_tot     ! total of dir, recip, nb_adjust, self
    double precision    :: elec_dir
    double precision    :: elec_recip
    double precision    :: elec_nb_adjust
    double precision    :: elec_self
    double precision    :: hbond
    double precision    :: bond
    double precision    :: angle
    double precision    :: dihedral
    double precision    :: vdw_14
    double precision    :: elec_14
    double precision    :: restraint
    double precision    :: angle_ub !CHARMM Urey-Bradley Energy
    double precision    :: imp      !CHARMM Improper Energy
    double precision    :: cmap     !CHARMM CMAP
    double precision    :: amd_boost ! AMD boost energy
    double precision    :: gamd_boost ! GAMD boost energy
    double precision    :: emap ! EMAP restraint energy
    double precision    :: efield
    double precision    :: nfe  ! restraint energy in NFE module
    double precision    :: gamd_ppi  ! temp for envronment
    double precision    :: gamd_bond_ppi ! temp for bonded
    double precision    :: phmd ! virtual particles in continuous CpHMD
  end type pme_pot_ene_rec

  integer, parameter    :: pme_pot_ene_rec_size = 27

  type(pme_pot_ene_rec), parameter      :: null_pme_pot_ene_rec = &
    pme_pot_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0, &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                    0.d0,0.d0)
                         

  ! Virials, with breakdown, from pme.  This is intended to be the
  ! structure for storing virials; it currently is not exported but could be.

  type pme_virial_rec
    sequence
    double precision    :: molecular(3, 3)
    double precision    :: atomic(3, 3)
    double precision    :: elec_direct(3, 3)
    double precision    :: elec_nb_adjust(3, 3)
    double precision    :: elec_recip(3, 3)
    double precision    :: elec_recip_vdw_corr(3, 3)
    double precision    :: elec_recip_self(3, 3)
    double precision    :: elec_14(3, 3)
    double precision    :: ep_frame(3, 3)       ! for extra points...
    double precision    :: eedvir               ! used in Darden's error est.
  end type pme_virial_rec

  integer, parameter    :: pme_virial_rec_size = 82

  type(pme_virial_rec), parameter      :: null_pme_virial_rec = &
    pme_virial_rec(9*0.d0, 9*0.d0, 9*0.d0, &
                   9*0.d0, 9*0.d0, 9*0.d0, &
                   9*0.d0, 9*0.d0, 9*0.d0, &
                   0.d0)

!we need to be able to write out softcore energies for afe. this seemed 
!preferrable to passing in individual variables as far as clarity was 
!concerned given the large number of terms
  type afe_gpu_sc_ene_rec
    sequence
    double precision    :: dvdl !ti dvdl
    double precision    :: bond_R1
    double precision    :: bond_R2
    double precision    :: angle_R1
    double precision    :: angle_R2
    double precision    :: dihedral_R1
    double precision    :: dihedral_R2
    double precision    :: sc_res_dist_R1
    double precision    :: sc_res_dist_R2
    double precision    :: sc_res_ang_R1
    double precision    :: sc_res_ang_R2
    double precision    :: sc_res_tors_R1
    double precision    :: sc_res_tors_R2
    double precision    :: vdw_dir_R1
    double precision    :: vdw_dir_R2
    double precision    :: elec_dir_R1
    double precision    :: elec_dir_R2
    double precision    :: vdw_14_R1
    double precision    :: vdw_14_R2
    double precision    :: elec_14_R1
    double precision    :: elec_14_R2
    double precision    :: vdw_der_R1
    double precision    :: vdw_der_R2
    double precision    :: elec_der_R1
    double precision    :: elec_der_R2
  end type afe_gpu_sc_ene_rec

  integer, parameter    :: afe_gpu_sc_ene_rec_size = 25 

  type(afe_gpu_sc_ene_rec), parameter   :: null_afe_gpu_sc_ene_rec =  &
    afe_gpu_sc_ene_rec(0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,  &
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,&
                    0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)

end module energy_records_mod
