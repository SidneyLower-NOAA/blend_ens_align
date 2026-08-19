MODULE align_mod
  IMPLICIT NONE
  SAVE  ! Ensures variables retain their values across subroutine calls

  ! --------------------------------------------------------
  ! Parameters
  ! --------------------------------------------------------
  INTEGER, PARAMETER :: nvarshf = 1
  INTEGER, PARAMETER :: mx_shfpass = 5
  INTEGER, PARAMETER :: mx_members = 25

  ! --------------------------------------------------------
  ! mpi_vars
  ! --------------------------------------------------------
  INTEGER :: nprocs
  INTEGER :: myproc

  ! --------------------------------------------------------
  ! grid data
  ! --------------------------------------------------------
  INTEGER :: ibgn, iend, jbgn, jend

  ! --------------------------------------------------------
  ! align_parm
  ! --------------------------------------------------------
  INTEGER :: nshfpass, minkdat,applyshft
  INTEGER :: nbaksmth, nshfsmth, noutsmth
  REAL    :: hrzlap, minkdratio, slnratio0h, slnratio48h, slen, reflmin

  ! --------------------------------------------------------
  ! align_zones
  ! --------------------------------------------------------
  INTEGER :: iborder(mx_shfpass), jborder(mx_shfpass)
  INTEGER :: izsize(mx_shfpass), jzsize(mx_shfpass)
  INTEGER :: nizone(mx_shfpass), njzone(mx_shfpass)
  INTEGER :: loopstep(mx_shfpass), procspg(mx_shfpass)

  ! --------------------------------------------------------
  ! shift_wgts
  ! --------------------------------------------------------
  REAL    :: wgtvar, threshvar
  INTEGER :: thresh_flag

  ! --------------------------------------------------------
  ! lpm_parm
  ! --------------------------------------------------------
  INTEGER :: patch_nx, patch_ny, ovx, ovy, gauss_sigma
  REAL    :: filt_min

END MODULE align_mod