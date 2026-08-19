! MAIN PROGRAM: ensalign.f90
! Program to create ensemble mean after spatially aligning results
!
! Author: Keith Brewster
! CAPS/University of Oklahoma
! March, 2017
!
! Modification:
! MPI Version
!    Keith Brewster, CAPS, October 2017
!
!
! Modification:
!    ChangJae Lee, KMA, October 2022
!    - Use eccodes to read Grib file
!    - Write shifted means (also include PM and LPM) and shift vectors
!      into NETCDF file
!
! Modification:
!    Sidney Lower, NWS/OMD/WPGD August 2026
!    - refactor program to subroutine for f2py
!    - use xarray for NetCDF & GRIB2 IO
    
    
    
SUBROUTINE ensalign(ensfcst,py_ny,py_nx,py_dy,py_dx,nmembers,ifhr, &  ! <-- in
                    py_nbaksmth,py_nshfsmth,py_applyshft, &
                    py_slnratio0h,py_slnratio48h, py_nshfpass, &
                    py_noutsmth,py_minkdat,py_minkdratio, &  ! <-- shift_const
                    py_hrzlap,py_iborder,py_izsize,py_jborder, &
                    py_jzsize,py_loopstep,py_procspg, & ! <-- shift_zone
                    py_wgtvar, py_thresh_flag, py_threshvar, & ! <-- shift_weight
                    py_patch_nx, py_patch_ny, py_ovx, py_ovy, &
                    py_gauss_sigma, py_filt_min, & ! <-- lpm_const
                    xshiftmn, yshiftmn, ensfcst_shf, ensshfpm, ensshflpm) ! <-- out

    
    use mpi
    use align_mod
    IMPLICIT NONE

    ! parameters
    INTEGER, PARAMETER :: root = 0
    INTEGER, PARAMETER :: stdout = 6
    INTEGER, PARAMETER :: max_elem_send = 100000
    REAL,    PARAMETER :: smt_param = 0.5
    REAL,    PARAMETER :: rmisg_data = -999.
    INTEGER, PARAMETER :: posdef=1
    LOGICAL, PARAMETER :: lposdef=.TRUE.
    

    ! input
    INTEGER, INTENT(IN) :: py_nx,py_ny,py_dy,py_dx
    INTEGER, INTENT(IN) :: nmembers
    REAL,    INTENT(IN) :: ensfcst(py_nx,py_ny,nmembers)
    INTEGER, INTENT(IN) :: ifhr ! lead time?

    ! config / namelist
    INTEGER, INTENT(IN) :: py_nshfpass,py_noutsmth,py_nbaksmth,py_nshfsmth,py_applyshft
    REAL,    INTENT(IN) :: py_minkdat,py_minkdratio
    REAL,    INTENT(IN) :: py_hrzlap, py_slnratio0h,py_slnratio48h
    INTEGER, INTENT(IN) :: py_iborder(py_nshfpass),py_jborder(py_nshfpass),py_izsize(py_nshfpass)
    INTEGER, INTENT(IN) :: py_jzsize(py_nshfpass),py_loopstep(py_nshfpass), py_procspg(py_nshfpass)
    REAL,    INTENT(IN) :: py_wgtvar, py_threshvar
    INTEGER, INTENT(IN) :: py_thresh_flag
    INTEGER, INTENT(IN) :: py_patch_nx,py_patch_ny,py_ovx,py_ovy,py_gauss_sigma,py_filt_min

    ! output
    REAL, INTENT(OUT) :: xshiftmn(py_nx,py_ny,nmembers,py_nshfpass)
    REAL, INTENT(OUT) :: yshiftmn(py_nx,py_ny,nmembers,py_nshfpass)
    REAL, INTENT(OUT) :: ensfcst_shf(py_nx,py_ny,nmembers,py_nshfpass)
    REAL, INTENT(OUT) :: ensshflpm(py_nx,py_ny,py_nshfpass),ensshfpm(py_nx,py_ny,py_nshfpass)

    ! internal
    INTEGER :: ipass,mxzone
    INTEGER :: i,j,k,k1,k2,ny,nx
    INTEGER :: nelem2d,nelem3d
    INTEGER :: ierr
    INTEGER :: ilap,jlap,istep,jstep
    REAL :: dx,dy
    REAL :: rninv,fnorm,fcstmin,fcstmax,favgmin,favgmax
    REAL :: time,vmin,vmax
    REAL :: cput0,cput1,cput2,cput3,cput4,cput5,cput6,cput7

    REAL, ALLOCATABLE :: xs(:)
    REAL, ALLOCATABLE :: ys(:)
    REAL, ALLOCATABLE :: ensmean(:,:)
    REAL, ALLOCATABLE :: enspm(:,:)
    REAL, ALLOCATABLE :: enslpm(:,:)
    REAL, ALLOCATABLE :: ensmax(:,:)
    REAL, ALLOCATABLE :: ensfcst_smt(:,:,:,:)
    REAL, ALLOCATABLE:: enspm_smt(:,:,:)
    REAL, ALLOCATABLE :: ensshfmn(:,:,:)
    REAL, ALLOCATABLE :: ensshfmax(:,:,:)
    REAL, ALLOCATABLE :: xshift(:,:)
    REAL, ALLOCATABLE :: yshift(:,:)
    INTEGER, ALLOCATABLE :: ibkshift(:,:,:,:), jbkshift(:,:,:,:)
    REAL, ALLOCATABLE :: xsum(:,:)
    REAL, ALLOCATABLE :: ysum(:,:)
    REAL, ALLOCATABLE :: wgtsum(:,:)
    REAL, ALLOCATABLE :: tem2d1(:,:)
    REAL, ALLOCATABLE :: tem2d2(:,:)
    REAL, ALLOCATABLE :: tem2d3(:,:) 
    REAL, ALLOCATABLE :: recv_buf(:,:)
    INTEGER, ALLOCATABLE :: zero2d(:,:) 
    CHARACTER(LEN=1), ALLOCATABLE :: mpi_work_buf(:)
    REAL, ALLOCATABLE :: istart(:)
    REAL, ALLOCATABLE :: jstart(:)
    REAL, ALLOCATABLE :: ifinish(:)
    REAL, ALLOCATABLE :: jfinish(:)
    
    REAL, ALLOCATABLE :: dxfld(:)
    REAL, ALLOCATABLE :: dyfld(:)
    REAL, ALLOCATABLE :: rdxfld(:)
    REAL, ALLOCATABLE :: rdyfld(:)
    REAL, ALLOCATABLE :: slopey(:,:)
    REAL, ALLOCATABLE :: alphay(:,:)
    REAL, ALLOCATABLE :: betay(:,:)
    REAL, ALLOCATABLE :: gs_weight(:)
    
    LOGICAL :: is_init


    nbaksmth = 0
    iborder  = 0
    jborder  = 0
    izsize   = 0
    jzsize   = 0
    loopstep = 0
    procspg  = 0
    nizone   = 0
    njzone   = 0
    slen     = 0.0
    reflmin  = 0.0
    ibgn     = 1
    jbgn     = 1

    ! -----------------------------------------------------------------------
    ! Populate Module Variables from Python Inputs
    ! -----------------------------------------------------------------------
    ny         = py_ny
    nx         = py_nx
    dx         = py_dx
    dy         = py_dy
    iend       = nx
    jend       = ny
    
    nshfpass   = py_nshfpass
    nbaksmth   = py_nbaksmth
    noutsmth   = py_noutsmth
    nshfsmth   = py_nshfsmth
    minkdat    = py_minkdat
    minkdratio = py_minkdratio
    slnratio0h = py_slnratio0h
    slnratio48h= py_slnratio48h
    applyshft  = py_applyshft
    hrzlap     = py_hrzlap
    
    patch_nx    = py_patch_nx
    patch_ny    = py_patch_ny
    ovx         = py_ovx
    ovy         = py_ovy
    gauss_sigma = py_gauss_sigma
    filt_min    = py_filt_min

    wgtvar      = py_wgtvar
    threshvar   = py_threshvar
    thresh_flag = py_thresh_flag
    
    ! Map array inputs to module arrays for active passes
    iborder(1:nshfpass)  = py_iborder(1:nshfpass)
    jborder(1:nshfpass)  = py_jborder(1:nshfpass)
    izsize(1:nshfpass)   = py_izsize(1:nshfpass)
    jzsize(1:nshfpass)   = py_jzsize(1:nshfpass)
    loopstep(1:nshfpass) = py_loopstep(1:nshfpass)
    procspg(1:nshfpass)  = py_procspg(1:nshfpass)

    ! Check if MPI is already running
    CALL MPI_INITIALIZED(is_init, ierr)
    
    IF (.NOT. is_init) THEN
        CALL MPI_INIT(ierr)
    END IF
    
    CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nprocs, ierr)
    CALL MPI_COMM_RANK (MPI_COMM_WORLD, myproc, ierr)

    
    IF(myproc == root) THEN
        WRITE(*, '(//a)') ' ----------------------------------------'
        WRITE(*, '(a)') ' --    ENTERING ENSEMBLE ALIGNMENT     --'
        WRITE(*, '(a)') ' ----------------------------------------'
        WRITE(6,'(//a,i5)') ' Number of processors: ',nprocs
        WRITE(6,'(a,i5)') ' Number of ensemble members: ',nmembers
        CALL CPU_TIME(cput1)
    END IF
    
    !!! begin ---------------------------------------
    
    mxzone=0
    DO ipass=1,nshfpass
        ilap=MAX(IFIX((izsize(ipass)*hrzlap)+0.5),1)
        jlap=MAX(IFIX((jzsize(ipass)*hrzlap)+0.5),1)
        istep=izsize(ipass)-ilap
        jstep=jzsize(ipass)-jlap
        nizone(ipass)=MAX((nx-(2*iborder(ipass)))/istep,1)
        njzone(ipass)=MAX((ny-(2*jborder(ipass)))/jstep,1)
        mxzone=MAX(mxzone,(nizone(ipass)*njzone(ipass)))
        IF (myproc == root) THEN
            WRITE(6,'(4(a,i0))') ' nizone(', ipass, ') : ', nizone(ipass), &
                ', njzone(', ipass, '): ', njzone(ipass)
        END IF
    END DO

    ALLOCATE(xs(nx))
    ALLOCATE(ys(ny))
    ALLOCATE(ensmax(nx,ny))
    ALLOCATE(ensmean(nx,ny))
    ALLOCATE(enspm(nx,ny))
    ALLOCATE(enslpm(nx,ny))
    ALLOCATE(xshift(nx,ny))
    ALLOCATE(yshift(nx,ny))
    ALLOCATE(ibkshift(nx,ny,nmembers,nshfpass))
    ALLOCATE(jbkshift(nx,ny,nmembers,nshfpass))
    ALLOCATE(xsum(nx,ny))
    ALLOCATE(ysum(nx,ny))
    ALLOCATE(wgtsum(nx,ny))
    ALLOCATE(ensshfmn(nx,ny,nshfpass),ensshfmax(nx,ny,nshfpass))
    ALLOCATE(gs_weight(nx*ny))
    
    ALLOCATE(ensfcst_smt(nx,ny,nmembers,nshfpass))
    ALLOCATE(enspm_smt(nx,ny,nshfpass))
    
    ALLOCATE(istart(mxzone))
    ALLOCATE(jstart(mxzone))
    ALLOCATE(ifinish(mxzone))
    ALLOCATE(jfinish(mxzone))
    
    ALLOCATE(tem2d1(nx,ny))
    ALLOCATE(tem2d2(nx,ny))
    ALLOCATE(tem2d3(nx,ny))
    ALLOCATE(recv_buf(nx,ny))
    ALLOCATE(dxfld(nx))
    ALLOCATE(dyfld(ny))
    ALLOCATE(rdxfld(nx))
    ALLOCATE(rdyfld(ny))
    ALLOCATE(slopey(nx,ny))
    ALLOCATE(alphay(nx,ny))
    ALLOCATE(betay(nx,ny))
    ALLOCATE(zero2d(nx,ny))


    
    dx = 2500.
    dy = 2500.
    
    ibgn = 1
    jbgn = 1
    iend = nx
    jend = ny
    
    
    IF(myproc == root) THEN
        CALL gaussian_weight(nx,ny,gauss_sigma,gs_weight)
    END IF
    
    IF (gauss_sigma > 0) THEN
        IF ( nprocs > 1 ) THEN
          nelem2d=nx*ny
          IF( nelem2d < max_elem_send) THEN
            CALL MPI_BCAST(gs_weight,nelem2d,MPI_REAL,root, &
                       MPI_COMM_WORLD,ierr)
          ELSE
            CALL LGARRAY_BCASTR(gs_weight,nelem2d,max_elem_send, &
                            myproc,root,MPI_COMM_WORLD,ierr)
          END IF
        END IF
    END IF
    
    DO i=1,nx
        xs(i)=(float(i-2)+0.5)*dx
    END DO
    DO j=1,ny
        ys(j)=(float(j-2)+0.5)*dy
    END DO
    
    CALL setdxdy(nx,ny,ibgn,iend,jbgn,jend, &
           xs,ys,dxfld,dyfld,rdxfld,rdyfld)

    IF( myproc == root ) THEN
        ensmean(:,:) = -999.
        ensmean(ibgn:iend,jbgn:jend) = 0.
        DO k=1,nmembers
          fcstmax=-999.
          fcstmin=999.
          DO j=jbgn,jend
            DO i=ibgn,iend
              fcstmin=min(fcstmin,ensfcst(i,j,k))
              fcstmax=max(fcstmax,ensfcst(i,j,k))
              ensmean(i,j)=ensmean(i,j)+ensfcst(i,j,k)
              ensmax(i,j)=max(ensmax(i,j),ensfcst(i,j,k))
            END DO
          END DO
        END DO
        rninv=1.0/float(nmembers)
        fcstmax=-999.
        fcstmin=999.
        DO j=jbgn,jend
          DO i=ibgn,iend
            ensmean(i,j)=rninv*ensmean(i,j)
            IF(lposdef) ensmean(i,j)=max(0.,ensmean(i,j))
            favgmin=min(favgmin,ensmean(i,j))
            favgmax=max(favgmax,ensmean(i,j))
          END DO
        END DO
        
        CALL pm_mean(nx,ny,nmembers,RESHAPE(ensfcst(:,:,1:nmembers),(/nx*ny*nmembers/)), &
                  ensmean,enspm(:,:))
        
        enspm_smt(:,:,1) = enspm(:,:)
        DO i=1,nbaksmth
        CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,smt_param, &
                  enspm(:,:),tem2d1,enspm_smt(:,:,1))
        END DO
    END IF
    
    !! Broadcast PM mean for alignopt = 2 (align to ensemble mean location)
    IF( nprocs > 1 ) THEN
      nelem2d=nx*ny
    
      IF( nelem2d < max_elem_send) THEN
        CALL MPI_BCAST(enspm_smt(:,:,1),nelem2d,MPI_REAL,root, &
                 MPI_COMM_WORLD,ierr)
      ELSE
        CALL LGARRAY_BCASTR(enspm_smt(:,:,1),nelem2d,max_elem_send, &
                    myproc,root,MPI_COMM_WORLD,ierr)
      END IF
    END IF

    IF( myproc == root ) THEN
        CALL CPU_TIME(cput2)
        WRITE(6,'(a,f10.2,a)') ' Set-up,initial ens mean: ',(cput2-cput1),' seconds'
    END IF
    
    !-----------------------------------------------------------------------
    ! Begin alignment
    ! Outer loop: number of shift iterations
    ! Inner loop: n ens members
    !-----------------------------------------------------------------------
    
    ensfcst_shf(:,:,:,:) = 0.
    ibkshift(:,:,:,:) = 0
    jbkshift(:,:,:,:) = 0
    fnorm=1.0/float(nmembers)
    zero2d(:,:) = 0
    
    DO ipass=1,nshfpass
        IF( myproc == root ) THEN
            CALL CPU_TIME(cput3)
            WRITE(6,'(//a,i4)') ' Alignment pass: ', ipass
        END IF
        DO k1=1,nmembers
            xshiftmn(:,:,k1,ipass) = 0.
            yshiftmn(:,:,k1,ipass) = 0.
        
            time = real(ifhr)
        
            IF( myproc == root ) THEN
               WRITE(6,'(//a,i4)') ' Processing ensemble member: ',k1
            END IF
        
            xshift(:,:) = 0.
            yshift(:,:) = 0.
            CALL rshift2dgrd(nx,ny,nvarshf,ipass,mx_shfpass,mxzone,    &
                   max_elem_send,posdef,time,xs,ys,              &
                   ensfcst_smt(:,:,k1,ipass),enspm_smt(:,:,ipass),         &
                   ensfcst_shf(:,:,k1,ipass),                              &
                   istart,jstart,ifinish,jfinish,                          &
                   ibkshift(:,:,k1,ipass),jbkshift(:,:,k1,ipass),          &
                   zero2d(:,:),zero2d(:,:),                                &
                   xshift,yshift,xsum,ysum,wgtsum,                         &
                   dxfld,dyfld,rdxfld,rdyfld,                              &
                   slopey,alphay,betay,                                    &
                   tem2d3,recv_buf)
                 
        !-----------------------------------------------------------------------
        ! Calculate mean shift for ensemble member k1. 
        ! Note that zero shift is implied for k1,k1, so denominator is nmembers.
        ! + Save backshift for next pass
        !-----------------------------------------------------------------------
        
            IF( myproc == root ) THEN
               xshiftmn(:,:,k1,ipass)=xshift(:,:)
               yshiftmn(:,:,k1,ipass)=yshift(:,:)
            
              IF (ipass .ne. nshfpass) THEN
                ibkshift(:,:,k1,ipass+1)=ibkshift(:,:,k1,ipass)+NINT(xshiftmn(:,:,k1,ipass))
                jbkshift(:,:,k1,ipass+1)=jbkshift(:,:,k1,ipass)+NINT(yshiftmn(:,:,k1,ipass))
              END IF
            
              IF (ipass .gt. 1) THEN
                xshiftmn(:,:,k1,ipass)=xshiftmn(:,:,k1,ipass)+xshiftmn(:,:,k1,ipass-1)
                yshiftmn(:,:,k1,ipass)=yshiftmn(:,:,k1,ipass)+yshiftmn(:,:,k1,ipass-1)
              END IF
            
              DO i=1,nshfsmth
                CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                            xshiftmn(:,:,k1,ipass),tem2d3,xshiftmn(:,:,k1,ipass))
                CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                            yshiftmn(:,:,k1,ipass),tem2d3,yshiftmn(:,:,k1,ipass))
              END DO
            END IF
        
        !-----------------------------------------------------------------------
        ! Apply shift to ensemble using mean shift and sum 
        ! to calculate ensemble mean (and broadcast ensemble mean)
        !-----------------------------------------------------------------------
            IF( myproc == root ) THEN
                WRITE(6,'(a)') ' Applying shift vectors'
                CALL a2dmax0(ensfcst(:,:,k1),1,nx,ibgn,iend,           &
                             1,ny,jbgn,jend,vmin,vmax)
                WRITE(6,'(1x,2(a,f13.4))')                             &
                      '  Pre-shift min = ', vmin,',  max =',vmax
                CALL movegr(nx,ny,ensfcst(:,:,k1),tem2d3,                    &
                        xshiftmn(:,:,k1,ipass),yshiftmn(:,:,k1,ipass),     &
                        ensfcst_shf(:,:,k1,ipass),                         &
                        ibgn,iend,jbgn,jend,                               &
                        dxfld,dyfld,rdxfld,rdyfld,                         &
                        slopey,alphay,betay)
                
                !! do smooth to shifted fields
                DO i=1,noutsmth
                    CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,0.5,                     &
                                ensfcst_shf(:,:,k1,ipass),tem2d3,ensfcst_shf(:,:,k1,ipass))
                END DO
                
                !! ** IMPORTANT **
                !! rescale shifted field with its original PDF
                CALL pm_mean(nx,ny,1,RESHAPE(ensfcst(:,:,k1),(/nx*ny/)), &
                ensfcst_shf(:,:,k1,ipass),ensfcst_shf(:,:,k1,ipass))
                
                CALL a2dmax0(ensfcst_shf(:,:,k1,ipass),1,nx,ibgn,iend,       &
                       1,ny,jbgn,jend,vmin,vmax)
                WRITE(6,'(1x,2(a,f13.4))')                                   &
                  '  Post-shift min = ', vmin,',  max =',vmax
                
                ensshfmn(:,:,ipass)=ensshfmn(:,:,ipass)+ensfcst_shf(:,:,k1,ipass)
            END IF

        
            IF( nprocs > 1 ) THEN
              nelem2d=nx*ny
            
              IF (ipass .ne. nshfpass) THEN
                IF( nelem2d < max_elem_send) THEN
                  CALL MPI_BCAST(ibkshift(:,:,k1,ipass+1),nelem2d,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
                ELSE
                  CALL LGARRAY_INT_BCASTR(ibkshift(:,:,k1,ipass+1),nelem2d,max_elem_send, &
                              myproc,root,MPI_COMM_WORLD,ierr)
                END IF
            
                IF( nelem2d < max_elem_send) THEN
                  CALL MPI_BCAST(jbkshift(:,:,k1,ipass+1),nelem2d,MPI_INTEGER,root,MPI_COMM_WORLD,ierr)
                ELSE
                  CALL LGARRAY_INT_BCASTR(jbkshift(:,:,k1,ipass+1),nelem2d,max_elem_send, &
                              myproc,root,MPI_COMM_WORLD,ierr)
                END IF
            
              END IF
            
            END IF
        
        !-----------------------------------------------------------------------
        ! Calculate shifted ensemble mean and PM
        !-----------------------------------------------------------------------
          IF( myproc == root ) THEN
              WRITE(6,'(a)') ' Calculating aligned ensemble mean'
                fnorm=1.0/float(nmembers)
                ensshfmn(:,:,ipass)=fnorm*ensshfmn(:,:,ipass)
                IF(lposdef) THEN
                  DO j=jbgn,jend
                    DO i=ibgn,iend
                      ensshfmn(i,j,ipass)=max(0.,ensshfmn(i,j,ipass))
                      DO k=1,nmembers
                        ensshfmax(i,j,ipass)=max(ensshfmax(i,j,ipass),ensfcst_shf(i,j,k,ipass))
                      END DO
                    END DO
                  END DO
                END IF
            
                CALL pm_mean(nx,ny,nmembers,RESHAPE(ensfcst_shf(:,:,1:nmembers,ipass),(/nx*ny*nmembers/)), &
                        ensshfmn(:,:,ipass),ensshfpm(:,:,ipass))
            
                IF (ipass .ne. nshfpass) THEN
                  enspm_smt(:,:,ipass+1) = ensshfpm(:,:,ipass)
            
                  DO i=1,nbaksmth
                    CALL smooth2d(nx,ny,ibgn,iend,jbgn,jend,smt_param,          &
                              ensshfpm(:,:,ipass),tem2d1,enspm_smt(:,:,ipass+1))
                  END DO
                END IF
          END IF
        
          IF (ipass .ne. nshfpass) THEN
            IF( nprocs > 1 ) THEN
              nelem2d=nx*ny
        
              IF( nelem2d < max_elem_send) THEN
                CALL MPI_BCAST(enspm_smt(:,:,ipass+1),nelem2d,MPI_REAL,root, &
                         MPI_COMM_WORLD,ierr)
              ELSE
                CALL LGARRAY_BCASTR(enspm_smt(:,:,ipass+1),nelem2d,max_elem_send, &
                            myproc,root,MPI_COMM_WORLD,ierr)
              END IF
            END IF
          END IF
        END DO !k1
        IF( myproc == root ) THEN
        CALL CPU_TIME(cput4)
          WRITE(*,'(//a,i4,a,f10.2,a)') ' ** Pass ',ipass,': ',(cput4-cput3),' seconds **'
        END IF
    END DO  ! ipass


    !---------------------
    ! Calculate LPM and PM
    !------------------------
    
    IF( myproc == root ) THEN
        CALL CPU_TIME(cput5)
        WRITE(6,'(//a)') ' Calculating LPM for raw ensemble'
    END IF
    
    IF( nprocs > 1 ) THEN
        nelem3d=nx*ny*nmembers
        nelem2d=nx*ny

        IF( nelem3d < max_elem_send) THEN
          CALL MPI_BCAST(ensfcst(:,:,:),nelem3d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
        ELSE IF( nelem2d < max_elem_send) THEN
          DO k=1,nmembers
            CALL MPI_BCAST(ensfcst(:,:,k),nelem2d,MPI_REAL,root, &
                       MPI_COMM_WORLD,ierr)
          END DO
        ELSE
          DO k=1,nmembers
            CALL LGARRAY_BCASTR(ensfcst(:,:,k),nelem2d,max_elem_send, &
                            myproc,root,MPI_COMM_WORLD,ierr)
          END DO
        END IF
        
        IF( nelem2d < max_elem_send) THEN
          CALL MPI_BCAST(ensmean(:,:),nelem2d,MPI_REAL,root, &
                   MPI_COMM_WORLD,ierr)
          CALL MPI_BCAST(ensmax(:,:),nelem2d,MPI_REAL,root, &
                   MPI_COMM_WORLD,ierr)
        ELSE
          CALL LGARRAY_BCASTR(ensmean(:,:),nelem2d,max_elem_send, &
                      myproc,root,MPI_COMM_WORLD,ierr)
          CALL LGARRAY_BCASTR(ensmax(:,:),nelem2d,max_elem_send, &
                      myproc,root,MPI_COMM_WORLD,ierr)
        END IF
    END IF

    
    CALL lpm_mean_la(myproc,nprocs,nx,ny,ibgn,iend,jbgn,jend,      &
            nmembers,patch_nx,patch_ny,ovx,ovy,     &
            filt_min,float(gauss_sigma),gs_weight, &
            ensfcst(:,:,1:nmembers),ensmean,ensmax, &
            enslpm,myproc,nprocs)
    
    IF( myproc == root ) THEN
        CALL CPU_TIME(cput6)
        WRITE(*,'(a,f10.2,a)') ' Raw ensemble LPM ',(cput6-cput5),' seconds'
        CALL CPU_TIME(cput5)
    END IF

    IF( myproc == root ) THEN
        WRITE(6,'(//a)') ' Calculating LPM for aligned ensemble'
    END IF
    
    IF( nprocs > 1 ) THEN
        nelem3d=nx*ny*nmembers
        nelem2d=nx*ny
        DO ipass=1,nshfpass
          IF( nelem3d < max_elem_send) THEN
            CALL MPI_BCAST(ensfcst_shf(:,:,:,ipass),nelem3d,MPI_REAL,root,MPI_COMM_WORLD,ierr)
          ELSE IF( nelem2d < max_elem_send) THEN
            DO k=1,nmembers
              CALL MPI_BCAST(ensfcst_shf(:,:,k,ipass),nelem2d,MPI_REAL,root, &
                         MPI_COMM_WORLD,ierr)
            END DO
          ELSE
            DO k=1,nmembers
              CALL LGARRAY_BCASTR(ensfcst_shf(:,:,k,ipass),nelem2d,max_elem_send, &
                              myproc,root,MPI_COMM_WORLD,ierr)
            END DO
          END IF
        
          IF( nelem2d < max_elem_send) THEN
            CALL MPI_BCAST(ensshfmn(:,:,ipass),nelem2d,MPI_REAL,root, &
                     MPI_COMM_WORLD,ierr)
            CALL MPI_BCAST(ensshfmax(:,:,ipass),nelem2d,MPI_REAL,root, &
                     MPI_COMM_WORLD,ierr)
          ELSE
            CALL LGARRAY_BCASTR(ensshfmn(:,:,ipass),nelem2d,max_elem_send, &
                        myproc,root,MPI_COMM_WORLD,ierr)
            CALL LGARRAY_BCASTR(ensshfmax(:,:,ipass),nelem2d,max_elem_send, &
                        myproc,root,MPI_COMM_WORLD,ierr)
          END IF
        END DO
    END IF
    
    DO ipass=1,nshfpass
        IF( myproc == root ) THEN
            WRITE(6,'(a,i4)') ' Pass number:',ipass
        END IF
        CALL lpm_mean_la(myproc,nprocs,nx,ny,ibgn,iend,jbgn,jend, &
                nmembers,patch_nx,patch_ny,ovx,ovy,  &
                filt_min,float(gauss_sigma),gs_weight,        &
                ensfcst_shf(:,:,1:nmembers,ipass),ensshfmn(:,:,ipass),ensshfmax(:,:,ipass), &
                ensshflpm(:,:,ipass),myproc,nprocs)
        IF( myproc == root .AND. ipass==1) THEN
            WRITE(6,'(a)') ' '
        END IF
    END DO

    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    
    IF( myproc == root ) THEN
        CALL CPU_TIME(cput6)
        WRITE(*,'(a,f10.2,a)') ' Aligned ensemble LPM ',(cput6-cput5),' seconds'
        WRITE(*, '(//a)') ' ----------------------------------------'
        WRITE(*, '(a)') ' --      Sending back to python       -- '
        WRITE(*, '(a)') ' ----------------------------------------'
    END IF

    
RETURN
END SUBROUTINE ensalign
