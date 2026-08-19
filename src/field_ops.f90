!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE SETDXDY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE setdxdy(nx,ny,ibeg,iend,jbeg,jend,x1d,y1d,dxfld,dyfld,rdxfld,rdyfld)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the local delta-x, delta-y and their inverses.
!    Precalculating these variables speeds up later calculations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!
!    x1d     Array of x-coordinate grid locations (m)
!    y1d     Array of y-coordinate grid locations (m)
!
!  OUTPUT:
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, intent(in) :: nx
  INTEGER, intent(in) :: ny
  INTEGER, intent(in) :: ibeg
  INTEGER, intent(in) :: iend
  INTEGER, intent(in) :: jbeg
  INTEGER, intent(in) :: jend
  REAL, intent(in) :: x1d(nx)
  REAL, intent(in) :: y1d(ny)
  REAL, intent(out) :: dxfld(nx)
  REAL, intent(out) :: dyfld(ny)
  REAL, intent(out) :: rdxfld(nx)
  REAL, intent(out) :: rdyfld(ny)

!  Misc. local variables
  INTEGER :: i,j,istop,jstop

  istop=MIN((iend-1),(nx-1))
  DO i=ibeg,istop
    dxfld(i)=(x1d(i+1)-x1d(i))
    rdxfld(i)=1./(x1d(i+1)-x1d(i))
  END DO
  jstop=MIN((jend-1),(ny-1))
  DO j=jbeg,jstop
    dyfld(j)=(y1d(j+1)-y1d(j))
    rdyfld(j)=1./(y1d(j+1)-y1d(j))
  END DO
  RETURN
END SUBROUTINE setdxdy

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE MOVEGR                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE movegr(nx,ny, var,wrk, xshf,yshf, varout,                    &
           ibgn,iend,jbgn,jend,                                         &
           dxfld,dyfld,rdxfld,rdyfld,                                   &
           slopey,alphay,betay)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Using the shift vectors, xshf and yshf, translate the variables
!  in array var horizontally.
!
!
!-----------------------------------------------------------------------
!
!  INPUT :
!    nx,ny,nz Array dimensions for forecast field.
!
!  OUTPUT :
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny
  REAL, INTENT(IN) :: var(nx,ny)
  REAL, INTENT(OUT) :: wrk(nx,ny)
  REAL, INTENT(IN) :: xshf(nx,ny)
  REAL, INTENT(IN) :: yshf(nx,ny)
  REAL, INTENT(OUT) :: varout(nx,ny)
  INTEGER, INTENT(IN) :: ibgn,iend
  INTEGER, INTENT(IN) :: jbgn,jend
  REAL, INTENT(IN) :: dxfld(nx)
  REAL, INTENT(IN) :: dyfld(ny)
  REAL, INTENT(IN) :: rdxfld(nx)
  REAL, INTENT(IN) :: rdyfld(ny)
  REAL, INTENT(OUT) :: slopey(nx,ny)
  REAL, INTENT(OUT) :: alphay(nx,ny)
  REAL, INTENT(OUT) :: betay(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,idev,jdev,ii,jj
  REAL :: xdev,ydev,delx,dely
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Compute y-derivative terms
!
!-----------------------------------------------------------------------
!
  CALL setdrvy(nx,ny,1,                                                 &
               ibgn,iend,jbgn,jend,1,1,                                 &
               dyfld,rdyfld,var,                                        &
               slopey,alphay,betay)

  wrk(:,:) = var(:,:)
  DO j=jbgn,jend
    DO i=ibgn,iend
      idev=nint(xshf(i,j)-0.5)
      jdev=nint(yshf(i,j)-0.5)
      xdev=xshf(i,j)-FLOAT(idev)
      ydev=yshf(i,j)-FLOAT(jdev)
!    print *, '  i,j = ',i,j
!    print *, '  xshift,idev,xdev= ',u(i,j),idev,xdev
      jj=j+jdev
      ii=i+idev
!
!  Periodic x boundary conditions
!
!    IF(ii.ge.nx) ii=ii-ixlen
!    IF(ii.lt.1)  ii=ii+ixlen
!
!  Mirror boundary conditions   Mirror at j=3 and j=ny-2
!
!    IF(jj.ge.ny) THEN
!      jj=twnym2-jj
!      ydev=1.-ydev
!    ELSE IF(jj.lt.1) THEN
!      jj=5-jj
!      ydev=1.-ydev
!    END IF
!
!  For now assume zero gradiant boundaries
!
      ii = min(ii,iend-1)
      ii = max(ii,ibgn)
      jj = min(jj,jend)
      jj = max(jj,jbgn)

!    c1=xdev
!    c2=ydev
!    c3=1.-xdev
!    c4=1.-ydev
!    wrk(i,j)=
!    +        c3*(c4*var(  ii,jj)+c2*var(  ii,jj+1))+
!    +        c1*(c4*var(ii+1,jj)+c2*var(ii+1,jj+1))
!
      delx=xdev*dxfld(ii)
      dely=ydev*dyfld(jj)
!     print *, ' i,j,ii,jj:',i,j,ii,jj
      wrk(i,j)=(1.-delx*rdxfld(ii))*                                    &
               (var(ii  ,jj)+slopey(ii  ,jj)*dely)+                     &
               (delx*rdxfld(ii))*                                       &
               (var(ii+1,jj)+slopey(ii+1,jj)*dely)

    END DO
  END DO
!
!  Transfer shifted and original array into output array
!
  varout(:,:) = wrk(:,:)
  RETURN
END SUBROUTINE movegr
!

SUBROUTINE smooth2d(nx,ny,ibgn,iend,jbgn,jend,s,zin,zwork,zout)
!
!  Performs symmetrical two-dimensional smoothing of input field
!  zin which is output as zout, the smoothed field.  A work array
!  zwork, dimension (nx,ny) is required.
!
!  K. Brewster, October, 1991
!
  IMPLICIT NONE
!
!  Arguments
!
  INTEGER, INTENT(IN) :: nx,ny
  INTEGER, INTENT(IN) :: ibgn,iend,jbgn,jend
  REAL, INTENT(IN)    :: s
  REAL, INTENT(IN)    :: zin(nx,ny)
  REAL, INTENT(OUT)   :: zwork(nx,ny)
  REAL, INTENT(OUT)   :: zout(nx,ny)
!
!  Misc internal variables
!
  INTEGER :: i,j
  REAL :: wcen,wsid
!
  zwork(:,:)=zin(:,:)
  wcen=1.-s
  wsid=s*0.5
  DO j=jbgn,jend
    DO i=ibgn+1,iend-1
      zwork(i,j)=zin(i  ,j)*wcen +                                      &
                 zin(i+1,j)*wsid +                                      &
                 zin(i-1,j)*wsid
    END DO
  END DO
!
!
!
  zout(:,:)=zwork(:,:)
  DO j=jbgn+1, jend-1
    DO i=ibgn, iend
      zout(i,j)=zwork(i  ,j)*wcen +                                     &
                zwork(i,j+1)*wsid +                                     &
                zwork(i,j-1)*wsid
    END DO
  END DO
  RETURN
END SUBROUTINE smooth2d


SUBROUTINE vmaxmin(varray,nx,ny,ibgn,iend,jbgn,jend,vmin,vmax)
  IMPLICIT NONE
  INTEGER :: nx,ny
  REAL :: varray(nx,ny)
  INTEGER :: ibgn,iend
  INTEGER :: jbgn,jend
  REAL :: vmin,vmax

  INTEGER :: i,j

  vmin=varray(ibgn,jbgn)
  vmax=varray(ibgn,jbgn)
  DO j=jbgn,jend
    DO i=ibgn,iend
      vmin=min(varray(i,j),vmin) 
      vmax=max(varray(i,j),vmax) 
!     IF(varray(i,j) > 1000.) THEN
!       WRITE(6,'(a,f12.2,a,2i5)') ' Huge value of: ',varray(i,j),' at i,j: ',i,j
!     END IF
    END DO
  END DO
  RETURN
END SUBROUTINE vmaxmin

SUBROUTINE a2dmax0(array,idim0,idim1,ibgn,iend,jdim0,jdim1,jbgn,jend,vmin,vmax)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: idim0,idim1,ibgn,iend
  INTEGER, INTENT(IN) :: jdim0,jdim1,jbgn,jend
  REAL, INTENT(IN) :: array(idim0:idim1,jdim0:jdim1)
  REAL, INTENT(OUT) :: vmin
  REAL, INTENT(OUT) :: vmax

  INTEGER :: i,j

  vmin=1.0E30
  vmax=-1.0E30
  DO j=jbgn,jend
    DO i=ibgn,iend
      vmin=MIN(vmin,array(i,j))
      vmax=MAX(vmax,array(i,j))
    END DO
  END DO
  RETURN
END SUBROUTINE a2dmax0


SUBROUTINE setdrvy(nx,ny,nz,                                            &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           dyfld,rdyfld,var,                                            &
           slopey,alphay,betay)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Calculate the coefficients of interpolating polynomials
!    in the y-direction.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, November, 1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of model grid points in the x-direction (east/west)
!    ny       Number of model grid points in the y-direction (north/south)
!    nz       Number of model grid points in the vertical
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    var      variable to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
  REAL :: dyfld(ny)
  REAL :: rdyfld(ny)
  REAL :: var(nx,ny,nz)
  REAL :: slopey(nx,ny,nz)
  REAL :: alphay(nx,ny,nz)
  REAL :: betay(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: jstart,jstop
  REAL :: rtwody
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  jstart=MAX(jbeg,2)
  jstop=MIN((jend-1),(ny-2))
  DO k=kbeg,kend
    DO j=jstart,jstop
      DO i=ibeg,iend
        slopey(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j)
        rtwody=1./(dyfld(j-1)+dyfld(j))
        alphay(i,j,k)=((var(i,j+1,k)-var(i,j,k))*rdyfld(j) +            &
                 (var(i,j-1,k)-var(i,j,k))*rdyfld(j-1))*rtwody
        betay(i,j,k)=(var(i,j+1,k)-var(i,j,k))*rdyfld(j) -              &
                   dyfld(j)*alphay(i,j,k)
      END DO
    END DO
  END DO
  RETURN
END SUBROUTINE setdrvy
