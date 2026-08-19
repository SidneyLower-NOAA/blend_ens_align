SUBROUTINE LGARRAY_BCASTR( array, nelem, max_elem_send, &
                             myproc, srcproc, comm_channel, ierr)
!
! Broadcasts data to other processors from a large array using MPI_BCAST.
! Solves size limit problem of MPI_BCAST by sending the data in
! smaller chunks.
!
! array: REAL array to be broadcasted
! nelem: Number of total elements in array (can be nx*ny in calling program)
! max_elem_send: Limit to MPI_BCAST (Intel is about 500,000)
! myproc: Processor number of calling processor
! srcproc: Source processor number for broadcast
! comm_channel: MPI Communications channel, e.g., MPI_COMM_WORLD
! ierr: Output status
!
! Keith Brewster, CAPS/Univ of Oklahoma
! March, 2018
!
  use mpi
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nelem
  REAL, INTENT(INOUT) :: array(nelem)
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: myproc
  INTEGER, INTENT(IN) :: srcproc
  INTEGER, INTENT(IN) :: comm_channel
  INTEGER, INTENT(OUT) :: ierr
!
! Misc local variables
!
  INTEGER :: k,ncall,nelemcall,idxbgn,idxend,nelembcast

  ierr = 0
  ncall=nelem/max_elem_send
  IF( (ncall*max_elem_send) < nelem ) ncall=ncall+1
! IF( myproc == srcproc ) WRITE(6,'(a,i8,a,i8,a,i5)') &
!   ' LGARRAY_BCASTR: nelm=',nelem,' max elem:',max_elem_send,' ncall:',ncall 
  nelemcall=nelem/ncall
  IF( (nelemcall*ncall) < nelem ) nelemcall=nelemcall+1
! IF( myproc == srcproc ) THEN
!   WRITE(6,'(a,i8)') ' LGARRAY_BCASTR: nelem=',nelem
!   WRITE(6,'(a,i5,a,i8)') &
!   ' LGARRAY_BCASTR: ncall=',ncall,' nelem/call:',nelemcall 
! END IF
! WRITE(6,'(a,i5,a)') 'Processor ',myproc,' here 1.'
! FLUSH(6)
  DO k=1,ncall
    IF(k == 1) THEN
      idxbgn=1
    ELSE
      idxbgn=idxend+1
    END IF
    idxend=min((idxbgn+nelemcall-1),nelem)
    nelembcast=(idxend-idxbgn)+1
!   IF( myproc == srcproc ) THEN
!     WRITE(6,'(a,i5,a,i8,a,i8,a,i8)') &
!     '   k=',k,' idxbgn=',idxbgn,' idxend=',idxend,' nelembcast:',nelembcast
!     FLUSH(6)
!   END IF
    CALL MPI_BCAST(array(idxbgn),nelembcast,MPI_REAL,srcproc, &
                       comm_channel,ierr)
    !CALL MPI_BARRIER(comm_channel,ierr)
  END DO
  RETURN
END SUBROUTINE LGARRAY_BCASTR

SUBROUTINE LGARRAY_INT_BCASTR( iarray, nelem, max_elem_send, &
                             myproc, srcproc, comm_channel, ierr)

  use mpi
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nelem
  INTEGER, INTENT(INOUT) :: iarray(nelem)
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: myproc
  INTEGER, INTENT(IN) :: srcproc
  INTEGER, INTENT(IN) :: comm_channel
  INTEGER, INTENT(OUT) :: ierr
!
! Misc local variables
!
  INTEGER :: k,ncall,nelemcall,idxbgn,idxend,nelembcast

  ierr = 0
  ncall=nelem/max_elem_send
  IF( (ncall*max_elem_send) < nelem ) ncall=ncall+1
! IF( myproc == srcproc ) WRITE(6,'(a,i8,a,i8,a,i5)') &
!   ' LGARRAY_BCASTR: nelm=',nelem,' max elem:',max_elem_send,' ncall:',ncall
  nelemcall=nelem/ncall
  IF( (nelemcall*ncall) < nelem ) nelemcall=nelemcall+1
! IF( myproc == srcproc ) THEN
!   WRITE(6,'(a,i8)') ' LGARRAY_BCASTR: nelem=',nelem
!   WRITE(6,'(a,i5,a,i8)') &
!   ' LGARRAY_BCASTR: ncall=',ncall,' nelem/call:',nelemcall
! END IF
! WRITE(6,'(a,i5,a)') 'Processor ',myproc,' here 1.'
! FLUSH(6)
  DO k=1,ncall
    IF(k == 1) THEN
      idxbgn=1
    ELSE
      idxbgn=idxend+1
    END IF
    idxend=min((idxbgn+nelemcall-1),nelem)
    nelembcast=(idxend-idxbgn)+1
!   IF( myproc == srcproc ) THEN
!     WRITE(6,'(a,i5,a,i8,a,i8,a,i8)') &
!     '   k=',k,' idxbgn=',idxbgn,' idxend=',idxend,' nelembcast:',nelembcast
!     FLUSH(6)
!   END IF
    CALL MPI_BCAST(iarray(idxbgn),nelembcast,MPI_INTEGER,srcproc, &
                       comm_channel,ierr)
    !CALL MPI_BARRIER(comm_channel,ierr)
  END DO
  RETURN
END SUBROUTINE LGARRAY_INT_BCASTR

SUBROUTINE LGARRAY_REDUCER( array,recv_buf,nelem,max_elem_send, &
                            myproc,root,reduce_op,comm_channel,ierr)
!
! Reduces array of data to recv_buf for a large array using MPI_REDUCE.
! Solves size limit problem of MPI_REDUCE by sending the data in
! smaller chunks.
!
! array: REAL array to be reduced
! recv_buf: REAL array to hold reduced result in root processor
! nelem: Number of total elements in array (can be nx*ny in calling program)
! max_elem_send: Limit to MPI_BCAST (Intel is about 500,000)
! myproc: Processor number of calling processor
! root: Root processor number, where reduce result is collected
! reduce_op: Reduce operation (such as MPI_SUM)
! comm_channel: MPI Communications channel, e.g., MPI_COMM_WORLD
! ierr: Output status
!
! Keith Brewster, CAPS/Univ of Oklahoma
! March, 2018
!

  use mpi
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nelem
  REAL, INTENT(INOUT) :: array(nelem)
  REAL, INTENT(OUT) :: recv_buf(nelem)
  INTEGER, INTENT(IN) :: max_elem_send
  INTEGER, INTENT(IN) :: myproc
  INTEGER, INTENT(IN) :: root
  INTEGER, INTENT(IN) :: reduce_op
  INTEGER, INTENT(IN) :: comm_channel
  INTEGER, INTENT(OUT) :: ierr
!
! Misc local variables
!
  INTEGER :: k,ncall,nelemcall,idxbgn,idxend,nelemreduce
  
  ierr = 0
  ncall=nelem/max_elem_send
  IF( (ncall*max_elem_send) < nelem ) ncall=ncall+1
! IF( myproc == srcproc ) WRITE(6,'(a,i8,a,i8,a,i5)') &
!   ' LGARRAY_BCASTR: nelm=',nelem,' max elem:',max_elem_send,' ncall:',ncall 
  nelemcall=nelem/ncall
  IF( (nelemcall*ncall) < nelem ) nelemcall=nelemcall+1
! IF( myproc == srcproc ) THEN
!   WRITE(6,'(a,i8)') ' LGARRAY_BCASTR: nelem=',nelem
!   WRITE(6,'(a,i5,a,i8)') &
!   ' LGARRAY_BCASTR: ncall=',ncall,' nelem/call:',nelemcall 
! END IF
! WRITE(6,'(a,i5,a)') 'Processor ',myproc,' here 1.'
! FLUSH(6)
  DO k=1,ncall
    IF(k == 1) THEN
      idxbgn=1
    ELSE
      idxbgn=idxend+1
    END IF
    idxend=min((idxbgn+nelemcall-1),nelem)
    nelemreduce=(idxend-idxbgn)+1
!   IF( myproc == srcproc ) THEN
!     WRITE(6,'(a,i5,a,i8,a,i8,a,i8)') &
!     '   k=',k,' idxbgn=',idxbgn,' idxend=',idxend,' nelemreduce:',nelemreduce
!     FLUSH(6)
!   END IF
    CALL MPI_REDUCE(array(idxbgn),recv_buf(idxbgn),nelemreduce, &
                    MPI_REAL,reduce_op,root,comm_channel,ierr)
    !CALL MPI_BARRIER(comm_channel,ierr)
  END DO
  RETURN
END SUBROUTINE LGARRAY_REDUCER