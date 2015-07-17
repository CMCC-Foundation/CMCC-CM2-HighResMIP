MODULE mod_global_memory

  INTERFACE allocate_global_memory
    MODULE PROCEDURE Allocate_global_memory_r8,    &
                     Allocate_global_memory_i8,    &
!                     Allocate_global_memory_r4,    &
                     Allocate_global_memory_i4
  END INTERFACE Allocate_global_memory
  
CONTAINS

  SUBROUTINE Allocate_global_memory_r8(size,Pt)
  USE mpi_mod
  IMPLICIT NONE
    REAL(kind=8),POINTER :: Pt(:)
    INTEGER              :: size

    POINTER (Pbuffer,MPI_Buffer(size))
    REAL(kind=8) :: MPI_Buffer
    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
    INTEGER :: ierr
    
    BS=8*size
    CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
    CALL associate_buffer(MPI_Buffer,Pt)
  
  CONTAINS
  
    SUBROUTINE associate_buffer(MPI_buffer,Pt)
    IMPLICIT NONE
      REAL(kind=8),DIMENSION(:),target :: MPI_Buffer
      REAL(kind=8),POINTER             :: Pt(:)
      Pt=>MPI_buffer
    END SUBROUTINE associate_buffer
  
  END SUBROUTINE Allocate_global_memory_r8



  SUBROUTINE Allocate_global_memory_i8(size,Pt)
  USE mpi_mod
  IMPLICIT NONE
    INTEGER(kind=8),POINTER :: Pt(:)
    INTEGER              :: size

    POINTER (Pbuffer,MPI_Buffer(size))
    INTEGER(kind=8) :: MPI_Buffer
    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
    INTEGER :: ierr
    
    BS=8*size
    CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
    CALL associate_buffer(MPI_Buffer,Pt)
  
  CONTAINS
  
    SUBROUTINE associate_buffer(MPI_buffer,Pt)
    IMPLICIT NONE
      INTEGER(kind=8),DIMENSION(:),target :: MPI_Buffer
      INTEGER(kind=8),POINTER             :: Pt(:)
      Pt=>MPI_buffer
    END SUBROUTINE associate_buffer
  
  END SUBROUTINE Allocate_global_memory_i8
  
  
!  SUBROUTINE Allocate_global_memory_r4(size,Pt)
!  USE mpi_mod
!  IMPLICIT NONE
!    REAL(kind=4),POINTER :: Pt(:)
!    INTEGER              :: size
!
!    POINTER (Pbuffer,MPI_Buffer(size))
!    REAL(kind=4) :: MPI_Buffer
!    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
!    INTEGER :: ierr
!    
!    BS=4*size
!    CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
!    CALL associate_buffer(MPI_Buffer,Pt)
!  
!  CONTAINS
!  
!    SUBROUTINE associate_buffer(MPI_buffer,Pt)
!    IMPLICIT NONE
!      REAL(kind=4),DIMENSION(:),target :: MPI_Buffer
!      REAL(kind=4),POINTER             :: Pt(:)
!      Pt=>MPI_buffer
!    END SUBROUTINE associate_buffer
!  
!  END SUBROUTINE Allocate_global_memory_r4




  SUBROUTINE Allocate_global_memory_i4(size,Pt)
  USE mpi_mod
  IMPLICIT NONE
    INTEGER(kind=4),POINTER :: Pt(:)
    INTEGER              :: size

    POINTER (Pbuffer,MPI_Buffer(size))
    INTEGER(kind=4) :: MPI_Buffer
    INTEGER(KIND=MPI_ADDRESS_KIND) :: BS 
    INTEGER :: ierr
    
    BS=4*size
    CALL MPI_ALLOC_MEM(BS,MPI_INFO_NULL,Pbuffer,ierr)
    CALL associate_buffer(MPI_Buffer,Pt)
  
  CONTAINS
  
    SUBROUTINE associate_buffer(MPI_buffer,Pt)
    IMPLICIT NONE
      INTEGER(kind=4),DIMENSION(:),target :: MPI_Buffer
      INTEGER(kind=4),POINTER             :: Pt(:)
      Pt=>MPI_buffer
    END SUBROUTINE associate_buffer
  
  END SUBROUTINE Allocate_global_memory_i4
  
  
END MODULE mod_global_memory
