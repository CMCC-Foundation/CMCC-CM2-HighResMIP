PROGRAM server
  USE mod_ioserver_para
  USE mod_mpi_buffer_server
  USE mod_wait
  USE mod_global_memory
  USE ioipsl
  USE iomanager
  USE mod_interface_ioipsl
  USE mod_ioserver_namelist
  
  IMPLICIT NONE
#if defined key_mpp_mpi
  INCLUDE 'mpif.h'
  
  INTEGER :: i,j,ierr
  LOGICAL :: is_terminated
  
!  CALL SLEEP(60)
  PRINT *,'je suis un serveur'
  CALL read_namelist
  CALL init_parallel
  CALL init_mpi_buffer
  CALL init_wait
  CALL init_interface_ioipsl
  
  CALL iom__init(nb_client,mpi_size,mpi_rank)
  PRINT *,'je suis un serveur'
  

!  CALL MPI_BARRIER(iocomm,ierr)
  
  is_terminated=.FALSE.
  
  DO WHILE (.NOT. is_terminated)
    CALL Check_buffer
    CALL process_request(is_terminated)
!    CALL wait_us(10)
  ENDDO

!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL histclo
  
  CALL Finalize_parallel

  PRINT *,"THAT'S ALL FOLK"
#endif
  
END PROGRAM server
