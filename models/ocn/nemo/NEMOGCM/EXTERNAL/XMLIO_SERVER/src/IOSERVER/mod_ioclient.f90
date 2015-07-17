MODULE mod_ioclient


CONTAINS

  SUBROUTINE init_ioclient(NEW_COMM)
  USE mod_ioclient_para
  USE mod_mpi_buffer_client
  USE mod_wait
  USE mod_ioserver_namelist
  USE mod_event_client
  USE iomanager
  USE mpi_mod
  IMPLICIT NONE
    INTEGER, INTENT(INOUT),OPTIONAL :: NEW_COMM

    INTEGER :: Comm
    INTEGER :: nb_server
    INTEGER :: rank
    INTEGER :: ierr
    LOGICAL :: init
 
    CALL read_namelist

    IF (using_mpi) THEN
      IF (using_server) THEN
        CALL Init_parallel
        CALL Init_mpi_buffer
        CALL Init_wait
        IF (PRESENT(NEW_COMM)) THEN
          NEW_COMM=intracomm
        ENDIF
      ELSE
        CALL MPI_INITIALIZED(init,ierr)
        IF (init) THEN
          IF (.NOT. PRESENT(NEW_COMM)) THEN
             Comm=MPI_COMM_WORLD
          ELSE
            Comm=New_Comm
          ENDIF
        ELSE
          CALL MPI_INIT(ierr)
          Comm=MPI_COMM_WORLD
          
          IF (PRESENT(NEW_COMM)) THEN
            New_Comm=MPI_COMM_WORLD
          ENDIF
        ENDIF  
        CALL MPI_COMM_SIZE(Comm,nb_server,ierr)     
        CALL MPI_COMM_RANK(Comm,rank,ierr)
        CALL iom__init(1,nb_server,rank)
        CALL iom__set_current_rank(1)
      ENDIF
    ELSE
      CALL iom__init(1,1,0)
      CALL iom__set_current_rank(1)
    ENDIF
    
  END SUBROUTINE init_ioclient

END MODULE mod_ioclient
