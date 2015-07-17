  MODULE mod_ioclient_para
    INTEGER,PARAMETER :: color_client=1
    INTEGER,PARAMETER :: color_server=2 
    INTEGER,SAVE      :: iocomm
    INTEGER,SAVE      :: iosize
    INTEGER,SAVE      :: iorank
    INTEGER,SAVE      :: server_rank
    INTEGER,SAVE      :: intracomm 
  CONTAINS
  
  
  SUBROUTINE Init_parallel
  USE mpitrace
  USE mod_ioserver_namelist
#if defined key_oasis3 || defined key_oasis4
  USE mod_prism_get_comm  
#endif
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: NEW_COMM
    INTEGER :: ierr
    INTEGER :: global_rank
    INTEGER :: global_size
    INTEGER :: mpi_rank
    INTEGER :: mpi_size
    INTEGER :: nb_server_io  
    INTEGER,ALLOCATABLE :: proc_color(:) 
    INTEGER :: i
    INTEGER :: div,remain
    INTEGER :: group_color
    INTEGER :: Comm_client_server
    CHARACTER(LEN=6) :: oasis_server_id
    
    IF (using_oasis) THEN
      oasis_server_id=server_id
      PRINT *,'prism_get_intracomm'
#if defined key_oasis3 || defined key_oasis4
      CALL prism_get_intracomm(Comm_client_server,oasis_server_id,ierr)
#endif
    ELSE
      CALL MPI_INIT(ierr)
      Comm_client_server=MPI_COMM_WORLD
    ENDIF

    CALL MPI_COMM_RANK(Comm_client_server,global_rank,ierr)
    CALL MPI_COMM_SIZE(Comm_client_server,global_size,ierr)

    CALL MPI_COMM_SPLIT(Comm_client_server,color_client,global_rank,intracomm,ierr)
    CALL MPI_COMM_SIZE(intracomm,mpi_size,ierr)
    CALL MPI_COMM_RANK(intracomm,mpi_rank,ierr)

    nb_server_io=global_size-mpi_size
    div=mpi_size/nb_server_io
    remain=MOD(mpi_size,nb_server_io)
    
    IF (mpi_rank<remain*(div+1)) THEN
      group_color=mpi_rank/(div+1)
    ELSE
      group_color=(nb_server_io-1)-(mpi_size-1-mpi_rank)/div
    ENDIF

    CALL MPI_COMM_SPLIT(Comm_client_server,group_color,global_rank,iocomm,ierr)
    
    CALL MPI_COMM_SIZE(iocomm,iosize,ierr)
    CALL MPI_COMM_RANK(iocomm,iorank,ierr)

    ALLOCATE(proc_color(0:iosize-1))
    CALL MPI_ALLGATHER(color_client,1,MPI_INTEGER,proc_color,1,MPI_INTEGER,iocomm,ierr)
    
    DO i=0,iosize-1
      IF (proc_color(i)==color_server) THEN
        server_rank=i
        EXIT
      ENDIF
    ENDDO
    
    PRINT *,"Proces No",mpi_rank,"--> server",server_rank
  END SUBROUTINE Init_parallel
  
  SUBROUTINE Finalize_parallel
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: ierr
    
    CALL MPI_FINALIZE(ierr)

  END SUBROUTINE Finalize_parallel
  
  END MODULE mod_ioclient_para
