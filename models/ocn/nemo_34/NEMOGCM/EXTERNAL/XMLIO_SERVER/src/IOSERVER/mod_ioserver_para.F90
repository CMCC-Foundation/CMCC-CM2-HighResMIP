  MODULE mod_ioserver_para
    INTEGER,PARAMETER :: color_client=1
    INTEGER,PARAMETER :: color_server=2 

    INTEGER,SAVE      :: iocomm
    INTEGER,SAVE      :: iosize
    INTEGER,SAVE      :: iorank
    
    INTEGER,SAVE      :: mpi_rank
    INTEGER,SAVE      :: mpi_size
    INTEGER,SAVE      :: mpi_master
    LOGICAL,SAVE      :: is_mpi_master
  
    INTEGER,SAVE      :: nb_client
    INTEGER,ALLOCATABLE,SAVE      :: client_rank(:)
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
    INTEGER :: ierr
    INTEGER :: global_rank
    INTEGER :: global_size
    INTEGER,ALLOCATABLE :: proc_color(:)
    INTEGER :: i
    INTEGER :: group_color
    INTEGER :: Comm_client_server
    INTEGER :: comp_id
    CHARACTER(LEN=6) :: oasis_server_id, oasis_client_id

    IF (using_oasis) THEN
       oasis_server_id=server_id
       oasis_client_id=client_id
#if defined key_oasis3 || defined key_oasis4
       CALL prism_init_comp_proto (comp_id, oasis_server_id, ierr)
       CALL prism_get_intracomm(Comm_client_server,oasis_client_id,ierr)
#endif
    ELSE
      CALL MPI_INIT(ierr)
      Comm_client_server=MPI_COMM_WORLD
    ENDIF
        
    CALL MPI_COMM_RANK(Comm_client_server,global_rank,ierr)
    CALL MPI_COMM_SIZE(Comm_client_server,global_size,ierr)
       
    CALL MPI_COMM_SPLIT(Comm_client_server,color_server,global_rank,intracomm,ierr)
    CALL MPI_COMM_SIZE(intracomm,mpi_size,ierr)
    CALL MPI_COMM_RANK(intracomm,mpi_rank,ierr)
    
    group_color=mpi_rank
    PRINT *,'group_color',group_color

    CALL MPI_COMM_SPLIT(Comm_client_server,group_color,global_rank,iocomm,ierr)    
    
    CALL MPI_COMM_SIZE(iocomm,iosize,ierr)
    CALL MPI_COMM_RANK(iocomm,iorank,ierr)

    PRINT *,"io_size-> ",iosize,"iorank-> ",iorank
    
    ALLOCATE(proc_color(0:iosize-1))
    CALL MPI_ALLGATHER(color_server,1,MPI_INTEGER,proc_color,1,MPI_INTEGER,iocomm,ierr)
    print *,"proc_color -> ",proc_color
    
    ALLOCATE(client_rank(iosize-1))
    nb_client=0
    DO i=0,iosize-1
      IF (proc_color(i)==color_client) THEN
        nb_client=nb_client+1
        client_rank(nb_client)=i
      ENDIF
    ENDDO

    PRINT *,"Proces No",mpi_rank,"--> client ",client_rank
  END SUBROUTINE Init_parallel
  
  SUBROUTINE Finalize_parallel
  USE mod_ioserver_namelist
#if defined key_oasis3 || defined key_oasis4
  USE mod_prism_proto
#endif
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: ierr
    
    IF (using_oasis) THEN
#if defined key_oasis3 || defined key_oasis4
      CALL prism_terminate_proto(ierr)
#endif
    ELSE
      CALL MPI_FINALIZE(ierr)
    ENDIF
    
  END SUBROUTINE Finalize_parallel

END MODULE mod_ioserver_para
