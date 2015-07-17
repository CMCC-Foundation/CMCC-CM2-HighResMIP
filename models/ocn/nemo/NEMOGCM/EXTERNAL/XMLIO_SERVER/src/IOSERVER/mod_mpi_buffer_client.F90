MODULE mod_mpi_buffer_client
  USE mod_mpi_buffer_parameters
#if defined key_mpp_mpi
  USE mpi, ONLY : status_size=>MPI_STATUS_SIZE
#else
  INTEGER, PARAMETER :: status_size = 1
#endif

  
  INTEGER(KIND=8),POINTER :: MPI_Buffer(:)
  
  TYPE mpi_requests
    INTEGER :: request
    INTEGER :: status(status_size)
    INTEGER :: Pos
  END TYPE mpi_requests
  
  TYPE(mpi_requests),DIMENSION(max_request) :: pending_request
  
  INTEGER :: Buffer_pos
  INTEGER :: request_pos
  INTEGER :: buffer_begin
  INTEGER :: request_begin
  INTEGER :: nb_request_pending
  INTEGER,SAVE :: start_pos
  LOGICAL,SAVE :: ok_new_request
  LOGICAL,SAVE :: is_last_request

  INTEGER,SAVE      :: mpi_buffer_size
  INTEGER,SAVE      :: aggregated_request

CONTAINS

  SUBROUTINE Init_mpi_buffer
  USE mod_global_memory
  USE mod_pack
  USE mod_ioclient_para
  IMPLICIT NONE

  INTEGER  :: nb_client
  
    nb_client=iosize-1  
    mpi_buffer_size=global_mpi_buffer_size/nb_client
    
    CALL allocate_global_memory(mpi_buffer_size,MPI_Buffer)
    buffer_begin=1
    request_begin=1
    Buffer_pos=1
    nb_request_pending=0
    Request_pos=1
    ok_new_request=.TRUE.
    is_last_request=.FALSE.
    aggregated_request=0
    
    CALL set_pack_buffer(MPI_Buffer,buffer_begin)
    
  END SUBROUTINE Init_mpi_buffer
  

  SUBROUTINE create_request(request_id)
  USE mod_pack
  USE mod_ioclient_para
  USE mpitrace
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: request_id
    
    CALL VTb(VTprocess_event)
    IF (ok_new_request) THEN  
      Pending_request(Request_pos)%Pos = pack_pos
      ok_new_request=.FALSE.
    ENDIF
    start_pos=pack_pos
    pack_pos=pack_pos+1
!    PRINT *,"Pos in Buffer",Pending_request(Request_pos)%Pos,"pack_pos",pack_pos
    CALL pack_data(request_id)
  END SUBROUTINE create_request
  
  
  SUBROUTINE Finalize_request
  USE mod_pack
  USE mod_ioclient_para
  USE mpitrace
  use mod_wait
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: ierr
    INTEGER :: message_size
    INTEGER(KIND=8) :: request_size
    INTEGER :: buffer_free
    LOGICAL :: ok_out
    LOGICAL :: is_Buffer_full
    request_size=pack_pos-start_pos    
    pack_buffer(start_pos)=request_size
    message_size=pack_pos-Pending_request(Request_pos)%Pos
    aggregated_request=aggregated_request+1


!! ICI verifier que le buffer ne se recouvre pas ainsi que les requetes
    
    ok_out=.FALSE.
    is_buffer_full=.FALSE.
    DO WHILE (.NOT. ok_out)
      CALL check_request
    
      IF ( buffer_begin <= pack_pos) THEN
        Buffer_free=mpi_buffer_size-pack_pos+1
      ELSE
        Buffer_free=buffer_begin-pack_pos
      ENDIF
      
!      Print *,"message_size",message_size,"buffer_free",buffer_free
!      PRINT *,"Request_pos",request_pos
!      PRINT *,"Pos in Buffer",Pending_request(Request_pos)%Pos,"pack_pos",pack_pos
      IF ( nb_request_pending==1 .AND. ( (buffer_free < MPI_buffer_size * 0.4) .OR. is_last_request .OR.  &
                                         (aggregated_request>=max_request/2) )) THEN
        ok_out=.FALSE.
        CALL Wait_us(10)
        IF (.NOT. is_buffer_full) THEN
          CALL VTb(VTbuffer_full)
        ENDIF      
        is_Buffer_full=.TRUE.
        
      ELSE
        ok_out=.TRUE.
        IF (is_buffer_full) THEN
          CALL VTe(VTbuffer_full)
        ENDIF
        is_buffer_full=.FALSE.
      
      ENDIF
      
   ENDDO
   
   IF (nb_request_pending==0 .OR. (buffer_free < MPI_buffer_size* 0.4 ) .OR. (aggregated_request> max_request/2) ) THEN
      
     CALL MPI_ISSEND(MPI_Buffer(Pending_request(Request_pos)%Pos),message_size,MPI_INTEGER8,     &
                    server_rank,tag_iocomm,iocomm,Pending_request(Request_pos)%request,ierr )
    
!    PRINT *,"Requete envoy�e !!!!" 
!    PRINT *,"Message : ",MPI_Buffer(Pending_request(Request_pos)%Pos:Pending_request(Request_pos)%Pos+message_size-1)              
      IF ( Pack_Pos > MPI_buffer_size*0.6 ) THEN
        Pack_Pos=1
      ENDIF
    
      IF (Request_Pos==max_request) THEN
        Request_Pos=1
      ELSE 
        Request_Pos=Request_Pos+1
      ENDIF
      nb_request_pending=nb_request_pending+1
    
      ok_new_request=.TRUE.
      aggregated_request=0
    ENDIF
    CALL VTe(VTprocess_event)    
  END SUBROUTINE Finalize_request


  SUBROUTINE Check_request
  USE mpi_mod
  IMPLICIT NONE
  LOGICAL :: ok_out
  LOGICAL :: OK_complete
  INTEGER :: ierr
  
!    PRINT *, 'on entre dans Check_request'
!    PRINT *, 'nb_request_pending',nb_request_pending
    
    IF (nb_request_pending>0) THEN
      ok_out=.FALSE.
    ELSE
      ok_out=.TRUE.
    ENDIF
    
    DO WHILE (.NOT. ok_out)
    
!      PRINT *,'Testing_request...'
!      PRINT *,'request_begin',request_begin
      CALL MPI_TEST(Pending_request(request_begin)%request,ok_complete,Pending_request(request_begin)%status,ierr)
!      PRINT *,'Request has been tested...'
      IF (ok_complete) THEN
!        PRINT *,"Request_completed"
        IF (Request_begin==max_request) THEN
          Request_begin=1
        ELSE
          request_begin=request_begin+1
        ENDIF
        
        buffer_begin=Pending_request(request_begin)%Pos
        
        nb_request_pending=nb_request_pending-1
        
        IF (nb_request_pending==0) THEN 
          ok_out=.TRUE.
        ELSE
          ok_out=.FALSE.
        ENDIF
      ELSE
        ok_out=.TRUE.
      ENDIF
    
    ENDDO
!    PRINT *, 'on sort de Check_request'  
  END SUBROUTINE Check_Request

END MODULE mod_mpi_buffer_client
