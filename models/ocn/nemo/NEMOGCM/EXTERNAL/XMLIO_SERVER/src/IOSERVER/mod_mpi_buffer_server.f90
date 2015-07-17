MODULE mod_mpi_buffer_server
  USE mod_mpi_buffer_parameters

  TYPE mpi_requests
    INTEGER :: Pos
    INTEGER :: size
  END TYPE mpi_requests
  
  TYPE buffer
    INTEGER(KIND=8),POINTER  :: mpi_buffer(:)
    INTEGER :: pos
    INTEGER :: begin
    TYPE(mpi_requests),POINTER :: request(:)
    INTEGER :: request_pos
    INTEGER :: nb_request
    INTEGER :: request_begin
    INTEGER :: request_free
    LOGICAL :: is_posted
    LOGICAL :: is_terminated
    INTEGER :: last_request
  END TYPE buffer

  TYPE(buffer),ALLOCATABLE,SAVE :: buffers(:)
  
  INTEGER,SAVE :: current_rank
  INTEGER,SAVE :: mpi_buffer_size
  
CONTAINS
  
  SUBROUTINE Init_mpi_buffer
  USE mod_ioserver_para
  USE mod_global_memory
  IMPLICIT NONE
    INTEGER :: n
    
    ALLOCATE(buffers(nb_client))
    mpi_buffer_size=global_mpi_buffer_size/nb_client 
    DO n=1,nb_client
      CALL allocate_global_memory(mpi_buffer_size,buffers(n)%MPI_Buffer)
      buffers(n)%pos=1
      buffers(n)%begin=1
      ALLOCATE(buffers(n)%request(max_request))
      buffers(n)%nb_request=0
      buffers(n)%request_free=max_request
      buffers(n)%request_pos=1
      buffers(n)%request_begin=1
      buffers(n)%is_posted=.FALSE.
      buffers(n)%is_terminated=.FALSE.
    ENDDO

  END SUBROUTINE Init_mpi_buffer

  SUBROUTINE Check_buffer
  USE mod_ioserver_para
  USE mpitrace
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: n
    INTEGER :: status(MPI_STATUS_SIZE)
    LOGICAL :: ok_complete
    LOGICAL :: ok_out
    INTEGER :: ierr
    INTEGER :: message_size
    INTEGER :: Buffer_space 
    INTEGER :: request
    
    CALL VTb(VTcheck_buffer)
    DO n=1,nb_client
      ok_out=.FALSE.
      
      DO WHILE (.NOT. ok_out)
        
        IF (buffers(n)%is_posted) THEN
!          PRINT*,"MPI_TEST"
!          PRINT *,"MPI_TEST LAST_REQUEST (avant)",buffers(n)%last_request,"ok_complete ?",ok_complete
          CALL MPI_TEST(buffers(n)%last_request,ok_complete,status,ierr)
!          PRINT *,"MPI_TEST LAST_REQUEST (apres)",buffers(n)%last_request,"ok_complete ?",ok_complete

          IF (ok_complete) THEN
            CALL MPI_GET_COUNT(status,MPI_INTEGER8,message_size,ierr)
            CALL Fill_request(n,Buffers(n)%pos,message_size)
            buffers(n)%is_posted=.FALSE.
            buffers(n)%pos=buffers(n)%pos+message_size
         
!            PRINT *,"buffer_pos",buffers(n)%pos
!            PRINT *,"Message recu de taille", message_size
!              PRINT *,"Message : ",buffers(n)%mpi_buffer(buffers(n)%pos:buffers(n)%pos+message_size-1)
!              PRINT *,"Nouvelle taille buffer_pos",buffers(n)%pos
              ok_out=.FALSE.
          ELSE
              ok_out=.TRUE.
          ENDIF
       
        ELSE  ! (.NOT. is_posted)        
!          PRINT *,"MPI_IPROBE"
     
          CALL MPI_IPROBE(client_rank(n),tag_iocomm,iocomm,ok_complete,status,ierr)
!          PRINT *,"Ok_complete.. ?",ok_complete
          IF (ok_complete) THEN
            CALL MPI_GET_COUNT(status,MPI_INTEGER8,message_size,ierr)
!            PRINT *,"message_size",message_size
            IF (buffers(n)%pos >= buffers(n)%begin) THEN
              IF (buffers(n)%pos+message_size>mpi_buffer_size) THEN
                buffer_space=buffers(n)%begin-1
              ELSE
                buffer_space=mpi_buffer_size-buffers(n)%begin+1
              ENDIF
            ELSE
              buffer_space=buffers(n)%begin-buffers(n)%pos
            ENDIF
            
            IF (buffer_space < message_size .OR. buffers(n)%request_free<=max_request/2) THEN
              ok_out=.TRUE.
!              PRINT *,"BUFFER FULL !!!!"
!              PRINT *,"buffer_space",buffer_space,"request_free",buffers(n)%request_free
!              PRINT *,"buffer_pos",buffers(n)%pos,"buffer_begin",buffers(n)%begin
            ELSE
              
              IF (buffers(n)%pos+message_size>mpi_buffer_size) buffers(n)%pos=1
          
!              PRINT *,'reception du message'
              CALL MPI_IRECV(buffers(n)%mpi_buffer(buffers(n)%pos),message_size,MPI_INTEGER8,          &
                          client_rank(n),tag_iocomm,iocomm,buffers(n)%last_request,ierr)
                          
              buffers(n)%is_posted=.TRUE.
!              PRINT *,"MPI_IRECV LAST_REQUEST",buffers(n)%last_request
              CALL MPI_TEST(buffers(n)%last_request,ok_complete,status,ierr)
!              PRINT *,"MPI_TEST LAST_REQUEST",buffers(n)%last_request,"ok_complete ?",ok_complete
              
              IF (ok_complete) THEN
                buffers(n)%is_posted=.FALSE.
                CALL Fill_request(n,Buffers(n)%pos,message_size)

                buffers(n)%pos=buffers(n)%pos+message_size
                 
!                PRINT *,"buffer_pos",buffers(n)%pos
!                PRINT *,"Message recu de taille", message_size
!                PRINT *,"Message : ",buffers(n)%mpi_buffer(buffers(n)%pos:buffers(n)%pos+message_size-1)
!                PRINT *,"Nouvelle taille buffer_pos",buffers(n)%pos
                ok_out=.FALSE.
              ELSE
                ok_out=.TRUE.
              ENDIF
         
            ENDIF
          ELSE
            ok_out=.TRUE.
          ENDIF
        ENDIF  ! (is_posted)
        ok_out=.TRUE.
      ENDDO
    ENDDO

    CALL VTe(VTcheck_buffer)

  END SUBROUTINE Check_buffer
  
  SUBROUTINE Fill_request(n,pos,message_size)
  USE mpi_mod
  IMPLICIT NONE
    INTEGER :: n
    INTEGER :: pos
    INTEGER :: message_size
    INTEGER :: current_size
    INTEGER :: current_pos
    INTEGER :: total_size
    LOGICAL :: ok_out
    INTEGER :: ierr
!    PRINT *,'Fill_request!!!'
!    PRINT *,'Buffer_pos',pos
!    PRINT *,"First integer",Buffers(n)%mpi_buffer(pos)
    total_size=0
    current_pos=pos
    
    ok_out=.FALSE.
    DO WHILE (.NOT. ok_out)
      
      IF (Buffers(n)%request_free==0) THEN
        PRINT *,'Plus de requete disponible !!!!'
        CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      ENDIF
      
      current_size=buffers(n)%mpi_buffer(current_pos)
      Buffers(n)%request(Buffers(n)%request_pos)%pos=current_pos+1
      Buffers(n)%request(Buffers(n)%request_pos)%size=current_size
      total_size=total_size+current_size
      
      IF (Buffers(n)%request_pos==max_request) THEN 
        Buffers(n)%request_pos=1
      ELSE
        Buffers(n)%request_pos=Buffers(n)%request_pos+1
      ENDIF
      Buffers(n)%request_free=Buffers(n)%request_free-1
      Buffers(n)%nb_request=Buffers(n)%nb_request+1
      
      IF (total_size==message_size) ok_out=.TRUE.
      IF (total_size>message_size) THEN
        PRINT *,"Probleme : la taille du message ne coincide pas avec la taille de l'enveloppe"
        CALL MPI_ABORT(MPI_COMM_WORLD,-1,ierr)
      ENDIF
      
      current_pos=current_pos+current_size
    ENDDO
  
  END SUBROUTINE Fill_request 
    
  
  SUBROUTINE Process_request(is_terminated)
  USE mod_event_server
  USE mod_ioserver_para
  USE mod_pack
  USE mpitrace
  USE mod_wait
  IMPLICIT NONE
    LOGICAL  :: is_terminated
    INTEGER  :: n
    
    is_terminated=.FALSE.
    
    IF (ALL(buffers(:)%nb_request>0)) THEN
      DO n=1,nb_client
!        PRINT *,"nb_request ?",buffers(n)%nb_request
        IF ( (buffers(n)%nb_request > 0) .AND. (.NOT. buffers(n)%is_terminated)) THEN
!          PRINT *,"request_pos",buffers(n)%request_begin
!          PRINT *,"buffer_pos",buffers(n)%request(buffers(n)%request_begin)%pos
          CALL set_pack_buffer(buffers(n)%mpi_buffer,buffers(n)%request(buffers(n)%request_begin)%pos)
          CALL VTb(VTprocess_event)
          current_rank=n
          CALL process_event(n,buffers(n)%is_terminated)
          CALL VTe(VTprocess_event)
          buffers(n)%nb_request=buffers(n)%nb_request-1
          buffers(n)%request_free=buffers(n)%request_free+1
        
       
          IF (Buffers(n)%request_begin==max_request) THEN
            Buffers(n)%request_begin=1
          ELSE
            buffers(n)%request_begin=Buffers(n)%request_begin+1
          ENDIF
          
          IF (buffers(n)%nb_request > 0) THEN 
            buffers(n)%begin=buffers(n)%request(buffers(n)%request_begin)%pos
          ELSE
            buffers(n)%begin=buffers(n)%pos
          ENDIF
        ENDIF
      ENDDO 
    
!      PRINT *,"buffers_is_terminated",buffers(1:nb_client)%is_terminated
!      PRINT *,"buffers_nb_request",buffers(1:nb_client)%nb_request
      IF (ALL(buffers(1:nb_client)%is_terminated)) THEN
        IF (ALL(buffers(1:nb_client)%nb_request==0)) THEN
          is_terminated=.TRUE.
        ENDIF
      ENDIF
    ELSE
      CALL Wait_us(5)
    ENDIF     
  END SUBROUTINE process_request  
  
    
END MODULE mod_mpi_buffer_server
