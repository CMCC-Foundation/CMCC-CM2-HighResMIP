MODULE mod_ioserver_namelist

  USE nc4interface
  LOGICAL,SAVE :: using_server
  LOGICAL,SAVE :: using_oasis
  CHARACTER(LEN=100),SAVE :: client_id
  CHARACTER(LEN=100),SAVE :: server_id
  NAMELIST /coupling_param/ using_server,using_oasis,client_id,server_id  

  INTEGER,SAVE :: global_mpi_buffer_size            
  NAMELIST /mpi_param/ global_mpi_buffer_size  

  INTEGER,PARAMETER,PRIVATE :: starting_unit=356
  TYPE(snc4_ctl), SAVE :: snc4ioset
  INTEGER :: nn_nchunks_i, nn_nchunks_j, nn_nchunks_k
  LOGICAL :: ln_nc4zip=.false.
  NAMELIST /namnc4/ nn_nchunks_i, nn_nchunks_j, nn_nchunks_k, ln_nc4zip
  
CONTAINS

  SUBROUTINE set_default_namelist
  IMPLICIT NONE
  
    using_server=.FALSE.
    using_oasis=.FALSE.
    client_id='unknown'
    server_id='unknown'
    global_mpi_buffer_size=512
    snc4ioset%luse = .false.
    
  END SUBROUTINE set_default_namelist
  
    
  SUBROUTINE read_namelist
  IMPLICIT NONE
    LOGICAL :: opened
    INTEGER :: unit
    INTEGER :: ierr

    CALL set_default_namelist   
    unit=starting_unit
    opened=.TRUE.
    DO WHILE (opened) 
      unit=unit+1
      INQUIRE(unit,OPENED=opened)
    ENDDO
      
    OPEN(UNIT=unit,FILE='xmlio_server.def',STATUS='old',IOSTAT=ierr)
    
    
    IF (ierr==0) THEN
      READ(unit,nml=coupling_param)
      READ(unit,nml=mpi_param)
      READ(unit,nml=namnc4,ERR=666,END=666)
 666  global_mpi_buffer_size=global_mpi_buffer_size*1024*128  
      CLOSE(unit)
      IF (ln_nc4zip) THEN
        snc4ioset%ni  = nn_nchunks_i
        snc4ioset%nj  = nn_nchunks_j
        snc4ioset%nk  = nn_nchunks_k
        snc4ioset%luse = ln_nc4zip
      ENDIF
    ELSE
      PRINT *,'WARNING : mod_ioserver::read_namelist : file xmlio_server.def is absent', &
              ' ---> using default parameter'
    ENDIF 
    
    
    PRINT *,'namelist read --> ',using_server,using_oasis,client_id,server_id
    
  END SUBROUTINE read_namelist


END MODULE mod_ioserver_namelist 
