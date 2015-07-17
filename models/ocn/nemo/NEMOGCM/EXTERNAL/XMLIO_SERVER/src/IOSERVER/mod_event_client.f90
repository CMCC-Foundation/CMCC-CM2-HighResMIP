MODULE mod_event_client
  USE mod_pack, ONLY : pack_data, pack_field
  USE mod_mpi_buffer_client, ONLY : create_request, finalize_request,is_last_request
  USE mod_event_parameters 
  USE mod_ioserver_namelist
      
CONTAINS

  SUBROUTINE event__swap_context(id)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: id
    
    IF (using_server) THEN
      CALL create_request(event_id_swap_context)
      CALL pack_data(LEN(TRIM(id)))
      CALL pack_data(TRIM(id))
      CALL Finalize_request
    ELSE
      CALL iom__swap_context(TRIM(id))
    ENDIF
    
  END SUBROUTINE  event__swap_context
  
  
  SUBROUTINE event__parse_xml_file(filename)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: filename
    
    IF (using_server) THEN
      CALL create_request(event_id_parse_xml_file)
      CALL pack_data(LEN(TRIM(filename)))
      CALL pack_data(TRIM(filename))
      CALL Finalize_request
    ELSE
      CALL iom__parse_xml_file(filename)
    ENDIF
    
  END SUBROUTINE  event__parse_xml_file

  SUBROUTINE event__set_vert_axis(vert_name,vert_value)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: vert_name
    REAL,INTENT(IN)             :: vert_value(:)
    
    
    IF (using_server) THEN
      CALL create_request(event_id_set_vert_axis)
      CALL pack_data(LEN(TRIM(vert_name)))
      CALL pack_data(TRIM(vert_name))
      CALL pack_data(SIZE(vert_value))
      CALL pack_data(vert_value)
      CALL Finalize_request
    ELSE
      CALL iom__set_vert_axis(vert_name,vert_value)
    ENDIF
    
  END SUBROUTINE  event__set_vert_axis

  SUBROUTINE event__set_grid_dimension(name,ni_glo,nj_glo)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,INTENT(IN) :: ni_glo
    INTEGER,INTENT(IN) :: nj_glo
    
    IF (using_server) THEN
      CALL create_request(event_id_set_grid_dimension)
      CALL pack_data(LEN(TRIM(name)))
      CALL pack_data(TRIM(name))
      CALL pack_data(ni_glo)
      CALL pack_data(nj_glo)
      CALL Finalize_request
    ELSE
      CALL iom__set_grid_dimension(name,ni_glo,nj_glo)
    ENDIF
    
  END SUBROUTINE event__set_grid_dimension

  
  SUBROUTINE event__set_grid_domain(name,ni,nj,ibegin,jbegin,lon,lat)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,INTENT(IN) :: ni
    INTEGER,INTENT(IN) :: nj
    INTEGER,INTENT(IN) :: ibegin
    INTEGER,INTENT(IN) :: jbegin
    REAL,INTENT(IN)    :: lon(ni,nj)
    REAL,INTENT(IN)    :: lat(ni,nj)
    
    IF (using_server) THEN
      CALL create_request(event_id_set_grid_domain)
      CALL pack_data(LEN(TRIM(name)))
      CALL pack_data(TRIM(name))
      CALL pack_data(ni)
      CALL pack_data(nj)
      CALL pack_data(ibegin)
      CALL pack_data(jbegin)
      CALL pack_data(lon)
      CALL pack_data(lat)
      CALL Finalize_request
    ELSE
     CALL iom__set_grid_domain(name,ni,nj,ibegin,jbegin,lon,lat)
    ENDIF
      
  END SUBROUTINE event__set_grid_domain

  SUBROUTINE event__set_grid_type_nemo(name)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    
    IF (using_server) THEN
      CALL create_request(event_id_set_grid_type_nemo)
      CALL pack_data(LEN(TRIM(name)))
      CALL pack_data(TRIM(name))
      CALL Finalize_request
    ELSE
     CALL iom__set_grid_type_nemo(name)
    ENDIF
      
  END SUBROUTINE event__set_grid_type_nemo
  
  SUBROUTINE event__set_grid_type_lmdz(name,nbp,offset)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,INTENT(IN)          :: nbp
    INTEGER,INTENT(IN)          :: offset
        
    IF (using_server) THEN
      CALL create_request(event_id_set_grid_type_lmdz)
      CALL pack_data(LEN(TRIM(name)))
      CALL pack_data(TRIM(name))
      CALL pack_data(nbp)
      CALL pack_data(offset)
      CALL Finalize_request
    ELSE
     CALL iom__set_grid_type_lmdz(name,nbp,offset)
    ENDIF
      
  END SUBROUTINE event__set_grid_type_lmdz
  
  SUBROUTINE event__set_time_parameters(itau0,zjulian,zdt)
  USE iomanager
  IMPLICIT NONE
    INTEGER,INTENT(IN)  :: itau0
    REAL,   INTENT(IN)  :: zjulian
    REAL,   INTENT(IN)  :: zdt
        
    IF (using_server) THEN
      CALL create_request(event_id_set_time_parameters)
      CALL pack_data(itau0)
      CALL pack_data(zjulian)
      CALL pack_data(zdt)
      CALL Finalize_request
    ELSE
      CALL iom__set_time_parameters(itau0,zjulian,zdt)
    ENDIF
    
  END SUBROUTINE event__set_time_parameters


  SUBROUTINE event__set_calendar(str_calendar)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: str_calendar
        
    IF (using_server) THEN
      CALL create_request(event_id_set_calendar)
      CALL pack_data(LEN(TRIM(str_calendar)))
      CALL pack_data(TRIM(str_calendar))
      CALL Finalize_request
    ELSE
     CALL iom__set_calendar(str_calendar)
    ENDIF
      
  END SUBROUTINE event__set_calendar
  
  
  SUBROUTINE event__enable_field(varname)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)     :: varname

    IF (using_server) THEN
      CALL create_request(event_id_enable_field)
      CALL pack_data(len(varname))
      CALL pack_data(varname)
      CALL Finalize_request 
    ELSE
      CALL iom__enable_field(varname)
    ENDIF
    
  END SUBROUTINE event__enable_field

  SUBROUTINE event__disable_field(varname)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)     :: varname

    IF (using_server) THEN
      CALL create_request(event_id_disable_field)
      CALL pack_data(len(varname))
      CALL pack_data(varname)
      CALL Finalize_request 
    ELSE
      CALL iom__disable_field(varname)
    ENDIF
    
  END SUBROUTINE event__disable_field
  
  SUBROUTINE event__write_field1d(varname,var)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)     :: varname
    REAL, DIMENSION(:),INTENT(IN) :: var

    IF (using_server) THEN
      CALL create_request(event_id_write_field1d)
      CALL pack_data(len(varname))
      CALL pack_data(size(var,1))
      CALL pack_data(varname)
      CALL pack_field(var)
      CALL Finalize_request 
    ELSE
      CALL iom__write_field1d(varname,var)
    ENDIF
    
  END SUBROUTINE event__write_field1d

  SUBROUTINE event__write_field2d(varname,var)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)     :: varname
    REAL, DIMENSION(:,:),INTENT(IN) :: var

    IF (using_server) THEN
      CALL create_request(event_id_write_field2d)
      CALL pack_data(len(varname))
      CALL pack_data(size(var,1))
      CALL pack_data(size(var,2))
      CALL pack_data(varname)
      CALL pack_field(var)
      CALL Finalize_request 
    ELSE
      CALL iom__write_field2d(varname,var)
    ENDIF
    
  END SUBROUTINE event__write_field2d
  
  SUBROUTINE event__write_field3d(varname,var)
  USE iomanager
  IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN)     :: varname
    REAL, DIMENSION(:,:,:),INTENT(IN) :: var

    IF (using_server) THEN
      CALL create_request(event_id_write_field3d)
      CALL pack_data(len(varname))
      CALL pack_data(size(var,1))
      CALL pack_data(size(var,2))
      CALL pack_data(size(var,3))
      CALL pack_data(varname)
      CALL pack_field(var)
      CALL Finalize_request 
    ELSE
      CALL iom__write_field3d(varname,var)
    ENDIF
  END SUBROUTINE event__write_field3d
  

  SUBROUTINE event__set_timestep(timestep)
  USE iomanager
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: timestep
    
    IF (using_server) THEN
      CALL create_request(event_id_set_timestep)
      CALL pack_data(timestep)
      CALL Finalize_request
    ELSE
      CALL iom__set_timestep(timestep)
    ENDIF

  END SUBROUTINE event__set_timestep
    
    
  SUBROUTINE event__close_io_definition
  USE iomanager
  IMPLICIT NONE
    
    IF (using_server) THEN
      CALL create_request(event_id_close_io_definition)
      CALL Finalize_request
    ELSE
      CALL iom__close_io_definition
    ENDIF

  END SUBROUTINE event__close_io_definition

  
  SUBROUTINE event__set_attribut(id,attrib)
  USE iomanager
  USE mod_attribut
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: id
    TYPE(attribut),INTENT(IN)   :: attrib
 
    IF (using_server) THEN
      CALL create_request(event_id_set_attribut)
      CALL pack_data(len(id))
      CALL pack_data(id)
      CALL pack_data(attrib)
      CALL Finalize_request
    ELSE
      CALL iom__set_attribut(id,attrib)
    ENDIF
 
  END SUBROUTINE event__set_attribut  

  SUBROUTINE event__stop_ioserver
  USE iomanager
  IMPLICIT NONE
 
    IF (using_server) THEN
      CALL create_request(event_id_stop_ioserver)
      is_last_request=.TRUE.
      CALL Finalize_request
    ELSE
      CALL iom__finalize
    ENDIF

  END SUBROUTINE event__stop_ioserver
    
  
END MODULE mod_event_client
