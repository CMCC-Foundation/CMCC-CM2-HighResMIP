MODULE mod_event_server
  USE mod_pack, ONLY : unpack_data, unpack_field
  USE mod_event_parameters
  USE iomanager
  
CONTAINS

  SUBROUTINE Process_event(current_rank,is_terminated)
  IMPLICIT NONE
    INTEGER :: event_id
    INTEGER, INTENT(IN) :: current_rank
    LOGICAL,INTENT(OUT) :: is_terminated
     
    CALL iom__set_current_rank(current_rank)
     
    is_terminated=.FALSE.
    
    CALL unpack_data(event_id)
    
    SELECT CASE (event_id)
    
      CASE (event_id_swap_context)
        CALL event__swap_context
        
      CASE (event_id_parse_xml_file)
        CALL event__parse_xml_file 

      CASE (event_id_set_vert_axis)
        CALL event__set_vert_axis 
     
      CASE (event_id_set_grid_dimension)
        CALL event__set_grid_dimension
      
      CASE (event_id_set_grid_domain)
        CALL event__set_grid_domain

      CASE (event_id_set_grid_type_nemo)
        CALL event__set_grid_type_nemo

      CASE (event_id_set_grid_type_lmdz)
        CALL event__set_grid_type_lmdz

      CASE (event_id_set_time_parameters)
        CALL event__set_time_parameters

      CASE (event_id_close_io_definition)
        CALL event__close_io_definition  

      CASE (event_id_set_timestep)
        CALL event__set_timestep

      CASE (event_id_set_calendar)
        CALL event__set_calendar

      CASE (event_id_enable_field)
        CALL event__enable_field
      
      CASE (event_id_disable_field)
        CALL event__disable_field

      CASE (event_id_write_Field1d)
        CALL event__write_Field1d

      CASE (event_id_write_Field2d)
        CALL event__write_Field2d

      CASE (event_id_write_Field3d)
        CALL event__write_Field3d

      CASE (event_id_set_attribut)
        CALL event__set_attribut

      CASE (event_id_stop_ioserver)
        is_terminated=.TRUE.  
        PRINT *,"TERMINATE_EVENT RECEIVED"
        

      CASE DEFAULT 
        STOP 'UNDEFINED EVENT'
      
     END SELECT
     
   END SUBROUTINE Process_event
   
  SUBROUTINE event__swap_context
  IMPLICIT NONE
    INTEGER :: id_size
   
    CALL unpack_data(id_size)
    CALL sub_internal(id_size)
     
  CONTAINS
    
    SUBROUTINE sub_internal(id_size)
      INTEGER :: id_size
      CHARACTER(LEN=id_size) :: id     
     
       CALL unpack_data(id)
     
       CALL iom__swap_context(id)
       
     END SUBROUTINE sub_internal
  
  END SUBROUTINE event__swap_context
  

  SUBROUTINE event__parse_xml_file
  IMPLICIT NONE
    INTEGER :: name_size
   
    CALL unpack_data(name_size)
    CALL sub_internal(name_size)
     
  CONTAINS
    
    SUBROUTINE sub_internal(name_size)
      INTEGER :: name_size
      CHARACTER(LEN=name_size) :: filename     
     
       CALL unpack_data(filename)
     
       CALL iom__parse_xml_file(filename)
     END SUBROUTINE sub_internal
  
  END SUBROUTINE event__parse_xml_file
  
  
  SUBROUTINE event__set_grid_dimension
  IMPLICIT NONE    
    INTEGER :: name_size
    INTEGER :: ni_glo
    INTEGER :: nj_glo
   
    CALL unpack_data(name_size)
    CALL sub_internal(name_size)
     
  CONTAINS
    
    SUBROUTINE sub_internal(name_size)
      INTEGER :: name_size
      CHARACTER(LEN=name_size) :: name     
     
       CALL unpack_data(name)
       CALL unpack_data(ni_glo)
       CALL unpack_data(nj_glo)
     
       CALL iom__set_grid_dimension(name,ni_glo,nj_glo)
     END SUBROUTINE sub_internal
   
   END SUBROUTINE event__set_grid_dimension


  SUBROUTINE event__set_grid_domain 
  IMPLICIT NONE    
    INTEGER :: name_size
    INTEGER :: ni
    INTEGER :: nj
    INTEGER :: ibegin
    INTEGER :: jbegin
    REAL,ALLOCATABLE :: lon(:,:)
    REAL,ALLOCATABLE :: lat(:,:)

    CALL unpack_data(name_size)
    CALL sub_internal(name_size)
     
  CONTAINS
    
    SUBROUTINE sub_internal(name_size)
      INTEGER :: name_size
      CHARACTER(LEN=name_size) :: name     
     
       CALL unpack_data(name)
    
       CALL unpack_data(ni)
       CALL unpack_data(nj)
       CALL unpack_data(ibegin)
       CALL unpack_data(jbegin)
    
       ALLOCATE(lon(ni,nj))
       ALLOCATE(lat(ni,nj))
       CALL unpack_data(lon)
       CALL unpack_data(lat)
    
       CALL iom__set_grid_domain(name,ni,nj,ibegin,jbegin,lon,lat)

     END SUBROUTINE sub_internal
   
  END SUBROUTINE event__set_grid_domain   


  SUBROUTINE event__set_grid_type_nemo 
  IMPLICIT NONE    
    INTEGER :: name_size

    CALL unpack_data(name_size)
    CALL sub_internal(name_size)
     
  CONTAINS
    
    SUBROUTINE sub_internal(name_size)
      INTEGER :: name_size
      CHARACTER(LEN=name_size) :: name     
     
       CALL unpack_data(name)
       CALL iom__set_grid_type_nemo(name)

     END SUBROUTINE sub_internal
   
  END SUBROUTINE event__set_grid_type_nemo   

  SUBROUTINE event__set_grid_type_lmdz 
  IMPLICIT NONE    
    INTEGER :: name_size

    CALL unpack_data(name_size)
    CALL sub_internal(name_size)
     
  CONTAINS
    
    SUBROUTINE sub_internal(name_size)
      INTEGER :: name_size
      CHARACTER(LEN=name_size) :: name     
      INTEGER                  :: nbp
      INTEGER                  :: offset
      
       CALL unpack_data(name)
       CALL unpack_data(nbp)
       CALL unpack_data(offset)
       CALL iom__set_grid_type_lmdz(name,nbp,offset)

     END SUBROUTINE sub_internal
   
  END SUBROUTINE event__set_grid_type_lmdz   

  SUBROUTINE event__set_vert_axis
  IMPLICIT NONE
    INTEGER :: name_size
    INTEGER :: vert_size 
    REAL,ALLOCATABLE :: vert_value(:)
    
    CALL unpack_data(name_size)
    CALL sub_internal(name_size)

  CONTAINS
    
    SUBROUTINE sub_internal(name_size)
      INTEGER :: name_size
      CHARACTER(LEN=name_size) :: name
    
      CALL unpack_data(name)
      CALL unpack_data(vert_size)
      ALLOCATE(vert_value(vert_size))
      CALL unpack_data(vert_value) 
      
      CALL iom__set_vert_axis(name,vert_value)
      
    END SUBROUTINE sub_internal
  END SUBROUTINE event__set_vert_axis

  SUBROUTINE event__set_time_parameters
  IMPLICIT NONE
    INTEGER   :: itau0
    REAL      :: zjulian
    REAL      :: zdt
      
    CALL unpack_data(itau0)
    CALL unpack_data(zjulian)
    CALL unpack_data(zdt)
    
    CALL iom__set_time_parameters(itau0,zjulian,zdt)
     
  END SUBROUTINE event__set_time_parameters
  

  SUBROUTINE event__enable_field
  IMPLICIT NONE
    INTEGER :: lenc
      
    CALL unpack_data(lenc)
    CALL sub_internal(lenc)
  
  CONTAINS
    SUBROUTINE sub_internal(lenc)
    IMPLICIT NONE
      INTEGER :: lenc
      CHARACTER(len=lenc) :: varname
      
      CALL unpack_data(varname)
      
      CALL iom__enable_field(varname)

    END SUBROUTINE sub_internal
  END SUBROUTINE event__enable_field
  

  SUBROUTINE event__disable_field
  IMPLICIT NONE
    INTEGER :: lenc
      
    CALL unpack_data(lenc)
    CALL sub_internal(lenc)
  
  CONTAINS
    SUBROUTINE sub_internal(lenc)
    IMPLICIT NONE
      INTEGER :: lenc
      CHARACTER(len=lenc) :: varname
      
      CALL unpack_data(varname)
      
      CALL iom__disable_field(varname)

    END SUBROUTINE sub_internal
    
  END SUBROUTINE event__disable_field
  
      
  SUBROUTINE event__write_field1D
  IMPLICIT NONE
    INTEGER :: lenc
    INTEGER :: dim1
      
    CALL unpack_data(lenc)
    CALL unpack_data(dim1)
    CALL sub_internal(lenc,dim1)
  
  CONTAINS
    SUBROUTINE sub_internal(lenc,dim1)
    IMPLICIT NONE
      INTEGER :: lenc
      INTEGER :: dim1
      CHARACTER(len=lenc) :: varname
      REAL                :: var(dim1)
      
      CALL unpack_data(varname)
      CALL unpack_field(var)
      
      CALL iom__write_Field1d(varname,var)

    END SUBROUTINE sub_internal
  END SUBROUTINE event__write_field1d

  SUBROUTINE event__write_field2D
  IMPLICIT NONE
    INTEGER :: lenc
    INTEGER :: dim1
    INTEGER :: dim2
      
    CALL unpack_data(lenc)
    CALL unpack_data(dim1)
    CALL unpack_data(dim2)
    CALL sub_internal(lenc,dim1,dim2)
  
  CONTAINS
    SUBROUTINE sub_internal(lenc,dim1,dim2)
    IMPLICIT NONE
      INTEGER :: lenc
      INTEGER :: dim1
      INTEGER :: dim2
      CHARACTER(len=lenc) :: varname
      REAL                :: var(dim1,dim2)
      
      CALL unpack_data(varname)
      CALL unpack_field(var)
      
      CALL iom__write_Field2d(varname,var)

    END SUBROUTINE sub_internal
  END SUBROUTINE event__write_field2d
    
    
  SUBROUTINE event__write_field3d
  IMPLICIT NONE
    INTEGER :: lenc
    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3
    
    CALL unpack_data(lenc)
    CALL unpack_data(dim1)
    CALL unpack_data(dim2)
    CALL unpack_data(dim3)
    CALL sub_internal(lenc,dim1,dim2,dim3)
 
  CONTAINS
 
    SUBROUTINE sub_internal(lenc,dim1,dim2,dim3)
    IMPLICIT NONE
      INTEGER :: lenc
      INTEGER :: dim1
      INTEGER :: dim2
      INTEGER :: dim3
 
      CHARACTER(len=lenc) :: varname
      REAL                :: var(dim1,dim2,dim3)
        
      CALL unpack_data(varname)
      CALL unpack_field(var)
       
      CALL iom__write_field3d(varname,var)
         
    END SUBROUTINE sub_internal

  END SUBROUTINE event__write_field3d
    
    
  SUBROUTINE event__set_timestep
  IMPLICIT NONE
    INTEGER :: timestep
     
    CALL unpack_data(timestep)
    CALL iom__set_timestep(timestep)
    
  END SUBROUTINE event__set_timestep
    

  SUBROUTINE event__set_calendar
  IMPLICIT NONE
    INTEGER :: lenc
      
    CALL unpack_data(lenc)
    CALL sub_internal(lenc)
  
  CONTAINS
    SUBROUTINE sub_internal(lenc)
    IMPLICIT NONE
      INTEGER :: lenc
      CHARACTER(len=lenc) :: str_calendar
      
      CALL unpack_data(str_calendar)
      
      CALL iom__set_calendar(str_calendar)

    END SUBROUTINE sub_internal
  END SUBROUTINE event__set_calendar
  
    
  SUBROUTINE event__close_io_definition
  IMPLICIT NONE
    
    CALL iom__close_io_definition
    
  END SUBROUTINE event__close_io_definition
  
  SUBROUTINE event__set_attribut
   USE mod_attribut
   IMPLICIT NONE
     TYPE(attribut) :: attrib
     INTEGER        :: len_id
     
     CALL unpack_data(len_id)
     CALL sub_internal
   CONTAINS
      
     SUBROUTINE sub_internal
       CHARACTER(LEN=len_id) :: id
       
       CALL unpack_data(id)
       CALL unpack_data(attrib)
       CALL iom__set_attribut(id,attrib)
       CALL attr_deallocate(attrib)
     END SUBROUTINE sub_internal
  END SUBROUTINE event__set_attribut    

END MODULE mod_event_server   
