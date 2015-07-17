MODULE iomanager
  INTEGER,PRIVATE,SAVE :: nb_client
  INTEGER,PRIVATE,SAVE :: nb_server
  INTEGER,PRIVATE,SAVE :: server_rank
  INTEGER,PRIVATE,SAVE :: current_rank
  
CONTAINS

  SUBROUTINE iom__init(nb_client_,nb_server_,server_rank_)
  IMPLICIT NONE
   INTEGER,INTENT(IN) :: nb_client_
   INTEGER,INTENT(IN) :: nb_server_
   INTEGER,INTENT(IN) :: server_rank_
   
     nb_client=nb_client_ 
     nb_server=nb_server_
     server_rank=server_rank_
     
  END SUBROUTINE iom__init
 
  SUBROUTINE iom__parse_xml_file(filename)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*) :: filename
    
    IF (current_rank==nb_client) CALL xmlio__init(filename)
     
  END SUBROUTINE iom__parse_xml_file
   
  SUBROUTINE iom__swap_context(id)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*) :: id
    
    IF (current_rank==nb_client) CALL context__swap(id)
     
  END SUBROUTINE iom__swap_context

  SUBROUTINE iom__set_current_rank(rank)
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: rank
    
    current_rank=rank
  
  END SUBROUTINE iom__set_current_rank
  
  SUBROUTINE iom__set_vert_axis(name,vert_value)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    REAL,INTENT(IN)             :: vert_value(:)
    TYPE(axis), POINTER :: pt_axis
    
    IF (current_rank==nb_client) THEN
      CALL axis__get(name, pt_axis)
      CALL axis__set(pt_axis, a_size=size(vert_value), values=vert_value)
    ENDIF
 
  END SUBROUTINE iom__set_vert_axis

 
  SUBROUTINE iom__set_grid_dimension(name,ni_glo,nj_glo)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,INTENT(IN) :: ni_glo
    INTEGER,INTENT(IN) :: nj_glo
    TYPE(grid), POINTER :: pt_grid
  
    IF (current_rank==nb_client) THEN
      CALL grid__get(name,pt_grid)
      CALL grid__set_dimension(pt_grid,ni_glo,nj_glo)
    ENDIF

  END SUBROUTINE iom__set_grid_dimension
  
  SUBROUTINE iom__set_grid_domain(name,ni,nj,ibegin,jbegin,lon,lat)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,INTENT(IN) :: ni
    INTEGER,INTENT(IN) :: nj
    INTEGER,INTENT(IN) :: ibegin
    INTEGER,INTENT(IN) :: jbegin
    REAL,INTENT(IN)    :: lon(ni,nj)
    REAL,INTENT(IN)    :: lat(ni,nj)

    TYPE(grid),   POINTER :: pt_grid
    TYPE(domain), POINTER :: pt_domain

      CALL grid__get(name,pt_grid)
      CALL grid__get_new_subdomain(pt_grid,current_rank,pt_domain)
      CALL domain__set(pt_domain,current_rank,ni,nj,ibegin,jbegin,lon,lat)

  END SUBROUTINE iom__set_grid_domain

  SUBROUTINE iom__set_grid_type_nemo(name)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name

    TYPE(grid),   POINTER :: pt_grid
    TYPE(domain), POINTER :: pt_domain
      
      CALL grid__get(name,pt_grid)
      CALL grid__get_subdomain(pt_grid,current_rank,pt_domain)
      CALL domain__set_type_box(pt_domain)

  END SUBROUTINE iom__set_grid_type_nemo

  SUBROUTINE iom__set_grid_type_lmdz(name,nbp,offset)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: name
    INTEGER,INTENT(IN)          :: nbp
    INTEGER,INTENT(IN)          :: offset

    TYPE(grid),   POINTER :: pt_grid
    TYPE(domain), POINTER :: pt_domain
    LOGICAL,ALLOCATABLE   :: mask(:,:)
      
      CALL grid__get(name,pt_grid)
      CALL grid__get_subdomain(pt_grid,current_rank,pt_domain)
      ALLOCATE(mask(pt_domain%ni,pt_domain%nj))
      mask(:,:)=.TRUE.
      mask(1:offset,1)=.FALSE.
      mask(MOD(offset+nbp-1,pt_domain%ni)+2:pt_domain%ni,pt_domain%nj)=.FALSE.
      CALL domain__set_type_box(pt_domain,mask)

  END SUBROUTINE iom__set_grid_type_lmdz
      
  SUBROUTINE iom__set_time_parameters(itau0,zjulian,zdt)
  USE mod_interface_ioipsl
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: itau0
    REAL,INTENT(IN)     :: zjulian
    REAL,INTENT(IN)     :: zdt
      
    IF (current_rank==nb_client) THEN
      CALL set_time_parameters(itau0, zjulian, zdt)
    ENDIF
  END SUBROUTINE iom__set_time_parameters
  
  SUBROUTINE iom__close_io_definition
  USE mod_interface_ioipsl
  IMPLICIT NONE
  
    IF (current_rank==nb_client) CALL Create_file_definition(nb_server,server_rank)
  
  END SUBROUTINE iom__close_io_definition
      

  SUBROUTINE iom__set_timestep(timestep)
  USE mod_interface_ioipsl
  IMPLICIT NONE
    INTEGER, INTENT(IN) :: timestep
    
    IF (current_rank==nb_client) CALL set_timestep(timestep)
    
  END SUBROUTINE iom__set_timestep

  SUBROUTINE iom__set_calendar(str_calendar)
  USE mod_interface_ioipsl
  IMPLICIT NONE
    CHARACTER(LEN=*) :: str_calendar
    
    IF (current_rank==nb_client) CALL set_calendar(str_calendar)
    
  END SUBROUTINE iom__set_calendar

  SUBROUTINE iom__enable_field(varname)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*) :: varname
    TYPE(field),POINTER :: pt_field
    
    CALL field__get(TRIM(varname), pt_field)
    CALL field__set(pt_field,enabled=.TRUE.)
    
  END SUBROUTINE iom__enable_field
  
  SUBROUTINE iom__disable_field(varname)
  USE xmlio
  IMPLICIT NONE
    CHARACTER(LEN=*) :: varname
    TYPE(field),POINTER :: pt_field
    
    CALL field__get(TRIM(varname), pt_field)
    CALL field__set(pt_field,enabled=.FALSE.)
    
  END SUBROUTINE iom__disable_field
   

  SUBROUTINE iom__write_field1d(varname,var)
  USE xmlio
  USE field_bufferize
  USE mod_interface_ioipsl
  IMPLICIT NONE
    CHARACTER(LEN=*) :: varname
    REAL             :: var(:)

    TYPE(field),POINTER  :: pt_field
    TYPE(domain),POINTER :: subdomain
    TYPE(domain),POINTER :: local_domain
    INTEGER :: ni,nj,ibegin,jbegin
    INTEGER :: id_rank
    
    CALL field__get(varname,pt_field)
    id_rank=pt_field%grid%ranks(current_rank)
    local_domain=>pt_field%grid%domain
    
    subdomain=>pt_field%grid%subdomain%at(id_rank)%pt
    ni=subdomain%ni
    nj=subdomain%nj
    ibegin=subdomain%ibegin-local_domain%ibegin+1
    jbegin=subdomain%jbegin-local_domain%jbegin+1 

    IF (subdomain%type==box) THEN

      WRITE(message,*) 'Field must have 2 or 3 dimensions for box domain. There, it have only one.' 
      CALL error("iom__write_field2d")

    ELSE IF (subdomain%type==orange) THEN
      IF (size(var,1)/=subdomain%nbp) THEN 
         WRITE(message,*) 'Field dimensions are not compliant with the associated grid ',  &
                          'field nbp :',size(var,1),'   grid nbp :',subdomain%nbp
         CALL error("iom__write_field2d")
      ENDIF
      IF (pt_field%axis%name /= 'none') THEN
        WRITE(message,*) "Missing axis dimension for this field"
        CALL error("iom__write_field2d")
      ENDIF
    ENDIF
    
    
    IF (current_rank==1) THEN
      CALL init_field_bufferize(local_domain%ni,local_domain%nj,1)
    ENDIF
    
    CALL bufferize_field(ni,ibegin,nj,jbegin,1,1,subdomain%nbp,var,         &
                         subdomain%i_index,subdomain%j_index,subdomain%mask)
    
    IF (current_rank==nb_client) THEN
      ni=local_domain%ni
      nj=local_domain%nj
      CALL write_ioipsl_2d(varname,Field_buffer(1:ni,1:nj,1))
    ENDIF
    
  END SUBROUTINE iom__write_field1d

   
  SUBROUTINE iom__write_field2d(varname,var)
  USE xmlio
  USE field_bufferize
  USE mod_interface_ioipsl
  IMPLICIT NONE
    CHARACTER(LEN=*) :: varname
    REAL             :: var(:,:)

    TYPE(field),POINTER  :: pt_field
    TYPE(domain),POINTER :: subdomain
    TYPE(domain),POINTER :: local_domain
    INTEGER :: ni,nj,nk,ibegin,jbegin
    INTEGER :: id_rank
    
    CALL field__get(varname,pt_field)
    id_rank=pt_field%grid%ranks(current_rank)
    local_domain=>pt_field%grid%domain
    
    subdomain=>pt_field%grid%subdomain%at(id_rank)%pt
    ni=subdomain%ni
    nj=subdomain%nj
    ibegin=subdomain%ibegin-local_domain%ibegin+1
    jbegin=subdomain%jbegin-local_domain%jbegin+1 
    nk=pt_field%axis%size

    IF (subdomain%type==box) THEN
      IF (size(var,1)/=ni .OR. size(var,2)/=nj) THEN
        WRITE(message,*) 'Field dimensions are not compliant with the associated grid ',  &
                          'field dim :',size(var,1),',',size(var,2),'   grid dim',ni,',',nj
        CALL error("iom__write_field2d")
      ENDIF
      
      IF (pt_field%axis%name/="none") THEN
        WRITE(message,*) "Missing axis dimension for this field"
        CALL error("iom__write_field2d")
      ENDIF
    
    ELSE IF (subdomain%type==orange) THEN
      IF (size(var,1)/=subdomain%nbp) THEN 
         WRITE(message,*) 'Field dimensions are not compliant with the associated grid ',  &
                          'field nbp :',size(var,1),'   grid nbp : ',subdomain%nbp
         CALL error("iom__write_field2d")
      ENDIF
      IF (nk /= size(var,2)) THEN
         WRITE(message,*) 'Field dimensions are not compliant with the associated axis', &
                          'field dim :',size(var,2),'   axis dim :', nk
         CALL error("iom__write_field2d")
      ENDIF
    ENDIF
    
    
    IF (current_rank==1) THEN
      CALL init_field_bufferize(local_domain%ni,local_domain%nj,nk)
    ENDIF
    
    CALL bufferize_field(ni,ibegin,nj,jbegin,nk,1,subdomain%nbp,var,         &
                         subdomain%i_index,subdomain%j_index,subdomain%mask)
    
    IF (current_rank==nb_client) THEN
      ni=local_domain%ni
      nj=local_domain%nj
      
      IF (pt_field%axis%name=="none") THEN
        CALL write_ioipsl_3d(varname,Field_buffer(1:ni,1:nj,1:nk))
      ELSE 
        CALL write_ioipsl_2d(varname,Field_buffer(1:ni,1:nj,1))
      ENDIF
    
    ENDIF
    
  END SUBROUTINE iom__write_field2d

  SUBROUTINE iom__write_field3d(varname,var)
  USE xmlio
  USE field_bufferize
  USE mod_interface_ioipsl
  IMPLICIT NONE
    CHARACTER(LEN=*) :: varname
    REAL             :: var(:,:,:)

    TYPE(field),POINTER  :: pt_field
    TYPE(domain),POINTER :: subdomain
    TYPE(domain),POINTER :: local_domain
    INTEGER :: ni,nj,nk,ibegin,jbegin
    INTEGER :: id_rank

    CALL field__get(varname,pt_field)
    id_rank=pt_field%grid%ranks(current_rank)
    local_domain=>pt_field%grid%domain
    
    subdomain=>pt_field%grid%subdomain%at(id_rank)%pt
    ni=subdomain%ni
    nj=subdomain%nj
    ibegin=subdomain%ibegin-local_domain%ibegin+1
    jbegin=subdomain%jbegin-local_domain%jbegin+1 
    nk=pt_field%axis%size

    IF (subdomain%type==box) THEN
      IF (size(var,1)/=ni .OR. size(var,2)/=nj) THEN
        WRITE(message,*) 'Field dimensions are not compliant with the associated grid ',  &
                          'field dim :',size(var,1),',',size(var,2),'   grid dim',ni,',',nj
        CALL error("iom__write_field3d")
      ENDIF
      
      IF (pt_field%axis%name=='none' .OR. nk/=size(var,3)) THEN
        WRITE(message,*) 'Field dimension is not compliant with the associated axis : ', &
                         'field dim : ',size(var,3),'  axis dim : ',nk
        CALL error("iom__write_field3d")
      ENDIF

    ELSE IF (subdomain%type==orange) THEN

      WRITE(message,*) 'Field have too much dimensions for box domain. There, it has 3.' 
      CALL error("iom__write_field3d")
      
    ENDIF
    
    
    IF (current_rank==1) THEN
      CALL init_field_bufferize(local_domain%ni,local_domain%nj,nk)
    ENDIF
    
    CALL bufferize_field(ni,ibegin,nj,jbegin,nk,1,subdomain%nbp,var,         &
                         subdomain%i_index,subdomain%j_index,subdomain%mask)
    
    IF (current_rank==nb_client) THEN
      ni=local_domain%ni
      nj=local_domain%nj
      
      CALL write_ioipsl_3d(varname,Field_buffer(1:ni,1:nj,1:nk))
    
    ENDIF
    
  END SUBROUTINE iom__write_field3d
 
  SUBROUTINE iom__Finalize
  USE mod_interface_ioipsl
  IMPLICIT NONE
  
    IF (current_rank==nb_client) CALL ioipsl_finalize
    
  END SUBROUTINE iom__Finalize


  SUBROUTINE iom__set_attribut(id,attrib)
  USE mod_attribut
  USE mod_object
  USE mod_field
  USE mod_field_group
  USE mod_file
  USE mod_file_group
  USE mod_axis
  USE mod_axis_group
  USE mod_grid
  USE mod_grid_group
  USE mod_zoom
  IMPLICIT NONE
    CHARACTER(LEN=*) :: id
    TYPE(attribut)   :: attrib  
    LOGICAL          :: success
    IF (current_rank==nb_client) THEN
    
      SELECT CASE(attrib%object)
        CASE(field_object)
          CALL field_group__set_attribut(id,attrib,success)
          IF (.NOT. success) CALL field__set_attribut(id,attrib,success)         
        CASE(file_object)
          CALL file_group__set_attribut(id,attrib,success)
          IF (.NOT. success) CALL file__set_attribut(id,attrib,success)         
        CASE(axis_object)
          CALL axis_group__set_attribut(id,attrib,success)
          IF (.NOT. success) CALL axis__set_attribut(id,attrib,success)         
        CASE(grid_object)
          CALL grid_group__set_attribut(id,attrib,success)
          IF (.NOT. success) CALL grid__set_attribut(id,attrib,success)         
        CASE(zoom_object)
          CALL zoom__set_attribut(id,attrib,success)
      END SELECT
    ENDIF
    
  END SUBROUTINE iom__set_attribut
   
END MODULE iomanager 
 
