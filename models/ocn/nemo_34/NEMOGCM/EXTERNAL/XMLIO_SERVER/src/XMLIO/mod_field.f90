MODULE mod_field

  USE mod_xmlio_parameters
  USE mod_sorted_list
  USE mod_axis
  USE mod_grid
  USE mod_zoom
   
  IMPLICIT NONE

  TYPE, PUBLIC :: field
    CHARACTER(len=str_len)         :: id
    LOGICAL                        :: has_id
    CHARACTER(len=str_len)         :: name
    LOGICAL                        :: has_name
    CHARACTER(len=str_len)         :: description
    LOGICAL                        :: has_description
    CHARACTER(len=str_len)         :: unit
    LOGICAL                        :: has_unit
    CHARACTER(len=str_len)         :: operation
    LOGICAL                        :: has_operation
    INTEGER                        :: freq_op
    LOGICAL                        :: has_freq_op
    CHARACTER(len=str_len)         :: axis_ref
    LOGICAL                        :: has_axis_ref
    CHARACTER(len=str_len)         :: grid_ref
    LOGICAL                        :: has_grid_ref
    CHARACTER(len=str_len)         :: zoom_ref
    LOGICAL                        :: has_zoom_ref
    INTEGER                        :: level
    LOGICAL                        :: has_level
    INTEGER                        :: prec
    LOGICAL                        :: has_prec
    CHARACTER(len=str_len)         :: field_ref
    LOGICAL                        :: has_field_ref
    TYPE(field),POINTER            :: field_base
    LOGICAL                        :: has_field_base
    LOGICAL                        :: enabled
    LOGICAL                        :: has_enabled
    LOGICAL                        :: solved_field_ref
    TYPE(axis), POINTER            :: axis
    LOGICAL                        :: has_axis
    TYPE(grid),POINTER             :: grid
    LOGICAL                        :: has_grid
    TYPE(zoom),POINTER             :: zoom
    LOGICAL                        :: has_zoom
    INTEGER                        :: internal(internal_field)
    
    
  END TYPE field
  
  INTERFACE field__set_attribut
    MODULE PROCEDURE field__set_attribut_id,field__set_attribut_pt
  END INTERFACE
  
  INCLUDE 'vector_field_def.inc'
  
  TYPE(vector_field),POINTER,SAVE            :: field_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids
  

CONTAINS
  INCLUDE 'vector_field_contains.inc'

  SUBROUTINE field__swap_context(saved_field_ids,saved_ids)
  IMPLICIT NONE
    TYPE(vector_field),POINTER :: saved_field_ids
    TYPE(sorted_list),POINTER  :: saved_ids
    
    field_Ids=>saved_field_ids
    Ids=>saved_Ids

  END SUBROUTINE field__swap_context

  SUBROUTINE field__init
  IMPLICIT NONE
    
    CALL vector_field__new(field_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE field__init
  
  SUBROUTINE field__get(Id,Pt_field)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(field),POINTER             :: Pt_field

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_field=>field_ids%at(Pos)%Pt
    ELSE
      Pt_field=>NULL()
    ENDIF
    
  END SUBROUTINE field__get
    
  
  SUBROUTINE field__new(pt_field,Id)
   USE string_function
   IMPLICIT NONE
   TYPE(field), POINTER          :: pt_field
   CHARACTER(LEN=*),OPTIONAL     :: Id
   
   INTEGER              :: Pos

   Pt_field%has_id=.FALSE.
   pt_field%has_name = .FALSE.
   pt_field%has_description = .FALSE.
   pt_field%has_unit = .FALSE.
   pt_field%has_operation = .FALSE.
   pt_field%has_freq_op = .FALSE.
   pt_field%has_axis_ref = .FALSE.
   pt_field%has_grid_ref = .FALSE.
   pt_field%has_zoom_ref = .FALSE.
   pt_field%has_prec = .FALSE.
   pt_field%has_level = .FALSE. 
   pt_field%has_field_ref = .FALSE.
   pt_field%has_field_base = .FALSE.
   pt_field%has_enabled = .FALSE.
   Pt_field%solved_field_ref=.FALSE.
   Pt_field%has_axis=.FALSE.
   Pt_field%has_grid=.FALSE.
   Pt_field%has_zoom=.FALSE.
      
   IF (PRESENT(Id)) THEN
     Pt_field%id=TRIM(ADJUSTL(Id))
     Pt_field%has_id=.TRUE.
     CALL vector_field__set_new(field_Ids,Pt_field,Pos)
     CALL sorted_list__Add(Ids,hash(id),Pos)
   ENDIF
   
 END SUBROUTINE field__new


 SUBROUTINE field__set(p_field, name, ref, description, unit, operation, freq_op, axis_ref, grid_ref, zoom_ref, prec, level, &
                       enabled)

    TYPE(field), pointer :: p_field
    CHARACTER(len=*), OPTIONAL :: name
    CHARACTER(len=*), OPTIONAL :: ref
    CHARACTER(len=*), OPTIONAL  :: description
    CHARACTER(len=*), OPTIONAL  :: unit
    CHARACTER(len=*), OPTIONAL :: operation
    INTEGER, OPTIONAL  :: freq_op
    CHARACTER(len=*),OPTIONAL :: axis_ref
    CHARACTER(len=*),OPTIONAL :: grid_ref
    CHARACTER(len=*),OPTIONAL :: zoom_ref
    INTEGER, OPTIONAL :: prec
    INTEGER, OPTIONAL :: level
    LOGICAL, OPTIONAL :: enabled

    IF (PRESENT(name)) THEN
        p_field%name=TRIM(ADJUSTL(name))
        p_field%has_name = .TRUE.
    ENDIF
    IF (PRESENT(ref)) THEN
        p_field%field_ref=TRIM(ADJUSTL(ref))
        p_field%has_field_ref = .TRUE.
    ENDIF
    IF (PRESENT(description)) THEN
        p_field%description=TRIM(ADJUSTL(description))
        p_field%has_description = .TRUE.
    ENDIF
 
    IF (PRESENT(unit)) then
        p_field%unit=TRIM(ADJUSTL(unit))
        p_field%has_unit = .TRUE.
    ENDIF
    IF (PRESENT(operation)) THEN
        p_field%operation=TRIM(ADJUSTL(operation))
        p_field%has_operation = .TRUE.
    ENDIF
    IF (PRESENT(freq_op)) THEN
        p_field%freq_op=freq_op
        p_field%has_freq_op = .TRUE.
    ENDIF
    IF (PRESENT(axis_ref)) THEN
        p_field%axis_ref=TRIM(ADJUSTL(axis_ref))
        p_field%has_axis_ref = .TRUE.
    ENDIF
    IF (PRESENT(grid_ref)) THEN
        p_field%grid_ref=TRIM(ADJUSTL(grid_ref))
        p_field%has_grid_ref = .TRUE.
    ENDIF

    IF (PRESENT(zoom_ref)) THEN
        p_field%zoom_ref=TRIM(ADJUSTL(zoom_ref))
        p_field%has_zoom_ref = .TRUE.
    ENDIF

    IF (PRESENT(prec)) then
        p_field%prec=prec
        p_field%has_prec = .TRUE.
    ENDIF
    
    IF (PRESENT(level)) then
        p_field%level=level
        p_field%has_level = .TRUE.
    ENDIF

    IF (PRESENT(enabled)) then
        p_field%enabled=enabled
        p_field%has_enabled = .TRUE.
    ENDIF

  END SUBROUTINE field__set

  SUBROUTINE field__set_attribut_id(id,attrib,Ok)
  USE mod_attribut
  USE mod_field_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)   :: id
    TYPE(attribut),INTENT(IN)     :: attrib
    LOGICAL,OPTIONAL,INTENT(out)  :: Ok
    
    TYPE(field),POINTER             :: Pt_field
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_field=>field_ids%at(Pos)%Pt
      CALL field__set_attribut(Pt_field,attrib)
      IF (PRESENT(OK)) ok=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'Field id :',id,'is undefined'
        CALL error('mod_field::field__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF  
  END SUBROUTINE field__set_attribut_id

    
  SUBROUTINE field__set_attribut_pt(pt_field,attrib)
  USE mod_attribut
  USE mod_field_attribut
  USE error_msg
  IMPLICIT NONE
    TYPE(field),POINTER             :: Pt_field
    TYPE(attribut),INTENT(IN) :: attrib

    SELECT CASE(attrib%name)
      CASE (field__name)
        IF (attrib%type==string0) CALL  field__set(pt_field,name=attrib%string0_ptr) ; RETURN
      CASE (field__field_ref)
        IF (attrib%type==string0) CALL  field__set(pt_field,ref=attrib%string0_ptr) ; RETURN
      CASE (field__description)
        IF (attrib%type==string0) CALL  field__set(pt_field,description=attrib%string0_ptr) ; RETURN
      CASE (field__unit)
        IF (attrib%type==string0) CALL  field__set(pt_field,unit=attrib%string0_ptr) ; RETURN
      CASE (field__operation)
        IF (attrib%type==string0) CALL  field__set(pt_field,operation=attrib%string0_ptr) ; RETURN
      CASE (field__freq_op)
        IF (attrib%type==integer0) CALL  field__set(pt_field,freq_op=attrib%integer0_ptr) ; RETURN
      CASE (field__axis_ref)
        IF (attrib%type==string0) CALL  field__set(pt_field,axis_ref=attrib%string0_ptr) ; RETURN
      CASE (field__grid_ref)
        IF (attrib%type==string0) CALL  field__set(pt_field,grid_ref=attrib%string0_ptr) ; RETURN
      CASE (field__zoom_ref)
        IF (attrib%type==string0) CALL  field__set(pt_field,zoom_ref=attrib%string0_ptr) ; RETURN
      CASE (field__prec)
        IF (attrib%type==integer0) CALL  field__set(pt_field,prec=attrib%integer0_ptr) ; RETURN
      CASE (field__level)
        IF (attrib%type==integer0) CALL  field__set(pt_field,level=attrib%integer0_ptr) ; RETURN
      CASE (field__enabled)
        IF (attrib%type==logical0) CALL  field__set(pt_field,enabled=attrib%logical0_ptr) ; RETURN
     END SELECT

     WRITE(message,*) 'field attribut ',attrib%name,' : type :',attrib%type,   &
                      ' : Attribute type is incompatible with the provided value'
     CALL error('mod_field::field__set_attribut')
    
  END SUBROUTINE field__set_attribut_pt    

  SUBROUTINE field__print(pt_field)

    TYPE(field), POINTER :: pt_field
    
    PRINT *,"--- FIELD ---"

    IF (pt_field%has_id) THEN
        PRINT *, 'id : ',TRIM(pt_field%id)
    ELSE
        PRINT *, 'id undefined '
    ENDIF

    IF (pt_field%has_name) THEN
        PRINT *, 'name : ',TRIM(pt_field%name)
    ELSE
        PRINT *, 'name undefined '
    ENDIF
    IF (pt_field%has_description) THEN
        PRINT *, 'description : ',TRIM(pt_field%description)
    ELSE
        PRINT *, 'description undefined '
    ENDIF
    IF (pt_field%has_unit) THEN
        PRINT *, 'unit : ',TRIM(pt_field%unit)
    ELSE
        PRINT *, 'unit undefined '
    ENDIF
    IF (pt_field%has_operation) THEN
        PRINT *, 'operation ',TRIM(pt_field%operation)
    ELSE
        PRINT *, 'operation undefined '
    ENDIF
    IF (pt_field%has_freq_op) THEN
        PRINT *, 'freq_op ',pt_field%freq_op
    ELSE
        PRINT *, 'freq_op undefined '
    ENDIF

    IF (pt_field%has_axis_ref) THEN
        PRINT *, 'axis_ref : ',TRIM(pt_field%axis_ref)
    ELSE
        PRINT *, 'axis_ref undefined '
    ENDIF

    IF (pt_field%has_grid_ref) THEN
        PRINT *, 'grid_ref : ',TRIM(pt_field%grid_ref)
    ELSE
        PRINT *, 'grid_ref undefined '
    ENDIF

    IF (pt_field%has_zoom_ref) THEN
        PRINT *, 'zoom_ref : ',TRIM(pt_field%zoom_ref)
    ELSE
        PRINT *, 'zoom_ref undefined '
    ENDIF
    
    IF (pt_field%has_field_ref) THEN
        PRINT *, 'field_ref : ',TRIM(pt_field%field_ref)
    ELSE
        PRINT *, 'field_ref undefined '
    ENDIF

!    call vert_axis__print(pt_field%p_vert_axis)
!
!    IF (pt_field%is_vert_axis_ref_def) THEN
!        PRINT *, 'pt_field%vert_axis_ref ',TRIM(pt_field%vert_axis_ref)
!    ELSE
!        PRINT *, 'pt_field%vert_axis_ref undefined '
!    ENDIF
    IF (pt_field%has_prec) THEN
        PRINT *, 'prec ',pt_field%prec
    ELSE
        PRINT *, 'prec undefined '
    ENDIF
    IF (pt_field%has_level) then
        PRINT *, 'level ',pt_field%level
    ELSE
        PRINT *, 'level undefined '
    ENDIF
    IF (pt_field%has_field_base) THEN
        PRINT *, 'field_base :',TRIM(Pt_field%field_base%id)
    ELSE
        PRINT *, 'field_base indefini'
    ENDIF

    IF (pt_field%has_enabled) THEN
        PRINT *, 'enabled : ',pt_field%enabled
    ELSE
        PRINT *, 'enabled indefini'
    ENDIF
 
    PRINT *,"------------"
    
  END SUBROUTINE field__print

!  SUBROUTINE field__resolve_ref_vert_axis(p_field)
!
!    TYPE(field), POINTER :: p_field
!    CHARACTER(len=str_len) :: name
!
!    IF (p_field%is_vert_axis_ref_def) THEN
!        name=p_field%vert_axis_ref
!        IF (vert_axis_def__is_exist(name)) THEN
!            CALL vert_axis_def__get(name,p_field%p_vert_axis)
!            p_field%is_vert_axis_def = .TRUE.
!        ENDIF
!    ENDIF
!
!  END SUBROUTINE field__resolve_ref_vert_axis

  SUBROUTINE field__apply_default(pt_field_default, pt_field_in, pt_field_out)

    TYPE(field), POINTER :: pt_field_default, pt_field_in, pt_field_out

    IF (pt_field_in%has_name) THEN
        pt_field_out%name=pt_field_in%name
        pt_field_out%has_name=.TRUE.
    ELSE IF ( pt_field_default%has_name) THEN
        pt_field_out%name=pt_field_default%name
        pt_field_out%has_name=.TRUE.
    ELSE
        pt_field_out%has_name=.FALSE.
    ENDIF
        
    IF (pt_field_in%has_description) THEN
        pt_field_out%description=pt_field_in%description
        pt_field_out%has_description=.TRUE.
    ELSE IF ( pt_field_default%has_description ) THEN
        pt_field_out%description=pt_field_default%description
        pt_field_out%has_description=.TRUE.
    ELSE
        pt_field_out%has_description=.FALSE.
    ENDIF

    IF (pt_field_in%has_unit) THEN
        pt_field_out%unit=pt_field_in%unit
        pt_field_out%has_unit=.TRUE.
    ELSE IF ( pt_field_default%has_unit ) THEN
        pt_field_out%unit=pt_field_default%unit
        pt_field_out%has_unit=.TRUE.
    ELSE
        pt_field_out%has_unit=.FALSE.
    ENDIF

    IF (pt_field_in%has_operation) THEN
        pt_field_out%operation=pt_field_in%operation
        pt_field_out%has_operation=.TRUE.
    ELSE IF ( pt_field_default%has_operation ) THEN
        pt_field_out%operation=pt_field_default%operation
        pt_field_out%has_operation=.TRUE.
    ELSE
        pt_field_out%has_operation=.FALSE.
    ENDIF

    IF (pt_field_in%has_freq_op) THEN
        pt_field_out%freq_op=pt_field_in%freq_op
        pt_field_out%has_freq_op=.TRUE.
    ELSE IF ( pt_field_default%has_freq_op ) THEN
        pt_field_out%freq_op=pt_field_default%freq_op
        pt_field_out%has_freq_op=.TRUE.
    ELSE
        pt_field_out%has_freq_op=.FALSE.
    ENDIF

!    IF (pt_field_in%has_axis) THEN
!        pt_field_out%p_axis => pt_field_in%p_axis
!        pt_field_out%has_axis=.TRUE.
!    ELSE IF ( pt_field_default%has_axis ) THEN
!        pt_field_out%p_axis => pt_field_default%p_axis
!        pt_field_out%has_axis=.TRUE.
!    ELSE
!        pt_field_out%has_axis=.FALSE.
!    ENDIF
    
    IF (pt_field_in%has_axis_ref) THEN
        pt_field_out%axis_ref=pt_field_in%axis_ref
        pt_field_out%has_axis_ref=.TRUE.
    ELSE IF ( pt_field_default%has_axis_ref ) THEN
        pt_field_out%axis_ref=pt_field_default%axis_ref
        pt_field_out%has_axis_ref=.TRUE.
    ELSE
        pt_field_out%has_axis_ref=.FALSE.
    ENDIF

    IF (pt_field_in%has_grid_ref) THEN
        pt_field_out%grid_ref=pt_field_in%grid_ref
        pt_field_out%has_grid_ref=.TRUE.
    ELSE IF ( pt_field_default%has_grid_ref ) THEN
        pt_field_out%grid_ref=pt_field_default%grid_ref
        pt_field_out%has_grid_ref=.TRUE.
    ELSE
        pt_field_out%has_grid_ref=.FALSE.
    ENDIF

    IF (pt_field_in%has_zoom_ref) THEN
        pt_field_out%zoom_ref=pt_field_in%zoom_ref
        pt_field_out%has_zoom_ref=.TRUE.
    ELSE IF ( pt_field_default%has_zoom_ref ) THEN
        pt_field_out%zoom_ref=pt_field_default%zoom_ref
        pt_field_out%has_zoom_ref=.TRUE.
    ELSE
        pt_field_out%has_zoom_ref=.FALSE.
    ENDIF

    IF (pt_field_in%has_prec) THEN
        pt_field_out%prec=pt_field_in%prec
        pt_field_out%has_prec=.TRUE.
    ELSE IF ( pt_field_default%has_prec ) THEN
        pt_field_out%prec=pt_field_default%prec
        pt_field_out%has_prec=.TRUE.
    ELSE
        pt_field_out%has_prec=.FALSE.
    ENDIF

    IF (pt_field_in%has_level) THEN
        pt_field_out%level=pt_field_in%level
        pt_field_out%has_level=.TRUE.
    ELSE IF ( pt_field_default%has_level ) THEN
        pt_field_out%level=pt_field_default%level
        pt_field_out%has_level=.TRUE.
    ELSE
        pt_field_out%has_level=.FALSE.
    ENDIF

    IF (pt_field_in%has_enabled) THEN
        pt_field_out%enabled=pt_field_in%enabled
        pt_field_out%has_enabled=.TRUE.
    ELSE IF ( pt_field_default%has_enabled ) THEN
        pt_field_out%enabled=pt_field_default%enabled
        pt_field_out%has_enabled=.TRUE.
    ELSE
        pt_field_out%has_enabled=.FALSE.
    ENDIF
    
  END SUBROUTINE field__apply_default

!  FUNCTION field__is_vert_axis_attached(p_field, vert_axis_name)
!
!    LOGICAL :: field__is_vert_axis_attached
!    TYPE(field), POINTER :: p_field
!    CHARACTER(len=*), INTENT(IN) :: vert_axis_name
!
!    field__is_vert_axis_attached = .false.
!    IF (p_field%is_vert_axis_def) THEN
!        IF (vert_axis_name == p_field%p_vert_axis%name) field__is_vert_axis_attached = .TRUE.
!    ENDIF
!
!  END FUNCTION field__is_vert_axis_attached

  RECURSIVE SUBROUTINE field__solve_field_ref(pt_field)
  USE error_msg
  IMPLICIT NONE
    TYPE(field), POINTER :: pt_field
    
    TYPE(field), POINTER :: field_ref
    
    IF (.NOT. pt_field%solved_field_ref) THEN
      
      IF (pt_field%has_field_ref) THEN
      
        CALL field__get(pt_field%field_ref,field_ref)
      
        IF (.NOT. ASSOCIATED(field_ref)) THEN
          WRITE (message,*) "The field : id = ",pt_field%id,"  name = ",Pt_field%name,   &
                        " has a unknown reference to field : id =",pt_field%field_ref
          CALL error("field__solve_field_ref")
        ENDIF
      
        CALL field__get_field_base(field_ref,pt_field%field_base) 
        Pt_field%has_field_base=.TRUE.
        
        CALL field__apply_default(field_ref,pt_field,pt_field)
        
          
      ELSE 

        IF (pt_field%has_id) THEN
          pt_field%field_base=>pt_field
        ENDIF
      
      ENDIF
      
      IF (.NOT. pt_field%has_name) THEN
        IF (pt_field%has_id) THEN
          pt_field%name=pt_field%id
          pt_field%has_name=.TRUE.
        ENDIF
      ENDIF
    
      Pt_field%solved_field_ref=.TRUE.
    
    ENDIF
    
  END SUBROUTINE field__solve_field_ref
  

  RECURSIVE SUBROUTINE field__get_field_base(pt_field,pt_field_base)
  IMPLICIT NONE
    TYPE(field), POINTER :: pt_field
    TYPE(field), POINTER :: pt_field_base
    
     
    IF (.NOT. Pt_field%solved_field_ref) THEN
      CALL field__solve_field_ref(Pt_field)
    ENDIF
     
    IF (pt_field%has_field_base) THEN
      pt_field_base=>pt_field%field_base
    ELSE
      pt_field_base=>pt_field
    ENDIF
   
 END SUBROUTINE field__get_field_base

 SUBROUTINE field__solve_axis_ref(pt_field)
 USE error_msg
 IMPLICIT NONE
   TYPE(field), POINTER :: pt_field
   
   IF (pt_field%has_axis_ref) THEN
     CALL axis__get(pt_field%axis_ref,pt_field%axis)
     IF (ASSOCIATED(pt_field%axis)) THEN
       pt_field%has_axis=.TRUE.
     ELSE
       WRITE (message,*) "The field : id = ",pt_field%id,"  name = ",Pt_field%name,   &
                        " has a unknown reference to axis : id =",pt_field%axis_ref
       CALL error("mod_field::field__solve_axis_ref")
     ENDIF
   ENDIF
   
 END SUBROUTINE field__solve_axis_ref
   
 SUBROUTINE field__solve_grid_ref(pt_field)
 USE error_msg
 IMPLICIT NONE
   TYPE(field), POINTER :: pt_field
   
   IF (pt_field%has_grid_ref) THEN
     CALL grid__get(pt_field%grid_ref,pt_field%grid)
     IF (ASSOCIATED(pt_field%grid)) THEN
       pt_field%has_grid=.TRUE.
     ELSE
       WRITE (message,*) "The field : id = ",pt_field%id,"  name = ",Pt_field%name,   &
                        " has a unknown reference to grid : id =",pt_field%grid_ref
       CALL error("mod_field::field__solve_grid_ref")
     ENDIF
   ENDIF
   
 END SUBROUTINE field__solve_grid_ref
   
 SUBROUTINE field__solve_zoom_ref(pt_field)
 USE error_msg
 IMPLICIT NONE
   TYPE(field), POINTER :: pt_field
   
   IF (.NOT. pt_field%has_zoom_ref) THEN
     IF (pt_field%has_grid_ref) THEN
       pt_field%has_zoom_ref=.TRUE.
       pt_field%zoom_ref=pt_field%grid_ref
     ENDIF
   ENDIF
    
   IF (pt_field%has_zoom_ref) THEN
     CALL zoom__get(pt_field%zoom_ref,pt_field%zoom)
     IF (ASSOCIATED(pt_field%zoom)) THEN
       pt_field%has_zoom=.TRUE.
     ELSE
       WRITE (message,*) "The field : id = ",pt_field%id,"  name = ",Pt_field%name,   &
                        " has a unknown reference to zoom : id =",pt_field%zoom_ref
       CALL error("mod_field::field__solve_zoom_ref")
     ENDIF
   ENDIF
   
 END SUBROUTINE field__solve_zoom_ref

   
END MODULE mod_field
