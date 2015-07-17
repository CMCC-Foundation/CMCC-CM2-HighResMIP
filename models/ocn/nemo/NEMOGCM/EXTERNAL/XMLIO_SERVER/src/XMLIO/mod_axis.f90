MODULE mod_axis

  USE mod_xmlio_parameters
  USE mod_sorted_list

  
  IMPLICIT NONE

  TYPE, PUBLIC :: axis
    CHARACTER(len=str_len)      :: id
    LOGICAL                     :: has_id
    CHARACTER(len=str_len)      :: name
    LOGICAL                     :: has_name
    INTEGER                     :: size
    LOGICAL                     :: has_size
    CHARACTER(len=str_len)      :: description
    LOGICAL                     :: has_description
    CHARACTER(len=str_len)      :: unit
    LOGICAL                     :: has_unit
    LOGICAL                     :: positive
    LOGICAL                     :: has_positive
    REAL, DIMENSION(:), POINTER :: values
    LOGICAL                     :: has_values

  END TYPE axis

  INCLUDE 'vector_axis_def.inc'
  
  TYPE(vector_axis),POINTER,SAVE             :: axis_Ids
  TYPE(sorted_list),POINTER,SAVE,PRIVATE     :: Ids 

  INTERFACE axis__set_attribut
    MODULE PROCEDURE axis__set_attribut_id,axis__set_attribut_pt
  END INTERFACE

CONTAINS
  INCLUDE 'vector_axis_contains.inc'

  SUBROUTINE axis__swap_context(saved_axis_Ids,saved_Ids)
  IMPLICIT NONE
    TYPE(vector_axis),POINTER          :: saved_axis_Ids
    TYPE(sorted_list),POINTER          :: saved_Ids 
    
    axis_ids=>saved_axis_ids
    ids=>saved_ids
    
  END SUBROUTINE axis__swap_context

  SUBROUTINE axis__init
  IMPLICIT NONE
    
    CALL vector_axis__new(axis_Ids)
    CALL sorted_list__new(Ids)
   
  END SUBROUTINE axis__init
  
  SUBROUTINE axis__get(Id,Pt_axis)
  USE string_function
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)     :: Id
    TYPE(axis),POINTER              :: Pt_axis

    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_axis=>axis_ids%at(Pos)%Pt
    ELSE
      Pt_axis=>NULL()
    ENDIF
    
  END SUBROUTINE axis__get
  
  SUBROUTINE axis__new(pt_axis,Id)
  USE string_function
  IMPLICIT NONE
   TYPE(axis), POINTER           :: pt_axis
   CHARACTER(LEN=*),OPTIONAL     :: Id
   INTEGER                       :: Pos
   
   pt_axis%has_id          = .FALSE.
   pt_axis%has_name        = .FALSE.
   pt_axis%has_size        = .FALSE.
   pt_axis%has_description = .FALSE.
   pt_axis%has_unit        = .FALSE.
   pt_axis%has_values      = .FALSE.
   pt_axis%has_positive    = .FALSE. 
     
   IF (PRESENT(Id)) THEN
     Pt_axis%id=TRIM(ADJUSTL(Id))
     Pt_axis%has_id=.TRUE.
     CALL vector_axis__set_new(axis_Ids,Pt_axis,Pos)
     CALL sorted_list__Add(Ids,hash(id),Pos)
   ENDIF

 END SUBROUTINE axis__new

  SUBROUTINE axis__set(pt_axis, name, description, unit, a_size, values, positive)
  IMPLICIT NONE
    TYPE(axis), POINTER         :: pt_axis
    CHARACTER(len=*)  ,OPTIONAL :: name
    CHARACTER(len=*)  ,OPTIONAL :: description
    CHARACTER(len=*)  ,OPTIONAL :: unit
    INTEGER           ,OPTIONAL :: a_size
    REAL, DIMENSION(:),OPTIONAL :: values
    LOGICAL           ,OPTIONAL :: positive

    IF (PRESENT(name)) THEN
        pt_axis%name=TRIM(ADJUSTL(name))
        pt_axis%has_name = .TRUE.
    ENDIF

    IF (PRESENT(description)) THEN
        pt_axis%description=TRIM(ADJUSTL(description))
        pt_axis%has_description = .TRUE.
    ENDIF
 
    IF (PRESENT(unit)) then
        pt_axis%unit=TRIM(ADJUSTL(unit))
        pt_axis%has_unit = .TRUE.
    ENDIF

    IF (PRESENT(a_size)) then
        pt_axis%size=a_size
        pt_axis%has_size = .TRUE.
    ENDIF
    
    IF (PRESENT(values)) then
        IF (pt_axis%has_values) DEALLOCATE(pt_axis%values)  
        ALLOCATE(pt_axis%values(size(values)))
        pt_axis%values(:)=values(:)
        pt_axis%has_values = .TRUE.
    ENDIF

    IF (PRESENT(positive)) then
        pt_axis%positive=positive
        pt_axis%has_positive = .TRUE.
    ENDIF

  END SUBROUTINE axis__set

  SUBROUTINE axis__set_attribut_id(id,attrib,ok)
  USE mod_attribut
  USE error_msg
  IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)   :: id
    TYPE(attribut),INTENT(IN)     :: attrib
    LOGICAL,OPTIONAL,INTENT(OUT)  :: ok
    
    TYPE(axis),POINTER              :: Pt_axis
    INTEGER                         :: Pos
    LOGICAL                         :: success
    
    CALL sorted_list__find(Ids,hash(Id),Pos,success)
    IF (success) THEN
      Pt_axis=>axis_ids%at(Pos)%Pt
      CALL axis__set_attribut_pt(Pt_axis,attrib)
      IF (PRESENT(OK)) OK=.TRUE.
    ELSE
      IF (.NOT.PRESENT(OK)) THEN
        WRITE(message,*) 'axis id :',id,'is undefined'
        CALL error('mod_axis::axis__set_attribut')
      ELSE
        OK=.FALSE.
      ENDIF
    ENDIF  
  
  END SUBROUTINE axis__set_attribut_id
      
  SUBROUTINE axis__set_attribut_pt(Pt_axis,attrib)
  USE mod_attribut
  USE mod_axis_attribut
  USE error_msg
  IMPLICIT NONE
    TYPE(axis),POINTER        :: Pt_axis
    TYPE(attribut),INTENT(IN) :: attrib
    
    SELECT CASE(attrib%name)
      CASE (axis__name)
        IF (attrib%type==string0) CALL  axis__set(pt_axis,name=attrib%string0_ptr) ; RETURN
      CASE (axis__description)
        IF (attrib%type==string0) CALL  axis__set(pt_axis,description=attrib%string0_ptr) ; RETURN
      CASE (axis__unit)
        IF (attrib%type==string0) CALL  axis__set(pt_axis,unit=attrib%string0_ptr) ; RETURN
      CASE (axis__size)
        IF (attrib%type==integer0) CALL  axis__set(pt_axis,a_size=attrib%integer0_ptr) ; RETURN
      CASE (axis__values)
        IF (attrib%type==real1) CALL  axis__set(pt_axis,values=attrib%real1_ptr) ; RETURN
      CASE (axis__positive)
        IF (attrib%type==logical0) CALL  axis__set(pt_axis,positive=attrib%logical0_ptr) ; RETURN
       END SELECT

     WRITE(message,*) 'axis attribut ',attrib%name,' : type :',attrib%type,   &
                      ' : Attribute type is incompatible with the provided value'
     CALL error('mod_axis::axis__set_attribut')
    
  END SUBROUTINE axis__set_attribut_pt
  
  SUBROUTINE axis__print(pt_axis)
  IMPLICIT NONE
    TYPE(axis), POINTER         :: pt_axis

    PRINT *,"---- AXIS ----"
    IF (pt_axis%has_id) THEN
      PRINT *,"id = ",TRIM(pt_axis%id)
    ELSE
      PRINT *,"id undefined"
    ENDIF
    
    IF (pt_axis%has_name) THEN
      PRINT *,"name = ",TRIM(pt_axis%name)
    ELSE
      PRINT *,"name undefined"
    ENDIF
    
    IF (pt_axis%has_description) THEN
      PRINT *,"description = ",TRIM(pt_axis%description)
    ELSE
      PRINT *,"description undefined"
    ENDIF
 
    IF (pt_axis%has_unit) THEN
      PRINT *,"unit = ",TRIM(pt_axis%unit)
    ELSE
      PRINT *,"unit undefined"
    ENDIF

    IF (pt_axis%has_size) THEN
      PRINT *,"size = ",pt_axis%size
    ELSE
      PRINT *,"size undefined"
    ENDIF

    IF (pt_axis%has_values) THEN
      PRINT *,"values = ",pt_axis%values
    ELSE
      PRINT *,"values undefined"
    ENDIF

    IF (pt_axis%has_positive) THEN
      PRINT *,"positive = ",pt_axis%positive
    ELSE
      PRINT *,"positive undefined"
    ENDIF

  END SUBROUTINE axis__print


  SUBROUTINE axis__apply_default(pt_axis_default, pt_axis_in, pt_axis_out)

    TYPE(axis), POINTER :: pt_axis_default, pt_axis_in, pt_axis_out

    IF (pt_axis_in%has_name) THEN
        pt_axis_out%name=pt_axis_in%name
        pt_axis_out%has_name=.TRUE.
    ELSE IF ( pt_axis_default%has_name) THEN
        pt_axis_out%name=pt_axis_default%name
        pt_axis_out%has_name=.TRUE.
    ELSE
        pt_axis_out%has_name=.FALSE.
    ENDIF
        
    IF (pt_axis_in%has_description) THEN
        pt_axis_out%description=pt_axis_in%description
        pt_axis_out%has_description=.TRUE.
    ELSE IF ( pt_axis_default%has_description ) THEN
        pt_axis_out%description=pt_axis_default%description
        pt_axis_out%has_description=.TRUE.
    ELSE
        pt_axis_out%has_description=.FALSE.
    ENDIF

    IF (pt_axis_in%has_unit) THEN
        pt_axis_out%unit=pt_axis_in%unit
        pt_axis_out%has_unit=.TRUE.
    ELSE IF ( pt_axis_default%has_unit ) THEN
        pt_axis_out%unit=pt_axis_default%unit
        pt_axis_out%has_unit=.TRUE.
    ELSE
        pt_axis_out%has_unit=.FALSE.
    ENDIF

    IF (pt_axis_in%has_size) THEN
        pt_axis_out%size=pt_axis_in%size
        pt_axis_out%has_size=.TRUE.
    ELSE IF ( pt_axis_default%has_size ) THEN
        pt_axis_out%size=pt_axis_default%size
        pt_axis_out%has_size=.TRUE.
    ELSE
        pt_axis_out%has_size=.FALSE.
    ENDIF

    IF (pt_axis_in%has_values) THEN
        pt_axis_out%values(:)=pt_axis_in%values(:)
        pt_axis_out%has_values=.TRUE.
    ELSE IF ( pt_axis_default%has_values ) THEN
        pt_axis_out%values(:)=pt_axis_default%values(:)
        pt_axis_out%has_values=.TRUE.
    ELSE
        pt_axis_out%has_values=.FALSE.
    ENDIF

    IF (pt_axis_in%has_positive) THEN
        pt_axis_out%positive=pt_axis_in%positive
        pt_axis_out%has_positive=.TRUE.
    ELSE IF ( pt_axis_default%has_positive ) THEN
        pt_axis_out%positive=pt_axis_default%positive
        pt_axis_out%has_positive=.TRUE.
    ELSE
        pt_axis_out%has_positive=.FALSE.
    ENDIF
    
  END SUBROUTINE axis__apply_default

  SUBROUTINE axis__check(pt_axis)
  USE error_msg
  IMPLICIT NONE
    TYPE(axis), POINTER :: pt_axis
      
    IF (.NOT. pt_axis%has_name) THEN
      IF (pt_axis%has_id) THEN
        pt_axis%name=TRIM(pt_axis%id)
      ELSE
        WRITE(message,*) "Axis has no name and no id" 
        CALL error("mod_axis::axis__check")
      ENDIF 
    ENDIF
 
 END SUBROUTINE axis__Check

END MODULE mod_axis
