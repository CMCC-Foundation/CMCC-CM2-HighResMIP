MODULE mod_interface_ioipsl

 
  INTEGER,PARAMETER     :: id_file=1


CONTAINS

  SUBROUTINE init_interface_ioipsl
  USE xmlio
  IMPLICIT NONE
  
    
  END SUBROUTINE init_interface_ioipsl

  SUBROUTINE set_calendar(str_calendar)
  USE ioipsl
    CHARACTER(LEN=*) :: str_calendar
    
    CALL ioconf_calendar(str_calendar)
    
  END SUBROUTINE set_calendar
  
  
  SUBROUTINE set_time_parameters(ini_timestep0,zjulian0,timestep0)
  USE xmlio
  IMPLICIT NONE
    INTEGER :: ini_timestep0
    REAL :: zjulian0, timestep0

    timestep_value=timestep0
    initial_timestep=ini_timestep0
    initial_date=zjulian0

  END SUBROUTINE set_time_parameters

  SUBROUTINE Create_file_definition(nb_server,server_rank)
  USE ioipsl
  USE xmlio
  USE mod_ioserver_namelist
  IMPLICIT NONE
    INTEGER,INTENT(IN)  :: nb_server
    INTEGER,INTENT(IN)  :: server_rank
    
    TYPE(file_dep),POINTER :: pt_file_dep
    TYPE(file),POINTER     :: pt_file
    TYPE(field),POINTER    :: pt_field
    TYPE(grid),POINTER     :: pt_grid
    TYPE(zoom),POINTER     :: pt_zoom
    TYPE(axis),POINTER     :: pt_axis
    TYPE(domain),POINTER   :: pt_domain

    TYPE(sorted_list),POINTER :: axis_id
    LOGICAL :: found                  
    INTEGER :: ioipsl_axis_id
    INTEGER :: ioipsl_file_id
    INTEGER :: ioipsl_hori_id
    INTEGER :: ioipsl_domain_id
    INTEGER :: i,j
    CHARACTER(LEN=20) :: direction
    CHARACTER(LEN=255) :: full_name
    CALL xmlio__close_definition 
    
    ALLOCATE(axis_id)
    
    DO i=1,file_enabled%size
      
      pt_file_dep=>file_enabled%at(i)%pt

      IF (pt_file_dep%fields%size>0) THEN
        CALL sorted_list__new(axis_id)

        pt_file=>pt_file_dep%file
      
        pt_grid=>pt_file_dep%grids%at(1)%pt
        pt_domain=>pt_grid%domain
        pt_zoom=>pt_file_dep%zooms%at(1)%pt
!        print *,TRIM(pt_file%name),' ',TRIM(pt_zoom%id)
!        print*,'Global --->',pt_zoom%ni_glo,pt_zoom%nj_glo,pt_zoom%ibegin_glo,pt_zoom%jbegin_glo
!        print*,'Local  --->',pt_zoom%ni_loc,pt_zoom%nj_loc,pt_zoom%ibegin_loc,pt_zoom%jbegin_loc
      
        IF (pt_zoom%ni_loc*pt_zoom%nj_loc > 0) THEN
        
          full_name=TRIM(pt_file%name)
          IF (pt_file%has_name_suffix) full_name=TRIM(full_name)//TRIM(pt_file%name_suffix) 
          IF ( (pt_zoom%ni_loc == pt_zoom%ni_glo) .AND. (pt_zoom%nj_loc == pt_zoom%nj_glo) ) THEN 

            CALL histbeg(TRIM(full_name),pt_domain%ni,pt_domain%lon, pt_domain%nj, pt_domain%lat, & 
                       pt_zoom%ibegin_loc, pt_zoom%ni_loc,pt_zoom%jbegin_loc,pt_zoom%nj_loc,         &
                       initial_timestep, initial_date, timestep_value,                               &
                       ioipsl_hori_id, ioipsl_file_id, snc4chunks=snc4ioset)
           ELSE                                              

            CALL set_ioipsl_domain_id(pt_grid,nb_server,server_rank,ioipsl_domain_id)
            CALL histbeg(TRIM(full_name),pt_domain%ni,pt_domain%lon, pt_domain%nj, pt_domain%lat,  &
                       pt_zoom%ibegin_loc, pt_zoom%ni_loc,pt_zoom%jbegin_loc,pt_zoom%nj_loc,          &
                       initial_timestep, initial_date, timestep_value,                                &
                       ioipsl_hori_id, ioipsl_file_id,domain_id=ioipsl_domain_id, snc4chunks=snc4ioset)                                              
          
           ENDIF
        
      
          DO j=1,pt_file_dep%axis%size
            pt_axis=>pt_file_dep%axis%at(j)%pt
            CALL sorted_list__find(axis_id,hash(Pt_axis%name),ioipsl_axis_id,found)
            IF (.NOT. found) THEN
              IF (TRIM(pt_axis%name) /= "none") THEN
              
                IF (pt_axis%has_positive) THEN 
                  IF (pt_axis%positive) THEN
                    direction="up"
                  ELSE
                    direction="down"
                  ENDIF
                ELSE
                  direction='unknown'
                ENDIF

                CALL histvert(ioipsl_file_id, TRIM(pt_axis%name),TRIM(pt_axis%description),    &
                             TRIM(pt_axis%unit), pt_axis%size,pt_axis%values, ioipsl_axis_id,  &
                           pdirect=direction)
                CALL sorted_list__add(axis_id,hash(Pt_axis%name),ioipsl_axis_id)
              ENDIF
            ENDIF
         ENDDO
        
          DO j=1,pt_file_dep%fields%size
            pt_field=>pt_file_dep%fields%at(j)%pt
            IF (pt_field%axis%name=="none") THEN
              pt_field%internal(id_file)=ioipsl_file_id
              CALL histdef(ioipsl_file_id, TRIM(pt_field%name), pt_field%description,            &
                        &  pt_field%unit, pt_field%grid%domain%ni, pt_field%grid%domain%nj,      &
                        &  ioipsl_hori_id, 1, 1, 1, -99, 32, pt_field%operation,                 &
                        &  real(pt_field%freq_op), real(pt_file%output_freq) )
            ELSE
              pt_field%internal(id_file)=ioipsl_file_id
              CALL sorted_list__find(axis_id,hash(Pt_field%axis%name),ioipsl_axis_id,found)
              CALL histdef(ioipsl_file_id, TRIM(pt_field%name), pt_field%description,          &
                         & pt_field%unit, pt_field%grid%domain%ni, pt_field%grid%domain%nj,    &
                         & ioipsl_hori_id, pt_field%axis%size, 1, pt_field%axis%size,          &
                         & ioipsl_axis_id, 32, pt_field%operation, real(pt_field%freq_op),     &
                         & real(pt_file%output_freq) )
            ENDIF
          ENDDO
          CALL histend(ioipsl_file_id, snc4chunks=snc4ioset)
        ENDIF
        CALL sorted_list__delete(axis_id)
      ENDIF
    ENDDO
    
    DEALLOCATE(axis_id)
     
  END SUBROUTINE Create_file_definition  



   SUBROUTINE write_ioipsl_2d(varname,var)
   USE ioipsl
   USE xmlio
   IMPLICIT NONE
     CHARACTER(len=*),INTENT(IN) :: varname
     REAL            ,INTENT(IN) :: var(:,:)
     
     TYPE(field_dep),POINTER :: pt_field_base
     TYPE(field)    ,POINTER :: pt_field
     INTEGER :: nindex(size(var))
     INTEGER :: ioipsl_file_id
     INTEGER :: pos
     LOGICAL :: found
     INTEGER :: i
      
     CALL sorted_list__find(sorted_id,hash(varname),pos,found)

     IF (found) THEN
       pt_field_base=>field_id%at(pos)%pt
     
       DO i=1,pt_field_base%field_out%size
         pt_field=>pt_field_base%field_out%at(i)%pt%field
         IF ( pt_field%zoom%ni_loc * pt_field%zoom%nj_loc > 0) THEN           
           ioipsl_file_id=pt_field%internal(id_file)
           CALL histwrite(ioipsl_file_id, TRIM(pt_field%name), timestep_number, var, size(var), nindex)
         ENDIF
       ENDDO
     ENDIF
     
   END SUBROUTINE write_ioipsl_2d
       
     
   SUBROUTINE write_ioipsl_3d(varname,var)
   USE ioipsl
   USE xmlio
   IMPLICIT NONE
     CHARACTER(len=*),INTENT(IN) :: varname
     REAL            ,INTENT(IN) :: var(:,:,:)
     
     TYPE(field_dep),POINTER :: pt_field_base
     TYPE(field)    ,POINTER :: pt_field
     INTEGER :: nindex(size(var))
     INTEGER :: ioipsl_file_id
     INTEGER :: pos
     LOGICAL :: found
     INTEGER :: i
        
     CALL sorted_list__find(sorted_id,hash(varname),pos,found)

     IF (found) THEN
       pt_field_base=>field_id%at(pos)%pt
     
       DO i=1,pt_field_base%field_out%size
         pt_field=>pt_field_base%field_out%at(i)%pt%field
         IF ( pt_field%zoom%ni_loc * pt_field%zoom%nj_loc > 0) THEN           
           ioipsl_file_id=pt_field%internal(id_file)
           CALL histwrite(ioipsl_file_id, TRIM(pt_field%name), timestep_number, var, size(var), nindex)
         ENDIF
       ENDDO
     ENDIF
   END SUBROUTINE write_ioipsl_3d
           

  SUBROUTINE set_timestep(timestep_nb0)
  USE xmlio
  IMPLICIT NONE
    INTEGER,INTENT(IN) :: timestep_nb0
     
    timestep_number=timestep_nb0

  END SUBROUTINE set_timestep


  SUBROUTINE set_ioipsl_domain_id(pt_grid, nb_server,server_rank,domain_id)
  USE xmlio
  USE mod_ioserver_para
  USE ioipsl
  IMPLICIT NONE
    TYPE(grid), POINTER :: pt_grid
    INTEGER,INTENT(IN)  :: nb_server
    INTEGER,INTENT(IN)  :: server_rank
    INTEGER,INTENT(OUT) :: domain_id
    TYPE(domain), POINTER :: pt_domain
               
    INTEGER,DIMENSION(2) :: ddid
    INTEGER,DIMENSION(2) :: dsg
    INTEGER,DIMENSION(2) :: dsl
    INTEGER,DIMENSION(2) :: dpf
    INTEGER,DIMENSION(2) :: dpl
    INTEGER,DIMENSION(2) :: dhs
    INTEGER,DIMENSION(2) :: dhe 
    
    pt_domain=>pt_grid%domain

    ddid = (/ 1,2 /)
    dsg  = (/ pt_grid%ni, pt_grid%nj /)
    dsl  = (/ pt_domain%ni, pt_domain%nj /)
    dpf  = (/ pt_domain%ibegin,pt_domain%jbegin /)
    dpl  = (/ pt_domain%iend, pt_domain%jend /)
    dhs  = (/ 0,0 /)
    dhe  = (/ 0,0 /)
    
    call flio_dom_set(nb_server,server_rank,ddid,dsg,dsl,dpf,dpl,dhs,dhe, &
                      'BOX',domain_id)
  
  END SUBROUTINE set_ioipsl_domain_id


  SUBROUTINE ioipsl_finalize
  USE ioipsl
  IMPLICIT NONE

    CALL histclo
    
  END SUBROUTINE ioipsl_finalize
  
END MODULE mod_interface_ioipsl
