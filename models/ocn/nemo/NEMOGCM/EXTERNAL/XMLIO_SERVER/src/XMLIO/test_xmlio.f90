PROGRAM test_xmlio
  USE xmlio
  USE mod_interface_ioipsl
  IMPLICIT NONE
  
  TYPE(field),POINTER :: Pt_f
  TYPE(field_group),POINTER :: Pt_fg
  TYPE(file),POINTER :: Pt_file
  TYPE(grid),POINTER :: pt_grid 
  INTEGER,PARAMETER  :: ni=360
  INTEGER,PARAMETER  :: nj=180
  REAL,DIMENSION(ni,nj) :: lon
  REAL,DIMENSION(ni,nj) :: lat
  INTEGER               :: i,j
  TYPE(axis),POINTER :: pt_axis 
  
  DO i=1,ni
    DO j=1,nj
      lon(i,j)=i
      lat(i,j)=j
    ENDDO
  ENDDO  
  
  
  
  CALL xmlio__init("iodef.xml")

  CALL grid__get("grid_T",pt_grid)
  CALL domain__set(pt_grid%domain,0,ni,nj,1,1,lon,lat)

  CALL grid__get("grid_U",pt_grid)
  CALL domain__set(pt_grid%domain,0,ni,nj,1,1,lon,lat)

  CALL grid__get("grid_V",pt_grid)
  CALL domain__set(pt_grid%domain,0,ni,nj,1,1,lon,lat)

  CALL grid__get("grid_W",pt_grid)
  CALL domain__set(pt_grid%domain,0,ni,nj,1,1,lon,lat)
  
  CALL axis__get("deptht",pt_axis)
  CALL axis__set(pt_axis,a_size=5,values=(/1.,2.,3.,4.,5./))

  CALL axis__get("depthu",pt_axis)
  CALL axis__set(pt_axis,a_size=5,values=(/1.,2.,3.,4.,5./))

  CALL axis__get("depthv",pt_axis)
  CALL axis__set(pt_axis,a_size=5,values=(/1.,2.,3.,4.,5./))

  CALL axis__get("depthw",pt_axis)
  CALL axis__set(pt_axis,a_size=5,values=(/1.,2.,3.,4.,5./))


  CALL axis_group__apply_default(axis_definition)
  CALL grid_group__apply_default(grid_definition)
  CALL field_group__apply_default(field_definition)
  CALL file_group__apply_default(file_definition)
  
  CALL field_group__solve_ref(field_definition)
  CALL file_group__solve_field_ref(file_definition)

  CALL file_group__check(file_definition)
  CALL axis_group__check(axis_definition)

!
  CALL axis_group__print(axis_definition)
  CALL grid_group__print(grid_definition)
  CALL field_group__print(field_definition)
  CALL file_group__print(file_definition)
!  CALL field_group__get_new_field(field_definition,pt_f,"toto")
!  CALL field__set(Pt_f,name="toto",unit="kg")
!  CALL field_group__get_new_field(field_definition,pt_f,"tata")
!  CALL field__set(Pt_f,name="tata",unit="kg")
!  CALL field_group__get_new_field(field_definition,pt_f,"titi")
!  CALL field__set(Pt_f,name="titi",unit="kg")
!  CALL field_group__get_new_field(field_definition,pt_f,"xyz")
!  CALL field__set(Pt_f,name="xyz",unit="kg")
   CALL set_file_dependency
   CALL set_field_enabled
   CALL set_field_dependency
   
   CALL set_time_parameters(0,0.,3600.)
   
   CALL create_header
   
END PROGRAM test_xmlio
