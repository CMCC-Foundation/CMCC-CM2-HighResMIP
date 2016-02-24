MODULE parsing_xml
  USE flib_dom
  USE error_msg
  USE string_function
  
CONTAINS

  SUBROUTINE parsing_xml_file(name)
  IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: name
    TYPE(fnode), POINTER :: myDoc

    myDoc => parsefile(name)
    CALL parsing_root(myDoc)

    PRINT *, 'Le parsing est termine !!! '
  
  END SUBROUTINE parsing_xml_file
  
 
  SUBROUTINE parsing_root(root)
  IMPLICIT NONE
    TYPE(fnode), POINTER :: root
    TYPE(fnode), POINTER :: child_node
    TYPE(fnodeList), POINTER :: child_list

    INTEGER :: il
    CHARACTER(len=100) :: node_name
     
    IF (hasChildNodes(root)) THEN
      child_list => getChildnodes(root)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
          
        SELECT CASE (TRIM(node_name)) 
          
           CASE ('simulation') 
             CALL parsing_simulation(child_node)

           CASE DEFAULT
             IF (is_bad_node(node_name)) THEN
               WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing root'
               CALL Warning("mod_parse_xml:parsing_root")
             ENDIF
        END SELECT
      ENDDO
    ENDIF
             
  END SUBROUTINE parsing_root

  SUBROUTINE parsing_simulation(root)
  IMPLICIT NONE
    TYPE(fnode), POINTER :: root
    TYPE(fnode), POINTER :: child_node
    TYPE(fnodeList), POINTER :: child_list

    INTEGER :: il
    CHARACTER(len=100) :: node_name
     
     IF (hasChildNodes(root)) THEN
      child_list => getChildnodes(root)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
          
        SELECT CASE (TRIM(node_name)) 
          
           CASE ('context') 
             CALL parsing_context(child_node)

           CASE DEFAULT
             IF (is_bad_node(node_name)) THEN
               WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing simulation'
               CALL Warning("mod_parse_xml:parsing_simulationt")
             ENDIF
         END SELECT
       
       ENDDO
    ENDIF
             
  END SUBROUTINE parsing_simulation

  SUBROUTINE parsing_context(node)
  USE mod_context
  USE mod_axis_definition
  USE mod_grid_definition
  USE mod_field_definition
  USE mod_file_definition
  IMPLICIT NONE  
    TYPE(fnode), POINTER         :: node
    
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    TYPE(axis),POINTER       :: attribute
    LOGICAL                  :: is_root
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
    
    IF (is_attribute_exist(node,"id")) THEN
      value=getAttribute(node,"id")
      CALL context__create(TRIM(value))  
      CALL context__swap(TRIM(value))
    ENDIF
    
    IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
          
           CASE ('axis_definition') 
             CALL parsing_axis_group(child_node,axis_definition,root=.TRUE.)
               
           CASE ('grid_definition')
             CALL parsing_grid_group(child_node,grid_definition,root=.TRUE.)
               
           CASE ('field_definition')
             CALL parsing_field_group(child_node,field_definition,root=.TRUE.)
               
           CASE ('file_definition')
             CALL parsing_file_group(child_node,file_definition,root=.TRUE.)
               
           CASE DEFAULT
             IF (is_bad_node(node_name)) THEN
               WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing context'
               CALL Warning("mod_parse_xml:parsing_context")
             ENDIF   
        END SELECT
      ENDDO

    ENDIF

  END SUBROUTINE parsing_context


  RECURSIVE SUBROUTINE parsing_axis_group(node,parent,root)
  USE mod_axis_definition
  USE mod_axis_group
  USE mod_axis
  IMPLICIT NONE  
    TYPE(fnode), POINTER         :: node
    TYPE(axis_group),POINTER     :: parent
    LOGICAL,INTENT(IN),OPTIONAL :: root
    
    TYPE(axis_group),POINTER :: pt_axis_group
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    TYPE(axis),POINTER       :: attribute
    LOGICAL                  :: is_root
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
    
    is_root=.FALSE.
    IF (PRESENT(root)) is_root=root
  
    IF (is_root) THEN
      pt_axis_group=>parent
    ELSE  
      IF (is_attribute_exist(node,"id")) THEN
        value=getAttribute(node,"id")
        CALL axis_group__get_new_group(parent,pt_axis_group,TRIM(value))
      ELSE
        CALL axis_group__get_new_group(parent,pt_axis_group)
      ENDIF
    ENDIF
      
    CALL axis_group__get_default_attrib(pt_axis_group,attribute)
    CALL parsing_axis_attribute(node,attribute)
    
    IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
        
          CASE ('group') 
            CALL parsing_axis_group(child_node,pt_axis_group)
             
          CASE ('axis')
            CALL parsing_axis(child_node,pt_axis_group)

          CASE DEFAULT
            IF (is_bad_node(node_name)) THEN
              WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing group_axis'
              CALL Warning("mod_parse_xml:parsing_group_axis")
            ENDIF
        END SELECT
      ENDDO
    ENDIF

  END SUBROUTINE parsing_axis_group

  SUBROUTINE parsing_axis(node,parent)
  USE mod_axis_group
  USE mod_axis
  IMPLICIT NONE  
    TYPE(fnode), POINTER     :: node
    TYPE(axis_group),POINTER :: parent
    
    TYPE(axis),POINTER :: pt_axis
    TYPE(axis),POINTER       :: attribute
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
  
      
    IF (is_attribute_exist(node,"id")) THEN
      value=getAttribute(node,"id")
      CALL axis_group__get_new_axis(parent,pt_axis,TRIM(value))
    ELSE
      CALL axis_group__get_new_axis(parent,pt_axis)
    ENDIF
      
    CALL parsing_axis_attribute(node,pt_axis)
  
  END SUBROUTINE parsing_axis

  
  SUBROUTINE parsing_axis_attribute(node,pt_axis)
  USE mod_axis
  IMPLICIT NONE
    TYPE(fnode), POINTER     :: node
    TYPE(axis),POINTER :: pt_axis

    CHARACTER(len=100) :: value
   
    IF (is_attribute_exist(node,"name")) THEN
      value =  getAttribute(node,"name")
      CALL axis__set(pt_axis,name=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"description")) THEN
      value =  getAttribute(node,"description")
      CALL axis__set(pt_axis,description=TRIM(value))
    ENDIF
        
    IF (is_attribute_exist(node,"unit")) THEN
      value =  getAttribute(node,"unit")
      CALL axis__set(pt_axis,unit=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"size")) THEN
      value =  getAttribute(node,"size")
      CALL axis__set(pt_axis,a_size=string_to_integer(value))
    ENDIF

    IF (is_attribute_exist(node,"positive")) THEN
      value =  getAttribute(node,"positive")
      CALL axis__set(pt_axis,positive=string_to_logical(value))
    ENDIF
   
  END SUBROUTINE parsing_axis_attribute 








  RECURSIVE SUBROUTINE parsing_grid_group(node,parent,root)
  USE mod_grid_definition
  USE mod_grid_group
  USE mod_grid
  IMPLICIT NONE  
    TYPE(fnode), POINTER     :: node
    TYPE(grid_group),POINTER :: parent
    LOGICAL,INTENT(IN),OPTIONAL :: root

    TYPE(grid_group),POINTER :: pt_grid_group
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    TYPE(grid),POINTER       :: attribute
    LOGICAL                  :: is_root
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
  
    is_root=.FALSE.
    IF (PRESENT(root)) is_root=root
  
    IF (is_root) THEN
      pt_grid_group=>parent
    ELSE  
      IF (is_attribute_exist(node,"id")) THEN
        value=getAttribute(node,"id")
        CALL grid_group__get_new_group(parent,pt_grid_group,TRIM(value))
      ELSE
        CALL grid_group__get_new_group(parent,pt_grid_group)
      ENDIF
    ENDIF
      
    CALL grid_group__get_default_attrib(pt_grid_group,attribute)
    CALL parsing_grid_attribute(node,attribute)
    
    IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
        
          CASE ('group') 
            CALL parsing_grid_group(child_node,pt_grid_group)
             
          CASE ('grid')
            CALL parsing_grid(child_node,pt_grid_group)

          CASE DEFAULT
            IF (is_bad_node(node_name)) THEN
              WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing group_grid'
              CALL Warning("mod_parse_xml:parsing_group_grid")
            ENDIF
        END SELECT
      ENDDO
    ENDIF

  END SUBROUTINE parsing_grid_group

  SUBROUTINE parsing_grid(node,parent)
  USE mod_grid_group
  USE mod_grid
  IMPLICIT NONE  
    TYPE(fnode), POINTER     :: node
    TYPE(grid_group),POINTER :: parent
    
    TYPE(grid),POINTER       :: pt_grid
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
  
      
    IF (is_attribute_exist(node,"id")) THEN
      value=getAttribute(node,"id")
      CALL grid_group__get_new_grid(parent,pt_grid,TRIM(value))
    ELSE
      CALL grid_group__get_new_grid(parent,pt_grid)
    ENDIF
      
    CALL parsing_grid_attribute(node,pt_grid)
    

   IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
        
          CASE ('zoom') 
            CALL parsing_zoom(child_node,pt_grid)
             
          CASE DEFAULT
            IF (is_bad_node(node_name)) THEN
              WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing grid'
              CALL Warning("mod_parse_xml:parsing_grid")
            ENDIF
        END SELECT
      ENDDO
    ENDIF
  
  END SUBROUTINE parsing_grid

  SUBROUTINE parsing_grid_attribute(node,pt_grid)
  USE mod_grid
  IMPLICIT NONE
    TYPE(fnode), POINTER     :: node
    TYPE(grid),POINTER :: pt_grid

    CHARACTER(len=100) :: value
   
    IF (is_attribute_exist(node,"name")) THEN
      value =  getAttribute(node,"name")
      CALL grid__set(pt_grid,name=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"description")) THEN
      value =  getAttribute(node,"description")
      CALL grid__set(pt_grid,description=TRIM(value))
    ENDIF
        
  END SUBROUTINE parsing_grid_attribute 

  SUBROUTINE parsing_zoom(node,parent)
  USE mod_zoom
  USE mod_grid 
  IMPLICIT NONE  
    TYPE(fnode), POINTER     :: node
    TYPE(grid),POINTER       :: parent
    
    TYPE(zoom),POINTER       :: pt_zoom
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
  
      
    IF (is_attribute_exist(node,"id")) THEN
      value=getAttribute(node,"id")
      CALL grid__get_new_zoom(parent,pt_zoom,TRIM(value))
    ELSE
      CALL grid__get_new_zoom(parent,pt_zoom)
    ENDIF
      
    CALL parsing_zoom_attribute(node,pt_zoom)
    
  END SUBROUTINE parsing_zoom


  SUBROUTINE parsing_zoom_attribute(node,pt_zoom)
  USE mod_zoom
  IMPLICIT NONE
    TYPE(fnode), POINTER     :: node
    TYPE(zoom),POINTER       :: pt_zoom

    CHARACTER(len=100) :: value
   
    IF (is_attribute_exist(node,"name")) THEN
      value =  getAttribute(node,"name")
      CALL zoom__set(pt_zoom,name=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"description")) THEN
      value =  getAttribute(node,"description")
      CALL zoom__set(pt_zoom,description=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"ni")) THEN
      value =  getAttribute(node,"ni")
      CALL zoom__set(pt_zoom,ni_glo=string_to_integer(value))
    ENDIF

    IF (is_attribute_exist(node,"nj")) THEN
      value =  getAttribute(node,"nj")
      CALL zoom__set(pt_zoom,nj_glo=string_to_integer(value))
    ENDIF

    IF (is_attribute_exist(node,"ibegin")) THEN
      value =  getAttribute(node,"ibegin")
      CALL zoom__set(pt_zoom,ibegin_glo=string_to_integer(value))
    ENDIF

    IF (is_attribute_exist(node,"jbegin")) THEN
      value =  getAttribute(node,"jbegin")
      CALL zoom__set(pt_zoom,jbegin_glo=string_to_integer(value))
    ENDIF
        
  END SUBROUTINE parsing_zoom_attribute 



  RECURSIVE SUBROUTINE parsing_field_group(node,parent,root)
  USE mod_field_definition
  USE mod_field_group
  USE mod_field
  IMPLICIT NONE  
    TYPE(fnode), POINTER          :: node
    TYPE(field_group),POINTER     :: parent
    LOGICAL,INTENT(IN),OPTIONAL   :: root
    
    TYPE(field_group),POINTER :: pt_field_group
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    TYPE(field),POINTER       :: attribute
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
    LOGICAL            :: is_root
    
    is_root=.FALSE.
    IF (PRESENT(root)) is_root=root
    
    IF (is_root) THEN
      pt_field_group=>parent
    ELSE  
      IF (is_attribute_exist(node,"id")) THEN
        value=getAttribute(node,"id")
        CALL field_group__get_new_group(parent,pt_field_group,TRIM(value))
      ELSE
        CALL field_group__get_new_group(parent,pt_field_group)
      ENDIF
    ENDIF
      
    CALL field_group__get_default_attrib(pt_field_group,attribute)
    CALL parsing_field_attribute(node,attribute)
    
    IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
        
          CASE ('group') 
            CALL parsing_field_group(child_node,pt_field_group)
             
          CASE ('field')
            CALL parsing_field(child_node,pt_field_group)

          CASE DEFAULT
            IF (is_bad_node(node_name)) THEN
              WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing group_field'
              CALL Warning("mod_parse_xml:parsing_group_field")
            ENDIF
        END SELECT
      ENDDO
    ENDIF

  END SUBROUTINE parsing_field_group

  SUBROUTINE parsing_field(node,parent)
  USE mod_field_group
  USE mod_field
  IMPLICIT NONE  
    TYPE(fnode), POINTER     :: node
    TYPE(field_group),POINTER :: parent
    
    TYPE(field),POINTER :: pt_field
    TYPE(field),POINTER       :: attribute
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
  
      
    IF (is_attribute_exist(node,"id")) THEN
      value=getAttribute(node,"id")
      CALL field_group__get_new_field(parent,pt_field,TRIM(value))
    ELSE
      CALL field_group__get_new_field(parent,pt_field)
    ENDIF
      
    CALL parsing_field_attribute(node,pt_field)
  
  END SUBROUTINE parsing_field

  
  SUBROUTINE parsing_field_attribute(node,pt_field)
  USE mod_field
  IMPLICIT NONE
    TYPE(fnode), POINTER     :: node
    TYPE(field),POINTER :: pt_field

    CHARACTER(len=100) :: value
   
    IF (is_attribute_exist(node,"name")) THEN
      value =  getAttribute(node,"name")
      CALL field__set(pt_field,name=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"description")) THEN
      value =  getAttribute(node,"description")
      CALL field__set(pt_field,description=TRIM(value))
    ENDIF
        
    IF (is_attribute_exist(node,"unit")) THEN
      value =  getAttribute(node,"unit")
      CALL field__set(pt_field,unit=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"operation")) THEN
      value =  getAttribute(node,"operation")
      CALL field__set(pt_field,operation=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"freq_op")) THEN
      value =  getAttribute(node,"freq_op")
      CALL field__set(pt_field,freq_op=string_to_integer(value))
    ENDIF
   
    IF (is_attribute_exist(node,"axis_ref")) THEN
      value =  getAttribute(node,"axis_ref")
      CALL field__set(pt_field,axis_ref=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"grid_ref")) THEN
      value =  getAttribute(node,"grid_ref")
      CALL field__set(pt_field,grid_ref=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"zoom_ref")) THEN
      value =  getAttribute(node,"zoom_ref")
      CALL field__set(pt_field,zoom_ref=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"level")) THEN
      value =  getAttribute(node,"level")
      CALL field__set(pt_field,level=string_to_integer(value))
    ENDIF
    
    IF (is_attribute_exist(node,"prec")) THEN
      value =  getAttribute(node,"prec")
      CALL field__set(pt_field,prec=string_to_integer(value))
    ENDIF

    IF (is_attribute_exist(node,"ref")) THEN
      value =  getAttribute(node,"ref")
      CALL field__set(pt_field,ref=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"enabled")) THEN
      value =  getAttribute(node,"enabled")
      CALL field__set(pt_field,enabled=string_to_logical(value))
    ENDIF
  

  END SUBROUTINE parsing_field_attribute 


  RECURSIVE SUBROUTINE parsing_file_group(node,parent,root)
  USE mod_file_definition
  USE mod_file_group
  USE mod_file
  IMPLICIT NONE  
    TYPE(fnode), POINTER          :: node
    TYPE(file_group),POINTER      :: parent
    LOGICAL,INTENT(IN),OPTIONAL   :: root

    TYPE(file_group),POINTER :: pt_file_group
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    TYPE(file),POINTER       :: attribute
    LOGICAL                  :: is_root
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
    
    is_root=.FALSE.
    IF (PRESENT(root)) is_root=root
  
    IF (is_root) THEN
      pt_file_group=>parent
    ELSE  
      IF (is_attribute_exist(node,"id")) THEN
        value=getAttribute(node,"id")
        CALL file_group__get_new_group(parent,pt_file_group,TRIM(value))
      ELSE
        CALL file_group__get_new_group(parent,pt_file_group)
      ENDIF
    ENDIF
      
    CALL file_group__get_default_attrib(pt_file_group,attribute)
    CALL parsing_file_attribute(node,attribute)
    
    IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
        
          CASE ('group') 
            CALL parsing_file_group(child_node,pt_file_group)
             
          CASE ('file')
            CALL parsing_file(child_node,pt_file_group)

          CASE DEFAULT
            IF (is_bad_node(node_name)) THEN
              WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing group_file'
              CALL Warning("mod_parse_xml:parsing_group_file")
            ENDIF
        END SELECT
      ENDDO
    ENDIF

  END SUBROUTINE parsing_file_group

  SUBROUTINE parsing_file(node,parent)
  USE mod_file_group
  USE mod_file
  IMPLICIT NONE  
    TYPE(fnode), POINTER     :: node
    TYPE(file_group),POINTER :: parent
    
    TYPE(file),POINTER :: pt_file
    TYPE(file),POINTER       :: attribute
    TYPE(fnode), POINTER     :: child_node
    TYPE(fnodeList), POINTER :: child_list
    INTEGER :: il
    CHARACTER(len=100) :: node_name
    CHARACTER(len=100) :: value
  
      
    IF (is_attribute_exist(node,"id")) THEN
      value=getAttribute(node,"id")
      CALL file_group__get_new_file(parent,pt_file,TRIM(value))
    ELSE
      CALL file_group__get_new_file(parent,pt_file)
    ENDIF
      
    CALL parsing_file_attribute(node,pt_file)
  
    IF (hasChildNodes(node)) THEN
      child_list => getChildnodes(node)

      DO il=0,getLength(child_list)-1
        child_node => item(child_list,il)
        node_name=getNodename(child_node)
        
        SELECT CASE (TRIM(node_name)) 
        
          CASE ('group') 
            CALL parsing_field_group(child_node,pt_file%field_list)
             
          CASE ('field')
            CALL parsing_field(child_node,pt_file%field_list)

          CASE DEFAULT
            IF (is_bad_node(node_name)) THEN
              WRITE(message,*) 'Unknown node <<',TRIM(node_name),'>> while parsing group_file'
              CALL Warning("mod_parse_xml:parsing_group_file")
            ENDIF
        END SELECT
      ENDDO
    ENDIF
  END SUBROUTINE parsing_file

  
  SUBROUTINE parsing_file_attribute(node,pt_file)
  USE mod_file
  IMPLICIT NONE
    TYPE(fnode), POINTER     :: node
    TYPE(file),POINTER :: pt_file

    CHARACTER(len=100) :: value
   
    IF (is_attribute_exist(node,"name")) THEN
      value =  getAttribute(node,"name")
      CALL file__set(pt_file,name=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"name_suffix")) THEN
      value =  getAttribute(node,"name_suffix")
      CALL file__set(pt_file,name_suffix=TRIM(value))
    ENDIF

    IF (is_attribute_exist(node,"description")) THEN
      value =  getAttribute(node,"description")
      CALL file__set(pt_file,description=TRIM(value))
    ENDIF
        
    IF (is_attribute_exist(node,"output_freq")) THEN
      value =  getAttribute(node,"output_freq")
      CALL file__set(pt_file,output_freq=string_to_integer(value))
    ENDIF
   
    IF (is_attribute_exist(node,"output_level")) THEN
      value =  getAttribute(node,"output_level")
      CALL file__set(pt_file,output_level=string_to_integer(value))
    ENDIF

    IF (is_attribute_exist(node,"enabled")) THEN
      value =  getAttribute(node,"enabled")
      CALL file__set(pt_file,enabled=string_to_logical(value))
    ENDIF

  END SUBROUTINE parsing_file_attribute
  

   
  
  FUNCTION is_attribute_exist(node, name)

    LOGICAL :: is_attribute_exist
    TYPE(fnode), POINTER :: node
    CHARACTER(len=*) :: name
    CHARACTER(len=100) :: value

    value=""
    
    is_attribute_exist= .false.
    value=getAttribute(node, TRIM(name))
    IF (value .NE. "") is_attribute_exist= .true.

  END FUNCTION is_attribute_exist

  FUNCTION is_bad_node(node_name)
  IMPLICIT NONE
    CHARACTER(len=*),INTENT(IN) :: node_name
    LOGICAL                     :: is_bad_node
    
    IF (TRIM(node_name)=='#text' .OR. TRIM(node_name)=='#comment') THEN
      is_bad_node=.FALSE.
    ELSE
      is_bad_node=.TRUE.
    ENDIF
      
  END FUNCTION is_bad_node
  
END MODULE parsing_xml
