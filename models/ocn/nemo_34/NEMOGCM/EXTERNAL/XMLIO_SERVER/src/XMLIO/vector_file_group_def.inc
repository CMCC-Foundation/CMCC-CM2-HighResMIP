
  TYPE Pt_file_group__ 
    TYPE (file_group),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_file_group__

  TYPE vector_file_group
  
    TYPE (Pt_file_group__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_file_group
  
