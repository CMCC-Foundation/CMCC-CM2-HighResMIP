
  TYPE Pt_file_dep__ 
    TYPE (file_dep),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_file_dep__

  TYPE vector_file_dep
  
    TYPE (Pt_file_dep__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_file_dep
  
