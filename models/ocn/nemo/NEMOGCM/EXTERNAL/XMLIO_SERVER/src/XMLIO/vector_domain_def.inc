
  TYPE Pt_domain__ 
    TYPE (domain),POINTER :: Pt
    LOGICAL               :: owned
  END TYPE Pt_domain__

  TYPE vector_domain
  
    TYPE (Pt_domain__), POINTER, DIMENSION(:) :: at => NULL()
    INTEGER :: size
    
    INTEGER :: size_max
    REAL    :: grow_factor
    INTEGER :: chunk
  
  END TYPE vector_domain
  
