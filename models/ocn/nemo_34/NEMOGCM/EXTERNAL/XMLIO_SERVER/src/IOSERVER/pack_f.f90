SUBROUTINE packf_r(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    REAL,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    REAL,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size

    tot_size=var_size*REAL_SIZE

    DO i=1,var_size
      buffer(i)=var(i)
    ENDDO

    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1

END SUBROUTINE packf_r


SUBROUTINE packf_field(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    REAL,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    REAL,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*REAL_SIZE
        
    DO i=1,var_size
      buffer(i)=var(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE packf_field

SUBROUTINE packf_i(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    INTEGER,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    INTEGER,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*INTEGER_SIZE
        
    DO i=1,var_size
      buffer(i)=var(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE packf_i

SUBROUTINE packf_l(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    LOGICAL,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    LOGICAL,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*LOGICAL_SIZE
        
    DO i=1,var_size
      buffer(i)=var(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE packf_l

SUBROUTINE packf_c(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER                  :: var_size
    CHARACTER(LEN=var_size)  :: buffer
    INTEGER                  :: position
    CHARACTER(LEN=var_size)  :: var
   
    INTEGER :: tot_size
    
    tot_size=var_size*CHARACTER_SIZE
        
    buffer(1:var_size)=var(1:var_size)
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE packf_c






SUBROUTINE unpackf_r(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    REAL,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    REAL,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*REAL_SIZE
        
    DO i=1,var_size
      var(i)=buffer(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE unpackf_r


SUBROUTINE unpackf_field(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    REAL,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    REAL,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*REAL_SIZE
        
    DO i=1,var_size
      var(i)=buffer(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE unpackf_field

SUBROUTINE unpackf_i(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    INTEGER,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    INTEGER,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*INTEGER_SIZE
        
    DO i=1,var_size
      var(i)=buffer(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE unpackf_i

SUBROUTINE unpackf_l(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER           :: var_size
    LOGICAL,DIMENSION(var_size) :: buffer
    INTEGER           :: position
    LOGICAL,DIMENSION(var_size) :: var
   
    INTEGER :: i
    INTEGER :: tot_size
    
    tot_size=var_size*LOGICAL_SIZE
        
    DO i=1,var_size
      var(i)=buffer(i)
    ENDDO
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE unpackf_l

SUBROUTINE unpackf_c(buffer,position,var,var_size)
  USE mod_prec
  IMPLICIT NONE
    INTEGER                  :: var_size
    CHARACTER(LEN=var_size)  :: buffer
    INTEGER                  :: position
    CHARACTER(LEN=var_size)  :: var
    
    INTEGER :: tot_size
    
    tot_size=var_size*CHARACTER_SIZE
        
    var(1:var_size)=buffer(1:var_size)
          
    position=position+tot_size/BLOCK_SIZE
    IF (MOD(tot_size,BLOCK_SIZE)/=0) position=position+1
    
END SUBROUTINE unpackf_c
