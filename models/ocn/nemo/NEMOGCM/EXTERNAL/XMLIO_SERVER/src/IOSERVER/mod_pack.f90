  MODULE mod_pack
    INTEGER,PARAMETER :: integer_size=4
    INTEGER,PARAMETER :: pack_buffer_max_size=integer_size*1024*1024
    INTEGER(kind=8),POINTER,SAVE      :: pack_buffer(:)
    INTEGER,SAVE      :: pack_pos
    
    INTERFACE pack_data
      MODULE PROCEDURE pack_r,pack_r1,pack_r2,pack_r3,pack_r4,                 &
                       pack_i,pack_i1,pack_i2,pack_i3,pack_i4,                 &
                       pack_l,pack_l1,pack_l2,pack_l3,pack_l4,                 &
                       pack_c,pack_c1,pack_c2,pack_c3,pack_c4,                 &
                       pack_attr
    END INTERFACE pack_data

    INTERFACE unpack_data
      MODULE PROCEDURE unpack_r,unpack_r1,unpack_r2,unpack_r3,unpack_r4,       &
                       unpack_i,unpack_i1,unpack_i2,unpack_i3,unpack_i4,       &
                       unpack_l,unpack_l1,unpack_l2,unpack_l3,unpack_l4,       &
                       unpack_c,unpack_c1,unpack_c2,unpack_c3,unpack_c4,       &
                       unpack_attr
    END INTERFACE unpack_data

    INTERFACE pack_field
      MODULE PROCEDURE pack_field1,pack_field2,pack_field3,pack_field4        
    END INTERFACE pack_field

    INTERFACE unpack_field
      MODULE PROCEDURE unpack_field1,unpack_field2,unpack_field3,unpack_field4 
    END INTERFACE unpack_field
           
  CONTAINS
  
    SUBROUTINE set_pack_buffer(pack_buffer0,pack_pos0)
    IMPLICIT NONE
    INTEGER(kind=8),POINTER  :: pack_buffer0(:)
    INTEGER                  :: pack_pos0
    
      pack_buffer=>pack_buffer0
      pack_pos=pack_pos0

    END SUBROUTINE set_pack_buffer
    
    
    SUBROUTINE pack_r(arg)
      IMPLICIT NONE
        REAL :: arg
         
        CALL packf_r(pack_buffer(pack_pos),pack_pos,arg,1)
     
     END SUBROUTINE pack_r
     
     SUBROUTINE pack_r1(arg)
      IMPLICIT NONE
        REAL :: arg(:)
        
        CALL packf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_r1

    
     SUBROUTINE pack_r2(arg)
      IMPLICIT NONE
        REAL :: arg(:,:)
        
        CALL packf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_r2
     
     
     SUBROUTINE pack_r3(arg)
      IMPLICIT NONE
        REAL :: arg(:,:,:)
        
        CALL packf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_r3

     
     SUBROUTINE pack_r4(arg)
      IMPLICIT NONE
        REAL :: arg(:,:,:,:)
        
        CALL packf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_r4


    SUBROUTINE unpack_r(arg)
      IMPLICIT NONE
        REAL :: arg
         
        CALL unpackf_r(pack_buffer(pack_pos),pack_pos,arg,1)
     
     END SUBROUTINE unpack_r
     
     SUBROUTINE unpack_r1(arg)
      IMPLICIT NONE
        REAL :: arg(:)
        
        CALL unpackf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_r1

    
     SUBROUTINE unpack_r2(arg)
      IMPLICIT NONE
        REAL :: arg(:,:)
        
        CALL unpackf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_r2
     
     
     SUBROUTINE unpack_r3(arg)
      IMPLICIT NONE
        REAL :: arg(:,:,:)
        
        CALL unpackf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_r3

     
     SUBROUTINE unpack_r4(arg)
      IMPLICIT NONE
        REAL :: arg(:,:,:,:)
        
        CALL unpackf_r(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_r4




    SUBROUTINE pack_i(arg)
      IMPLICIT NONE
        INTEGER :: arg
         
        CALL packf_i(pack_buffer(pack_pos),pack_pos,arg,1)
     
     END SUBROUTINE pack_i
     
     SUBROUTINE pack_i1(arg)
      IMPLICIT NONE
        INTEGER :: arg(:)
        
        CALL packf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_i1

    
     SUBROUTINE pack_i2(arg)
      IMPLICIT NONE
        INTEGER :: arg(:,:)
        
        CALL packf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_i2
     
     
     SUBROUTINE pack_i3(arg)
      IMPLICIT NONE
        INTEGER :: arg(:,:,:)
        
        CALL packf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_i3

     
     SUBROUTINE pack_i4(arg)
      IMPLICIT NONE
        INTEGER :: arg(:,:,:,:)
        
        CALL packf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_i4
     

    SUBROUTINE unpack_i(arg)
      IMPLICIT NONE
        INTEGER :: arg
         
        CALL unpackf_i(pack_buffer(pack_pos),pack_pos,arg,1)
     
     END SUBROUTINE unpack_i
     
     SUBROUTINE unpack_i1(arg)
      IMPLICIT NONE
        INTEGER :: arg(:)
        
        CALL unpackf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_i1

    
     SUBROUTINE unpack_i2(arg)
      IMPLICIT NONE
        INTEGER :: arg(:,:)
        
        CALL unpackf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_i2
     
     
     SUBROUTINE unpack_i3(arg)
      IMPLICIT NONE
        INTEGER :: arg(:,:,:)
        
        CALL unpackf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_i3

     
     SUBROUTINE unpack_i4(arg)
      IMPLICIT NONE
        INTEGER :: arg(:,:,:,:)
        
        CALL unpackf_i(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_i4     







    SUBROUTINE pack_l(arg)
      IMPLICIT NONE
        LOGICAL :: arg
         
        CALL packf_l(pack_buffer(pack_pos),pack_pos,arg,1)
     
     END SUBROUTINE pack_l
     
     SUBROUTINE pack_l1(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:)
        
        CALL packf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_l1

    
     SUBROUTINE pack_l2(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:,:)
        
        CALL packf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_l2
     
     
     SUBROUTINE pack_l3(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:,:,:)
        
        CALL packf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_l3

     
     SUBROUTINE pack_l4(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:,:,:,:)
        
        CALL packf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE pack_l4
     

    SUBROUTINE unpack_l(arg)
      IMPLICIT NONE
        LOGICAL :: arg
         
        CALL unpackf_l(pack_buffer(pack_pos),pack_pos,arg,1)
     
     END SUBROUTINE unpack_l
     
     SUBROUTINE unpack_l1(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:)
        
        CALL unpackf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_l1

    
     SUBROUTINE unpack_l2(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:,:)
        
        CALL unpackf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_l2
     
     
     SUBROUTINE unpack_l3(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:,:,:)
        
        CALL unpackf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_l3

     
     SUBROUTINE unpack_l4(arg)
      IMPLICIT NONE
        LOGICAL :: arg(:,:,:,:)
        
        CALL unpackf_l(pack_buffer(pack_pos),pack_pos,arg,size(arg))
        
     END SUBROUTINE unpack_l4     




   
     
   SUBROUTINE pack_c(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg
         
        CALL packf_c(pack_buffer(pack_pos),pack_pos,arg,len(arg))
     
     END SUBROUTINE pack_c
     
     SUBROUTINE pack_c1(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:)
        
        CALL packf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1)))
        
     END SUBROUTINE pack_c1

    
     SUBROUTINE pack_c2(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:,:)
        
        CALL packf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1,1)))
        
     END SUBROUTINE pack_c2
     
     
     SUBROUTINE pack_c3(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:,:,:)
        
        CALL packf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1,1,1)))
        
     END SUBROUTINE pack_c3

     
     SUBROUTINE pack_c4(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:,:,:,:)
        
        CALL packf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1,1,1,1)))
        
     END SUBROUTINE pack_c4
     

    SUBROUTINE unpack_c(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg
         
        CALL unpackf_c(pack_buffer(pack_pos),pack_pos,arg,len(arg))
     
     END SUBROUTINE unpack_c
     
     SUBROUTINE unpack_c1(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:)
        
        CALL unpackf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1)))
        
     END SUBROUTINE unpack_c1

    
     SUBROUTINE unpack_c2(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:,:)
        
        CALL unpackf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1,1)))
        
     END SUBROUTINE unpack_c2
     
     
     SUBROUTINE unpack_c3(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:,:,:)
        
        CALL unpackf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1,1,1)))
        
     END SUBROUTINE unpack_c3

     
     SUBROUTINE unpack_c4(arg)
      IMPLICIT NONE
        CHARACTER(len=*) :: arg(:,:,:,:)
        
        CALL unpackf_c(pack_buffer(pack_pos),pack_pos,arg,size(arg)*len(arg(1,1,1,1)))
        
     END SUBROUTINE unpack_c4     
     




     SUBROUTINE pack_field1(arg)
     IMPLICIT NONE
       REAL :: arg(:)
     
       CALL packf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE pack_field1


     SUBROUTINE pack_field2(arg)
     IMPLICIT NONE
       REAL :: arg(:,:)
       
       CALL packf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE pack_field2
     
     SUBROUTINE pack_field3(arg)
     IMPLICIT NONE
       REAL :: arg(:,:,:)
     
       CALL packf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE pack_field3

     SUBROUTINE pack_field4(arg)
     IMPLICIT NONE
       REAL :: arg(:,:,:,:)
     
       CALL packf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE pack_field4



     SUBROUTINE unpack_field1(arg)
     IMPLICIT NONE
       REAL :: arg(:)
     
       CALL unpackf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE unpack_field1


     SUBROUTINE unpack_field2(arg)
     IMPLICIT NONE
       REAL :: arg(:,:)
     
       CALL unpackf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE unpack_field2
     
     SUBROUTINE unpack_field3(arg)
     IMPLICIT NONE
       REAL :: arg(:,:,:)
     
       CALL unpackf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE unpack_field3

     SUBROUTINE unpack_field4(arg)
     IMPLICIT NONE
       REAL :: arg(:,:,:,:)
     
       CALL unpackf_field(pack_buffer(pack_pos),pack_pos,arg,size(arg)) 
     END SUBROUTINE unpack_field4
      
     SUBROUTINE pack_attr(attrib)
     USE mod_attribut
     USE mod_stdtype
     IMPLICIT NONE
       TYPE(attribut) :: attrib

       CALL pack_data(attrib%object)
       CALL pack_data(attrib%name)
       CALL pack_data(attrib%type)
       CALL pack_data(attrib%dim)
       CALL pack_data(attrib%ndim)
       CALL pack_data(attrib%string_len)
              
       SELECT CASE(attrib%type)
         CASE (integer0)
           CALL pack_data(attrib%integer0_ptr)
         CASE (integer1)
           CALL pack_data(attrib%integer1_ptr)
         CASE (integer2)
           CALL pack_data(attrib%integer2_ptr)
         CASE (real0)
           CALL pack_data(attrib%real0_ptr)
         CASE (real1)
           CALL pack_data(attrib%real1_ptr)
         CASE (real2)
           CALL pack_data(attrib%real2_ptr)
         CASE (logical0)
           CALL pack_data(attrib%logical0_ptr)
         CASE (logical1)
           CALL pack_data(attrib%logical1_ptr)
         CASE (logical2)
           CALL pack_data(attrib%logical2_ptr)
         CASE (string0)
           CALL pack_string0(attrib%string0_ptr)
         CASE (string1)
           CALL pack_string1(attrib%string1_ptr)
         CASE (string2)
           CALL pack_data(attrib%string2_ptr)
       END SELECT

     CONTAINS

       SUBROUTINE pack_string0(str)
         CHARACTER(LEN=attrib%string_len) ::str
           CALL pack_data(str)
       END SUBROUTINE

       SUBROUTINE pack_string1(str)
         CHARACTER(LEN=attrib%string_len) ::str(:)
           CALL pack_data(str)
       END SUBROUTINE
       
       SUBROUTINE pack_string2(str)
         CHARACTER(LEN=attrib%string_len) ::str(:,:)
           CALL pack_data(str)
       END SUBROUTINE

     END SUBROUTINE pack_attr  

     SUBROUTINE unpack_attr(attrib)
     USE mod_attribut
     USE mod_stdtype
     IMPLICIT NONE
       TYPE(attribut) :: attrib

       CALL unpack_data(attrib%object)
       CALL unpack_data(attrib%name)
       CALL unpack_data(attrib%type)
       CALL unpack_data(attrib%dim)
       CALL unpack_data(attrib%ndim)
       CALL unpack_data(attrib%string_len)
              
       SELECT CASE(attrib%type)
         CASE (integer0)
!!$           ALLOCATE(attrib%integer0_ptr)
           CALL unpack_data(attrib%integer0_ptr)
         CASE (integer1)
           ALLOCATE(attrib%integer1_ptr(attrib%dim(1)))
           CALL unpack_data(attrib%integer1_ptr)
         CASE (integer2)
           ALLOCATE(attrib%integer2_ptr(attrib%dim(1),attrib%dim(2)))
           CALL unpack_data(attrib%integer2_ptr)
         CASE (real0)
!!$           ALLOCATE(attrib%real0_ptr)
           CALL unpack_data(attrib%real0_ptr)
         CASE (real1)
           ALLOCATE(attrib%real1_ptr(attrib%dim(1)))
           CALL unpack_data(attrib%real1_ptr)
         CASE (real2)
           ALLOCATE(attrib%real2_ptr(attrib%dim(1),attrib%dim(2)))
         CASE (logical0)
!!$           ALLOCATE(attrib%logical0_ptr)
           CALL unpack_data(attrib%logical0_ptr)
         CASE (logical1)
           ALLOCATE(attrib%logical1_ptr(attrib%dim(1)))
           CALL unpack_data(attrib%logical1_ptr)
         CASE (logical2)
           ALLOCATE(attrib%logical2_ptr(attrib%dim(1),attrib%dim(2)))
           CALL unpack_data(attrib%logical2_ptr)
         CASE (string0)
!!$           ALLOCATE(attrib%string0_ptr)
           CALL unpack_string0
         CASE (string1)
           ALLOCATE(attrib%string1_ptr(attrib%dim(1)))
           CALL unpack_string1
         CASE (string2)
           ALLOCATE(attrib%string2_ptr(attrib%dim(1),attrib%dim(2)))
           CALL unpack_string2
       END SELECT

     CONTAINS

       SUBROUTINE unpack_string0
         CHARACTER(LEN=attrib%string_len) ::str
           CALL unpack_data(str)
           attrib%string0_ptr=str
       END SUBROUTINE

       SUBROUTINE unpack_string1
         CHARACTER(LEN=attrib%string_len) ::str(attrib%dim(1))
           CALL unpack_data(str)
           attrib%string1_ptr=str
       END SUBROUTINE
       
       SUBROUTINE unpack_string2
         CHARACTER(LEN=attrib%string_len) ::str(attrib%dim(1),attrib%dim(2))
           CALL unpack_data(str)
           attrib%string2_ptr=str
       END SUBROUTINE

     END SUBROUTINE unpack_attr              

     
  END MODULE mod_pack
