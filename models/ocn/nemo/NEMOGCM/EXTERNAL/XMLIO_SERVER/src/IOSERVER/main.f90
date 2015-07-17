PROGRAM ioslave
USE mod_pack, ONLY : set_pack_buffer,pack_data,unpack_data
USE mod_wait
IMPLICIT NONE

  INTEGER,DIMENSION(3) :: a=(/1,2,3/)
  REAL                 :: b = 4.5
  REAL                 :: c = 6.7
  INTEGER,DIMENSION(3) :: aout
  REAL                 :: bout,cout
  LOGICAL,DIMENSION(5) :: d=(/.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE./)
  LOGICAL,DIMENSION(5) :: dout
  CHARACTER            :: e='C'
  CHARACTER            :: eout
  CHARACTER(len=11)    :: f(3)
  CHARACTER(len=11)    :: fout(3)
  INTEGER(kind=8),dimension(:),POINTER :: buffer
  REAL,dimension(:),POINTER :: buffer_field
  DOUBLE PRECISION :: t
  
  ALLOCATE(buffer(1024))
  ALLOCATE(buffer_field(1024))
  
  CALL set_pack_buffer(buffer,1)
  
  f(1)="COUCOU1"
  f(2)="COUCOU2"
  f(3)="COUCOU3"
  
  CALL pack_data(b)
  CALL pack_data(f)
  CALL pack_data(a)
  CALL pack_data(d)
  CALL pack_data(e)
  CALL pack_data(c)
  
  CALL set_pack_buffer(buffer,1)
  CALL unpack_data(bout)
  CALL unpack_data(fout)
  CALL unpack_data(aout)
  CALL unpack_data(dout)
  CALL unpack_data(eout)
  CALL unpack_data(cout)
  
  PRINT *,a,b,c,d,e,"  ",f
  PRINT *,aout,bout,cout,dout,eout,"  ",fout

  CALL Init_wait
  
  t=top()
  CALL wait_us(10)
  t=top()
  
  PRINT *,"ATTENTE 10 us",t

  t=top()
  CALL wait_us(250)
  t=top()
  
  PRINT *,"ATTENTE 250 us",t

  t=top()
  CALL wait_us(1235)
  t=top()
  
  PRINT *,"ATTENTE 1235 us",t


END PROGRAM ioslave

