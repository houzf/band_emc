  MODULE  LATTICEVECTOR
  use prec
  implicit none
  real(dp), dimension(3,3):: prim_lat, recip_lat
  real(dp), dimension(3,3):: inv_prim_lat, inv_recip_lat, tran_inv_prim_lat

  END MODULE


SUBROUTINE readposcar()
  use prec
  use constant
  use latticevector
  implicit none
  real(dp), dimension(3,3):: scaled_prim_lat
  real(dp)  ascale
  integer i, j, error_flag
  character(len=256) :: title
  character(len=256) :: filename
  logical :: lexist 
  real(dp):: a1m, a2m, a3m, a, b, c, alpha, beta, gama

  filename ="POSCAR"
  inquire(file=filename,exist=lexist)
  if ( .not.  lexist) then
       write(*,'(3A)') "Error:",trim(filename), "does not exist"
       STOP
  end if
  open(10,file=filename,status='old')
  read(10,*) title
  read(10,*) ascale
  do i = 1, 3
      read(10,*) (scaled_prim_lat(i,j), j = 1, 3)
      !! Unit converted from Ang to Bohr
      prim_lat(i,:) = ascale / b2a * scaled_prim_lat(i,:)
  end do
  close(10)
! From the relationship of reciprocal lattice vector and
! primitive lattice vector, we can obtain the reciprocal
! lattice vector.
! http://en.wikipedia.org/wiki/Reciprocal_lattice
! [b1,b2,b3]^T = 2 * pi [a1,a2,a3]^{-1}
  call findinv(prim_lat,inv_prim_lat, 3, error_flag)
  tran_inv_prim_lat=TRANSPOSE(inv_prim_lat)
  recip_lat= 2.0d0 * pi * tran_inv_prim_lat 

  call findinv(recip_lat,inv_recip_lat, 3, error_flag)

   open(7,file='lattice_info.out')
   write(7,*) 'Unit cell lattice vectors (in Bohr):' 
   do i = 1, 3
      write(7,'(3F15.8)'), (prim_lat(i,j), j=1,3) 
   end do
   write(7,*) 'Reciprocal lattice vectors containing 2*pi (in 1/Bohr):' 
   do i = 1, 3
      write(7,'(3F15.8)'), (recip_lat(i,j), j=1,3) 
   end do

   write(7,*), 'Inverse of reciprocal lattice vectors:' 
   do i = 1, 3
      write(7,'(3F15.8)'), (inv_recip_lat(i,j), j=1,3) 
   end do
   
   a1m = sqrt(DOT_PRODUCT(prim_lat(1,:),prim_lat(1,:)))
   a2m = sqrt(DOT_PRODUCT(prim_lat(2,:),prim_lat(2,:)))
   a3m = sqrt(DOT_PRODUCT(prim_lat(3,:),prim_lat(3,:)))
   alpha = acos(DOT_PRODUCT(prim_lat(2,:),prim_lat(3,:))/(a2m * a3m))
   beta  = acos(DOT_PRODUCT(prim_lat(1,:),prim_lat(3,:))/(a1m * a3m))
   gama  = acos(DOT_PRODUCT(prim_lat(1,:),prim_lat(2,:))/(a1m * a2m))
   write(7,*) 'Lengths of lattice vecotrs (in Angstrom):'
   write(7,'(3F15.8)') a1m * b2a, a2m * b2a, a3m *b2a
   write(7,*) 'Aangles between two of lattice vecotrs (in degree):'
   write(7,*) 'alpha(a2,a3), beta(a1,a3), gamma(a1,a2)'
   write(7,'(3F15.8)') alpha/pi*180, beta/pi*180, gama/pi*180
   close(7)

 END subroutine
