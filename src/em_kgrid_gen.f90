   program  main 
   use prec
   use LATTICEVECTOR 

   implicit none
   integer :: i, j, k, ik
   real(dp), ALLOCATABLE :: kptline(:,:), kptgrid(:,:)
   real(dp), dimension(3) :: kptc
   real(dp), dimension(2,3) :: kvcbm, kvcbm_cart
   character*1 :: icf
   real(dp) :: dk, w
   integer :: nkvc, nkg, nk
   integer, dimension(7,3):: idx

   call readposcar()

   open(40,file='inp.em',status='old')
   read(40,*) (kvcbm(1,i),i=1,3), (kvcbm(2,i),i=1,3)
   read(40,*) icf
   read(40,*) dk , nk
   close(40)
 

   if ( (abs(kvcbm(1,1) - kvcbm(2,1)) .LT. 0.001d0) .and.  &
&       (abs(kvcbm(1,2) - kvcbm(2,2)) .LT. 0.001d0) .and.  &
&       (abs(kvcbm(1,3) - kvcbm(2,3)) .LT. 0.001d0) ) then
     nkvc = 1
   else 
     nkvc = 2
   end if
  
   w=1.0d0
   nkg= (2*nk+1)**3
   allocate(kptline(2*nk+1,3))
   allocate(kptgrid(nkg,3))

   open(7,file='log.out')
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

   if ((icf .EQ. 'C') .OR. (icf .EQ. 'c'))  then
          write(7,*) 'Cartesian coordinate of k points in the input file'
           kvcbm_cart(:,:)=kvcbm(:,:)
           do i = 1, nkvc 
             write(7,'(3F15.8)') (kvcbm_cart(i,j),j=1,3)
           end do
       ELSE IF ((icf .EQ. 'F') .OR. (icf .EQ. 'f') .OR. (icf .EQ. 'R') &
           &  .OR. (icf .EQ. 'r')) then
           write(7,*) 'Fractional coordinate of k points in the input file'
           do i = 1, nkvc 
             write(7,'(3F15.8)') (kvcbm(i,j),j=1,3)
           end do
           write(7,*)  'Fractional coordinate converted into Cartesian coord.'
           do i = 1, nkvc
            kvcbm_cart(i,:)= MATMUL(kvcbm(i,:),recip_lat) 
            write(7,'(3F15.8)') (kvcbm_cart(i,j),j=1,3)
           end do
       ELSE
           print *, 'Unkonw coordination type of k point in the input file!'
           stop 
   end if
   close(7)
   
   idx(1,:)=reshape((/1,0,0/),shape(idx(1,:)))  
   idx(2,:)=reshape((/0,1,0/),shape(idx(2,:)))
   idx(3,:)=reshape((/0,0,1/),shape(idx(3,:)))
   idx(4,:)=reshape((/1,1,0/),shape(idx(4,:)))
   idx(5,:)=reshape((/1,0,1/),shape(idx(5,:)))
   idx(6,:)=reshape((/0,1,1/),shape(idx(6,:)))
   idx(7,:)=reshape((/1,1,1/),shape(idx(7,:)))

   open(50,file='KPOINTS.lines',status='UNKNOWN')
   write(50,'(A32)') 'K-line around a specific k point'
   write(50,'(I5)') (1 + 7 * 2 * nk) * nkvc 
   write(50,'(A10)')'Reciprocal' 
   
   open(60,file='KPOINTS.grid',status='UNKNOWN')
   write(60,'(A32)') 'K-grid around a specific k point'
   write(60,'(I5)')  nkg * nkvc 
   write(60,'(A10)')'Reciprocal' 


   do i = 1, nkvc 
       kptc= kvcbm_cart(i,:)
       if ( nkvc .GT. 1) then
          if (i .EQ. 1) then
             write(50,'(3F15.8,F10.3, A17)') MATMUL(kptc,inv_recip_lat), &
&               w, '!k point of VBM'
          else 
             write(50,'(3F15.8,F10.3, A17)') MATMUL(kptc,inv_recip_lat), &
&               w, '!k point of CBM'
          end if
       else
          write(50,'(3F15.8,F10.3, A26)') MATMUL(kptc,inv_recip_lat), &
&           w, '!k points of VBM and CBM'
       end if

       do j = 1, 7
          call gen_kline_around_kp(kptc, idx, nk, dk,kptline)
          do  k = 1, 2 * nk
             if (k .EQ. 1) then
                write(50,'(3F15.8,F10.3, A,3I1,A10)') &
&                   MATMUL(kptline(k,:),inv_recip_lat), &
&                   w, '  !', idx(j,1), idx(j,2), idx(j,3), 'direction'
             else
                write(50,'(3F15.8,F10.3)') MATMUL(kptline(i,:),inv_recip_lat), w
             end if
          end do
       end do

       call gen_kgrid_around_kp(kptc, nk, dk, kptgrid) 
       do ik=1, nkg 
         write(60,'(3F15.8,F10.4)') MATMUL(kptgrid(ik,:),inv_recip_lat), w
       end do
   end do

   close(50)
   close(60)
   deallocate(kptline,kptgrid)

 contains


    subroutine gen_kline_around_kp(kp, id, nk, dk, kpt)
    use prec
    implicit none
    integer:: i, j
    integer,  intent(in)::   id(3),nk
    real(dp), intent(in)::   kp(3),dk
    real(dp), intent(out)::  kpt(2*nk+1,3)
    j=0
    do  i= -nk, nk
       if ( i .NE. 0) then
          j= j+1
          kpt(j,1) =  kp(1) + i * dble(id(1)) * dk
          kpt(j,2) =  kp(2) + i * dble(id(2)) * dk
          kpt(j,3) =  kp(3) + i * dble(id(3)) * dk
      endif
    end do
    return
    end subroutine gen_kline_around_kp


    subroutine gen_kgrid_around_kp(kp, nk, dk, kptg)
    use prec
    implicit none
    integer:: i, j, k, ik
    integer, intent(in)::   nk
    real(dp), intent(in)::   kp(3),dk
    real(dp), intent(out)::  kptg((2*nk+1)**3,3)
    ik=0
    do i = -nk, nk
        do j = -nk, nk
          do k = -nk, nk
            ik= ik+1
            kptg(ik,1) = kp(1)+ dble(i) * dk 
            kptg(ik,2) = kp(2)+ dble(j) * dk 
            kptg(ik,3) = kp(3)+ dble(k) * dk 
          end do  
        end do
    end do 
    write(*,*) ik

    return
    end subroutine gen_kgrid_around_kp


   end program 
