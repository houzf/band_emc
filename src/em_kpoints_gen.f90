   program  main 
   use prec
   use LATTICEVECTOR 
   use emparam

   implicit none
   integer :: i, j, k, ik
   real(dp), ALLOCATABLE :: kptline(:,:), kptgrid(:,:)
   real(dp), dimension(3) :: kptc
   real(dp), dimension(2,3) :: kvcbm_cart
   real(dp) :: dk, w
   integer :: nkg, nkl
   integer, dimension(7,3):: idx
 ! Declare local variables
   INTEGER :: AllocateStatus,  DeAllocateStatus

   call readposcar()
   open(17,file='log.out')
   write(17,*) 'Unit cell lattice vectors (in Bohr):' 
   do i = 1, 3
      write(17,'(3F15.8)'), (prim_lat(i,j), j=1,3) 
   end do
   write(17,*) 'Reciprocal lattice vectors containing 2*pi (in 1/Bohr):' 
   do i = 1, 3
      write(17,'(3F15.8)'), (recip_lat(i,j), j=1,3) 
   end do

   write(17,*), 'Inverse of reciprocal lattice vectors:' 
   do i = 1, 3
      write(17,'(3F15.8)'), (inv_recip_lat(i,j), j=1,3) 
   end do

   call readinpem()

  
   w=1.0d0
   nkg= (2*np+1)**3
   nkl= (2*np+1)
   allocate(kptline(nkl,3), kptgrid(nkg,3), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

   if ((icf .EQ. 'C') .OR. (icf .EQ. 'c'))  then
          write(17,*) 'Cartesian coordinates of k points in the input file'
           kvcbm_cart=kvcbm
           do i = 1, nkvc 
             write(7,'(3F15.8)') (kvcbm_cart(i,j),j=1,3)
           end do
       ELSE IF ((icf .EQ. 'F') .OR. (icf .EQ. 'f') .OR. (icf .EQ. 'R') &
           &  .OR. (icf .EQ. 'r')) then
           write(17,*) 'Fractional coordinates of k points in the input file'
           do i = 1, nkvc 
             write(17,'(3F15.8)') (kvcbm(i,j),j=1,3)
           end do
           write(17,*)  'Fractional coordinates converted into Cartesian coord.'
           do i = 1, nkvc
            kvcbm_cart(i,:)= MATMUL(kvcbm(i,:),recip_lat) 
            write(17,'(3F15.8)') (kvcbm_cart(i,j),j=1,3)
           end do
       ELSE
           print *, 'Unkonw coordination type of k point in the input file!'
           stop 
   end if
   close(17)
   
   idx(1,:)=reshape((/1,0,0/),shape(idx(1,:)))  
   idx(2,:)=reshape((/0,1,0/),shape(idx(2,:)))
   idx(3,:)=reshape((/0,0,1/),shape(idx(3,:)))
   idx(4,:)=reshape((/1,1,0/),shape(idx(4,:)))
   idx(5,:)=reshape((/1,0,1/),shape(idx(5,:)))
   idx(6,:)=reshape((/0,1,1/),shape(idx(6,:)))
   idx(7,:)=reshape((/1,1,1/),shape(idx(7,:)))

   open(50,file='KPOINTS.lines',status='UNKNOWN')
   write(50,'(A32)') 'K-line around a specific k point'
   write(50,'(I5)') (1 + 7 * 2 * np) * nkvc 
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
          call gen_kline_around_kp(kptc, idx(j,:),np, delta, kptline)
          do  k = 1, 2 * np
             if (k .EQ. 1) then
                write(50,'(3F15.8,F10.3, A,3I1,A10)') &
&                   MATMUL(kptline(k,:),inv_recip_lat), &
&                   w, '  !', idx(j,1), idx(j,2), idx(j,3), 'direction'
             else
                write(50,'(3F15.8,F10.3)') MATMUL(kptline(k,:),inv_recip_lat), w
             end if
          end do
       end do

       call gen_kgrid_around_kp(kptc, np, delta, kptgrid) 
       do ik=1, nkg 
         write(60,'(3F15.8,F10.4)') MATMUL(kptgrid(ik,:),inv_recip_lat), w
       end do
   end do

   close(50)
   close(60)
   deallocate(kptline,kptgrid, STAT = DeAllocateStatus)
   IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"


   end program 
