program main 
 use prec
 use constant
 use latticevector
 use eigenvalue
 use fitting
 use emparam

  implicit none 
  integer  i,j,k,ii,jj, error_flag, ivm, im
  integer  icv
  real(dp), ALLOCATABLE :: kptc(:,:)

  real(dp), allocatable :: eighvb(:,:,:) 
  real(dp), allocatable :: eiglcb(:,:,:) 

  real(dp), allocatable :: vex(:),vey(:,:) 
  real(dp), allocatable :: coef(:,:)      
  real(dp), dimension(7,4):: em
  integer, dimension(7,3):: idx
  real(dp) normkbuf 
  integer degree, nkpt1m 
 ! Declare local variables
   INTEGER :: AllocateStatus,  DeAllocateStatus


  call readposcar()
  call readeigenval()

  if ( MOD(nelectrons,2) .GT. 0) then
     write(*,*) 'The calculated system contains odd electrons.'
     write(*,*) 'It will be no band gap. Please check your VASP calculations!'
     STOP
  end if

  allocate(kptc(nkpoints,3), eighvb(ispin,3,nkpoints), &
&       eiglcb(ispin,1,nkpoints), STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  do k=1,nkpoints
        do j=1,nbands
           if (j .EQ. (nelectrons/2-2)) then
               eighvb(:,1,k)=eigenvalues(:,j,k)  
           end if
           if (j .EQ. (nelectrons/2 -1)) then
               eighvb(:,2,k)=eigenvalues(:,j,k)  
           end if
           if (j .EQ. nelectrons/2) then
               eighvb(:,3,k)=eigenvalues(:,j,k)  
           end if
           if (j .EQ. (nelectrons/2 + 1)) then
               eiglcb(:,1,k)=eigenvalues(:,j,k)  
           end if
        end do
  end do

  open(25,file='k-eig.dat')
  write(25,*) nkpoints, ispin
  do k=1, nkpoints
     kptc(k,:) = MATMUL(kpoints(k,:),recip_lat)
     write(25,'(12F15.8)') (kptc(k,ii), ii = 1, 3), &
&       (eighvb(i,1,k)*ev2h, i=1,ispin),          &
&       (eighvb(i,2,k)*ev2h, i=1,ispin) ,         &
&       (eighvb(i,3,k)*ev2h, i=1,ispin) ,         &
&       (eiglcb(i,1,k)*ev2h, i=1,ispin) 
     
  end do
  close(25)

  call readinpem()
 
  degree=2
   allocate(coef(degree+1,4), vex(2*np+1), &
&           vey(2*np+1,4),STAT = AllocateStatus)
   IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"


  open(13,file='fitted-coff.dat')
  write(13,*)'Cofficients c, b, a in E(k) = c + b *|K-K0| + a *|K-K0|^2'
  write(13,'(A11,3F12.5)')'K0 at VBM:', (kvcbm(1,i), i=1,3) 
  write(13,'(A11,3F12.5)')'K0 at CBM:', (kvcbm(2,i), i=1,3) 
  write(13,*)'HVB: highest valence band; LCB: lowest conduction band'
  write(13,*)'(c, b, a) in order of HVB-2, HVB-1, HVB, LCB'

  open(14,file='efftive-mass.out')
  write(14,*)'Cofficients c, b, a in E(k) = c + b *|K-K0| + a *|K-K0|^2'
  write(14,'(A11,3F12.5)')'K0 at VBM:', (kvcbm(1,i), i=1,3) 
  write(14,'(A11,3F12.5)')'K0 at CBM:', (kvcbm(2,i), i=1,3) 
  write(14,*)'Effective masses of hole: m = -m_e * 1/(2 * a)'
  write(14,*)'Effective masses of electron: m = m_e * 1/(2 * a)'
  write(14,*)'m_e: electron rest mass'
  write(14,*)'HVB: highest valence band; LCB: lowest conduction band'
  write(14,*)'Effective masses of hole and electron'
  write(14,'(A18,4A10)')'*****************','****HVB-2*', '****HVB-1**',&
&            '*****HVB**', '*****LCB**'

  open(15,file='kd_eig4fit.dat')

   idx(1,:)=reshape((/1,0,0/),shape(idx(1,:)))  
   idx(2,:)=reshape((/0,1,0/),shape(idx(2,:)))
   idx(3,:)=reshape((/0,0,1/),shape(idx(3,:)))
   idx(4,:)=reshape((/1,1,0/),shape(idx(4,:)))
   idx(5,:)=reshape((/1,0,1/),shape(idx(5,:)))
   idx(6,:)=reshape((/0,1,1/),shape(idx(6,:)))
   idx(7,:)=reshape((/1,1,1/),shape(idx(7,:)))

  do i = 1, ispin
    if ( (i .EQ. 1) .AND. (ispin .EQ. 1)) then
       write(13,*) 'Non-spin-polarized results'
       write(14,*) 'Non-spin-polarized results'
    else  
       write(13,*) 'Spin channel:', i
       write(14,*) 'Spin channel:', i
    end if
         
    nkpt1m= 7 * 2 * np + 1 

    do icv = 1, nkvc 
    do j= 1, 7
      vex(1)=0.0d0
      do ivm= 1, 3
          vey(1,ivm)=eighvb(i,ivm,1+(icv-1)* nkpt1m)*ev2h
      end do
      vey(1,4)=eiglcb(i,1,1+(icv-1)* nkpt1m)*ev2h
      write(13,'(A1,3I1,A12)') '[',idx(j,1),idx(j,2),idx(j,3),'] direction:'

      jj=1
      if ( (j*2*np+1)+(icv-1)*nkpt1m .LE. nkpt1m*icv) then
          do  k = 2+ (j-1)*2*np+(icv-1)*nkpt1m,(j*2*np+1)+(icv-1)*nkpt1m
            jj=jj+1
            normkbuf= DOT_PRODUCT(kptc(k,:)-kptc(1+(icv-1)*nkpt1m,:),  &
&                          kptc(k,:)-kptc(1+(icv-1)*nkpt1m,:)) 
            if ( jj .LE. np+1) then
               vex(jj)=-sqrt(normkbuf) 
            else 
               vex(jj)=sqrt(normkbuf) 
            end if
            do ivm = 1, 3
               vey(jj,ivm)= eighvb(i,ivm,k)*ev2h
            end do
            vey(jj,4)= eiglcb(i,1,k)*ev2h
            if ( jj .EQ. np+2) then
              write(15,'(5F12.5)') vex(1), (vey(1,im),im=1,4)
            end if
            write(15,'(5F12.5)') vex(jj), (vey(jj,im),im=1,4)
          end do
          do im = 1, 4
             coef(:,im)=polyfit(vex,vey(:,im),degree)
          end do
          write(13,'(12F12.5)') (coef(:,im), im = 1, 4) 
          if (nkvc .GT. 1) then
              if ( icv .EQ. 1) then
                 do im = 1, 3
                    em(j,im)= -0.5d0/coef(degree+1,im)
                 end do
!               write(14,'(3F10.3,A10)') (em(j,im), im= 1,3), '---'
              else
                 em(j,4) = 0.5d0/coef(degree+1,4)
!                write(14,'(3A10,F10.3)') '---', '---','---', em(j,4) 
              end if
          else
              do im = 1, 3
                   em(j,im)= -0.5d0/coef(degree+1,im)
              end do
              em(j,4) = 0.5d0/coef(degree+1,4)
!             write(14,'(4F10.3)') (em(j,im), im= 1, 4)
          end if
        end if
        write(15,*) 
      end do
     end do
!! print out the effective masses of hole and electrons for each directions      
     do  j = 1, 7
          write(14,'(A2,3I1,A12,4F10.3)') '[',idx(j,1),idx(j,2),idx(j,3), &
                '] direction:', (em(j,im), im= 1, 4)
!         write(14,'(4F10.3)') (em(j,im), im= 1, 4)
     end do
         
   end do

   close(13)
   close(14)
   close(15)
   

    deallocate(eighvb,eiglcb,kptc,coef, vex, vey, &
 &         STAT = DeAllocateStatus) 
    IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"

   call  Dealloc_eigenvalues_kpoints()

end program
