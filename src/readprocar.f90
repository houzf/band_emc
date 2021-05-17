    SUBROUTINE readprocar()
      use prec
      use eigenvalue 
      use strings
      implicit none
      integer i, j, k, n, m
      character(len=100):: title, corb,filename
      character(len=10):: cbuf
      character(len=200):: cbuflong 
      character(len=200):: cbufkx,cbufky,cbufkz 
      real(dp):: kx, ky, kz
      integer:: ik, iion, ii, istat, iband,isp, iorb,ibuf
      real(dp) :: wbuf, occ,sumocc
      logical :: lexist 
 ! Declare local variables
      INTEGER :: AllocateStatus,  DeAllocateStatus
 ! parse special cases
      integer, parameter:: StrMax=40, Nmax=100
      character(len=StrMax), dimension(Nmax):: args, argy
      character(len=1):: delims
      
      filename='PROCAR'
      inquire(file=filename,exist=lexist)
      if ( .not.  lexist) then
           write(*,'(3A)') "Error:",trim(filename), " does not exist"
           STOP
      else 
          call read_head_of_procar() 
      end if

      write(*,'(A10,I6)')  adjustl('ISPIN ='), ispin
      write(*,'(A10,I6)')  adjustl('NKPOINTS ='), nkpoints
      write(*,'(A10,I6)')  adjustl('NBANDS ='),nbands
      write(*,'(A10,I6)')  adjustl('NATOMS ='),nions
      write(*,'(A10,I6)')  adjustl('NORBITS ='),norb

      allocate(eigenvalues(ispin,nbands,nkpoints), kpoints(nkpoints,3),  &
 &         lmpdos(ispin,nkpoints,nbands,nions+1,norb+1),STAT = AllocateStatus )
      IF (AllocateStatus /= 0) STOP "*** Not enough memory *** readprocar()"

      open(unit=10,file=filename, status='old')
      read(10,'(a)',iostat=istat) title 

      sumocc=0.0d0
      do isp = 1, ispin
         read(10,'(a)') cbuflong
!        write(*,*) trim(cbuflong)
         do ik = 1, nkpoints
           read(10,*) 
!!!  VASP  sphpro.F     
!!!    3201 FORMAT(/' k-point ',I4,' :',3X,3F11.8,'     weight = ',F10.8/)
!!    occasionally may meet error because kx,ky,kz are jointed togeter. 
!          read(10,*,iostat=istat) cbuf, ibuf, cbuf, (kpoints(ik,j), j=1,3), cbuf, cbuf, wbuf 
           read(10,*, iostat=istat) cbuf, ibuf, cbuf, cbufkx, cbufky, cbufkz 
           if (istat /=0) then
               write(*,*) 'Error in reading the',ik,'th k-point in PROCAR.'
               STOP
           end if
          delims='-'
          call parse(cbufkx, delims, args, n)
          call parse(cbufky, delims, argy, m)
          if (cbufkx(1:1) == '-' ) then
      ! kx is negative 
              select case(n)
                 case (4)
      !  kx, ky and kz are jointed.       
                  call value(args(2), kx, istat)
                  call value(args(3), ky, istat)
                  call value(args(4), kz, istat)
                  kpoints(ik,1) = -kx 
                  kpoints(ik,2) = -ky 
                  kpoints(ik,3) = -kz 
                 case (3)
      !  kx and ky are jointed.       
                    call value(args(2), kx, istat)
                    call value(args(3), ky, istat)
                    call value(cbufky, kz, istat)
                    kpoints(ik,1) = -kx 
                    kpoints(ik,2) = -ky 
                    kpoints(ik,3) =  kz 
                 case (2)
                    call value(cbufkx, kx, istat)
                    kpoints(ik,1) = kx 
                    if (m == 2) then
       ! ky and kz are jointed        
                     call value(argy(1), ky, istat)
                     call value(argy(2), kz, istat)
                     kpoints(ik,2) =  ky 
                     kpoints(ik,3) = -kz 
                    else
       ! kx, ky and kz are separated.        
                     call value(cbufky, ky, istat)
                     call value(cbufkz, kz, istat)
                     kpoints(ik,2) =  ky 
                     kpoints(ik,3) =  kz 
                    end if
              end select
          else
      ! kx is positive            
             select case (n) 
                case (3) 
      !  kx, ky and kz are jointed.       
                  call value(args(1), kx, istat)
                  call value(args(2), ky, istat)
                  call value(args(3), kz, istat)
                  kpoints(ik,1) =  kx 
                  kpoints(ik,2) = -ky 
                  kpoints(ik,3) = -kz 
                case (2)
      !  kx and ky are jointed.       
                  call value(args(1), kx, istat)
                  call value(args(2), ky, istat)
                  call value(cbufky, kz, istat)
                  kpoints(ik,1) =  kx 
                  kpoints(ik,2) = -ky 
                  kpoints(ik,3) =  kz 
                case (1)
                  call value(cbufkx, kx, istat)
                  kpoints(ik,1) =   kx 
       ! ky and kz are jointed        
                  if ( m == 2) then
                     call value(argy(1), ky, istat)
                     call value(argy(2), kz, istat)
                     kpoints(ik,2) =  ky 
                     kpoints(ik,3) = -kz 
                  else 
       ! kx, ky and kz are separated.        
                     call value(cbufky, ky, istat)
                     call value(cbufkz, kz, istat)
                     kpoints(ik,2) =  ky 
                     kpoints(ik,3) =  kz 
                  end if
                end select

          end if
!   for debug          
!        write(*,'(3F12.6)') kpoints(ik,1), kpoints(ik,2),  kpoints(ik,3)
           read(10,*)
           do iband = 1, nbands
              read(10,*) cbuf, ibuf, cbuf, cbuf, eigenvalues(isp,iband,ik), cbuf, cbuf, occ 
              read(10,*)
              sumocc= sumocc + occ
              read(10,'(a)',iostat=istat) corb 

              if (nions .GT. 1) then
                 do ii = 1, nions+1
                     read(10,*) cbuf, (lmpdos(isp,ik,iband,ii,iorb), iorb=1,norb+1) 
                 end do
              else
                 read(10,*) cbuf, (lmpdos(isp,ik,iband,1,iorb), iorb=1,norb+1) 
                  lmpdos(isp,ik,iband,2,:)=lmpdos(isp,ik,iband,1,:)
              end if
              if (procartype .EQ. 12) then
                 read(10,'(A)') cbuf
                 do ii = 1, 2*nions
                   read(10,'(A)') cbuf
                 end do
              end if
              read(10,*)
           end do
        end do
     end do
     nelectrons=INT(NINT(sumocc)/nkpoints)
     write(*,'(A10,I6)') adjustl('NELECT ='), nelectrons
 END SUBROUTINE readprocar 
      
