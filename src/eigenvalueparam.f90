    MODULE  EIGENVALUE 
    use prec
    integer :: nkpoints, nbands, nions, norb,ispin, nelectrons, procartype
    real(dp), allocatable, dimension(:,:) :: kpoints 
    real(dp), allocatable, dimension(:,:,:) :: eigenvalues
    real(dp), allocatable, dimension(:,:,:,:,:) ::lmpdos
    END MODULE

    SUBROUTINE  read_head_of_procar()
     use prec
     use strings
     use eigenvalue 
     implicit none
     integer i, j, k,n, istat
     character(len=100):: title, filename
     character(len=10):: cbuf
     character(len=800):: cbuflong,corb
     character(len=1):: delims
     integer:: ik, iion, iband, isp
     integer, parameter:: StrMax=30, Nmax=100
     character(len=StrMax), dimension(Nmax):: args
     real(dp), dimension(3):: kbuf
     real(dp) :: eigbuf, wbuf
     logical :: lexist 
     
     filename='PROCAR'
     inquire(file=filename,exist=lexist)
     if ( .not.  lexist) then
          write(*,'(3A)') "Error:",trim(filename), " does not exist"
          STOP
     end if
     open(10,file=filename,status='old')
     read(10,'(a)',iostat=istat) title 
     delims=' '

     call parse(title, delims, args, n)
     
     procartype = 11 
     do j=1, n
        if ( trim(args(j)) .EQ. 'phase' ) then
           procartype=12
        end if
     end do

     read(10,*) cbuf, cbuf, cbuf, nkpoints, cbuf, cbuf, cbuf, nbands, cbuf,  &
&               cbuf, cbuf, nions 

     read(10,*) 
!!!  VASP  sphpro.F     
!!!    3201 FORMAT(/' k-point ',I4,' :',3X,3F11.8,'     weight = ',F10.8/)
!!    occasionally may meet error because kx,ky,kz are jointed togeter. 
!!    so skip reading the k points here.
!!   read(10,*) cbuf, ik, cbuf, (kbuf(j), j=1,3), cbuf, cbuf, wbuf 
     read(10,*) cbuf, ik, cbuf  
     read(10,*)
     read(10,*) cbuf, iband, cbuf, cbuf, eigbuf, cbuf, cbuf, wbuf
     read(10,*)

     read(10,'(a)',iostat=istat) corb 
     call parse(corb, delims, args, n)
     do j = 1, n
       if ( trim(args(j)) .EQ. 'tot') then
               norb=j-2
       end if
     end do

     ispin = 1

     do while (.true.)
          read(10,'(a)', iostat=istat) cbuflong 
          if (istat /=0) exit
          call parse(cbuflong, delims, args, n)
          do j = 1, n
             if ( trim(args(j)) .EQ. 'k-points:') then
               ispin = ispin + 1 
             end if
          end do
     end do

     close(10)
    END SUBROUTINE read_head_of_procar 


  SUBROUTINE readeigenval()
   use prec
   use constant
   use eigenvalue

   implicit none
   integer::  nblocks
   real(dp), dimension(3):: lat_norm
   real(dp) ::  omega, potim, buf, w
   character (len=256) car
   character (len=256) title 
   integer i, j, k, ii, istat 
   character (len=256) iiband 
   character(len=256) :: filename
   logical :: lexist 
  
  ! Declare local variables
   INTEGER :: AllocateStatus

!---------read EIGENVAL
!
  filename ="EIGENVAL"
  inquire(file=filename,exist=lexist)
  if ( .not.  lexist) then
       write(*,'(3A)') "Error:",trim(filename), " does not exist"
       STOP
  end if
  open(20,file=filename,status='old')
  read(20,*) nions, nions, nblocks, ispin
  read(20,*) omega, (lat_norm(i),i=1,3), potim
  read(20,*) buf 
  read(20,*) car 
  read(20,*) title 
!!  to be consistent with main.F of VASP
!!  WRITE(22,'(3I5)') NINT(INFO%NELECT),KPOINTS%NKPTS,WDES%NB_TOT
! read(20,'(3I5)',iostat=istat) nelectrons,nkpoints, nbands
  read(20,*,iostat=istat) nelectrons,nkpoints, nbands
  if (istat /=0) then
     write(*,*) 'Error in reading the 6th line of EIGENVAL'
     STOP
  end if


  allocate(eigenvalues(ispin,nbands,nkpoints), kpoints(nkpoints,3),  &
&           STAT = AllocateStatus)
  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

  do k=1,nkpoints
     read(20,*) 
!!!  VASP main.F    
!!!!       WRITE(22,'(4E15.7)') WDES%VKPT(1,NK),WDES%VKPT(2,NK),WDES%VKPT(3,NK),KPOINTS%WTKPT(NK)
!!! debug  here
!    read(20,'(4E15.7)') (kpoints(k,ii),ii = 1, 3), w
     read(20,*,iostat=istat) (kpoints(k,ii),ii = 1, 3), w
     if (istat /=0) then
         write(*,*) 'Error in reading coordinates of the', k, 'th points in EIGENVAL.'
         STOP
     end if
        do j=1, nbands
           read(20,*) iiband,(eigenvalues(i,j,k), i=1,ispin)
        end do
  end do
  close(20)
 END SUBROUTINE readeigenval


 SUBROUTINE Dealloc_eigenvalues_kpoints_lmpdos
 USE eigenvalue 

 IMPLICIT NONE

 ! Declare local variables
 INTEGER :: DeAllocateStatus

 ! Deallocate storage for array A 
 DEALLOCATE(eigenvalues,kpoints, lmpdos, STAT = DeAllocateStatus)
 IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"

 END SUBROUTINE Dealloc_eigenvalues_kpoints_lmpdos 


 SUBROUTINE Dealloc_eigenvalues_kpoints
 USE eigenvalue 

 IMPLICIT NONE

 ! Declare local variables
 INTEGER :: DeAllocateStatus

 ! Deallocate storage for array A 
 DEALLOCATE(eigenvalues,kpoints, STAT = DeAllocateStatus)
 IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"

 END SUBROUTINE Dealloc_eigenvalues_kpoints 
