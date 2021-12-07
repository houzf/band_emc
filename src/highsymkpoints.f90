   MODULE highsymmetrykpoints
     use prec
     implicit none
     integer n_high_sym_kpoints, nkpaths
     real(dp), allocatable, dimension(:,:) :: high_sym_kpoints
     integer,  allocatable, dimension(:) :: segments
     integer   nkpoints_gen
     character(len=1) :: icfk
     character(len=1), allocatable, dimension(:) :: lab_high_sym_kpoints 
   END MODULE

   SUBROUTINE  readsyml()
     use prec
     use highsymmetrykpoints 
     use latticevector
     implicit none
     integer i, j, k, ii, jj
     real(dp), allocatable, dimension(:)::  hskdistances
     real(dp), dimension(3) :: kbuf
     logical lsyml
     character(len=256) :: filename,filename1
 ! Declare local variables
     INTEGER :: AllocateStatus,  DeAllocateStatus

     filename='syml'
     inquire(file=filename,exist=lsyml)
     if ( .not.  lsyml) then
       write(*,'(3A)') "Error: ",trim(filename), " does not exist"
       STOP
     end if
     open(5,file=filename,status='old')
     read(5,*)   n_high_sym_kpoints
     nkpaths = n_high_sym_kpoints -1
     allocate(high_sym_kpoints(n_high_sym_kpoints,3),    &
 &            segments(n_high_sym_kpoints-1) ,             &
 &            hskdistances(n_high_sym_kpoints),  &
 &            lab_high_sym_kpoints(n_high_sym_kpoints),  &
 &            STAT = AllocateStatus)
     IF (AllocateStatus /= 0) STOP "*** Not enough memory *** readsyml()"

     read(5,*) (segments(i), i = 1, n_high_sym_kpoints-1)
     do i = 1, n_high_sym_kpoints
        read(5,*) lab_high_sym_kpoints(i), (high_sym_kpoints(i,j), j= 1,3)
     end do

     nkpoints_gen = 1
     do i = 1, n_high_sym_kpoints - 1 
         nkpoints_gen = nkpoints_gen + segments(i)
     end do

     close(5)

     icfk = 'F'
     call get_kdistances(n_high_sym_kpoints,high_sym_kpoints,icfk,hskdistances) 
     call write_highk(n_high_sym_kpoints,hskdistances)

     call gnuplot_template(n_high_sym_kpoints,lab_high_sym_kpoints, &
 &               hskdistances) 
 ! Deallocate storage for array A 
       DEALLOCATE(hskdistances, STAT = DeAllocateStatus) 
       IF (DeAllocateStatus /= 0) STOP &
 &              "*** Trouble deallocating *** gnuplot_template()"

     END SUBROUTINE readsyml   

     SUBROUTINE gen_kpoints4bs()
       use prec
       use highsymmetrykpoints
       implicit none
       integer i, j, k, ii
       real(dp)  weight
       real(dp), allocatable, dimension(:,:) :: delkpt
       character(len=256) :: filename
       logical lexist
 ! Declare local variables
       INTEGER :: AllocateStatus,  DeAllocateStatus

       weight = 1.0d0
       filename='KPOINTS'
       inquire(file=filename,exist=lexist)
       if ( .not.  lexist) then
          write(*,'(3A)') trim(filename), ' does not exist,', &
 &            ' it will be prepared for band structure calculations.'        
       else
         write(*,'(4A)') trim(filename), ' exist,', &
 &            ' it will be overwritted for band structure calculations,', &
 &            ' if syml exists.' 
       end if

       open(22,file=filename)

       call readsyml()

       allocate(delkpt(n_high_sym_kpoints-1,3),STAT = AllocateStatus)
       IF (AllocateStatus /= 0) STOP "*** Not enough memory *** gen_kpoints4bs()"

       write(22,*) 'K-points along high symmetry lines'
       write(22,*) nkpoints_gen 
       write(22,*) 'Reciprocal'
       write(22,'(3F12.6,F8.3,1X,A,1X,2A)') (high_sym_kpoints(1,ii), ii = 1,3), weight, &
 &                    '!', lab_high_sym_kpoints(1)
       do i = 1, n_high_sym_kpoints -1
          delkpt(i,:) = (high_sym_kpoints(i+1,:) - high_sym_kpoints(i,:)) &
 &                      /dble(segments(i))
          do j = 1, segments(i)
             if  (j == segments(i) )  then
                write(22,'(3F12.6,F8.3,1X,A,1X,2A)')  (high_sym_kpoints(i,ii) + delkpt(i,ii)  &
 &                   * dble(j), ii=1, 3), weight, '!', lab_high_sym_kpoints(i+1) 
             else
                write(22,'(3F12.6,F8.3)')  (high_sym_kpoints(i,ii) + delkpt(i,ii)  &
 &                   * dble(j), ii=1, 3), weight 
             end if
          end do
      end do
 ! Deallocate storage for array A 
      DEALLOCATE(delkpt, STAT = DeAllocateStatus) 
      IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"
    END SUBROUTINE gen_kpoints4bs



    SUBROUTINE  readkpoints()
     use prec
     use highsymmetrykpoints 

     implicit none
     integer i, j, k, ii, jj, intersection, istat1, istat2,istat3
     logical lkpoints
     character(len=1):: kmode
     character(len=256) :: title, filename, cbuf 
     real(dp), allocatable, dimension(:)::  hskdistances
     real(dp), dimension(3):: kbuf
 ! Declare local variables
     INTEGER :: AllocateStatus,  DeAllocateStatus

     filename='KPOINTS'
     inquire(file=filename,exist=lkpoints)
     if ( .not.  lkpoints) then
       write(*,'(3A)') "Error: ",trim(filename), " does not exist"
       call gen_kpoints4bs()
       STOP
     end if
     open(50,file=filename,status='old')
     read(50,'(A)') title
     read(50,*) intersection
     read(50,*) kmode
     if ( (kmode .EQ. 'L') .OR. (kmode .EQ. 'l') ) then
          read(50,*) icfk
          nkpaths = 0
          do while (.true.)
               read(50,*, iostat=istat1) (kbuf(j), j =1,3)
               read(50,*, iostat=istat2) (kbuf(j), j =1,3)
               read(50,*, iostat=istat3)  
               if ((istat1 /=0) .OR. (istat2 /=0) .OR.(istat3 /=0)) then
                  if (istat3 .LT. 0 )then
                    nkpaths = nkpaths + 1
                  end if
                  exit
               else  
                  nkpaths = nkpaths + 1
               end if
          end do
          n_high_sym_kpoints = 2*nkpaths 
          allocate(high_sym_kpoints(2*nkpaths,3),    &
 &            segments(nkpaths) ,             &
 &            hskdistances(2*nkpaths) ,             &
 &            lab_high_sym_kpoints(2*nkpaths),  &
 &            STAT = AllocateStatus)
          IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
          rewind(50)
          ! skip first four lines
          read(50,*)
          read(50,*)
          read(50,*)
          read(50,*)
          do i = 1, nkpaths
               read(50,*) (high_sym_kpoints(2*i-1,j), j =1,3), cbuf, &
&                          lab_high_sym_kpoints(2*i-1)               
               read(50,*) (high_sym_kpoints(2*i,j), j =1,3), cbuf, &
&                          lab_high_sym_kpoints(2*i)               
               if (istat3 .GE. 0) then
                  read(50,*)  
               end if
               segments(i) = intersection 
          end do
          close(50)
          
           call get_kdistances(n_high_sym_kpoints,high_sym_kpoints, icfk, hskdistances) 
           call write_highk(n_high_sym_kpoints,hskdistances)

           call gnuplot_template(n_high_sym_kpoints,lab_high_sym_kpoints, &
 &                   hskdistances) 
        ! Deallocate storage for array A 
           DEALLOCATE(hskdistances, STAT = DeAllocateStatus) 
           IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating *** 2"
       else
          write(*,*) 'No line mode in KPOINTS.'
          close(50)
          call gen_kpoints4bs()
!         STOP
     end if

     END SUBROUTINE readkpoints

 SUBROUTINE Dealloc_high_sym_kpoints()
 USE  highsymmetrykpoints 

 IMPLICIT NONE

 ! Declare local variables
 INTEGER :: DeAllocateStatus

 ! Deallocate storage for array A 
 DEALLOCATE(segments, high_sym_kpoints,   &
&         lab_high_sym_kpoints,  STAT = DeAllocateStatus)
 IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"

 END SUBROUTINE Dealloc_high_sym_kpoints
