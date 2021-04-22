    program fs_bxsf 
      use prec
      use eigenvalue
      use latticevector
      use strings
      use constant
      implicit none

      integer:: i, j, k, isp, ik, ix, iy, iz, nx, ny, nz, ib
      integer:: ib4fs_start, ib4fs_end, nb4fs
      integer, dimension(3) :: nmesh
      real(dp) :: efermi
      character(len=256):: foutcar,finput, filename,foutname
      character(len=6) :: index_spin, str_nelect
      logical :: lfoutcar, lfinput
      real(dp), allocatable, dimension(:,:,:,:,:):: eig4fs
      integer :: ncmdarg, n, ios
      character(len=32), allocatable, dimension(:) :: cmdarg
      character(len=1):: delims
      integer, parameter:: StrMax=10, Nmax=30
      character(len=StrMax), dimension(Nmax):: args

 ! Declare local variables
      INTEGER :: AllocateStatus,  DeAllocateStatus
     

     ncmdarg = command_argument_count() 
     if ( (ncmdarg .LE. 2) .OR. (ncmdarg .GE. 5) ) then
         call print_help()
         STOP
       else 
         allocate(cmdarg(ncmdarg),STAT=AllocateStatus)
         if (AllocateStatus .ne. 0 ) then
               write(*,*) 'Fails to allocate space for cmdarg.'
               STOP
         end if

         do i = 1, ncmdarg
             call get_command_argument(i, cmdarg(i))
             call compact(cmdarg(i))
         end do
         do i = 1 ,3
            call value(cmdarg(i), nmesh(i), ios)
         end do
         nx = nmesh(1)  
         ny = nmesh(2)  
         nz = nmesh(3)  

         select CASE (ncmdarg)
           CASE (3)
             call fs_kpoints_gen(nx, ny, nz)
             STOP
           CASE (4)
             call readeigenval()
             delims='='
             call parse(cmdarg(4), delims, args, n)
             if ( n .eq. 2) then 
                  call  value(args(2), efermi, ios)
             else 
                 foutcar=cmdarg(4)
                 inquire(file=trim(foutcar),exist=lfoutcar)
                 if ( .not. lfoutcar) then
                    write(*,'(5A)') 'Error: ',trim(foutcar),' does not exist.'
                    print '(2a)', 'Please enter the value (in eV) of Fermi energy by',&
&                       ' using option: -ef=xxx.xx'
                    call print_help()
                    STOP
                 else 
                     call getefermi(foutcar,efermi)
                 end if
           end if
           end select
      end if
      
     
      if ( nkpoints .NE. (nx  * ny * nz)  ) then
          print '(2a,I4,A3,I4,A3,I4)', 'The EIGENVAL file is not the one calculated by', &
&              ' a periodic grid of ', nx, ' * ', ny, ' * ', nz
          write(*,*) 'nkpoints= ',nkpoints, 'nx * ny * nz= ',nx * ny * nz
          print '(2a)', 'Please check your KPOINTS file used for the Fermi ', &
&          'surface calcalculations!'
          print '(a)', 'Please check your EIGENVAL file.'
          STOP
      end if    

      call readposcar()

      if  (nelectrons .LE. 1 ) then
             STOP 
      else if ((nelectrons .GE. 2) .AND. (nelectrons .LE. 7) )  then
              ib4fs_start = nelectrons/2      
              ib4fs_end = nelectrons/2 + 1      
      else if ((nelectrons .GE. 8) .AND. (nelectrons .LE. 13) )  then
              ib4fs_start = nelectrons/2 - 2     
              ib4fs_end = nelectrons/2 + 2      
      else
              ib4fs_start = nelectrons/2 - 3     
              ib4fs_end = nelectrons/2 + 3      
      end if
      nb4fs = ib4fs_end - ib4fs_start + 1 
      allocate(eig4fs(ispin,nb4fs,nx+1,ny+1,nz+1), STAT=AllocateStatus)
      if (AllocateStatus .ne. 0 ) then
             write(*,*)  'Fails to allocate space for eig4fs.'
             STOP
      end if

      call writenum(nelectrons,str_nelect,'I6')
      call compact(str_nelect)

      do isp = 1, ispin
         filename='fs' 
         call writenum(isp,index_spin,'I6')
         call compact(index_spin)
         foutname=trim(filename)//trim('-spin-')//trim(index_spin)//trim('.bxsf')
         open(12,file=foutname)
         write(12,*) 'BEGIN_INFO'
         write(12,*) '#'
         write(12,*) '# this is a sample Band-XCRYSDEN-Structure-File' 
         write(12,*) '# generated from VASP calculations and ' 
         write(12,*) '# aimed for Visualization of Fermi Surface'
         write(12,*) '#'
         write(12,*) '# Case:      just an example'
         write(12,*) '#'
         write(12,*) '# Launch as: xcrysden --bxsf example.bxsf'
         write(12,*) '#'
         write(12,'(A15,F12.5)')'Fermi Energy:', efermi 
         write(12,*) 'END_INFO'
         write(12,*)
         write(12,*) 'BEGIN_BLOCK_BANDGRID_3D'
         write(12,*) trim('#_Number_of_electrons_')//trim(str_nelect)
         write(12,*) 'BEGIN_BANDGRID_3D_simple_example'
         write(12,*) nb4fs 
         write(12,'(3I6)')  nx+1, ny+1, nz+1
         write(12,'(3F10.3)') 0.0, 0.0 , 0.0 
         do i = 1, 3
             write(12,'(3F15.8)')  (recip_lat(i,j)/(b2a), j=1,3) 
         end do
         do  ib = 1,  nb4fs 
              write(12,*) 'BAND: ', ib4fs_start + ib - 1   
              ik = 0
              do ix = 1, nx
                 do iy = 1, ny
                   do  iz = 1, nz
                       ik = ik +1
                      eig4fs(isp,ib, ix, iy, iz) =  &
 &                    eigenvalues(isp, (ib4fs_start+ib-1), ik) 
                   end do
                 end do
              end do
              eig4fs(isp,ib, nx+1, :, :)=eig4fs(isp,ib, 1, :, :)
              eig4fs(isp,ib, :, ny+1, :)=eig4fs(isp,ib, :, 1, :)
              eig4fs(isp,ib, :, :, nz+1)=eig4fs(isp,ib, :, :, 1)
              do ix = 1, nx+1
                  do iy = 1, ny+1
                       write(12,'(1000F12.5)') (eig4fs(isp,ib, ix, iy, iz), iz= 1, nz+1) 
                  end do
                  if ( ix .NE. (nx+1) ) then
                    write(12,*)
                  end if
               end do
!              write(12,'(8F12.5)') (((eig4fs(isp, ib, ix, iy, iz), &
! &                iz = 1, nz), iy=1, ny), ix =1 , nx)             
         end do
         write(12,*) '  END_BANDGRID_3D'
         write(12,*) 'END_BLOCK_BANDGRID_3D'
         close(12)
     end do



     if ( (ncmdarg .GE. 3) .AND. (ncmdarg .LE. 4) ) then
         deallocate(cmdarg, eig4fs, STAT=DeAllocateStatus)
         IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating 1***"
     end if

     call  Dealloc_eigenvalues_kpoints()

contains

   subroutine fs_kpoints_gen(nx, ny, nz) 
     use prec
     implicit none 

     integer, intent(in) ::  nx, ny, nz 
     integer :: nkpt, ix, iy, iz, ik, i
     real(dp) :: w
     real(dp), allocatable, dimension(:,:):: kpt 

     
     w = 1.0d0
     if ( (nx .LT. 1) .OR. (ny .LT. 1) .OR. (nz .LT. 1) ) then
           write(*,*) 'nx, ny, nz should be larger than or equal to 1.'
           STOP
     end if

     nkpt= nx * ny * nz
     if ( nkpt .GE. 10000 ) then
        print '(2a)', 'Warning: number of generated k-points ',&
&             'is larger than 10000.' 
        print '(2a)', 'This will lead to a problem in the 6th line ', &
&                  'of EIGENVAL:'
        print '(a)', 'i.e., number of electrons and k-points will be messed up.'
     end if

     allocate(kpt(nkpt,3), STAT=AllocateStatus)
     if (AllocateStatus .ne. 0 ) then
             write(*,*)  'Fails to allocate space for kpt.'
             STOP
      end if

     filename='KPOINTS'
     open(22,file=filename)
     write(22,'(A40,I4,A3, I4, A3, I4)') 'K-points generated by a periodic grid of ', &
&       nx, ' * ', ny, ' * ', nz
     write(22,*) nkpt 
     write(22,*) 'Reciprocal'
     ik = 0
     do  ix = 1, nx
         do iy = 1, ny
            do  iz = 1, nz
              ik = ik + 1  
              kpt(ik,1) = 1.0d0/REAL(nx) * (ix - 1) 
              kpt(ik,2) = 1.0d0/REAL(ny) * (iy - 1) 
              kpt(ik,3) = 1.0d0/REAL(nz) * (iz - 1) 
              write(22,'(3F12.6,F10.3)')  (kpt(ik,i), i=1,3), w 
            end do
          end do
     end do

     close(22)
     deallocate(kpt)
  end subroutine fs_kpoints_gen

  subroutine print_help()
    print '(a)', 'usage: fs_bxsf nx ny nz [options]'
    print '(a)', ''
    print '(a)', ''
    print '(2a)', 'nx, ny, nz: number of data-points in each direction', &
&               '  of reciprocal lattice vectors '
    print '(a)', 'fs_bxsf options: '
    print '(a)', ''
    print '(2a)', ' OUTCAR       find out the value of Fermi energy from ',&
&             ' a file named OUTCAR'
    print '(a)', ' -ef=xxx.xxx    xxx.xxx is the value of Fermi energy.'
    print '(a)', 'Example usages of fs_bxsf:'
    print '(a)', 'fs_bxsf  nx, ny, nz'
    print '(a)', 'fs_bxsf  nx, ny, nz  OUTCAR'
    print '(a)', 'fs_bxsf  nx, ny, nz  -ef=xxx.xxx'
  end subroutine print_help

end program fs_bxsf 
