    program bandstructureplot 
      use prec
      use eigenvalue
      use latticevector
      use strings
      implicit none

      integer i, j, k, iorb,isp
      real(dp), allocatable, dimension(:) :: kdistance 
      real(dp), allocatable, dimension(:,:,:):: eigenvalues_lcb
      real(dp), allocatable, dimension(:,:,:):: eigenvalues_hvb
      real(dp), allocatable, dimension(:) :: porbuf, dorbuf, forbuf
      real(dp), dimension(3):: kbuf
      real(dp) ::  evbm, ecbm, efermi
      integer :: ikvbm(1), ikcbm(1)
      character(len=1) :: icfkb
      character(len=256):: filename1,filename2,filename3,foutcardos,foutcarscf
      character(len=256):: finput
      character(len=256):: cbuf,fprefix,fsuffix, fsuffix1,fnbuf1, fsuffix2,fnbuf2,cbuf1,cbuf2
      character(len=256):: cbuf3,cbuf4,fprefix3,fsuffix3,fnbuf3 
      character(len=6) :: index_atom 
      logical :: lprocar,leigenval, lread_procar, lfoutcardos, lfoutcarscf
      logical :: lfinput
      integer :: iloop, nloop, mten, istart, iend

 ! Declare local variables
      INTEGER :: AllocateStatus,  DeAllocateStatus
     
      filename1='PROCAR'
      filename2="EIGENVAL"
      inquire(file=filename1,exist=lprocar)
      inquire(file=filename2,exist=leigenval)
      if ( (.not.  lprocar) .AND. (.not. leigenval) ) then
           write(*,'(5A)') 'Error:',trim(filename1),' and ', trim(filename2), ' do not exist.'
           STOP
      else if ( lprocar ) then
          call readprocar() 
          lread_procar=.TRUE.
      else 
          call readeigenval()
          lread_procar=.FALSE.
      end if

      call readposcar()

      allocate(kdistance(nkpoints),eigenvalues_hvb(ispin,1,nkpoints), &
 &            eigenvalues_lcb(ispin,1,nkpoints), STAT = AllocateStatus) 
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      icfkb = 'F'
      call get_kdistances(nkpoints, kpoints, icfkb, kdistance)

 !!! extract Fermi energy from OUTCAR of DOS or SCF calculations
      foutcardos='../dos/OUTCAR'
      foutcarscf='../scf/OUTCAR'
      inquire(file=foutcardos,exist=lfoutcardos)
      inquire(file=foutcarscf,exist=lfoutcarscf)
      if ( (.not.  lfoutcardos) .AND. (.not. lfoutcarscf) ) then
           write(*,'(5A)') 'Error: ',trim(foutcardos),' and ', trim(foutcarscf), ' do not exist.'
           write(*,*)'Please enter the location of OUTCAR of DOS or SCF calculations:'
           read(*,'(A)') finput
           call compact(finput)
           inquire(file=trim(finput),exist=lfinput)
           if ( .not. lfinput) then
              write(*,'(5A)') 'Error: ',trim(finput),' does not exist.'
              write(*,*)'Please enter the value (in eV) of Fermi energy:'
              write(*,*)'If you do konw it, please enter 0.0'
              read(*,*) efermi 
           else 
              call getefermi(finput,efermi)
           end if
      else if ( lfoutcardos ) then
          call getefermi(foutcardos,efermi) 
      else 
          call getefermi(foutcarscf,efermi) 
      end if
      write(*,'(A9,F10.5)')'EFermi = ', efermi


      filename3='bnd.dat'
      open(30,file=filename3)
      do j = 1, nbands
        do k = 1, nkpoints
           write(30,'(3F12.5)') kdistance(k), (eigenvalues(i,j,k)-efermi, i=1, ispin)
        end do
           write(30,*)
      end do
      close(30)

      if ( lprocar ) then 
         fprefix='fatbnd-atom-'
         fsuffix1='-all.dat'
         fsuffix2='-spdf.dat'
         do i =1, nions
 !           write(index_atom,'(I6)') i 
             call writenum(i,index_atom,'I6')
             call compact(index_atom)
             cbuf=trim(fprefix)//trim(index_atom)
             fnbuf1=trim(cbuf)//trim(fsuffix1)
             fnbuf2=trim(cbuf)//trim(fsuffix2)
             open(40,file=fnbuf1)
             open(50,file=fnbuf2)
             write(40,*) '#kdis   eig(isp)  s(isp)  py(isp) pz(isp) px(isp) dxy(isp)', &
&                      ' dyz(isp) dz2 (isp) dxz(isp) dx2(isp) ... tot(isp)'   
             write(50,*) '#kdis  eig(isp)  s(isp)  p(isp) d(isp) ... tot(isp)'
             do j = 1, nbands
               do k = 1, nkpoints
                  write(40,'(35F12.5)') kdistance(k),           &
 &                    (eigenvalues(isp,j,k)-efermi, isp=1, ispin),     &
 &                    ((lmpdos(isp,k,j,i,iorb),isp=1,ispin),iorb=1,norb+1)   
                  allocate(porbuf(ispin),dorbuf(ispin),forbuf(ispin), STAT = AllocateStatus)
                  IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
                  do isp=1, ispin
                      porbuf(isp)= lmpdos(isp,k,j,i,2) + lmpdos(isp,k,j,i,3) &
 &                              + lmpdos(isp,k,j,i,4) 
                      dorbuf(isp)= lmpdos(isp,k,j,i,5) + lmpdos(isp,k,j,i,6) &
 &                               + lmpdos(isp,k,j,i,7) + lmpdos(isp,k,j,i,8) & 
 &                               + lmpdos(isp,k,j,i,9) 
                      if ( (norb .GT. 9) .AND. (norb .LE. 16)) then
                        forbuf(isp)= lmpdos(isp,k,j,i,10) + lmpdos(isp,k,j,i,11) &
 &                                 + lmpdos(isp,k,j,i,12) + lmpdos(isp,k,j,i,13) & 
 &                                 + lmpdos(isp,k,j,i,14) + lmpdos(isp,k,j,i,15) & 
 &                                 + lmpdos(isp,k,j,i,16) 
                      end if
                  end do
                  if  (norb .LE. 9) then
                     write(50,'(11F12.5)') kdistance(k),           &
 &                       (eigenvalues(isp,j,k)-efermi, isp=1, ispin),     &
 &                       (lmpdos(isp,k,j,i,1),isp=1,ispin),        &
 &                       (porbuf(isp),isp=1,ispin),                &
 &                       (dorbuf(isp),isp=1,ispin),                &
 &                       (lmpdos(isp,k,j,i,norb+1),isp=1,ispin)        
                  end if
                  if ( (norb .GT. 9) .AND. (norb .LE. 16)) then
                     write(50,'(13F12.5)') kdistance(k),           &
 &                       (eigenvalues(isp,j,k)-efermi, isp=1, ispin),     &
 &                       (lmpdos(isp,k,j,i,1),isp=1,ispin),        &
 &                       (porbuf(isp),isp=1,ispin),                &
 &                       (dorbuf(isp),isp=1,ispin),                &
 &                       (forbuf(isp),isp=1,ispin),                &
 &                       (lmpdos(isp,k,j,i,norb+1),isp=1,ispin)        
                  end if
                  DEALLOCATE(porbuf,dorbuf,forbuf, &
 &                             STAT = DeAllocateStatus)
                  IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"
              end do
              write(40,*)
              write(50,*)
            end do
            close(40)
            close(50)
         end do    

         fprefix3='fatbnd-tot-atoms-'
         fsuffix3='.dat'
         nloop= INT(nions/10)
         do iloop = 0, nloop
            mten = nions - iloop * 10
            istart = 1+ iloop * 10
            iend = iloop * 10 + min(10, mten)
            if (mten .NE. 0) then 
                 if ( mten .EQ. 1 ) then
 !                  write(index_atom,'(I6)') istart
                    call writenum(istart,index_atom,'I6')
                    call compact(index_atom)
                    cbuf3=trim(fprefix3)//trim(index_atom)
                    fnbuf3=trim(cbuf3)//trim(fsuffix3)
                    open(60,file=fnbuf3)
                    write(60,*) '#kdis  eig(isp) (total(i,isp),isp=1,ispin),i=1,10)'
                    do j = 1, nbands
                      do k = 1, nkpoints
                        write(60,'(23F12.5)') kdistance(k),         &
 &                           (eigenvalues(isp,j,k)-efermi, isp=1, ispin),    &
 &                          ((lmpdos(isp,k,j,i,norb+1),isp=1,ispin),i=istart,istart)        
                      end do
                      write(60,*)
                    end do
                    close(60)
                 else
  !                 write(index_atom,'(I6)') istart
                    call writenum(istart,index_atom,'I6')
                    call compact(index_atom)
                    cbuf1=trim(fprefix3)//trim(index_atom)
  !                 write(index_atom,'(I6)') iend 
                    call writenum(iend,index_atom,'I6')
                    cbuf2='-'
                    cbuf3=trim(cbuf1)//trim(cbuf2)
                    call compact(index_atom)
                    cbuf4=trim(cbuf3)//trim(index_atom)
                    fnbuf3=trim(cbuf4)//trim(fsuffix3)
                    open(60,file=fnbuf3)
                    write(60,*) '#kdis  eig(isp) (total(i,isp),isp=1,ispin),i=1,10)'
                    do j = 1, nbands
                      do k = 1, nkpoints
                        write(60,'(23F12.5)') kdistance(k),         &
 &                           (eigenvalues(isp,j,k)-efermi, isp=1, ispin),    &
 &                       ((lmpdos(isp,k,j,i,norb+1),isp=1,ispin),i=istart,iend)        
                      end do
                      write(60,*)
                    end do
                    close(60)
                 end if
            end if
          end do
      end if

      eigenvalues_hvb(:,1,:)=eigenvalues(:,nelectrons/2,:)
      eigenvalues_lcb(:,1,:)=eigenvalues(:,nelectrons/2+1,:)

      if ( MOD(nelectrons,2) .GT. 0) then
         write(*,*) 'The calculated system contains odd electrons.'
         write(*,*) 'It will be no band gap. Please check your VASP calculations!'
         write(*,*) 'Effective mass could not be calculated in next step.'
       else
         write(*,*) 'Input file for estimation of effective mass will be prepared.'
      end if

      open(60,file='inp.em')
 !    Only for the first spin channel     
      do i = 1, 1 
          evbm=MAXVAL(eigenvalues_hvb(i,1,:))
          ecbm=MINVAL(eigenvalues_lcb(i,1,:))
          ikvbm=MAXLOC(eigenvalues_hvb(i,1,:))
          ikcbm=MINLOC(eigenvalues_lcb(i,1,:))
          write(60,'(6F9.5, A26, 2F10.4)') kpoints(ikvbm(1),:),   &
 &                                      kpoints(ikcbm(1),:),   & 
 &                  '!k points of VBM and CBM', evbm, ecbm
      end do
      write(60,'(A1, A20)') 'F', ' !Coordination type' 
      write(60,'(F7.4,I3,A65)') 0.002, 7,  &
 &        ' !delta_k, n (which gives 2*n+1 k-points along each direction)'
      close(60)

      DEALLOCATE(eigenvalues_hvb,eigenvalues_lcb, kdistance, &
 &                 STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"

      if ( lread_procar ) then 
        call Dealloc_eigenvalues_kpoints_lmpdos()
      else
        call  Dealloc_eigenvalues_kpoints()
      end if
    end program bandstructureplot 
