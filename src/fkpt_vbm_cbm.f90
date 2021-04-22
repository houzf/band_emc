    program fkpt_vbm_cbm
      use prec
      use eigenvalue
      implicit none

      integer i, j, k
      real(dp), allocatable, dimension(:,:,:):: eigenvalues_lcb
      real(dp), allocatable, dimension(:,:,:):: eigenvalues_hvb
      real(dp) evbm, ecbm
      integer  ikvbm(1), ikcbm(1)
 ! Declare local variables
      INTEGER :: AllocateStatus,  DeAllocateStatus

      call readeigenval()

      if ( MOD(nelectrons,2) .GT. 0) then
         write(*,*) 'The calculated system contains odd electrons.'
         write(*,*) 'It will be no band gap. Please check your VASP calculations!'
         write(*,*) 'Effective mass could not be calculated in next step.'
         STOP
       else
         write(*,*) 'Input file for estimation of effective mass will be prepared.'
      end if
        
      allocate(eigenvalues_hvb(ispin,1,nkpoints), &
 &            eigenvalues_lcb(ispin,1,nkpoints), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

      eigenvalues_hvb(:,1,:)=eigenvalues(:,nelectrons/2,:)
      eigenvalues_lcb(:,1,:)=eigenvalues(:,nelectrons/2+1,:)

      open(20,file='inp.em')
!    treat only the first spin channel      
      do i = 1, 1 
          evbm=MAXVAL(eigenvalues_hvb(i,1,:))
          ecbm=MINVAL(eigenvalues_lcb(i,1,:))
          ikvbm=MAXLOC(eigenvalues_hvb(i,1,:))
          ikcbm=MINLOC(eigenvalues_lcb(i,1,:))
          write(20,'(6F9.5, A26, 2F10.4)') kpoints(ikvbm(1),:),   &
 &                                      kpoints(ikcbm(1),:),   & 
 &                  '!k points of VBM and CBM', evbm, ecbm
      end do
      write(20,'(A)') 'F'
      write(20,'(F7.4,I3)') 0.002, 7
      close(20)

      DEALLOCATE(eigenvalues_hvb,eigenvalues_lcb, STAT = DeAllocateStatus)
      IF (DeAllocateStatus /= 0) STOP "*** Trouble deallocating ***"

      call  Dealloc_eigenvalues_kpoints()
    end program fkpt_vbm_cbm
