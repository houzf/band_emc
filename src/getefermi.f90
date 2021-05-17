  SUBROUTINE  getefermi(filename,efermi)
      use prec
      use strings
      implicit none
      real(dp),intent(out):: efermi
      character(len=256),intent(in):: filename
      integer i, j,n, istat
      character(len=256):: cbuflong, sbuf
      integer, parameter:: StrMax=10, Nmax=20
      character(len=StrMax), dimension(Nmax):: args
      integer:: ik, iion, iband, isp
      logical :: lexist 
      character(len=1):: delims
     
      inquire(file=filename,exist=lexist)
      if ( .not.  lexist) then
           write(*,'(3A)') "Error:",trim(filename), " does not exist"
           STOP
      end if

      delims=' '
      open(10,file=trim(filename),status='old')
      do while (.true.)
         read(10,'(a)',iostat=istat) cbuflong 
         if (istat /=0) exit
         call parse(cbuflong, delims, args, n)
         if ((trim(args(1)) .EQ. 'BZINTS:').AND.(trim(args(2)).EQ.'Fermi')) then
             sbuf = args(4)//args(5) 
         end if
       end do

      delims=';'
      call parse(sbuf,delims,args,n)
      call value(args(1),efermi,istat)
      return
       

     close(10)
   END SUBROUTINE getefermi
