      implicit none
      real*8 f,noise,dum
      integer n,i
      character*60 infile,outfile

      infile='curves.txt'
      outfile='AdLIGODwyer.dat'

      open (10,file=infile)
      open (11,file=outfile)

      do i=1,3000
         read (10,*) f,noise,dum,dum,dum,dum,dum,dum,dum,dum,dum
         write (11,'(2(e15.8,1x))') f,noise
      end do
      
      close(10)
      close(11)

      end
