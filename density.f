c     Aaproximación SPH para la densidad

      subroutine density(mspace,ntype,npairs,pairs,nfilas,w)
      implicit none
      include 'param.inc'

      integer i,j
      integer ntype(2),npairs,pairs(npairs,ntype(1)),nfilas(ntype(1))
      double precision mspace(25,nmax),w(npairs,ntype(1)),rho,r,dx(dim),
     +     selfdens,wp(ntype(1)),rhoa(ntype(1))

      do i=1,dim
         dx(i)=0.e0
      enddo
      
      r = 0.
      do i=1,ntype(1)
         call Kernel(r,dx,mspace(13,i),selfdens,dx)
         mspace(9,i) = selfdens * mspace(8,i)
c         rhoa(i) = mspace(9,i)
c         write(*,*)'density',i,mspace(13,i),mspace(9,i),mspace(12,i)
      enddo

      do i=1,ntype(1)
c         wp(i) = 0.0
         do j=1,nfilas(i)
c            wp(i)=wp(i) + (mspace(8,pairs(j,i))/rhoa(pairs(j,i)))*w(j,i)
            mspace(9,i) = mspace(9,i) + mspace(8,pairs(j,i))*w(j,i)
c            write(*,*)'d',i,mspace(9,i),mspace(12,i),mspace(13,1),w(j,i)
         enddo
c         mspace(9,i)=mspace(9,i)/wp(i)
c         write(*,*)'d',i,mspace(9,i),mspace(12,i),w(j-1,i),wp(i)
      enddo
      
      end
