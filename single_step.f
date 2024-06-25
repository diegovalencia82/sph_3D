
c      subroutine single_step(id, itimestep, dt, ntotal, hsml, mass, x, 
c     &     vx, u, s,rho, p, t, dx, dvx, du, ds, drho,itype, av)

      subroutine single_step(mspace,ntype,itimestep,dt,av)
c----------------------------------------------------------------------


c Subroutine to determine the right hand side of a differential
c equation in a single step for performing time integration

c In this routine and its subroutines the SPH algorithms are performed.

c     ntotal-- total particle number ues                             [in]
      
c     dt--- Time step used in the time integration                   [in]
c     time -- time of the snapshot
c     itimestep -- number of the time step.
c     mspace(1,i) = id -- label of each particle                   [out]      
c     mspace(2 to 4,i) = x-- coordinates of particles              [in/out]
c     mspace(5 to 7,i) = vx-- velocities of particles              [in/out]
c     mspace(8,i) = mass-- mass of particles                       [in]
c     mspace(9,i) = rho-- dnesities of particles                   [in/out]
c     mspace(10,i) = p-- pressure of particles                     [in/out]
c     mspace(11,i) = u-- internal energy of particles              [in/out]
c     mspace(12,i) = itype-- types of particles                    [in]
c     mspace(13,i) = hsml-- smoothing lengths of particles         [in/out]
c     mspace(14,i) = c-- sound velocity of particle                [out]
c     mspace(15,i) = s-- entropy of particles                      [out]
c     mspace(16,i) = e-- total energy of particles                 [out]
c     mspace(17 to 19,i) = dx  : dx = vx = dx/dt                   [out]
c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [out]
c     mspace(23,i) = du        : du = du/dt                        [out]
c     mspace(24,i) = ds        : ds = ds/dt                        [out]
c     mspace(25,i) = drho      : drho = drh,o/dt                   [out]

c     npairs    : maximum number of pairs interaction             [out]
c     nfilas    : number of interaction for each fluid particle   [out]
c     mrij      : matrix of ri-rj for all fluid particles for each interaction [out]
c     mxiij      : matrix of xi-xj for all fluid particles for each interaction [out]         
c     mvij      : matrix of vi-vj for all fluid particles for each interaction [in]
c     mvxij     : matrix of vxi-vxj for all fluid particles for each interaction [in]      

c     t         : Temperature                                      [in/out]
c tdsdt     : Production of viscous entropy t*ds/dt               [out]
c av        : Monaghan average velocity                           [out]

      
      
c rdomain   : smoothing length                                     [in]
      
      implicit none
      include 'param.inc'

      integer i,j,itimestep, ntotal,npairs,ntype(2),nfluid,nvirt
c     parameter ( npairs = (kappa0+2)**dim)
      parameter ( npairs = (kappa0+5)**dim)
      double precision mspace(25,nmax),av
      double precision t(maxtimestep),dt,rdomain
      integer pairs(npairs,ntype(1)),nfilas(ntype(1))
      integer pairsv(npairs,ntype(1)),nfilasv(ntype(1))
      double precision w(npairs,ntype(1)),dwdx(3,npairs,ntype(1)),
     +     mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1)),
     +     mrijv(npairs,ntype(1)),mxijv(3,npairs,ntype(1))
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))
      double precision rrr,xi,yi,dx(dim),r,dwdx0(dim),ww

      rdomain = kappa0 * h0 * 1.00000007 !equivalent to hsml

      do i=1,nmax
         mspace(13,i) = h0      !rdomain
      enddo

      if (mod(itimestep,save_step).eq.0) then
         write(*,*)'-------------------------------------'
         write(*,*)'Max number of n pairs = ',npairs
         write(*,*)'kappa0 = ',kappa0,'  h0 = ',h0
         write(*,*)'rdomain = ',rdomain
         write(*,*)'-------------------------------------'
      endif

      write(*,*)'ss111'
      
      call  input(mspace,ntotal,nfluid,nvirt,1)

      write(*,*)'ss222'
      
      call neighboring_search(rdomain,mspace,ntype,npairs,pairs,nfilas,
     +     mrij,mxij,mvij,mvxij)

      write(*,*)'ss333'
      
c      call neighboring_searchv(rdomain,mspace,ntype,npairs,pairsv,
c     +     nfilasv,mrijv,mxijv)

      call wijdwij(mspace,pairs,mrij,mxij,npairs,nfilas,ntype,w,dwdx)

      write(*,*)'ss444'
      
      call density(mspace,ntype,npairs,pairs,nfilas,w)

      write(*,*)'ss555'
      
c      do i=1,nmax
c         do j=1,nfilas(i)
c            write(*,*)i,j,pairs(j,i),w(j,i),
c     +           mspace(2,i),mspace(4,i),
c     +           mspace(2,pairs(j,i)),mspace(4,pairs(j,i)),mrij(j,i)
c         enddo
c      enddo
      
      call presioni(mspace,ntype)

      write(*,*)'ss666'
      
      call momento(mspace,ntype,npairs,pairs,nfilas,w,dwdx,mrij,mxij,
     +     mvij,mvxij,pairsv,nfilasv,mrijv,mxijv)

      write(*,*)'ss777'
      
c      do i=1,ntype(1)
c         do j=1,nfilas(i)
c            dx(1) = mspace(2,i)-mspace(2,pairs(j,i))
c            dx(2) = mspace(4,i)-mspace(4,pairs(j,i))
c            rrr = sqrt(dx(1)**2 + dx(2)**2)
c            call kernel(rrr,dx,mspace(13,i),ww,dwdx0)
c            if(i.eq.450)write(*,*)i,j,pairs(j,i),mspace(2,pairs(j,i)),
c     +           mspace(4,pairs(j,i)),w(j,i),ww,dwdx(1,j,i),dwdx0(1),
c     +           dwdx(3,j,i),dwdx0(2),rrr
c     enddo
c        write(*,*),mspace(20,i),mspace(21,i),mspace(22,i)
c      enddo
      
      end


