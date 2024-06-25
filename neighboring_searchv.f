
      subroutine neighboring_searchv(rdomain,mspace,ntype,npairsv,pairsv
     +     ,nfilasv,mrijv,mxijv)

c----------------------------------------------------------------------
c     Subroutine to determine the neighboring of each particle

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
c     mspace(14,i) = c-- sound velocity of particle                [in]
c     mspace(15,i) = s-- entropy of particles                      [in]
c     mspace(16,i) = e-- total energy of particles                 [in]
c     mspace(17 to 19,i) = dx  : dx = vx = dx/dt                   [in]
c     mspace(20 to 22,i) = dvx = dvx/dt, force per unit mass       [in]
c     mspace(23,i) = du        : du = du/dt                        [in]
c     mspace(24,i) = ds        : ds = ds/dt                        [in]
c     mspace(25,i) = drho      : drho = drh,o/dt                   [in]


c     t         : Temperature                                      [in]
c     tdsdt     : Production of viscous entropy t*ds/dt            [in]
c     av        : Monaghan average velocity                        [in]
c     rdomain   : smoothing length                                 [in]      

c     npairs    : maximum number of pairs interaction             [out]
c     nfilas    : number of interaction for each fluid particle   [out]
c     mrij      : matrix of ri-rj for all fluid particles for each interaction [out]
c     mxiij      : matrix of xi-xj for all fluid particles for each interaction [out]         

      implicit none
      include 'param.inc'

      integer ntype(2),npairsv
      integer pairsv(npairsv,ntype(1)),nfilasv(ntype(1)),i,j,d
      double precision mspace(25,nmax),rdomain
      double precision rij,xij,yij,zij
      double precision mrijv(npairsv,ntype(1)),mxijv(3,npairsv,ntype(1))

      mrijv = 0.
      mxijv = 0.
c      mvij = 0.
c      mvxij = 0.

c      nfilasv = 0
      
      do i = 1,ntype(1)-ntype(2)
c         if(i.eq.109)write(*,*)i,j,mspace(2,i),mspace(4,i)
         nfilasv(i) = 0
         do j = ntype(1)-ntype(2)+1,nmax!ntype(1)+ntype(2)
            if(i.ne.j)then
               call radioij(mspace(2,i),mspace(3,i),mspace(4,i),
     +              mspace(2,j),mspace(3,j),mspace(4,j),rij,xij,yij,zij)
c               call velij(mspace(5,i),mspace(6,i),mspace(7,i),
c     +              mspace(5,j),mspace(6,j),mspace(7,j),
c     +              vij,vxij,vyij,vzij)
               if(rij.le.rdomain)then
                  nfilasv(i) = nfilasv(i) + 1
                  if(nfilasv(i).gt.npairsv)then
                     write(*,*)'neighboring_searchv STOP',i,j,nfilasv(i)
     +                    ,ntype
                     stop
                  endif
                  pairsv(nfilasv(i),i) = int(mspace(1,j))
                  mrijv(nfilasv(i),i) = rij
                  mxijv(1,nfilasv(i),i)  = xij
                  mxijv(2,nfilasv(i),i)  = yij
                  mxijv(3,nfilasv(i),i)  = zij
c                  write(*,*)i,j,mrijv(nfilasv(i),i),
c     +                 mxijv(1,nfilasv(i),i),mxijv(2,nfilasv(i),i),
c     +                 mxijv(3,nfilasv(i),i),mspace(12,i),mspace(12,j)
c                   if(i.eq.109)write(*,*)i,j,mspace(2,int(mspace(1,j))),
c     +                 mspace(4,int(mspace(1,j)))
               endif
            endif
         enddo
c         nfilas(i) = nfilas(i) - 1
         
      enddo
      
      end
