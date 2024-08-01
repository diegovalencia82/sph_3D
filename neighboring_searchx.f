
      subroutine neighboring_searchx(rdomain,mspace,ntype,npairs,pairs,
     +     nfilas,mrij,mxij,mvij,mvxij,itimestep)

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

      integer ntype(2),sel,npairs,itimestep
      integer pairs(npairs,ntype(1)),nfilas(ntype(1)),i,j,l
      integer ind_x(ntype(1)),ind_z(ntype(1)),init,ifin,window_size,k,h
      double precision mspace(25,nmax)
      double precision rij,xij,yij,zij
      double precision mrij(npairs,ntype(1)),mxij(3,npairs,ntype(1))
      double precision vij,vxij,vyij,vzij,rdomain
      double precision mvij(npairs,ntype(1)),mvxij(3,npairs,ntype(1))
      double precision x(ntype(1)),z(ntype(1)),r(ntype(1)),dlt

      dlt = rdomain*1.2
      mrij = 0.
      mxij = 0.
      mvij = 0.
      mvxij = 0.
      nfilas = 0
c      ind_x = 0
      
c      write(*,*)'Ordena'
      do i = 1,ntype(1)
         x(i) = mspace(2,i)
         ind_x(i) = int(mspace(1,i))
      enddo

! Llamar a la subrutina indexx para ordenar por r
c      if(itimestep.gt.0)call indexx(ntype(1),x,ind_x)
      
      init = 1
      ifin = ntype(1)
      window_size = int(90*ntype(1)/100.)

c      write(*,*)'Busca vecinos'
      do k = 1,ntype(1) - ntype(2)
         i = ind_x(k)
         if(mspace(12,i).eq.1.)then
c            call calcular_rango(i,ntype(1),window_size,init,ifin)
            nfilas(i) = 0
            do h = 1,ntype(1)!init,ifin!ntype(1)   !1,ntype(1)!
               j = ind_x(h)
               if(mspace(2,j).gt.mspace(2,i)-dlt.and.mspace(2,j)
     +              .lt.mspace(2,i)+dlt)then
               if(mspace(4,j).gt.mspace(4,i)-dlt.and.mspace(4,j)
     +                 .lt.mspace(4,i)+dlt)then
               if(mspace(3,j).gt.mspace(3,i)-dlt.and.mspace(3,j)
     +                 .lt.mspace(3,i)+dlt)then
               if(int(mspace(1,i)).ne.int(mspace(1,j)))then
                  call radioij(mspace(2,i),mspace(3,i),mspace(4,i),
     +              mspace(2,j),mspace(3,j),mspace(4,j),rij,xij,yij,zij)
                  call velij(mspace(5,i),mspace(6,i),mspace(7,i),
     +                 mspace(5,j),mspace(6,j),mspace(7,j),
     +                 vij,vxij,vyij,vzij)
                  if(rij.lt.rdomain)then
                     nfilas(i) = nfilas(i) + 1
                     if(nfilas(i).gt.npairs)then
                        write(*,*)'neighboring_search.f STOP, i=',i
     +                       ,'j=',j,'nfilas = ',nfilas(i)
                        stop
                     endif
                     pairs(nfilas(i),i) = int(mspace(1,j))
                     mrij(nfilas(i),i) = rij
                     mxij(1,nfilas(i),i)  = xij
                     mxij(2,nfilas(i),i)  = yij
                     mxij(3,nfilas(i),i)  = zij
                     mvij(nfilas(i),i) = vij
                     mvxij(1,nfilas(i),i)  = vxij
                     mvxij(2,nfilas(i),i)  = vyij
                     mvxij(3,nfilas(i),i)  = vzij
c             write(*,*)l,'00',i,pairs(nfilas(i),i),nfilas(i),mspace(1,j)
c                write(*,*)i,j,rij,xij,yij,zij,mspace(12,i),mspace(12,j),
c     +           mspace(1,int(mspace(1,i))),mspace(1,int(mspace(1,j)))
c     write(*,*)i,j,mrij(nfilas(i),i),mxij(1,nfilas(i),i),
c     +                 mxij(2,nfilas(i),i),mxij(3,nfilas(i),i),  
c     +                 mvij(nfilas(i),i),mvxij(1,nfilas(i),i),
c     +                 mvxij(2,nfilas(i),i),mvxij(3,nfilas(i),i)
c     if(i.lt.3)write(*,*)'11111',i,mspace(2,int(mspace(1,j))),
c     +                 mspace(4,int(mspace(1,j)))
                  endif
               endif
               endif
               endif
               endif
            enddo
         endif
      enddo

c      do k=1,ntype(1)-ntype(2)
c         i = ind_x(k)
cc         if(nfilas(k).eq.0)write(*,*)k,nfilas(k),i,mspace(1,i)
c         do h=1,nfilas(k)
c            j = ind_x(h)
c            l = pairs(h,k)
c            write(*,*)'ii',k,pairs(h,k),h,mspace(1,l),'qqq',nfilas(k)
c         enddo
c     enddo
      
      end
c     ----------------------------------------------------------------------

      subroutine radioij(xi,yi,zi,xj,yj,zj,rij,xij,yij,zij)

      implicit none
      include 'param.inc'

      double precision xi,yi,zi,xj,yj,zj,xij,yij,zij,rij

      
      if(dim.eq.2)then
         xij = xi-xj
         yij = 0.0d0
         zij = zi-zj
         rij = sqrt(xij*xij + zij*zij)
      endif

      if(dim.eq.3)then
         xij = xi-xj
         yij = yi-yj
         zij = zi-zj
         rij = sqrt(xij*xij + yij*yij + zij*zij)
      endif
      
      end 

c     ----------------------------------------------------------------------

      subroutine velij(vxi,vyi,vzi,vxj,vyj,vzj,vij,vxij,vyij,vzij)

      implicit none
      include 'param.inc'

      double precision vxi,vyi,vzi,vxj,vyj,vzj,vxij,vyij,vzij,vij

      
      if(dim.eq.2)then
         vxij = vxi-vxj
         vyij = 0.0
         vzij = vzi-vzj
         vij = sqrt(vxij*vxij + vzij*vzij)
      endif

      if(dim.eq.3)then
         vxij = vxi-vxj
         vyij = vyi-vyj
         vzij = vzi-vzj
         vij = sqrt(vxij*vxij + vyij*vyij + vzij*vzij)
      endif
      
      end       
      
c     --------------------------------------------------------------------


      subroutine calcular_rango(i, n, window_size, init, ifint)
      implicit none
      
      integer, intent(in) :: i, n, window_size
      integer, intent(out) :: init, ifint
      integer :: half_window
     
      half_window = (window_size - 1) / 2
      
      if (i <= half_window + 1) then
         init = 1
         ifint = window_size
      elseif (i >= n - half_window) then
         init = n - window_size + 1
         ifint = n
      else
         init = i - half_window
         ifint = i + half_window
      end if
      end subroutine calcular_rango
