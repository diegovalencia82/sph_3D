
c      subroutine time_integration(id, x,vx, mass, rho, p, u, c, s, e, 
c     &          itype, hsml, ntotal, dt )
      subroutine time_integration(mspace,ntotal,ntype,dt)
c----------------------------------------------------------------------


c     ntype-- total particle number for components                   [in]
      
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

c     av        : Monaghan average velocity      
      
c     maxtimestep-- maximum timesteps                           [input]
c     dt-- timestep                                             [input]

      
      implicit none 
      include 'param.inc'

      integer i,j,l,ntype(2),ntotal,itimestep,nf,rangg(2),itoutfile,
     +     nfluid
      double precision t(maxtimestep),dt,av
      double precision mspace(25,nmax),v_m(3,nmax)
      character infilebas*80,outfile(10000)*80
      
c     --  Make an array of outfiles using a base outfile
      itoutfile = 0
      infilebas = 'snapshot'
      rangg(1)=1
      rangg(2)=20000
      call array_infilebase(infilebas,rangg,outfile,1,nf)  
      
c---  Definition of variables derivates

      itimestep = initime
      itoutfile = inioutfile
c      t = 0.0

 10   if(itimestep.ne.1)then
         do i=1,ntype(1)-ntype(2)
            v_m(1,i) = mspace(5,i)
            v_m(2,i) = mspace(6,i)
            v_m(3,i) = mspace(7,i)
            mspace(5,i) = mspace(5,i) + (dt0/2.)*mspace(20,i)
            mspace(6,i) = mspace(6,i) + (dt0/2.)*mspace(21,i)
            mspace(7,i) = mspace(7,i) + (dt0/2.)*mspace(22,i)
c            mspace(2,i) = mspace(2,i) + dt0 * mspace(5,i)
c            mspace(4,i) = mspace(4,i) + dt0 * mspace(7,i)
         enddo
      endif
      
! 10   call single_step(mspace,ntype,itimestep,dt,av)
      call single_step(mspace,ntype,itimestep,dt,av)

!      Método de Euler ----------10      
!      do i=1,ntype(1)-ntype(2)
!         xy(1) = mspace(2,i)
!         xy(2) = mspace(4,i)
!         xy(3) = mspace(5,i)
!         xy(4) = mspace(7,i)

c     call derivadas(i,t(i),xy,dxydt,mspace)

!         xy(3) = xy(3) + dt0*mspace(20,i)
!         xy(1) = xy(1) + dt0*xy(3)
!         xy(4) = xy(4) + dt0*mspace(22,i)
!         xy(2) = xy(2) + dt0*xy(4)
         
c         call rk4_dd(xy,dxydt,4,t(i),dt0,xy,derivadas,mspace,i)

!         mspace(2,i) = xy(1)
!         mspace(4,i) = xy(2)
!         mspace(5,i) = xy(3)
!         mspace(7,i) = xy(4)
c         write(*,*)i,xy,mspace(20,i),mspace(22,i)
!      enddo
!      fin método de Euler
      
      if(itimestep.eq.1)then
         do i=1,ntype(1)-ntype(2)
            mspace(5,i) = mspace(5,i) + (dt0/2.)*mspace(20,i)
            mspace(6,i) = mspace(6,i) + (dt0/2.)*mspace(21,i)
            mspace(7,i) = mspace(7,i) + (dt0/2.)*mspace(22,i)
            mspace(2,i) = mspace(2,i) + dt0 * mspace(5,i)
            mspace(3,i) = mspace(3,i) + dt0 * mspace(6,i)
            mspace(4,i) = mspace(4,i) + dt0 * mspace(7,i)
         enddo
      else
         do i=1,ntype(1)-ntype(2)
            mspace(5,i) = v_m(1,i) + (dt0/2.)*mspace(20,i)
            mspace(6,i) = v_m(2,i) + (dt0/2.)*mspace(21,i)
            mspace(7,i) = v_m(3,i) + (dt0/2.)*mspace(22,i)
            mspace(2,i) = mspace(2,i) + dt0 * mspace(5,i)
            mspace(3,i) = mspace(3,i) + dt0 * mspace(6,i)
            mspace(4,i) = mspace(4,i) + dt0 * mspace(7,i)
         enddo
      endif
      
      itimestep = itimestep + 1
      t(itimestep) = itimestep * dt0

      if (mod(itimestep,save_step).eq.0) then
         itoutfile = itoutfile + 1
         write(*,*)'Saving file = ',outfile(itoutfile)
         write(*,*)'iteration',itimestep,'time = ',t(itimestep)
         
         l = len_trim(outfile(itoutfile))
         open(1,file=outfile(itoutfile))
         nfluid = ntype(1)-ntype(2)
         write(1,*)itimestep,t(itimestep),ntype(1),nfluid,ntype(2)
      
         do i=1,nmax  
c     write(1, *) i, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim),
c     +        mass (i), rho(i), p(i), u(i), itype(i), hsml(i)

c id, x, y ,z, vx, vy, vz, mass, rho, p, u, itype, hsml
            
          write(1,*)int(mspace(1,i)),real(mspace(2,i)),real(mspace(3,i))
     +          ,real(mspace(4,i))
     +          ,real(mspace(5,i)),real(mspace(6,i)),real(mspace(7,i))
     +          ,real(mspace(8,i)),real(mspace(9,i)),real(mspace(10,i))
     +          ,real(mspace(11,i)),int(mspace(12,i)),real(mspace(13,i))
         enddo
         
         write(*,*)'**********************************************'
      endif
      
      if(itimestep.LT.fintime)goto 10
! ===== END OF INTEGRATE======
      
      end

      
c==========================================================================

c
c THIS SUBROUTINE MAKE AN ARRAY OF INFILES USING A BASE OF INFILE.

      SUBROUTINE array_infilebase(infilebas,rangg,infile,step,nf)
      implicit none

      integer n,rangg(2),i,l,j,step,nf,k
      character infilebas*80,infile(10000)*80,ch1,ch2*2,ch3*3,ch4*4
      character ch5*5,ch6*6

      l = len_trim(infilebas)
      k = 0
      
      j=rangg(1)-1
      DO 10 i=1,10000
         j=j+1
c         if(mod(i,step).eq.1)then
            if(j.le.rangg(2))then
               k=k+1
               if(j.lt.10)write(ch1,'(I1)')j
               if(j.lt.10)infile(k)=infilebas(1:l)//'_000'//ch1
               
               if(j.ge.10.and.j.lt.100)write(ch2,'(I2)')j
               if(j.ge.10.and.j.lt.100)
     +              infile(k)=infilebas(1:l)//'_00'//ch2
               
               if(j.ge.100.and.j.lt.1000)write(ch3,'(I3)')j
               if(j.ge.100.and.j.lt.1000)
     +              infile(k)=infilebas(1:l)//'_0'//ch3
               
               if(j.ge.1000.and.j.lt.10000)write(ch4,'(I4)')j
               if(j.ge.1000.and.j.lt.10000)
     +              infile(k)=infilebas(1:l)//'_'//ch4
               
               if(j.ge.10000.and.j.lt.100000)write(ch5,'(I5)')j
               if(j.ge.10000.and.j.lt.100000)
     +              infile(k)=infilebas(1:l)//'_'//ch5
               
               if(j.ge.100000.and.j.lt.1000000)write(ch6,'(I6)')j
               if(j.ge.100000.and.j.lt.1000000)
     +              infile(k)=infilebas(1:l)//'_'//ch6
               
c     write(*,*)i,infile(i),j
            endif
c         endif
 10   CONTINUE
      
      nf = k
      
      RETURN
      END      
