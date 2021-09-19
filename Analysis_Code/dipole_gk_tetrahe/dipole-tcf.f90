       Program acf_velocity    
       integer t,t0, tt0,IT0,T0MAx,TMAx,dt,Ntel,Npart,frames
       integer i,j
       PARAMETER (T0MAx=5,TMax=50,Npart=2000)
       real*8 time0(T0MAx),vxt(TMAx)
       real*8 vx0(Npart,T0Max),vy0(Npart,T0Max),vz0(Npart,T0Max)
       real*8 r2t(Tmax),delt,norm
       integer ntime(TMAx)
       double precision vx(Npart),vy(Npart),vz(Npart)
       double precision x(Npart),y(Npart),z(Npart)

       !Format1="(2F10.5)"
!c       frames---- no of data points 
!c       IT0 ----- frequncy of to
!c       Tmax=T0Max*IT0 

        frames=200 !3000fs or 3 ps
        delt=0.01 
        !Npart=2000

        IT0=10
        t0=0
        !Tmax=T0Max*IT0

        Do i=1,TMax
          vxt(i)=0.0
          ntime(i)=0
        END DO 

        open(unit=10,file='out_dmeachmolec.dat',status='old')
        
        Do Ntel=1,frames
          read(10,*)
         !read(10,*)Npart  ! skeping 1st line of which are description
          DO i=1,Npart
            read(10,*) vx(i),vy(i),vz(i)
          end do

        IF (MOD(Ntel,IT0).EQ.0) THEN
        !        print*, "Ntel=", Ntel,IT0
!c         --new t=0
            t0 = t0 + 1           !Number of initial points t=0
            tt0 = MOD(t0-1,T0Max) + 1
            time0(tt0) = Ntel
        !    print*, Ntel,t0,tt0,time0(tt0)
            Do i=1,Npart
            vx0(i,tt0) = vx(i) !Store velocity for given t=0
            vy0(i,tt0) = vy(i)
            vz0(i,tt0) = vz(i)
            END DO 
                  
        END IF
        !print*,MIN(t0,T0Max)
        DO t = 1, MIN(t0,T0Max)    !update correlation for t=0
         dt = Ntel - time0(t) + 1  !actual time minus time t=0
      !   print*,"t&MIN(t0,T0Max)", t, MIN(t0,T0Max)
         IF (dt.LE.TMAx) THEN
          ntime(dt) = ntime(dt) + 1
       !   print*, dt, ntime(dt)
          Do i=1,Npart
          vxt(dt)= vxt(dt)+vx(i)*vx0(i,t)+vy(i)*vy0(i,t)+vz(i)*vz0(i,t)
          END DO 
         END IF
        END DO
       
       !write(*,*) 'Frame= ',Ntel
        END DO   !number of frames ended here
        close(10)

        open(unit=20,file='out_tcf-smor.dat',status='unknown')
     
        norm=vxt(1)/(Npart*ntime(1))
        DO t=1,TMAx
        !print*,ntime(t)
       !   vxt(t)=vxt(t)/vxt(1)
       !   write(20,*) (t-1)*delt,vxt(t)
          vxt(t)=vxt(t)/(Npart*ntime(t))/norm
         write(20,"(F7.3,F10.7)") (t-1)*delt,vxt(t)
        END DO 

        close(20)

        Stop
        end 
