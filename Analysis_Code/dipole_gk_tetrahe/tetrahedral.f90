     Program TD
        implicit none
        integer:: i,j,k,m, natoms,Nmol, nframe,iframe,tbin 
        integer:: AId,Atype,Atomtype,nbin,acount,nk,nn,ii,jj,nj
        real*8:: xlo,xhi,ylo,yhi,zlo,zhi
        real*8::cxi,cyi,czi,szk,xij,yij,zij,dij,xik,yik,zik,dik,sjk,cosjik
        real*8::v38,v13,deltD,Tdi,Dist, boxx,boxy,boxz, x,y,z
        integer,dimension(-100:0):: nDistTd
        integer,dimension(0:100):: pDistTd
        real(kind=8),dimension(8000):: cmx,cmy,cmz !allocation for 8000 molec
        integer,dimension(4,8000) ::Neigh
        
     nframe=501
     iframe=200
     natoms=6000
     AtomType=2  ! Ow type 
     delTd = 0.01d0 ! <- increment for histogram of Td
     tbin = int(1.0/delTd)

     nDistTd(:)=0; pDistTd(:)=0 
!     npdistTd(:)=0.0; ppdistTd(:)=0.0
     cmx(:)=0.0; cmy(:)=0.0; cmz(:)=0.0
     
     v13 = 1.0d0 / 3.0d0
     v38 = 0.375d0 !=3.0d0 / 8.0d0

!   nframe=0
!   open(unit=22,file='dumpnvt.lammpstrj',status='old')
!    do while(.not. EOF(22))
!     read(22,*)
!     read(22,*)
!     read(22,*)
!     read(22,*) natoms
!     do k=1,natoms+5
!        read(22,*)
!     end do
!     nframe=nframe+1
!   print*, "nframe= ", nframe
!  end do
!close(22)        

OPEN(UNIT=22,FILE="dump", status="unknown")
  do i = 1, nframe-iframe
     do k=1,natoms+9
        read(22,*)
     end do
  end do
print*,"skip conf ", nframe-iframe

  do jj = 1, iframe
    read(22,*)
    read(22,*)
    read(22,*)
    read(22,*) natoms
    read(22,*)
    read(22,*) xlo, xhi
    read(22,*) ylo, yhi
    read(22,*) zlo, zhi
    read(22,*)
     boxx=xhi-xlo
     boxy=yhi-ylo
     boxz=zhi-zlo
     acount = 0
    do i=1,natoms
       read(22,*) AId,Atype, x, y, z
       if (Atype == AtomType) then
         acount=acount+1
         cmx(acount) = x !*Lx
         cmy(acount) = y !*Ly
         cmz(acount) = z !*Lz
       end if
    end do
    print*, "Frame & acount=nmol= ", jj, acount  
    nmol=acount
    call nearestmole (acount,cmx,cmy,cmz,nmol,Neigh,boxx,boxy,boxz)

    do i = 1, acount !Nmol !for each molecule
      cxi =-cmx(i)
      cyi =-cmy(i)
      czi =-cmz(i)
      sjk = 0.0d0
        
      do nj = 1, 3
        j = Neigh(nj,i) 
        xij=cmx(j)+cxi 
        xij=xij-boxx*anint(xij/boxx)
        yij=cmy(j)+cyi
        yij=yij-boxy*anint(yij/boxy)
        zij=cmz(j)+czi
        zij=zij-boxz*anint(zij/boxz)
        dij=xij*xij+yij*yij+zij*zij
        do nk = nj+1, 4
          k=Neigh(nk,i)
          xik=cmx(k)+cxi
          xik=xik-boxx*anint(xik/boxx)
          yik=cmy(k)+cyi
          yik=yik-boxy*anint(yik/boxy)
          zik=cmz(k)+czi
          zik=zik-boxz*anint(zik/boxz)
          dik = xik*xik+yik*yik+zik*zik
          cosjik =(xij*xik+yij*yik+zij*zik)/dsqrt(dij*dik)
          sjk = sjk +(cosjik + v13)**2
 !       write (*,*) 'costheta=',dik,dij,i,k,j,xik,yik,zik
        end do
      end do
      Tdi = 1.0d0 - v38*sjk ! <- tetrahedrality of molecule i
      If (Tdi .lt. 0) then
         nbin = int(Tdi/delTd)
         nDistTd(nbin) = nDistTd(nbin) + 1.0d0
      Endif 
      If (Tdi .ge. 0) then
         nbin = int(Tdi/delTd)
         pDistTd(nbin) = pDistTd(nbin) + 1.0d0
      Endif

    end do  ! end of loop for molecule

   end do  ! Frame ends here jj

   close(22)

 OPEN(UNIT=33,FILE="out_distetra.dat", status="replace")   
       Dist=0.0 
       do ii=-tbin+1, 0
         Dist= real(nDistTd(ii))/real(acount*iframe)
         write(33,'(F10.4,F10.6)') (ii-0.5)*delTd, Dist
       end do 

       do ii=0,tbin-1
          Dist= real(pDistTd(ii))/real(acount*iframe)
          write(33,'(F10.4,F10.6)') (ii+0.5)*delTd, Dist
       end do

  close(33)

 end program TD


!The following lines are for the calculation of the nearest 4 molecules of i
! labi and labj are list vectors with the size of Npair

subroutine nearestmole(acount,cmx,cmy,cmz,nmol,Neigh,boxx,boxy,boxz)
      implicit none
      integer:: i,nn,np,j,k,kk,Nsmall,Ndata
      integer:: Npair,ij,ii,jj,Nmol,acount
      real*8:: rij,ddx,ddy,ddz,boxx,boxy,boxz
      integer,allocatable,dimension(:) ::List,labi,labj
      integer,dimension(4,8000)::Neigh
      real(kind=8),allocatable,dimension(:)::Val
      real(kind=8),dimension(8000,8000):: Dist
      real(kind=8),dimension(8000):: cmx,cmy,cmz
        
!The following lines are for the list vect
     np=acount*nmol
     allocate(labi(np),labj(np))
     labi(:)=0
     labj(:)=0
     allocate(List(nmol),Val(nmol))

     ! total no. of pair
     do ii=1,nmol
       do jj=1,nmol
         Dist(ii,jj)=0.0
       end do
     end do

     ij = 0
     do i = 1, acount-1
       do j = i+1, nmol
          ij = ij + 1
          labi(ij) = i
          labj(ij) = j
        end do
     end do
     Npair=ij

      do ij = 1,  Npair
        i = labi(ij)
        j = labj(ij)
        ddx=cmx(i) - cmx(j) !distance between O-O
        ddx=ddx - boxx * anint(ddx/boxx) !PBC in x dire
        ddy=cmy(i) - cmy(j)
        ddy = ddy - boxy * anint(ddy/boxy)
        ddz = cmz(i) - cmz(j)
        ddz = ddz - boxz * anint(ddz/boxz)
        rij = dsqrt(ddx*ddx+ddy*ddy+ddz*ddz) ! <- distance
        Dist(i,j) = rij
        Dist(j,i) = rij
      end do 
      
      Ndata = nmol
      Nsmall = 4
      do i = 1, acount
         do j = 1, nmol
           List(j) = j
           Val(j) = Dist(j,i)
         end do
        Val(i) = 1.0d10

        call Sort_SmallestN (nmol, Nsmall, List, Val)

        do kk = 1, Nsmall
        Neigh(kk,i) = List(kk) ! <- 4 nearest neighbor mols of i
        end do
       
       end do

      return 
      deallocate(List,Val)

      end subroutine nearestmole


!The following is a subroutine to search molecules with Nsmall smallest values

      subroutine Sort_SmallestN (Ndata, Nsmall, List, Val)
      implicit none ! real(kind(1d0))(a-h,o-z)
      integer::Nsmall,Ndata,k,i,tmp,j
      real*8::itmp
      logical,dimension (8000):: mask
      integer,dimension(Ndata):: List
      integer,dimension(1)::ik
      real(kind=8),dimension(Ndata):: Val
      mask = .true.
      do i = 1, Nsmall
        ik = Minloc(Val,mask) !will return the minimum distance
        itmp = Val(i)
        Val(i) = Val(ik(1))
        Val(ik(1)) = itmp
        tmp=List(i)
        List(i) =List(ik(1))
        List(ik(1)) =tmp
        mask(i) = .false.
     ! write(*,*) 'nearest for first=',List(i),Val(i),ik
        end do

      return
      end subroutine Sort_SmallestN

