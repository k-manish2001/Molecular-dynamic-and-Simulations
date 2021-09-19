! subroutine to compute charge dipole moments
!*************************************************************************************
  SUBROUTINE compute_charge_dipole(dipx_box,dipy_box,dipz_box, nmol,acount, x,y,z, boxx,boxy,boxz)
  IMPLICIT NONE
  INTEGER       :: i, j, k, count, nmol, nsite = 3, acount
  REAL (kind=8), DIMENSION(nmol) :: dipx_box, dipy_box, dipz_box
  REAL (kind=8), DIMENSION(acount) :: x, y, z
  REAL*8           :: xx, yy, zz, dx, dy, dz
  REAL*8          :: boxx, boxy, boxz, dipx, dipy, dipz
  REAL*8          :: chgOxy = -1.1128, chgHyd = 0.5564
        
     dipx_box (:) = 0.0
     dipy_box (:) = 0.0
     dipz_box (:) = 0.0

        count = 0 
        
        DO i = 1, acount, 3  ! OW first atom id  
           dipx = 0.0
           dipy = 0.0 
           dipz = 0.0
           count = count + 1
           dipx = dipx + x(i) * chgOxy
           dipy = dipy + y(i) * chgOxy
           dipz = dipz + z(i) * chgOxy
           Do j = 1, nsite-1
              dx = 0.0 
              dy = 0.0
              dz = 0.0
              dx = x(i+j) - x(i)   
              dy = y(i+j) - y(i)  
              dz = z(i+j) - z(i)            
              dx = dx - boxx * ANINT (dx/boxx)
              dy = dy - boxy * ANINT (dy/boxy)
              dz = dz - boxz * ANINT (dz/boxz)
              xx = x(i) + dx
              yy = y(i) + dy
              zz = z(i) + dz
              dipx = dipx + xx * chgHyd
              dipy = dipy + yy * chgHyd
              dipz = dipz + zz * chgHyd
            END DO
            dipx_box (count) = dipx_box (count) + dipx
            dipy_box (count) = dipy_box (count) + dipy
            dipz_box (count) = dipz_box (count) + dipz
           END DO
       ! print*, count, dipx_box (count),dipy_box (count),dipz_box (count)
        RETURN

 End Subroutine compute_charge_dipole

 program main
 Implicit none
 Integer  :: i,j,k, natoms,atomType,nmol, nframe,iframe
 Integer  :: nwater, acount, aId, mID, tId,m,Otype,Htype, satoms
 Real*8   :: tempX, tempY, tempZ
 Real*8   :: xlo,xhi,ylo,yhi,zlo,zhi, boxx,boxy,boxz
 real(kind=8),allocatable,dimension(:,:,:) :: dipdata 
 real(kind=8),allocatable,dimension(:)     :: dipdata_box, gk, dipMx, dipMy, dipMz 
 real(kind=8),allocatable,dimension(:)     :: x, y, z, dipx_box, dipy_box, dipz_box, corrdata  
 CHARACTER(LEN=30)  :: Format1,Format2,Format3

 Format1 = "(3F10.5)"
 Format2 ="(I10,3F10.5)"
 Format3 ="(I10,2F10.5)"

  nframe=501
  iframe=200
  natoms=6000
  satoms=0  ! if surface is present in id 1 to ...
  Otype=2 ! OW type  First id
  Htype=1 ! HW type
  nmol = 2000  ! no of H2O

 Open(100,file = 'out_dminst.dat', status='replace')
 Open(200,file = 'out_dmtotal.dat', status='replace')
 open(300,file = 'out_dmeachmolec.dat', status='replace')

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

 OPEN(UNIT=20,FILE="dump.lammpstrj", status="unknown")

 do i = 1, nframe-iframe
     do k=1,natoms+9
        read(20,*)
     end do
   end do
print*,"skip conf ", nframe-iframe 

 ALLOCATE (x(natoms), y(natoms), z(natoms))
 ALLOCATE (dipx_box(nmol), dipy_box(nmol), dipz_box(nmol))
 ALLOCATE (dipdata(iframe,nmol,3)) 
 ALLOCATE (dipdata_box(iframe),gk(iframe))
 ALLOCATE (dipMx(iframe),dipMy(iframe),dipMz(iframe))
  dipdata(:,:,:) = 0.0

  do k = 1, iframe
    read(20,*) 
    read(20,*) 
    read(20,*) 
    read(20,*) natoms
    read(20,*) 
    read(20,*) xlo, xhi
    read(20,*) ylo, yhi
    read(20,*) zlo, zhi
    read(20,*) 
    boxx=xhi-xlo
    boxy=yhi-ylo
    boxz=zhi-zlo
    
    x(:) = 0.0
    y(:) = 0.0
    z(:) = 0.0
    acount = 0
    do i = 1, natoms
      read(20,*) aId,tID,tempX,tempY,tempZ
      If (tID==Otype .or. tID==Htype) then
        x(aID) = tempX
        y(aID) = tempY
        z(aID) = tempZ
        acount = acount + 1
      End if 
    end do

  nwater = acount/3
  !nmol = nwater
  If(nwater .ne. nmol) then 
     write(*,*)'nwater ne nmol'
     exit
  End if
!print*, acount, nmol
 !write(*,*)'computing charge'
call compute_charge_dipole(dipx_box,dipy_box,dipz_box, nmol,acount, x,y,z, boxx,boxy,boxz)

! the dipole components are stored for all the snapshots 
 !   print*, "nmol=", nmol
    DO j = 1, nmol
      dipdata (k, j, 1) = dipx_box (j) 
      dipdata (k, j, 2) = dipy_box (j) 
      dipdata (k, j, 3) = dipz_box (j) 
    END DO
!    print*,k,j, dipdata (k, 2000, 1),dipdata (k, 2000, 2),dipdata (k, 2000, 3)

 END DO ! end of loop over trajectories


! Output file, to compute the autocorrelation 
 !  write(*,*)'printing dipole'

 do k = 1, iframe
   write (300,*) k
   do j=1, nmol
     write (300,Format1) dipdata(k,j,1), dipdata(k,j,2), dipdata(k,j,3)
    !write (300,*) j, dipdata(k,j,1), dipdata(k,j,2), dipdata(k,j,3)
   end do
 end do

  dipdata_box(:) = 0.0 
  dipMx(:) = 0.0
  dipMy(:) = 0.0
  dipMz(:) = 0.0
  gk(:) = 0.0
  write (100,*) iframe
  Do k = 1, iframe
        !write (100,*) k, nmol 
        Do j = 1, nmol
          dipdata_box(k) = dipdata_box(k) + sqrt( (dipdata(k, j, 1))**2 + (dipdata(k, j, 2))**2 + (dipdata(k, j, 3))**2 )
          dipMx(k) = dipMx(k) + dipdata(k, j, 1)
          dipMy(k) = dipMy(k) + dipdata(k, j, 2)
          dipMz(k) = dipMz(k) + dipdata(k, j, 3)
           !Write (100,*) j, dipdata (k, j, 1), dipdata (k, j, 2), dipdata (k, j, 3) 
        End Do
        Write (100,Format2) k, dipMx(k), dipMy(k), dipMz(k) 
        !kirkwood constant
        gk(k) = gk(k) + (dipMx(k)*dipMx(k) + dipMy(k)*dipMy(k) + dipMz(k)*dipMz(k))/(dipdata_box(k)*dipdata_box(k)/nmol)
        ! Calculation for the total system dipole
        Write (200,Format3) k, (dipdata_box(k))/(nmol), gk(k)
        !dipdata_box = sum(dipdata, DIM = 2) !sum of dipoles over all molecular species
  End Do 

 write(*,*) 'individual dipole and gk factor= ', (sum(dipdata_box)*16)/(3.336*nmol*iframe) , sum(gk)/iframe

 close(100)
 Close(200)
 close (300)

 DEALLOCATE (dipx_box, dipy_box, dipz_box, dipdata, x, y , z, dipdata_box)

End Program 


