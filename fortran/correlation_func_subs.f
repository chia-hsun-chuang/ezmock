c This file include the basic functions for calculating correlation function
c-------------------------------------------

c      SUBROUTINE Fit_nofz( zmin, zmax, 
c     & LRGsourcefile_re, nofzfile)
c      SUBROUTINE xywtoradec(x,y,w,ra,dec,z)
c      SUBROUTINE chooseZ(zmin,zmax,
c     & LRGsourcefile, LRGsourcefile_range)
c      SUBROUTINE radeczToxyz(ra,dec,z,x,y,w)
c      SUBROUTINE XX_corr(
c     & LRGsourcefile_xyz, nofrfile)
c      SUBROUTINE DR_corr(
c     & LRGsourcefile_xyz, RANsourcefile_xyz, nofrfile)
c      SUBROUTINE chooseRange(zmin,zmax,ramin,ramax,decmin,decmax,
c     & LRGsourcefile, LRGsourcefile_range)
c      SUBROUTINE toLambdaEta(
c     & radecfile, lambdaetafile)
c      SUBROUTINE convertToxyz(SigmaM, SigmaX, 
c     & LRGsourcefile_re, LRGsourcefile_xyz, nofzfile)
c      SUBROUTINE radecztoxyz_unit(ra,dec,z,x,y,w)  
c      subroutine show_corr_result(
c     & DDcount_file, DRcount_file, RRcount_file,
c     & DD_DR_corr_file, DD_RR_corr_file, DD_DR_RR_corr_file)
c      subroutine delete_file( filename)
c      subroutine switch_z_of_data(datafile,randfile)
c      subroutine show_corr_result_0216(
c     & DDcount_file, DRcount_file, RRcount_file,
c     & DD_DR_corr_file, DD_RR_corr_file, DD_DR_RR_corr_file,
c     & DR_RR_corr_file, DD_bin_file, DR_bin_file, RR_bin_file) 

c---------------------------------------------
      SUBROUTINE Fit_nofz( zmin, zmax, 
     & LRGsourcefile_re, nofzfile)
c use cubic spline to fit n(z), seperation of z = 1/binN
c -------------I/O------------
c LRGsourcefile_re is the filename of data used
c zmin, zmax are redshift range
c nofzfile is the number of galaxies per redshift not dz^3 
c -------------------------------------------------
      IMPLICIT NONE
      INTEGER count, i, zdatamin, zdatamax
c      INTEGER, parameter :: binN=5000 ! dz = 0.0002
c      INTEGER, parameter :: binN=1000 ! dz = 0.001
c      INTEGER, parameter :: binN = 50 ! dz = 0.02
      INTEGER, parameter :: binN = 100 ! dz = 0.01
      REAL    xdata(binN), nofzdata(binN),y2(binN), nofz(0:60000)
      REAL    zdata, radata, decdata
c... input
      REAL         :: zmin 
      REAL         :: zmax
      CHARACTER*100 :: LRGsourcefile_re
      CHARACTER*100 :: nofzfile 

      zdatamin=Int(zmin*binN)
      zdatamax=Int(zmax*binN)
      do 21 i=1,binN
      nofzdata(i)=0
   21 continue
      do 22 i=0,60000
      nofz(i)=0
   22 continue  
      do 23 i=1,binN
      xdata(i)=1.0/binN*(zdatamin+i-1)+0.5/binN
   23 continue      

c source data file - - - - - - - - - - - - - - - -- - - --------------> (1) 
      OPEN(UNIT=100,FILE = Trim(LRGsourcefile_re), STATUS= 'OLD')
      count=0 
c ---LOOP START---

  100 READ(UNIT=100,fmt=*,END=200) radata, decdata, zdata
      count=count+1
cc      PRINT *, count
      if ((zdata .lt. zmax) .and. (zdata .gt. zmin)) 
     & nofzdata(Int(zdata*binN)-zdatamin+1)=
     & nofzdata(Int(zdata*binN)-zdatamin+1)+1
      GOTO 100
c---END OF LOOP---
  200 call spline(xdata, nofzdata,zdatamax-zdatamin
     & ,1.0e30,1.0e30,y2)

c      call splint(xdata, nofzdata,y2,zdatamax-zdatamin,0.001*300,yy)
      do 24 i= Int(100000*zmin) ,Int(100000*zmax)      
      call splint(xdata, nofzdata,
     & y2,zdatamax-zdatamin,0.00001*i+0.000005, nofz(i))
     
   24 continue   
   
      OPEN(UNIT = 10, FILE = Trim(nofzfile),status='New')
      Write(10,fmt=*) (nofz(i),i=0,60000)
#ifdef debug
      Write(*,fmt=*) "finish count nofz"
#endif
      close(10)
      close(100)
      END SUBROUTINE !Fit_nofz()

c ---------- xywtoradec(x,y,w,ra,dec,z) ------------
      SUBROUTINE xywtoradec(x,y,w,ra,dec,r)
      IMPLICIT NONE 
      REAL    x,y,w,ra,dec,r, pi
      pi=3.1415926
      r=sqrt(x**2+y**2+w**2)
c avoid divide by zero (Arithmetic Exception)   
c When I set .lt. 0.00000001, it shows Arithmetic Exception, too.    
      if ((r .lt. 0.000001) .or. (abs(x) .lt. 0.000001) ) then
      ra=0
      dec=0
      r=0
      return
      END if
  
      dec=asin(w/r)/pi*180      
      if (x .ge. 0) then
        if (y .ge. 0) then
          ra=atan(y/x)/pi*180
        else
          ra=atan(y/x)/pi*180+360
        endif   
      else
        ra=atan(y/x)/pi*180+180
      endif       
      return
      END SUBROUTINE !xywtoradec()
c---------------------------------------

      SUBROUTINE chooseZ(zmin,zmax,
     & LRGsourcefile, LRGsourcefile_range)
c to choose galaxies in some redshift range
c -----------I/O----------------
c zim, zmax set the redshift range
c LRGsourcefile is source file
c LRGsourcefile_range is data in choosed range
c -------------------------------------
      IMPLICIT NONE 
c... input
      REAL         :: zmin
      REAL         :: zmax
      CHARACTER*100 :: LRGsourcefile
      CHARACTER*100 :: LRGsourcefile_range 
c
      REAL          :: ra, dec, z
      integer       :: count = 0, totalcount =0
      count = 0
      totalcount = 0
      Open(unit=10, file=Trim(LRGsourcefile), status="old")
!      Open(unit=20, file=TRIM(LRGsourcefile_range), Status='new')
!      Close(20)
!      Open(unit=20, file=TRIM(LRGsourcefile_range), position="append")
      Open(unit=20, file=TRIM(LRGsourcefile_range))
  100 Read(unit=10, fmt=*, end=200) ra, dec, z
      totalcount = totalcount + 1
      If( z .ge. zmin .and. z .le. zmax) then
        write(20,fmt=*) ra,dec,z
        count = count + 1
      endif
      goto 100
  200 close(10)
      close(20)
#ifdef debug
      Print *,"finish choosing z range" , count ,
     & "data out of ",totalcount, "from ",Trim(LRGsourcefile)
#endif
      END SUBROUTINE !chooseZ()   
c----------------------------------------------

      SUBROUTINE XX_corr(
     & LRGsourcefile_xyz, nofrfile, max_r)
c calculate DD or RR nofr
c  -------------  blueprint ------------
c LRGsourcefile_xyz is xyz data file
c nofrfile is nofr file(DD or RR)
c --------------------------------------
      IMPLICIT NONE
c... input
      CHARACTER*100     :: LRGsourcefile_xyz
      CHARACTER*100     :: nofrfile
      integer           :: max_r  

      INTEGER           :: data_num = 1
      REAL, allocatable :: x(:),y(:),w(:),nr(:),x2(:),y2(:),w2(:),nr2(:) 
      REAL              :: nofr(500), dr, drdr
      integer           :: count1, count2, i, j
      REAL              :: xx, yy, ww, nrr
c count data number
      OPEN(UNIT=100, FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  110 Read(UNIT=100, fmt=*, END=210) xx,yy,ww,nrr
      data_num = data_num+1
      goto 110
  210 close(100)

      allocate(x(data_num))
      allocate(y(data_num))
      allocate(w(data_num))
      allocate(nr(data_num))
      allocate(x2(data_num))
      allocate(y2(data_num))
      allocate(w2(data_num))
      allocate(nr2(data_num))

      count1=0
      count2=0
      i=0
      do i=1,max_r
        nofr(i)=0
      end do
c----read data in-----      
      close(100)
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  100 Read(UNIT=100,fmt=*,END=200) x(count1+1),y(count1+1),w(count1+1),
     &  nr(count1+1)
      count1=count1+1
      goto 100
  200 close(100)
      close(200)
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  300 Read(UNIT=200,fmt=*,END=400) x2(count2+1), y2(count2+1),
     &  w2(count2+1), nr2(count2+1)
      count2=count2+1
      goto 300
  400 close(200) 
      Do i=1,count1
        Do j=i+1,count2
            if((abs(x(i)-x2(j)) .lt. max_r) .and. 
     &         (abs(y(i)-y2(j)) .lt. max_r) .and.
     &         (abs(w(i)-w2(j)) .lt. max_r)) then
               drdr=(x(i)-x2(j))**2+(y(i)-y2(j))**2+(w(i)-w2(j))**2
               if(drdr .lt. max_r**2) then
                  dr=sqrt(drdr)
                  nofr(int(dr+1))=nofr(int(dr+1))+
     &               1.0/(1.0+40000.0*nr(i))/(1.0+40000.0*nr2(j))
               end if
            end if
        end do
      end do

      Open(unit = 300, File=TRIM(nofrfile), status='New')
      close(300)
      Open(unit = 300, File=TRIM(nofrfile), position='APPEND')

      do i = 1, max_r
        write(unit=300,fmt=*) nofr(i)/count1/count1
      end do

      close(300)
#ifdef debug
      Write(*,fmt=*)"Finish xx_corr of file ", LRGsourcefile_xyz,
     & ", sample num = ", count1
#endif
      deallocate(x)
      deallocate(y)
      deallocate(w)
      deallocate(nr)
      deallocate(x2)
      deallocate(y2)
      deallocate(w2)
      deallocate(nr2)

      END SUBROUTINE ! XX_corr() 
c-----------------------------------------

      SUBROUTINE DR_corr(
     & LRGsourcefile_xyz, RANsourcefile_xyz, nofrfile, max_r)
      IMPLICIT NONE
c... input
      CHARACTER*100 :: LRGsourcefile_xyz
      CHARACTER*100 :: RANsourcefile_xyz
      CHARACTER*100 :: nofrfile 
      INTEGER       :: max_r

      INTEGER           :: data_num = 1 , rand_num = 1
      REAL, allocatable :: x(:),y(:),w(:),nr(:),x2(:),y2(:),w2(:),nr2(:) 
      REAL              :: nofr(500), dr, drdr
      integer           :: count1, count2, i, j
      REAL              :: xx, yy, ww, nrr

c count data number
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_xyz),
     & STATUS= 'OLD')
  110 Read(UNIT=100,fmt=*,END=210) xx,yy,ww,nrr
      data_num = data_num+1
      goto 110
  210 close(100)
c count rand number
      OPEN(UNIT=100,FILE=TRIM(RANsourcefile_xyz),
     & STATUS= 'OLD')
  310 Read(UNIT=100,fmt=*,END=410) xx,yy,ww,nrr
      rand_num = rand_num+1
      goto 310
  410 close(100)

      allocate(x(data_num))
      allocate(y(data_num))
      allocate(w(data_num))
      allocate(nr(data_num))
      allocate(x2(rand_num))
      allocate(y2(rand_num))
      allocate(w2(rand_num))
      allocate(nr2(rand_num))
      
      count1=0
      count2=0

      do i=1,max_r
        nofr(i)=0
      end do

c----read data in-----      
      close(100)
      OPEN(UNIT=100,FILE=Trim(LRGsourcefile_xyz), STATUS= 'OLD')
  100 Read(UNIT=100,fmt=*,END=200) x(count1+1),y(count1+1),w(count1+1),
     &  nr(count1+1)
      count1=count1+1
      goto 100
  200 close(100)
      close(200)
      OPEN(UNIT=200,FILE=Trim(RANsourcefile_xyz), STATUS= 'OLD')
  300 Read(UNIT=200,fmt=*,END=400) x2(count2+1), y2(count2+1),
     &  w2(count2+1), nr2(count2+1)
      count2=count2+1
      goto 300
  400 close(200) 
c----------------------

      Do 20 i=1,count1
         Do 30 j=1,count2
            if((abs(x(i)-x2(j)) .lt. max_r) .and.
     &         (abs(y(i)-y2(j)) .lt. max_r) .and.
     &         (abs(w(i)-w2(j)) .lt. max_r)) then
               drdr=(x(i)-x2(j))**2+(y(i)-y2(j))**2+(w(i)-w2(j))**2
               if(drdr .lt. max_r**2) then
                dr=sqrt(drdr)
                nofr(int(dr+1))=nofr(int(dr+1))+1.0/(1.0+40000.0*nr(i))/
     &                        (1.0+40000.0*nr2(j))
               end if
            end if
   30    continue
   20 continue  
  
      Open(unit = 300, File=TRIM(nofrfile), status='New')
      close(300)
      Open(unit = 300, File=TRIM(nofrfile), position='APPEND')
      do i = 1, max_r
      write(unit=300,fmt=*) nofr(i)/count2/count1 
      end do

      close(300)
#ifdef debug
      Write(*,fmt=*)"Finish DR_corr, data sample num = ", count2, 
     & " ,and random sample num = ", count1
#endif
      deallocate(x)
      deallocate(y)
      deallocate(w)
      deallocate(nr)
      deallocate(x2)
      deallocate(y2)
      deallocate(w2)
      deallocate(nr2)

      END SUBROUTINE ! DR_corr()
 
c --------------------------------------------------------------

      SUBROUTINE chooseRange(zmin,zmax,ramin,ramax,decmin,decmax,
     & LRGsourcefile, LRGsourcefile_range)
c choose data from certain range
c -------------------------------------------------------------      

      IMPLICIT NONE

c...input
      REAL          :: zmin
      REAL          :: zmax
      REAL          :: ramin
      REAL          :: ramax
      REAL          :: decmin
      REAL          :: decmax
      CHARACTER*100  :: LRGsourcefile 
      CHARACTER*100  :: LRGsourcefile_range 

      REAL          :: ra,dec,z      
      integer       :: count

      count = 0    
      Open(unit=10,file=TRIM(LRGsourcefile),status="old")
      Open(unit=20,file=TRIM(LRGsourcefile_range),
     & Status='new')
      Close(20)
      Open(unit=20,file=TRIM(LRGsourcefile_range),
     & position="append")
  100 Read(unit=10,fmt=*, end=200) ra, dec, z
      If( z .ge. zmin .and. z .le. zmax .and. 
     & ra .ge. ramin .and. ra .le. ramax .and. 
     & dec .ge. decmin .and. dec .le. decmax) then
      write(20,fmt=*) ra,dec,z
      count = count+1
      endif
      goto 100
  200 close(10)
      close(20)

#ifdef debug
      Print *,"finish choosing ",count," galaxies from z=",
     &  zmin,"~",zmax,", ra=",ramin,"~",ramax,", dec=",decmin,"~",decmax
#endif
      END SUBROUTINE ! chooseRange()                                 
c --------------------------------------------------------------
      
      SUBROUTINE toLambdaEta(
     & radecfile, lambdaetafile)
c change data coordinate to lambda, eta from ra, dec
C In fact, not really convert to lambda, eta coordinate.
c However, the poles is correct.
      IMPLICIT NONE
c... input
      CHARACTER*100 :: radecfile
      CHARACTER*100 :: lambdaetafile
      REAL    pi, rad, ra, dec, z, lambda, eta, coseta
      
      pi=asin(1.D0)*2
      rad=pi/180.D0
      
      Open(unit=100, file= TRIM(radecfile), STATUS= 'OLD')
     
      Open(unit = 200, File = TRIM(lambdaetafile), status='New')
      close(200)
      Open(unit = 200, File = TRIM(lambdaetafile), position='APPEND')
      
  100 Read(UNIT=100,fmt=*,END=200) ra, dec, z
      lambda=asin(-cos( rad*(ra-95.0) )*cos(rad*dec))
      
      if((lambda/rad .gt. 89.999) .or. (lambda/rad .lt. -89.999)) then
      eta=0.D0
      goto 110
      endif
       
      eta=asin(sin(rad*dec)/cos(lambda))
c      eta=asin(sin(rad*dec)/cos(lambda))-32.5
      coseta=sin(rad*(ra-95.0))*cos(rad*dec)/cos(lambda)
      if (coseta .ge. 1)  coseta=1.0
      if (coseta .le. -1)  coseta=-1.0
      
      if (eta .ge. 0) then              
      eta=acos(coseta)
      else
      eta=2.D0*pi-acos(coseta)
      endif
      
  110 lambda=lambda/rad
      eta=eta/rad
      eta=eta+57.5
      if (eta .lt. 0) eta = eta + 360.0
      if (eta .gt. 360) eta = eta - 360
      write(unit=200,fmt=*) eta, lambda, z
      goto 100
  200 close(100)
      close(200)

#ifdef debug
      write(*,fmt=*) "finish convert to eta,lambda "
#endif

      END SUBROUTINE ! toLambdaEta               
c -------------------------------------      
	
      SUBROUTINE convertToxyz(SigmaM, SigmaX, 
     & LRGsourcefile_re, LRGsourcefile_xyz, nofzfile)
c convert ra, dec ,z data to x, y, w
c -------- I/O--------------
c SigmaM , SigmaX are cosmology parameters
c LRGsourcefile_re is the file to convert
c nofzfile is the number of galaxies per redshift not dz^3
c LRGsourcefile_xyz is the xyz file
c----------------------------------

      IMPLICIT NONE
c... input
      REAL    SigmaM
      REAL    SigmaX
      CHARACTER*100 :: LRGsourcefile_re
      CHARACTER*100 :: LRGsourcefile_xyz
      CHARACTER*100 :: nofzfile

c      integer, parameter :: r_seperation = 500000 !should be the same as the one in function ztor()
      REAL :: dz
      REAL    ra, dec, z, x, y, w,nofr, ztor
c      real    r(0:r_seperation)
      REAL    nofz(0:60000)
      INTEGER zz,i
c      common r

c      dz = 0.5/r_seperation
c      r(0)=0
      call initialize_r(SigmaM,SigmaX)
      Open(unit=300, file=TRIM(nofzfile), STATUS= 'OLD')
      Read(unit=300,fmt=*) (nofz(i),i=0,60000)
      
c      do i=1,r_seperation
cc dr = c/H*dz (Mpc/h)
c      r(i)=r(i-1)+dz*3000/Sqrt(SigmaM*(1+dz*i)**3+SigmaX)
c      end do
       
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_re), STATUS= 'OLD')
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), status='new')
      close(200)
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), position='append')
  100 READ(UNIT=100,fmt=*,END=200) ra, dec, z
      call radeczToxyz(ra,dec,z,x,y,w)
c -- old one, the factor seems wrong.
c      nofr=50.D0*nofz(Int(z*100000))*sqrt(SigmaM*(1+z)**3+SigmaX)
c     & /ztor(z)**2/3000
c -- n(r)=dN/dz * H/c /sky_area/r^2, sky_area=4*pi*9380/41253
c      nofr=nofz(Int(z*100000)) *100*sqrt(SigmaM*(1+z)**3+SigmaX)/300000
c     & /(4*3.14*9380/41253)/ztor(z)**2
c -- n(r)=dN/dz * dz/dr /sky_area/r^2
      dz=0.001
      nofr=nofz(Int(z*100000)) *dz/(ztor(z+dz)-ztor(z))
     & /(4*3.14*9380/41253)/ztor(z)**2

      write(200,fmt=*)x,y,w,nofr
      goto 100
  200 close(100)
      close(200)
      close(300)
#ifdef debug
      Write(*,fmt=*) "finish converting data to xyz"
#endif
      END SUBROUTINE !convertToxyz() 	

c-----------------------------------------------------------------------
	
      SUBROUTINE convertToxyz_adv(SigmaM, SigmaX,
     & z1, rofz1, z2, rofz2, 
     & LRGsourcefile_re, LRGsourcefile_xyz, nofzfile)
c convert ra, dec ,z data to x, y, w
c -------- I/O--------------
c SigmaM , SigmaX are cosmology parameters
c LRGsourcefile_re is the file to convert
c nofzfile is the number of galaxies per redshift not dz^3
c LRGsourcefile_xyz is the xyz file
c----------------------------------

      IMPLICIT NONE
c... input
      REAL    :: SigmaM
      REAL    :: SigmaX
      REAL    :: z1
      REAL    :: rofz1
      REAL    :: z2
      REAL    :: rofz2
      CHARACTER*100 :: LRGsourcefile_re
      CHARACTER*100 :: LRGsourcefile_xyz
      CHARACTER*100 :: nofzfile

c      integer, parameter :: r_seperation = 500000 !should be the same as the one in function ztor()
      REAL :: dz
      REAL    ra, dec, z, x, y, w,nofr, ztor
      REAL    nofz(0:60000)
      INTEGER zz,i

      If (SigmaM .gt. 0) then
          call initialize_r(SigmaM,SigmaX)
      else
          call initialize_r_spline(z1, rofz1, z2, rofz2)
      end if

      Open(unit=300, file=TRIM(nofzfile), STATUS= 'OLD')
      Read(unit=300,fmt=*) (nofz(i),i=0,60000)
      
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_re), STATUS= 'OLD')
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), status='new')
      close(200)
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), position='append')
  100 READ(UNIT=100,fmt=*,END=200) ra, dec, z
      call radeczToxyz(ra,dec,z,x,y,w)
c -- n(r)=dN/dz * dz/dr /sky_area/r^2
      dz=0.001
      nofr=nofz(Int(z*100000)) *dz/(ztor(z+dz)-ztor(z))
     & /(4*3.14*9380/41253)/ztor(z)**2

      write(200,fmt=*)x,y,w,nofr
      goto 100
  200 close(100)
      close(200)
      close(300)
#ifdef debug
      Write(*,fmt=*) "finish converting data to xyz"
#endif
      END SUBROUTINE !convertToxyz_adv() 	

c-----------------------------------------------------------------------

c ------radecToxyz()---------------
      SUBROUTINE radeczToxyz(ra,dec,z,x,y,w)
      IMPLICIT NONE 
      REAL    ra, dec, z, x, y, w, pi, rr, ztor
      pi=3.1415927
      rr=ztor(z)     
      ra=ra*pi/180.0
      dec=dec*pi/180.0   
      w=rr*sin(dec)  
      x=rr*cos(dec)*cos(ra)
      y=rr*cos(dec)*sin(ra)
      return
      END SUBROUTINE !radeczToxyz()

c ------ initialize_r ------------------
      subroutine initialize_r(sigmam, sigmax)
      implicit none
      real sigmam, sigmax
      integer, parameter :: r_seperation = 500000 !should be the same as the one in function ztor()
      real dz
      integer i
      real    r(0:r_seperation)
      common r
      dz = 0.5/r_seperation
      r(0)=0
      do i=1,r_seperation
c dr = c/H*dz (Mpc/h)
      r(i)=r(i-1)+dz*3000/Sqrt(SigmaM*(1+dz*i)**3+SigmaX)
      end do
      end subroutine initialize_r
c --------------------------------------------------
      subroutine initialize_r_spline(z1, rofz1, z2, rofz2)
      implicit none
      real :: z1, rofz1, z2, rofz2

      integer :: i
      integer,parameter :: num_node = 3
      integer,parameter :: r_seperation = 500000
      real :: z(num_node), rofz(num_node), y2(num_node)
      real :: r(0:r_seperation)
      common r

      z(1) = 0
      rofz(1) = 0
      z(2) = z1
      rofz(2) = rofz1
      z(3) = z2
      rofz(3) = rofz2
      call spline(z, rofz,num_node,1.0e30,1.0e30,y2)
#ifdef debug
      Write(*,*) "y2=", y2(1), y2(2), y2(3)
#endif
      do i = 0, r_seperation
          call splint(z, rofz, y2,num_node,1.0*i/1000000, r(i))
      end do

      end subroutine initialize_r_spline

C ------ ztor()---------------
      REAL    FUNCTION ztor(z)
      IMPLICIT NONE
      integer, parameter :: r_seperation = 500000 ! should be the same as sub convertTOxyz() 
      REAL    z, r(0:r_seperation), dz
      common r

      dz = 0.5/r_seperation
      ztor=r(int(z/dz))
      return
      END FUNCTION ! ztor90      
c---------------------------------      

      SUBROUTINE radecztoxyz_unit(ra,dec,z,x,y,w)  
c It's slide different from radecztoxyz() in BAOSUBROUTINEs.f
C rr = 1 in this SUBROUTINE

      IMPLICIT NONE 
      REAL    ra, dec, z, x, y, w, pi, rr, ztor
      REAL    raa, decc, cosdecc
      pi=3.1415927   
!      rr=1.d0 
     
      raa=ra*pi/180.0    
      decc=dec*pi/180.0  
      cosdecc=cos(decc)
      w=sin(decc)
      x=cosdecc*cos(raa)
      y=cosdecc*sin(raa)
      return
      END SUBROUTINE ! radecztoxyz_unit()
c-------------------------------------------
      subroutine show_corr_result(
     & DDcount_file, DRcount_file, RRcount_file,
     & DD_DR_corr_file, DD_RR_corr_file, DD_DR_RR_corr_file,
     & max_r, bin_size)

c calculate correlation functions with three estimators
c save results to files
      implicit none
c... input
      character*100 :: DDcount_file
      character*100 :: DRcount_file
      character*100 :: RRcount_file
      integer       :: max_r
      integer       :: bin_size
c... output
      character*100 :: DD_DR_corr_file
      character*100 :: DD_RR_corr_file
      character*100 :: DD_DR_RR_corr_file

      integer           :: max_distance
      REAL   , allocatable :: DDcount(:), DRcount(:)
      REAL   , allocatable :: RRcount(:)
      REAL               :: DD_DR_bin, DD_RR_bin, DD_DR_RR_bin
      integer           :: i
      max_distance = max_r
      Allocate(DDcount(max_distance))
      Allocate(DRcount(max_distance))
      Allocate(RRcount(max_distance))

      open(unit = 11, file = TRIM(DDcount_file), status='old')
      open(unit = 12, file = TRIM(DRcount_file), status='old')
      open(unit = 13, file = TRIM(RRcount_file), status='old')

      do i= 1, max_r
        read(unit = 11,fmt=*) DDcount(i)
        read(unit = 12, fmt=*) DRcount(i)
        read(unit = 13, fmt=*) RRcount(i)
      end do !i
      close(11)
      close(12)
      close(13)

      open(unit = 14, file = TRIM(DD_DR_corr_file), status="replace")
      open(unit = 15, file = TRIM(DD_RR_corr_file), status="replace")
      open(unit = 16, file = TRIM(DD_DR_RR_corr_file), status="replace")

      do i= 1, max_r, bin_size
        DD_DR_bin = Sum(DDcount(i:i+bin_size-1))/
     &               Sum(DRcount(i:i+bin_size-1))*2-1
        DD_RR_bin = Sum(DDcount(i:i+bin_size-1))/
     &               Sum(RRcount(i:i+bin_size-1))-1
        DD_DR_RR_bin = (Sum(DDcount(i:i+bin_size-1)) -
     &                  Sum(DRcount(i:i+bin_size-1))) /
     &                    Sum(RRcount(i:i+bin_size-1))+1
        write(14,fmt=*) i+1.0*bin_size/2-1, DD_DR_bin
        write(15,fmt=*) i+1.0*bin_size/2-1, DD_RR_bin
        write(16,fmt=*) i+1.0*bin_size/2-1, DD_DR_RR_bin
      end do

      write(14,fmt=*) ' '
      write(15,fmt=*) ' ' 
      write(16,fmt=*) ' '
      close(14)
      close(15)
      close(16)
      deallocate(DDcount)
      deallocate(DRcount)
      deallocate(RRcount)
      end subroutine ! show_corr_result
c---------------------------------------------
      
      subroutine delete_file(filename)
      implicit none
      character*100 :: filename
      open(unit = 11, file = TRIM(filename),
     & status="replace")
      close(11, status = "delete")
     
      end subroutine ! delete_file 
c-------------------------------------------

      SUBROUTINE XX_corr_np(
     & LRGsourcefile_xyz, nofrfile, max_r)
c calculate DD or RR nofr
c  -------------  blueprint ------------
c LRGsourcefile_xyz is xyz data file
c nofrfile is nofr file(DD or RR)
c --------------------------------------
      IMPLICIT NONE
c... input
      CHARACTER*100     :: LRGsourcefile_xyz
      CHARACTER*100     :: nofrfile
      integer           :: max_r  

      INTEGER           :: data_num = 1
      REAL, allocatable :: x1(:),y1(:),w1(:),np1(:)
      REAL, allocatable :: x2(:),y2(:),w2(:),np2(:) 
      REAL              :: nofr(500), dr, drdr
      integer           :: count1, count2, i, j
      REAL              :: xx, yy, ww, nrr, totalnp1, totalnp2
c count data number
      OPEN(UNIT=100, FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  110 Read(UNIT=100, fmt=*, END=210) xx,yy,ww,nrr
      data_num = data_num+1
      goto 110
  210 close(100)

      allocate(x1(data_num))
      allocate(y1(data_num))
      allocate(w1(data_num))
      allocate(np1(data_num))
      allocate(x2(data_num))
      allocate(y2(data_num))
      allocate(w2(data_num))
      allocate(np2(data_num))

      count1=0
      count2=0
      totalnp1 = 0
      totalnp2 = 0

      i=0
      do i=1,max_r
        nofr(i)=0
      end do
c----read data in-----      
      close(100)
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  100 Read(UNIT=100,fmt=*,END=200) x1(count1+1),y1(count1+1),
     & w1(count1+1), np1(count1+1)
      totalnp1 = totalnp1+np1(count1+1)
      count1=count1+1
      goto 100
  200 close(100)
      close(200)
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  300 Read(UNIT=200,fmt=*,END=400) x2(count2+1), y2(count2+1),
     &  w2(count2+1), np2(count2+1)
      totalnp2 = totalnp2 + np2(count2+1)
      count2=count2+1
      goto 300
  400 close(200) 
      Do i=1,count1
        Do j=i+1,count2
            if((abs(x1(i)-x2(j)) .lt. max_r) .and. 
     &         (abs(y1(i)-y2(j)) .lt. max_r) .and.
     &         (abs(w1(i)-w2(j)) .lt. max_r)) then
               drdr=(x1(i)-x2(j))**2+(y1(i)-y2(j))**2+(w1(i)-w2(j))**2
               if(drdr .lt. max_r**2) then
                  dr=sqrt(drdr)
                  nofr(int(dr+1))=nofr(int(dr+1))+np1(i)*np2(j)
               end if
            end if
        end do
      end do

      Open(unit = 300, File=TRIM(nofrfile), status='New')
      close(300)
      Open(unit = 300, File=TRIM(nofrfile), position='APPEND')

      do i = 1, max_r
        write(unit=300,fmt=*) nofr(i)/totalnp1/totalnp1
      end do

      close(300)
#ifdef debug
      Write(*,fmt=*)"Finish xx_corr of file ", LRGsourcefile_xyz,
     & ", sample num = ", count1
#endif
      deallocate(x1)
      deallocate(y1)
      deallocate(w1)
      deallocate(np1)
      deallocate(x2)
      deallocate(y2)
      deallocate(w2)
      deallocate(np2)

      END SUBROUTINE ! XX_corr_np() 
c-----------------------------------------

      SUBROUTINE DR_corr_np(
     & LRGsourcefile_xyz, RANsourcefile_xyz, nofrfile, max_r)
      IMPLICIT NONE
c... input
      CHARACTER*100 :: LRGsourcefile_xyz
      CHARACTER*100 :: RANsourcefile_xyz
      CHARACTER*100 :: nofrfile 
      INTEGER       :: max_r

      INTEGER           :: data_num = 1 , rand_num = 1
      REAL, allocatable :: x1(:),y1(:),w1(:),np1(:)
      REAL, allocatable :: x2(:),y2(:),w2(:),np2(:) 
      REAL              :: nofr(500), dr, drdr
      integer           :: count1, count2, i, j
      REAL              :: xx, yy, ww, nrr, totalnp1, totalnp2

c count data number
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_xyz),
     & STATUS= 'OLD')
  110 Read(UNIT=100,fmt=*,END=210) xx,yy,ww,nrr
      data_num = data_num+1
      goto 110
  210 close(100)
c count rand number
      OPEN(UNIT=100,FILE=TRIM(RANsourcefile_xyz),
     & STATUS= 'OLD')
  310 Read(UNIT=100,fmt=*,END=410) xx,yy,ww,nrr
      rand_num = rand_num+1
      goto 310
  410 close(100)

      allocate(x1(data_num))
      allocate(y1(data_num))
      allocate(w1(data_num))
      allocate(np1(data_num))
      allocate(x2(rand_num))
      allocate(y2(rand_num))
      allocate(w2(rand_num))
      allocate(np2(rand_num))
      
      count1=0
      count2=0
      totalnp1 = 0
      totalnp2 = 0
      do i=1,max_r
        nofr(i)=0
      end do

c----read data in-----      

      close(100)
      OPEN(UNIT=100,FILE=Trim(LRGsourcefile_xyz), STATUS= 'OLD')
  100 Read(UNIT=100,fmt=*,END=200) x1(count1+1),y1(count1+1),
     &  w1(count1+1), np1(count1+1)
      totalnp1 = totalnp1 + np1(count1+1)
      count1=count1+1
      goto 100
  200 close(100)
      close(200)
      OPEN(UNIT=200,FILE=Trim(RANsourcefile_xyz), STATUS= 'OLD')
  300 Read(UNIT=200,fmt=*,END=400) x2(count2+1), y2(count2+1),
     &  w2(count2+1), np2(count2+1)
      totalnp2=totalnp2+np2(count2+1)
      count2 = count2 + 1
      goto 300
  400 close(200) 
c----------------------

      Do 20 i=1,count1
         Do 30 j=1,count2
            if((abs(x1(i)-x2(j)) .lt. max_r) .and.
     &         (abs(y1(i)-y2(j)) .lt. max_r) .and.
     &         (abs(w1(i)-w2(j)) .lt. max_r)) then
               drdr=(x1(i)-x2(j))**2+(y1(i)-y2(j))**2+(w1(i)-w2(j))**2
               if(drdr .lt. max_r**2) then
                dr=sqrt(drdr)
                nofr(int(dr+1))=nofr(int(dr+1))+1.0*np1(i)*np2(j)
               end if
            end if
   30    continue
   20 continue  
  
      Open(unit = 300, File=TRIM(nofrfile), status='New')
      close(300)
      Open(unit = 300, File=TRIM(nofrfile), position='APPEND')
      do i = 1, max_r
      write(unit=300,fmt=*) nofr(i)/totalnp2/totalnp1
      end do

      close(300)
#ifdef debug
      Write(*,fmt=*)"Finish DR_corr, data sample num = ", count2, 
     & " ,and random sample num = ", count1
#endif
      deallocate(x1)
      deallocate(y1)
      deallocate(w1)
      deallocate(np1)
      deallocate(x2)
      deallocate(y2)
      deallocate(w2)
      deallocate(np2)

      END SUBROUTINE ! DR_corr_np()
 
c --------------------------------------------------------------
      subroutine convertxyz2radecz(xyzfile,radeczfile,boxsize)
      implicit none
      character*100 :: xyzfile
      character*100 :: radeczfile
      integer       :: boxsize
      real :: xx, yy, ww, nnp
      real,allocatable :: x(:), y(:), w(:), np(:), ra(:), dec(:), z(:)
      integer :: count, i, j, count_ra
      real,parameter :: omega_m=0.3, omega_v=0.7
      real,parameter :: c=300000,pi=3.1415927
      real :: dz, zz, dr,r, spline_z(2000), spline_r(0:2000), y2(2000)

      

! initial z to r
      spline_r(0) = 0
      dz = 0.00005
      do i = 1, 2000
        zz = dz*i*10
        dr = 0
        do j = 1, 10
          zz = zz + dz*j
          dr = dr + c*dz/100/sqrt(omega_m*(1.+zz)**3+omega_v) 
        end do
        spline_z(i) = dz*(-0.5+i)*10
        spline_r(i) = spline_r(i-1)+dr
      end do

      call spline(spline_r(1:2000), spline_z,2000,1.0e30,1.0e30,y2)
! count data
      open(unit=10, file=Trim(xyzfile), status='old')
      count = 0
 100  read(10,*, end=200) xx, yy ,ww ,nnp
      count = count +1
      goto 100
 200  close(10)

      allocate(x(count))
      allocate(y(count))
      allocate(w(count))
      allocate(np(count))
      allocate(ra(count))
      allocate(dec(count))
      allocate(z(count))


      open(unit=10, file=Trim(xyzfile), status='old')
      do i = 1, count
        read(10,*) x(i), y(i), w(i), np(i)
        x(i) = x(i) - boxsize/2
        y(i) = y(i) - boxsize/2
        w(i) = w(i) - boxsize/2
      end do
      close(10)

      count_ra = 0
      do i=1, count
         if(x(i) .gt. 0.001) then ! only need ra=100~270, z not very close to 0
            r = sqrt(x(i)**2+y(i)**2+w(i)**2)
            if(r .lt. boxsize/2) then
              call splint(spline_r(1:2000),spline_z,y2,2000,r,zz)
              if( zz .lt. 0.5) then !choose z < 0.5
                dec(count_ra+1) = ASIN(w(i)/r)/pi*180
                ra(count_ra+1)  = ATAN(y(i)/x(i))/pi*180 + 90 + 100 !now ra=100-280
                z(count_ra+1) = zz
                count_ra = count_ra +1
              end if
            end if
         end if
      end do

      open(unit =11, file=Trim(radeczfile))
      do i = 1, count_ra
        write(11,*) ra(i), dec(i), z(i), np(i)
      end do
      close(11)
      write(*,*) "finish obtain ", count_ra," ra dec z np data from ",
     &  count, " x y z np data"
      deallocate(x)
      deallocate(y)
      deallocate(w)
      deallocate(np)
      deallocate(ra)
      deallocate(dec)
      deallocate(z)

      end subroutine

c-------------------------------------------------------- 
c -------------------------------------      
	
      SUBROUTINE convertToxyz_nonr(SigmaM, SigmaX, 
     & LRGsourcefile_re, LRGsourcefile_xyz)
c convert ra, dec ,z data to x, y, w
c -------- I/O--------------
c SigmaM , SigmaX are cosmology parameters
c LRGsourcefile_re is the file to convert
c nofzfile is the number of galaxies per redshift not dz^3
c LRGsourcefile_xyz is the xyz file
c----------------------------------

      IMPLICIT NONE
c... input
      REAL    SigmaM
      REAL    SigmaX
      CHARACTER*100 :: LRGsourcefile_re
      CHARACTER*100 :: LRGsourcefile_xyz

      REAL :: dz
      REAL    ra, dec, z, x, y, w, ztor
      INTEGER zz,i

      call initialize_r(SigmaM,SigmaX)
      
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_re), STATUS= 'OLD')
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz))
  100 READ(UNIT=100,fmt=*,END=200) ra, dec, z
      call radeczToxyz(ra,dec,z,x,y,w)

      write(200,fmt=*)x,y,w
      goto 100
  200 close(100)
      close(200)

#ifdef debug
      Write(*,fmt=*) "finish converting data to xyz"
#endif
      END SUBROUTINE !convertToxyz_nonr() 	

c-----------------------------------------------------------------------

      SUBROUTINE XX_corr_nonp(
     & LRGsourcefile_xyz, nofrfile, max_r)
c calculate DD or RR nofr
c  -------------  blueprint ------------
c LRGsourcefile_xyz is xyz data file
c nofrfile is nofr file(DD or RR)
c --------------------------------------
      IMPLICIT NONE
c... input
      CHARACTER*100     :: LRGsourcefile_xyz
      CHARACTER*100     :: nofrfile
      integer           :: max_r  

      INTEGER           :: data_num = 1
      REAL, allocatable :: x1(:),y1(:),w1(:)
      REAL, allocatable :: x2(:),y2(:),w2(:)
      REAL              :: nofr(500), dr, drdr
      integer           :: count1, count2, i, j
      REAL              :: xx, yy, ww, totalnp1, totalnp2
c count data number
      OPEN(UNIT=100, FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  110 Read(UNIT=100, fmt=*, END=210) xx,yy,ww
      data_num = data_num+1
      goto 110
  210 close(100)

      allocate(x1(data_num))
      allocate(y1(data_num))
      allocate(w1(data_num))
      allocate(x2(data_num))
      allocate(y2(data_num))
      allocate(w2(data_num))

      count1=0
      count2=0
      totalnp1 = 0
      totalnp2 = 0

      i=0
      do i=1,max_r
        nofr(i)=0
      end do
c----read data in-----      
      close(100)
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  100 Read(UNIT=100,fmt=*,END=200) x1(count1+1),y1(count1+1),
     & w1(count1+1)
      totalnp1 = totalnp1+1
      count1=count1+1
      goto 100
  200 close(100)
      close(200)
      OPEN(UNIT=200,FILE=TRIM(LRGsourcefile_xyz), STATUS= 'OLD')
  300 Read(UNIT=200,fmt=*,END=400) x2(count2+1), y2(count2+1),
     &  w2(count2+1)
      totalnp2 = totalnp2 + 1
      count2=count2+1
      goto 300
  400 close(200) 
      Do i=1,count1
        Do j=i+1,count2
            if((abs(x1(i)-x2(j)) .lt. max_r) .and. 
     &         (abs(y1(i)-y2(j)) .lt. max_r) .and.
     &         (abs(w1(i)-w2(j)) .lt. max_r)) then
               drdr=(x1(i)-x2(j))**2+(y1(i)-y2(j))**2+(w1(i)-w2(j))**2
               if(drdr .lt. max_r**2) then
                  dr=sqrt(drdr)
                  nofr(int(dr+1))=nofr(int(dr+1))+1
               end if
            end if
        end do
      end do

      Open(unit = 300, File=TRIM(nofrfile))

      do i = 1, max_r
        write(unit=300,fmt=*) nofr(i)/totalnp1/totalnp1
      end do

      close(300)
#ifdef debug
      Write(*,fmt=*)"Finish xx_corr of file ", LRGsourcefile_xyz,
     & ", sample num = ", count1
#endif
      deallocate(x1)
      deallocate(y1)
      deallocate(w1)
      deallocate(x2)
      deallocate(y2)
      deallocate(w2)

      END SUBROUTINE ! XX_corr_nonp() 
c-----------------------------------------

      SUBROUTINE DR_corr_nonp(
     & LRGsourcefile_xyz, RANsourcefile_xyz, nofrfile, max_r)
      IMPLICIT NONE
c... input
      CHARACTER*100 :: LRGsourcefile_xyz
      CHARACTER*100 :: RANsourcefile_xyz
      CHARACTER*100 :: nofrfile 
      INTEGER       :: max_r

      INTEGER           :: data_num = 1 , rand_num = 1
      REAL, allocatable :: x1(:),y1(:),w1(:)
      REAL, allocatable :: x2(:),y2(:),w2(:) 
      REAL              :: nofr(500), dr, drdr
      integer           :: count1, count2, i, j
      REAL              :: xx, yy, ww, totalnp1, totalnp2

c count data number
      OPEN(UNIT=100,FILE=TRIM(LRGsourcefile_xyz),
     & STATUS= 'OLD')
  110 Read(UNIT=100,fmt=*,END=210) xx,yy,ww
      data_num = data_num+1
      goto 110
  210 close(100)
c count rand number
      OPEN(UNIT=100,FILE=TRIM(RANsourcefile_xyz),
     & STATUS= 'OLD')
  310 Read(UNIT=100,fmt=*,END=410) xx,yy,ww
      rand_num = rand_num+1
      goto 310
  410 close(100)

      allocate(x1(data_num))
      allocate(y1(data_num))
      allocate(w1(data_num))
      allocate(x2(rand_num))
      allocate(y2(rand_num))
      allocate(w2(rand_num))
      
      count1=0
      count2=0
      totalnp1 = 0
      totalnp2 = 0
      do i=1,max_r
        nofr(i)=0
      end do

c----read data in-----      

      close(100)
      OPEN(UNIT=100,FILE=Trim(LRGsourcefile_xyz), STATUS= 'OLD')
  100 Read(UNIT=100,fmt=*,END=200) x1(count1+1),y1(count1+1),
     &  w1(count1+1)
      totalnp1 = totalnp1 + 1
      count1=count1+1
      goto 100
  200 close(100)
      close(200)
      OPEN(UNIT=200,FILE=Trim(RANsourcefile_xyz), STATUS= 'OLD')
  300 Read(UNIT=200,fmt=*,END=400) x2(count2+1), y2(count2+1),
     &  w2(count2+1)
      totalnp2=totalnp2+1
      count2 = count2 + 1
      goto 300
  400 close(200) 
c----------------------

      Do 20 i=1,count1
         Do 30 j=1,count2
            if((abs(x1(i)-x2(j)) .lt. max_r) .and.
     &         (abs(y1(i)-y2(j)) .lt. max_r) .and.
     &         (abs(w1(i)-w2(j)) .lt. max_r)) then
               drdr=(x1(i)-x2(j))**2+(y1(i)-y2(j))**2+(w1(i)-w2(j))**2
               if(drdr .lt. max_r**2) then
                dr=sqrt(drdr)
                nofr(int(dr+1))=nofr(int(dr+1))+1.0
               end if
            end if
   30    continue
   20 continue  

      Open(unit = 300, File=TRIM(nofrfile))
      do i = 1, max_r
      write(unit=300,fmt=*) nofr(i)/totalnp2/totalnp1
      end do

      close(300)
#ifdef debug
      Write(*,fmt=*)"Finish DR_corr, data sample num = ", count2, 
     & " ,and random sample num = ", count1
#endif
      deallocate(x1)
      deallocate(y1)
      deallocate(w1)
      deallocate(x2)
      deallocate(y2)
      deallocate(w2)


      END SUBROUTINE ! DR_corr_nonp()
 
!-------------------------------------------------------
      subroutine identify_den(datafile, data_num,iseed) 
      implicit none
      character*100 :: datafile
      integer       :: data_num
      integer       :: iseed

      integer :: count, i, test_count
      real    :: ra, dec, z, defactor, ran3
      real,allocatable :: datara(:),datadec(:),dataz(:)

      count = 0
      open(10, file=Trim(datafile), status='old')
 100  read(10,*,end=200) ra, dec, z
        count =  count + 1
      goto 100
 200  close(10)

      allocate(datara(count))
      allocate(datadec(count))
      allocate(dataz(count))

      open(10, file=Trim(datafile), status='old')
      do i = 1, count
         read(10,*) datara(i),datadec(i),dataz(i)
      end do
      close(10)
!overwrite datafile
      test_count = 0
      open(10,file=Trim(datafile))
      defactor = 1.0*data_num/count
      do i =1, count
         if(ran3(iseed) .lt. defactor) then
            write(10,*) datara(i),datadec(i),dataz(i)
            test_count = test_count + 1
         end if
      end do
      close(10)

      deallocate(datara)
      deallocate(datadec)
      deallocate(dataz)
      end subroutine
