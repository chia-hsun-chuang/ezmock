c This file includes the subroutines dealing with vagc data and the mask especially. 
c -----------------------------------------	

c      SUBROUTINE mask_filter_data(
c     & splyfile, LRGsourcefile, LRGsourcefile_in, LRGsourcefile_out)
c      SUBROUTINE mask_filter_data_with_mg(
c     & splyfile, LRGsourcefile, LRGsourcefile_in, LRGsourcefile_out)
c      SUBROUTINE mask_filter_vagc2Drandom(total_ran,
c     & splyfile, vagc2Drandom_file,
c     & vagc2Drandom_file_in, vagc2Drandom_file_ou)
c      logical FUNCTION gptin(rp,cm,np,rpi)
c      SUBROUTINE generate_3D_random(
c     & rannum, zmin, zmax, use_fgot,
c     & vagcRand, nofzfile, randomgalaxy)
c      subroutine vagc_generate_2Drandfile(
c     & sply2file,twoDrandfile,rand_num, iseed)

c -----------------------------------------	    
      SUBROUTINE mask_filter_data(
     & splyfile, LRGsourcefile, LRGsourcefile_in, LRGsourcefile_out)
c--- This SUBROUTINE is for using mask/window to filter the data into
c--- two files(inside and outside the mask/window)
      IMPLICIT NONE
C...input
      Character*100 :: splyfile
      Character*100 :: LRGsourcefile
      Character*100 :: LRGsourcefile_in
      Character*100 :: LRGsourcefile_out

      REAL    datara(1000000), datadec(1000000), dataz(1000000)
      integer  totalpoly, polyNum, capsNum(700000), i,j,k, countdata
      integer  inpoint, outpoint
      REAL    weight, str, rpx(30,700000), rpy(30,700000)
      REAL    rpz(30,700000), cm(30,700000), rpi(3), rp(3,30)
      Character*8 str1
      logical gptin
      
C...read data in an array
      Open(unit=100, file= TRIM(LRGsourcefile), STATUS= 'OLD')
      countdata=0
  100 countdata=countdata+1
      Read(unit=100,fmt=*, end=200) datara(countdata), 
     & datadec(countdata), dataz(countdata)
       
      goto 100
  200 countdata=countdata-1
      if(countdata .gt. 1000000) then
        write(*,fmt=*) "# of data is ",countdata,
     &   " which is larger than 1000000."
        stop
      END if
      close(100)

C... read mask in arrays
      Open(unit=100, file= TRIM(splyfile), STATUS= 'OLD')
      Read(unit=100,fmt=*) totalpoly, str1
      Write(*,fmt=*) totalpoly ,str1
      if( totalpoly .gt. 700000) then
        write(*,fmt=*) "# of polygons is larger than 700000"
        stop
      END if
      do j=1,totalpoly
        	 
        Read(unit=100,fmt=*) polyNum,capsNum(j),weight,str

        if( capsNum(j) .gt. 30) then
        	write(*,fmt=*) "# of caps is larger than 50"
        	stop
        END if
        do i=1, capsNum(j)
          Read(unit=100,fmt=*) rpx(i,j), rpy(i,j) , rpz(i,j), cm(i,j)

        END do
      END do
      close(100)
C... open a new file to save the data been filtered.
      Open(unit=200, file= TRIM(LRGsourcefile_in), STATUS= 'new')
      close(200)
      Open(unit=200, file= TRIM(LRGsourcefile_in), position='append')
      Open(unit=300, file= TRIM(LRGsourcefile_out), STATUS= 'new')
      close(300)
      Open(unit=300, file= TRIM(LRGsourcefile_out), position='append')
C... loop for check data
      inpoint=0
      outpoint=0

      do i=1, countdata    
  
        do j=1, totalpoly
 
            call radecztoxyz_unit(datara(i),datadec(i),dataz(i),
     &       rpi(1),rpi(2),rpi(3))
      
            do k=1,capsnum(j)
              rp(1,k)=rpx(k,j)
              rp(2,k)=rpy(k,j)
              rp(3,k)=rpz(k,j)
            END do  

          if(gptin(rp(:,1:capsnum(j)),cm(1:capsnum(j),j),capsnum(j),
     &       rpi(:))) then

c... find the point inside the mask
            write(200,fmt=*) datara(i),datadec(i),dataz(i)            
            inpoint=inpoint+1
            goto 300
          END if
        END do
c... has checked all polygons so write into file        	
        write(300,fmt=*) datara(i),datadec(i),dataz(i)
        outpoint=outpoint+1	
  300   k=k  !meaningless line     	
      END do
      close(200)  
      close(300)

      Write(*,fmt=*) inpoint," points are inside the mask. "
      Write(*,fmt=*) "There are ", outpoint,
     & " points outside the mask." 

      END SUBROUTINE !mask_filter_data()
c------------------------------------------------
      SUBROUTINE mask_filter_data_with_mg(
     & splyfile, LRGsourcefile, LRGsourcefile_in, LRGsourcefile_out)
c--- This SUBROUTINE is for using mask/window to filter the data into
c--- two files(inside and outside the mask/window)
      IMPLICIT NONE
C...input
      Character*100 :: splyfile
      Character*100 :: LRGsourcefile
      Character*100 :: LRGsourcefile_in
      Character*100 :: LRGsourcefile_out

      REAL    datara(1000000), datadec(1000000)
      REAL    dataz(1000000), datamg(1000000)
      integer  totalpoly, polyNum, capsNum(700000), i,j,k, countdata
      integer  inpoint, outpoint
      REAL    weight, str, rpx(30,700000), rpy(30,700000)
      REAL    rpz(30,700000), cm(30,700000), rpi(3), rp(3,30)
      Character*8 str1
      logical gptin
      
C...read data in an array
      Open(unit=100, file= TRIM(LRGsourcefile), STATUS= 'OLD')
      countdata=0
  100 countdata=countdata+1
      Read(unit=100,fmt=*, end=200) datara(countdata), 
     & datadec(countdata), dataz(countdata), datamg(countdata)
       
      goto 100
  200 countdata=countdata-1
      if(countdata .gt. 1000000) then
        write(*,fmt=*) "# of data is ",countdata,
     &   " which is larger than 1000000."
        stop
      END if
      close(100)

C... read mask in arrays
      Open(unit=100, file= TRIM(splyfile), STATUS= 'OLD')
      Read(unit=100,fmt=*) totalpoly, str1
      Write(*,fmt=*) totalpoly ,str1
      if( totalpoly .gt. 700000) then
        write(*,fmt=*) "# of polygons is larger than 700000"
        stop
      END if
      do j=1,totalpoly
        	 
        Read(unit=100,fmt=*) polyNum,capsNum(j),weight,str

        if( capsNum(j) .gt. 30) then
        	write(*,fmt=*) "# of caps is larger than 50"
        	stop
        END if
        do i=1, capsNum(j)
          Read(unit=100,fmt=*) rpx(i,j), rpy(i,j) , rpz(i,j), cm(i,j)

        END do
      END do
      close(100)
C... open a new file to save the data been filtered.
      Open(unit=200, file= TRIM(LRGsourcefile_in), STATUS= 'new')
      close(200)
      Open(unit=200, file= TRIM(LRGsourcefile_in), position='append')
      Open(unit=300, file= TRIM(LRGsourcefile_out), STATUS= 'new')
      close(300)
      Open(unit=300, file= TRIM(LRGsourcefile_out), position='append')
C... loop for check data
      inpoint=0
      outpoint=0

      do i=1, countdata    
  
        do j=1, totalpoly
 
            call radecztoxyz_unit(datara(i),datadec(i),dataz(i),
     &       rpi(1),rpi(2),rpi(3))
      
            do k=1,capsnum(j)
              rp(1,k)=rpx(k,j)
              rp(2,k)=rpy(k,j)
              rp(3,k)=rpz(k,j)
            END do  

          if(gptin(rp(:,1:capsnum(j)),cm(1:capsnum(j),j),capsnum(j),
     &       rpi(:))) then

c... find the point inside the mask
            write(200,fmt=*) datara(i),datadec(i),dataz(i),datamg(i)            
            inpoint=inpoint+1
            goto 300
          END if
        END do
c... has checked all polygons so write into file        	
        write(300,fmt=*) datara(i),datadec(i),dataz(i),datamg(i)
        outpoint=outpoint+1
  300   k=k  !meaningless line     	
      END do
      close(200)  
      close(300)

      Write(*,fmt=*) inpoint," points are inside the mask. "
      Write(*,fmt=*) "There are ", outpoint,
     & " points outside the mask." 

      END SUBROUTINE !mask_filter_data_with_mg()
c------------------------------------------------
      SUBROUTINE mask_filter_vagc2Drandom(total_ran,
     & splyfile, vagc2Drandom_file,
     & vagc2Drandom_file_in, vagc2Drandom_file_ou)
c--- This SUBROUTINE is for using mask/window to filter the vagc-2D-random-data
c--- into two files(inside and outside the mask/window)
      IMPLICIT NONE

C...input
      Integer      :: total_ran
      Character*100 :: splyfile
      Character*100 :: vagc2Drandom_file
      Character*100 :: vagc2Drandom_file_in
      Character*100 :: vagc2Drandom_file_ou

      Integer ran_data_index
      REAL   , Allocatable :: datara(:), datadec(:)
      REAL   , Allocatable :: var3(:), var4(:)
      REAL   , Allocatable :: FGOT(:), var6(:)
      Integer  totalpoly, polyNum, capsNum(700000), i,j,k
      Integer  inpoint, outpoint
      REAL    weight, str, rpx(30,700000), rpy(30,700000)
      REAL    rpz(30,700000), cm(30,700000), rpi(3), rp(3,30)
      Character*8 str1
      logical gptin
      
C... allocate arrays
      Allocate(datara(total_ran))
      Allocate(datadec(total_ran))
      Allocate(var3(total_ran))
      Allocate(var4(total_ran))
      Allocate(FGOT(total_ran))
      Allocate(var6(total_ran))
      
C...read data in an array
      Open(unit=100, file= TRIM(vagc2Drandom_file), STATUS= 'OLD')
      ran_data_index=0
      DO ran_data_index=1,total_ran 
      Read(unit=100,fmt=*, end=200) datara(ran_data_index), 
     & datadec(ran_data_index),var3(ran_data_index),
     & var4(ran_data_index),FGOT(ran_data_index), 
     & var6(ran_data_index)
      END DO 
      
  200 IF (ran_data_index .lt. total_ran) then
        write(*,fmt=*)
     &  "random data in file is less than what you asked ", total_ran
      END IF
      close(100)

C... read mask in arrays
      Open(unit=100, file=TRIM(splyfile), STATUS= 'OLD')
      Read(unit=100,fmt=*) totalpoly, str1
      Write(*,fmt=*) totalpoly ,str1
      if( totalpoly .gt. 700000) then
        write(*,fmt=*) "# of polygons is larger than 700000"
        stop
      END if
      do j=1,totalpoly
        	 
        Read(unit=100,fmt=*) polyNum,capsNum(j),weight,str

        if( capsNum(j) .gt. 30) then
        	write(*,fmt=*) "# of caps is larger than 30"
        	stop
        END if
        do i=1, capsNum(j)
          Read(unit=100,fmt=*) rpx(i,j), rpy(i,j) , rpz(i,j), cm(i,j)

        END do
      END do
      close(100)

C... open a new file to save the data been filtered.
      Open(unit=200, file=TRIM(vagc2Drandom_file_in), STATUS= 'new')
      close(200)
      Open(unit=200, file=TRIM(vagc2Drandom_file_in), position='append')
      Open(unit=300, file=TRIM(vagc2Drandom_file_ou), STATUS= 'new')
      close(300)
      Open(unit=300, file=TRIM(vagc2Drandom_file_ou), position='append')

C... loop for check data
      inpoint=0
      outpoint=0

      do i=1, total_ran    
  
        do j=1, totalpoly
 
            call radecztoxyz_unit(datara(i),datadec(i), 1.0,
     &       rpi(1),rpi(2),rpi(3))
      
            do k=1,capsnum(j)
              rp(1,k)=rpx(k,j)
              rp(2,k)=rpy(k,j)
              rp(3,k)=rpz(k,j)
            END do  !k

          if(gptin(rp(:,1:capsnum(j)),cm(1:capsnum(j),j),capsnum(j),
     &       rpi(:))) then

c... find the point inside the mask
 
            write(200,fmt=*) datara(i), datadec(i), var3(i),
     &                   var4(i), FGOT(i), var6(i)            
            inpoint=inpoint+1
            goto 300
          END if
        END do !j

c... has checked all polygons so write into file        	
        write(300,fmt=*) datara(i),datadec(i), var3(i),
     &               var4(i), FGOT(i), var6(i)
        outpoint=outpoint+1	
  300   k=k  !meaningless line     	
      END do !i
      close(200)  
      close(300)

      Write(*,fmt=*) inpoint," points are inside the mask. "
      Write(*,fmt=*) "There are ", outpoint," points outside the mask."
 
      DEALLOCATE(datara)
      DEALLOCATE(datadec)
      DEALLOCATE(var3)
      DEALLOCATE(var4)
      DEALLOCATE(FGOT)
      DEALLOCATE(var6)

      END SUBROUTINE !mask_filter_vagc2Drandom()
      
c-----------------------------------------------------------------------      
      logical FUNCTION gptin(rp,cm,np,rpi)
      integer np
      REAL    rp(3,np),cm(np),rpi(3)
c
c        intrinsics
      intrinsic abs
c        externals
      integer gzeroar
c        local (automatic) variables
      integer j
      REAL    cmij,cmj
c *
c * Determine whether unit direction rpi lies within region bounded by
c *    1 - r.rp(j) <= cm(j)  (if cm(j).ge.0)
c *    1 - r.rp(j) > -cm(j)  (if cm(j).lt.0)
c * for j=1,np where rp(j) are unit directions.
c *
c  Input: rp(3,j),j=1,np
c         cm(j),j=1,np
c         np
c         rpi(3)
c Output: gptin = .true. if point lies within region
c                 .false. if outside.
c

      gptin=.false.
c        check for point outside because one circle is null 
!! I skip this because I cannot find gzeroar
cc      if (gzeroar(cm,np).eq.0) goto 410
c        check each boundary
      do 140 j=1,np
c        null boundary means no constraint
        if (cm(j).ge.2.d0) goto 140
        cmj=abs(cm(j))
c        1-cos of angle between point and rp(j) direction
        cmij=((rpi(1)-rp(1,j))**2+(rpi(2)-rp(2,j))**2
     *    +(rpi(3)-rp(3,j))**2)/2.d0
c        check if point is outside rp(j) boundary
        if (cm(j).ge.0.d0) then
          if (cmij.gt.cmj) goto 410
        elseif (cm(j).lt.0.d0) then
          if (cmij.le.cmj) goto 410
        endif
  140 continue
c        point survived all assails
      gptin=.true.
c        done
  410 continue
      return
      END FUNCTION !gptin

c-----------------------------------------------------------------------


      SUBROUTINE generate_3D_random(
     & rannum, zmin, zmax, use_fgot,
     & vagcRand, nofzfile, randomgalaxy)
c--- This SUBROUTINE generate 3D random data from vagc2D random data file

c ------- blueprint -----
c read nofzfile
c decide maximum of nofz
c loop to read random data
c    pick random set 0<x<max, zmin<y<zmax
c    if x < nofz(y), write to file 
c    else  pick again
c END loop
c ------------------------
      IMPLICIT NONE
c... input
      REAL          :: zmin
      REAL          :: zmax
      Integer       :: rannum 
      logical       :: use_fgot
      Character*100 :: vagcRand
      Character*100 :: nofzfile
      Character*100 :: randomgalaxy

      Integer       :: i, iseed,ixz
      REAL          :: ran3
      REAL          :: nofz(0:60000),nofzmax, xz, ynofz, ranFgot
      REAL          :: radata, decdata, var3, var4, fgot, var6

c... set iseed  
      iseed = -1001

c read nofzfile into an array       
      Open(unit=100, file=TRIM(nofzfile), STATUS= 'OLD')
      Read(unit=100,fmt=*) (nofz(i),i=0,60000)
c decide maximum of nofz      
      nofzmax=0
      do i=0,59999
      if(nofz(i) .gt. nofzmax) nofzmax=nofz(i)
      END do
cc test
      write(*,fmt=*) "max of nofz is ",nofzmax 
      
c loop to read data
      Open(unit=200, file=TRIM(vagcRand), STATUS= 'OLD')  
      Open(unit =300, File=TRIM(randomgalaxy), status='New')
      close(300)
      Open(unit =300, File=TRIM(randomgalaxy), position='APPEND')
      do i=1,rannum
   50   READ(UNIT=200,fmt=*) radata, decdata, var3, var4, fgot, var6
        ranFgot=ran3(iseed)
        if (.NOT. use_fgot) fgot = 1.0 
        if (ranFgot .gt. fgot) goto 50 
  100   ynofz=ran3(iseed)*nofzmax
        xz=ran3(iseed)*(zmax-zmin)+zmin
        ixz=Int(100000*xz)
          if (ynofz .lt. nofz(ixz)) then
      	    write(300,fmt=*) radata, decdata, xz
          else
c      	    Write(*,*) "pick another"
      	    goto 100
          END if
      END do
      close(100)
      close(200)
      close(300)

      END SUBROUTINE !generat_3D_random()

c ----------------------------------------------

      SUBROUTINE generate_3D_random_skip(
     & skip_num,
     & rannum, zmin, zmax, use_fgot,
     & vagcRand, nofzfile, randomgalaxy)
c--- This SUBROUTINE generate 3D random data from vagc2D random data file

c ------- blueprint -----
c read nofzfile
c decide maximum of nofz
c loop to read random data
c    pick random set 0<x<max, zmin<y<zmax
c    if x < nofz(y), write to file 
c    else  pick again
c END loop
c ------------------------
      IMPLICIT NONE
c... input
      integer      :: skip_num
      REAL         :: zmin
      REAL         :: zmax
      Integer      :: rannum 
      logical      :: use_fgot
      Character*100 :: vagcRand
      Character*100 :: nofzfile
      Character*100 :: randomgalaxy

      Integer      :: i, iseed,ixz
      REAL         :: ran3
      REAL         :: nofz(0:60000),nofzmax, xz, ynofz, ranFgot
      REAL         :: radata, decdata, var3, var4, fgot, var6

c... set iseed  
      iseed = -1001

c read nofzfile into an array       
      Open(unit=100, file=TRIM(nofzfile), STATUS= 'OLD')
      Read(unit=100,fmt=*) (nofz(i),i=0,60000)
c decide maximum of nofz      
      nofzmax=0
      do i=0,59999
      if(nofz(i) .gt. nofzmax) nofzmax=nofz(i)
      END do
cc test
      write(*,fmt=*) "max of nofz is ",nofzmax 
      
c loop to read data
      Open(unit=200, file=TRIM(vagcRand), STATUS= 'OLD')  
      Open(unit =300, File =TRIM(randomgalaxy), status='New')
      close(300)
      Open(unit =300, File =TRIM(randomgalaxy), position='APPEND')

      do i= 1, skip_num
        READ(UNIT=200,fmt=*,end=500) 
     &      radata, decdata, var3, var4, fgot, var6
      end do

      do i= 1, rannum
   50   READ(UNIT=200,fmt=*,end=500) 
     &      radata, decdata, var3, var4, fgot, var6
        ranFgot=ran3(iseed)
        if (.NOT. use_fgot) fgot = 1.0 
        if (ranFgot .gt. fgot) goto 50 
  100   ynofz=ran3(iseed)*nofzmax
        xz=ran3(iseed)*(zmax-zmin)+zmin
        ixz=Int(100000*xz)
          if (ynofz .lt. nofz(ixz)) then
      	    write(300,fmt=*) radata, decdata, xz
          else
c      	    Write(*,*) "pick another"
      	    goto 100
          END if
      END do
      close(100)
      close(200)
      close(300)
      
      If (.false.) then
 500	 write(*,*) "vagc rand file is not enough for usage"
      end if
      END SUBROUTINE !generat_3D_random_skip()

c --------------------------------------------------

      subroutine vagc_generate_2Drandfile(
     & sply2file,twoDrandfile,rand_num, iseed)
c to generate 2D random catalog with vagc mask
      implicit none
c input
      character*100       :: sply2file
      character*100       :: twoDrandfile
      integer, intent(in) :: rand_num
      integer, intent(in) :: iseed

c
      integer, parameter  :: maxcaps = 50
      Integer             :: totalpoly, polyNum, i,j,k
      Integer,allocatable :: capsNum(:)
      real, allocatable   ::  fgot(:), mmax(:)
      Integer             :: inpoint = 0
      REAL                :: rpi(3), rp(3,maxcaps)
      REAL, allocatable   :: rpx(:,:), rpy(:,:), rpz(:,:), cm(:,:)
      Character*8         :: str1
      logical             :: gptin
      real                :: ran3

      real                      :: x,y,w,ra,dec,z,rr,fgot_test
      real, dimension(rand_num) :: randra, randdec, var3, var4, 
     &                             randfgot, randmmax

c read mask in arrays
      Open(unit=100, file=TRIM(sply2file), STATUS= 'OLD')
      Read(unit=100,fmt=*) totalpoly, str1
      Write(*,fmt=*) totalpoly ,str1

      allocate(capsNum(totalpoly))
      allocate(fgot(totalpoly))
      allocate(mmax(totalpoly))
      allocate(rpx(maxcaps,totalpoly))
      allocate(rpy(maxcaps,totalpoly))
      allocate(rpz(maxcaps,totalpoly))
      allocate(cm(maxcaps,totalpoly))

      do j=1,totalpoly
        	 
        Read(unit=100,fmt=*) polyNum,capsNum(j),fgot(j),mmax(j)

        if( capsNum(j) .gt. 50) then
        	write(*,fmt=*) "# of caps is larger than 50"
        	stop
        END if
        do i=1, capsNum(j)
          Read(unit=100,fmt=*) rpx(i,j), rpy(i,j) , rpz(i,j), cm(i,j)

        END do
      END do
      close(100)

c generate random data

      DO while (inpoint .lt. rand_num)

        x = ran3(iseed)*2 - 1
        y = ran3(iseed)*2 - 1
        w = ran3(iseed)*2 - 1
        rr = x**2+y**2+w**2
        if (rr .lt. 1) then

          call xywtoradec(x,y,w,ra,dec,z)
          call radecztoxyz_unit(ra,dec,1.0,rpi(1),rpi(2),rpi(3))

          do j=1, totalpoly
 
             do k=1,capsnum(j)
               rp(1,k)=rpx(k,j)
               rp(2,k)=rpy(k,j)
               rp(3,k)=rpz(k,j)
             END do  !k

             if(gptin(rp(:,1:capsnum(j)),cm(1:capsnum(j),j),capsnum(j),
     &          rpi(:))) then

c... find the point inside the mask

c               fgot_test = ran3(iseed)
c               if(fgot_test .le. fgot(j)) then 
                 randra(inpoint+1)   = ra
                 randdec(inpoint+1)  = dec 
                 var3(inpoint+1)     = 1
                 var4(inpoint+1)     = 1 
                 randFGOT(inpoint+1) = fgot(j)
                 randMMAX(inpoint+1) = mmax(j)            
                 inpoint             = inpoint+1
c               end if
             END if
          END do !j
        end if
      end do !while

c check file not exist
      open(unit=200, file= Trim(twoDrandfile), status='new')
      close(200)      
     
      open(unit=200, file= Trim(twoDrandfile), position='append')      
      do i = 1, rand_num
        write(200,*) randra(i), randdec(i), var3(i), var4(i),
     &               randfgot(i), randmmax(i)
      end do
      close(200)

c check in mask or not, filter with fgot
c convert to ra dec
c save as the vagc 2D rand format: ra dec v3 v4 fgot mmax

      deallocate(capsNum)
      deallocate(rpx)
      deallocate(rpy)
      deallocate(rpz)
      deallocate(cm)

      end subroutine ! vagc_generate_2Drandfile

c ---------------------------------------------------------

      SUBROUTINE generate_3D_random_mmax(
     & rannum, zmin, zmax, use_fgot, mmax_min, mmax_max,
     & vagcRand, nofzfile_prefix, randomgalaxy)
c--- This SUBROUTINE generate 3D random data from vagc2D random data file
c--- with nofzfiles of different mmax
      IMPLICIT NONE
c... input
      REAL          :: zmin
      REAL          :: zmax
      Integer       :: rannum 
      logical       :: use_fgot
      REAL          :: mmax_min
      REAL          :: mmax_max
      Character*100 :: vagcRand
      Character*100 :: nofzfile_prefix
      Character*100 :: randomgalaxy
c
      Integer       :: i,j, iseed,ixz
      REAL          :: ran3
      REAL          :: nofz(30,0:60000),nofzmax(30), xz, ynofz, ranFgot
      REAL          :: radata, decdata, var3, var4, fgot, mmax
      Integer       :: total_mmax, int_mmax_min, int_mmax_max
      character*4   :: mmax_string
      character*5   :: nofzfile_suffix = '.nofz'
      real          :: integral_nofz(30), integral_nofz_max
      real          :: norm_mmax(30)
      integer       :: mmax_index
c... set iseed  
      iseed = -1001

      int_mmax_min = int((mmax_min+0.005)*100)
      int_mmax_max = int((mmax_max+0.005)*100)
      total_mmax = int_mmax_max - int_mmax_min + 1

c read nofzfile into an array       
      Do j = 1, total_mmax
        Write(mmax_string,'(I0)') int_mmax_min+j-1
!test
c      write(*,*) "test0"
        Open(unit=100, file=
     &   TRIM(nofzfile_prefix)//TRIM(mmax_string)//nofzfile_suffix, 
     &   STATUS= 'OLD')
        Read(unit=100,fmt=*) (nofz(j,i),i=0,60000)
        close(100)
      END DO

c decide maximum of nofz and integral of nofz     
      do j = 1, total_mmax
        nofzmax(j) = 0
        integral_nofz(j) = 0
        do i=0,59999
         if(nofz(j,i) .gt. nofzmax(j)) nofzmax(j)=nofz(j,i)
c
         integral_nofz(j)=integral_nofz(j)+nofz(j,i)
        END do
      END do

c calculate normalization factor due to nofz
      integral_nofz_max = 0
      Do j =1, total_mmax
         IF(integral_nofz(j) .gt. integral_nofz_max) then
           integral_nofz_max = integral_nofz(j)
         END IF
      END DO
      Do j = 1, total_mmax
         norm_mmax(j)=integral_nofz(j)/integral_nofz_max
      END DO

c loop to read data
      Open(unit=200, file=TRIM(vagcRand), STATUS= 'OLD')  
      Open(unit =300, File=TRIM(randomgalaxy), status='New')
      close(300)
      Open(unit =300, File=TRIM(randomgalaxy), position='APPEND')
      do i=1,rannum
   50   READ(UNIT=200,fmt=*) radata, decdata, var3, var4, fgot, mmax
c decide which mmax
        mmax_index = int((mmax - mmax_min + 0.005)*100) + 1
        ranFgot=ran3(iseed)
        if (.NOT. use_fgot) fgot = 1.0 
        if (ranFgot .gt. fgot*norm_mmax(mmax_index)) goto 50 
  100   ynofz=ran3(iseed)*nofzmax(mmax_index)
        xz=ran3(iseed)*(zmax-zmin)+zmin
        ixz=Int(100000*xz)
          if (ynofz .lt. nofz(mmax_index,ixz)) then
      	    write(300,fmt=*) radata, decdata, xz
          else
c      	    Write(*,*) "pick another"
      	    goto 100
          END if
      END do
c      close(100)
      close(200)
      close(300)

      END SUBROUTINE !generat_3D_random_mmax()

c ----------------------------------------------
