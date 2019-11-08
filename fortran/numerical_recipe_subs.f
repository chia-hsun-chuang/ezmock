c This file includes the subroutines from numerical recipe and 
c the other random number generators.
c --------------------------------------
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
c see numerical recipe
      IMPLICIT NONE 
      INTEGER n,NMAX
      REAL    yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=20000)
      INTEGER i,k
      REAL    p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+
     *1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*
     *u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END SUBROUTINE !spline()
c ---------------------------------------------
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
c see numerical recipe
      IMPLICIT NONE 
      INTEGER n
      REAL    x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL    a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**
     *2)/6.
      return
      END SUBROUTINE ! splint()
c---------------------------------------------
      SUBROUTINE splie2(x1a,x2a,ya,m,n,y2a)
      INTEGER m,n,NN
      REAL x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
CU    USES spline
      INTEGER j,k
      REAL y2tmp(NN),ytmp(NN)
      do 13 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
11      continue
        call spline(x2a,ytmp,n,1.e30,1.e30,y2tmp)
        do 12 k=1,n
          y2a(j,k)=y2tmp(k)
12      continue
13    continue
      return
      END
c------------------------------------------
c from numerical recipe
      SUBROUTINE splin2(x1a,x2a,ya,y2a,m,n,x1,x2,y)
      INTEGER m,n,NN
      REAL x1,x2,y,x1a(m),x2a(n),y2a(m,n),ya(m,n)
      PARAMETER (NN=100)
CU    USES spline,splint
      INTEGER j,k
      REAL y2tmp(NN),ytmp(NN),yytmp(NN)
      do 12 j=1,m
        do 11 k=1,n
          ytmp(k)=ya(j,k)
          y2tmp(k)=y2a(j,k)
11      continue
        call splint(x2a,ytmp,y2tmp,n,x2,yytmp(j))
12    continue
      call spline(x1a,yytmp,m,1.e30,1.e30,y2tmp)
      call splint(x1a,yytmp,y2tmp,m,x1,y)
      return
      END

c-------------- ran3() ---------------------------------- 
      FUNCTION ran3(idum)
c see numerical recipe
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL    MBIG,MSEED,MZ
      REAL    ran3,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL    mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran3=mj*FAC
      return
      END FUNCTION ! ran3()
c------------------------------------------

       SUBROUTINE RMARIN(IJ,KL)
       IMPLICIT NONE
! This is the initialization routine for the random number generator RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
!The random number sequences created by these two seeds are of sufficient
! length to complete an entire calculation with. For example, if sveral
! different groups are working on different parts of the same calculation,
! each group could be assigned its own IJ seed. This would leave each group
! with 30000 choices for the second seed. That is to say, this random
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.
!
! Use IJ = 1802 & KL = 9373 to test the random number generator. The
! SUBROUTINE RANMAR should be used to generate 20000 random numbers.
! Then display the next six random numbers generated multiplied by 4096*4096
! If the random number generator is working properly, the random numbers
!    should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0
       REAL    U(97), C, CD, CM
       integer I97, J97, IJ, KL
       integer i,j,k,l,ii,jj,m
       REAL    s,t
!      INTEGER IRM(103)

       common /RASET1/ U, C, CD, CM, I97, J97
       if(IJ.lt.0.or.IJ.gt.31328.or.KL.lt.0.or.KL.gt.30081) then
           print '(A)', ' The first random number seed must'
       print '(A)', ' have a value  between 0 and 31328'
           print '(A)',' The second seed must have a value'
       print '(A)', 'between 0 and   30081'
             stop
       endif
       I = mod(IJ/177, 177) + 2
       J = mod(IJ    , 177) + 2
       K = mod(KL/169, 178) + 1
       L = mod(KL,     169)
       do 2 II = 1, 97
          S = 0.0
          T = 0.5
          do 3 JJ = 1, 24
             M = mod(mod(I*J, 179)*K, 179)
             I = J
             J = K
             K = M
             L = mod(53*L+1, 169)
             if (mod(L*M, 64) .ge. 32) then
                S = S + T
             endif
             T = 0.5 * T
    3        continue
          U(II) = S
    2     continue
       C = 362436.0 / 16777216.0
       CD = 7654321.0 / 16777216.0
       CM = 16777213.0 /16777216.0
       I97 = 97
       J97 = 33

       return
       END SUBROUTINE !RMARIN()
c-------------------------------------
       REAL    FUNCTION RANMAR()
       IMPLICIT NONE
! This is the random number generator proposed by George Marsaglia in
! Florida State University Report: FSU-SCRI-87-50
! It was slightly modified by F. James to produce an array of pseudorandom
! numbers.
       REAL    U(97), C, CD, CM
       integer I97, J97
       REAL    UNI
       common /RASET1/ U, C, CD, CM, I97, J97
!      INTEGER IVEC
          UNI = U(I97) - U(J97)
          if( UNI .lt. 0.0 ) UNI = UNI + 1.0
          U(I97) = UNI
          I97 = I97 - 1
          if(I97 .eq. 0) I97 = 97
          J97 = J97 - 1
          if(J97 .eq. 0) J97 = 97
          C = C - CD
          if( C .lt. 0.0 ) C = C + CM
          UNI = UNI - C
          if( UNI .lt. 0.0 ) UNI = UNI + 1.0 ! bug?
          RANMAR = UNI
       return
       END FUNCTION !RANMAR()
c -----------------------------------------	
