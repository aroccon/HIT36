      module RANDu
c..
      implicit none

      real*8  :: rseed

c..
c      double precision random
c..
      contains
c----------------------------------------------------
      function random(idum3)
c..
c..   initialize with idum<0; Then keep idum>0 unchanged for the same
c..   sequence
c..
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      double precision random,am,eps,rnmx,idum3
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     +   ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,
     +   ir2=3791,ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.0-eps)
      integer idum2,j,k,iv(ntab),iy
c..
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
c..
      idum=nint(idum3/987)
      if (idum.le.0) then 
       idum=max(-idum,1)
       idum2=idum
       do 11 j=ntab+8,1,-1
        k=idum/iq1
        idum=ia1*(idum-k*iq1)-k*ir1
        if(idum.lt.0) idum=idum+im1
        if(j.le.ntab) iv(j)=idum
 11    continue
      endif 
c.. start here when not initializing
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      random=min(am*iy,rnmx)
c..
      return
      end function random
c--------------------------------------------------
      end module RANDu
