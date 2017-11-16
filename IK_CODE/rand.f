      real  function urand(iy)
      integer iy
c
c      urand is a uniform random number generator based  on  theory  and
c  suggestions  given  in  d.e. knuth (1969),  vol  2.   the integer  iy
c  should be initialized to an arbitrary integer prior to the first call
c  to urand.  the calling program should  not  alter  the  value  of  iy
c  between  subsequent calls to urand.  values of urand will be returned
c  in the interval (0,1).
c
      integer  ia,ic,itwo,m2,m,mic
      double precision  halfm
      double precision  s
      double precision  datan,dsqrt
      data m2/0/,itwo/2/
      if (m2 .ne. 0) go to 20
c
c  if first entry, compute machine integer word length
c
      m = 1
   10 m2 = m
      m = itwo*m2
      if (m .gt. m2) go to 10
      halfm = m2
c
c  compute multiplier and increment for linear congruential method
c
      ia = 8*idint(halfm*datan(1.d0)/8.d0) + 5
      ic = 2*idint(halfm*(0.5d0-dsqrt(3.d0)/6.d0)) + 1
      mic = (m2 - ic) + m2
c
c  s is the scale factor for converting to floating point
c
      s = 0.5/halfm
c
c  compute next random number
c
   20 iy = iy*ia
c
c  the following statement is for computers which do not allow
c  integer overflow on addition
c
      if (iy .gt. mic) iy = (iy - m2) - m2
c
      iy = iy  + ic
c
c  the following statement is for computers where the
c double precision  word length for addition is greater than for multiplication
c
      if (iy/2 .gt. m2) iy = (iy - m2) - m2
c
c  the following statement is for computers where integer
c  overflow affects the sign bit
c
      if (iy .lt. 0) iy = (iy + m2) + m2
      urand = float(iy)*s
      return
      end


