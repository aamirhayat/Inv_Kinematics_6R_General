      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(lda,2),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     this routine has been modified by dinesh manocha to use complete 
c     pivoting instead of partial pivoting.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a as well as ipvt .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n,2)
c                an integer vector of pivot indices. Each entry of the 
c                vector has 2 coordinates, corresponding to the row 
c                and column number of the pivot element.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,i,j,k,kp1,l,row,col,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
c
c        l = idamax(n-k+1,a(k,k),1) + k - 1
c        ipvt(k) = l
c
c
c        here we do complete pivoting as compared to partial pivoting in 
c        the original code from netlib.
c
c        the main idea is to find the maximum in each column using the 
c        implementation from each column. 
         
         row = idamax(n-k+1,a(k,k),1) + k - 1
         col = k
         do 5 i = k+1,n,1
             l = idamax(n-k+1,a(k,i),1) + k - 1
             if (abs(a(row,col)) .ge. abs(a(l,i))) go to 5     
             row = l
             col = i
 5       continue
         ipvt(k,1) = row
         ipvt(k,2) = col

c
c        zero pivot implies this column already triangularized
c
         if (a(row,col) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
         if (col .eq. k) go to 8
c
c           perform the column interchange
c
c
c         write(6,*) "Interchanging Columns"

          do 7 i = 1,n,1
              t = a(i,k)
              a(i,k) = a(i,col)
              a(i,col) = t
 7        continue
 8        continue
               t = a(row,k)
               a(row,k) = a(k,k)
               a(k,k) = t
  10        continue
c
c           compute multipliers
c
            t = -1.0e0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           ROW elimination with column indexing
c
            do 30 j = kp1, n
               t = a(row,j)
               if (row .eq. k) go to 20
                  a(row,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n,1) = n
      ipvt(n,2) = n
      if (a(n,n) .eq. 0.0e0) info = n
      return
      end


      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(lda,2),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a and ipvt.
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n,2)
c                the pivot vector (for complete pivoting) from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,k1,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k,1)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k,1)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
c
c     Now we adjust for the columns that were being interchanged during 
c     the pivoting process, and that reflects in the order of the variables
c     that are being produced.
c
      if (nm1 .lt. 1) go to 120
      do 110 k = 1, nm1
          k1 = n - k
          l = ipvt(k1,2)
          t = b(l)
          if (l .eq. k1) go to 110
          b(l) = b(k1)
          b(k1) = t
 110  continue

 120  return
      end
