      SUBROUTINE DGEEVX( BALANC, COMPVL, COMPVR, SENSE, N, A, LDA, WR,
     $                   WI, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO,
     $                   RSHIFT, CSHIFT)
*
*  -- LAPACK driver routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     August 6, 1991
*
*     .. Scalar Arguments ..
      CHARACTER          BALANC, COMPVL, COMPVR, SENSE
      INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   ABNRM,RSHIFT,CSHIFT
*
* RSHIFT and CSHIFT are used for shifts of eigenvalues.
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ),
     $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  For an N by N real nonsymmetric matrix A, compute
*
*     the eigenvalues (WR and WI)
*     the left and/or right eigenvectors (VL and VR)
*     a balancing transformation to improve the conditioning of the
*        eigenvalues and eigenvectors (ILO, IHI, SCALE, and ABNRM)
*     reciprocal condition numbers for the eigenvalues (RCONDE)
*     reciprocal condition numbers for the right eigenvectors (RCONDV)
*
*  The last four outputs are optional:
*
*     COMPVL determines whether to compute left eigenvectors VL
*     COMPVR determines whether to compute right eigenvectors VR
*     BALANC determines how to balance the matrix (output in ILO, IHI,
*        SCALE, and ABNRM)
*     SENSE determines whether to compute reciprocal condition numbers
*        RCONDE and RCONDV
*
*  Balancing a matrix means permuting the rows and columns to make it
*  more nearly upper triangular, and computing a diagonal similarity
*  D * A * D**(-1), D a diagonal matrix, to make its rows and columns
*  closer in norm and the condition numbers of its eigenvalues and
*  eigenvectors smaller. These two steps, permuting and diagonal
*  scaling, may be applied independently as determined by BALANC.
*  The one-norm of the balanced matrix (the maximum of the sum of
*  absolute values of entries of any column) is returned in ABNRM.
*
*  The reciprocal condition numbers correspond to the balanced matrix.
*  Permuting rows and columns will not change the condition numbers
*  (in exact arithmetic) but diagonal scaling will.
*
*  The reciprocal of the condition number of an eigenvalue lambda
*  is defined as
*
*          RCONDE(lambda) = |v'*u| / norm(u) * norm(v)
*
*  where u and v are the right and left eigenvectors of T
*  corresponding to lambda (v' denotes the conjugate transpose
*  of v), and norm(u) denotes the Euclidean norm.
*  These reciprocal condition numbers always lie between zero
*  (very badly conditioned) and one (very well conditioned).
*  These reciprocal condition numbers are returned in the array RCONDE.
*  An approximate error bound for a computed eigenvalue
*  WR(i)+sqrt(-1)*WI(i) is given by
*
*                      EPS * ABRNM / RCONDE(i)
*
*  where EPS = DLAMCH( 'P' ) is the machine precision.
*
*  The reciprocal condition number of the right eigenvector u
*  corresponding to lambda is defined as follows. Suppose
*
*              T = [ lambda  c  ]
*                  [   0    T22 ]
*
*  Then the reciprocal condition number is
*
*          RCONDV(lambda,T22) = sigma-min(T22 - lambda*I)
*
*  where sigma-min denotes the smallest singular value.
*  We approximate the smallest singular value by the reciprocal of
*  an estimate of the one-norm of the inverse of T22 - lambda*I.
*  When RCONDV is small, small changes in the matrix can cause large
*  changes in u. These reciprocal condition numbers are returned in
*  the array RCONDV.
*  If BALANC = 'N' or 'P', an approximate error bound for a computed
*  right eigenvector VR(i) is given by
*
*                      EPS * ABRNM / RCONDV(i)
*
*  where EPS = DLAMCH( 'P' ) is the machine precision.  When
*  BALANC = 'S' or 'B', the interpretation of RCONDV(i) is more complex.
*
*  See ``On the Conditioning of the Nonsymmetric Eigenproblem: Theory
*  and Software'' by Z. Bai, J. Demmel and A. McKenney, LAPACK Working
*  Note 13, for a detailed description of these condition numbers.
*
*  Arguments
*  =========
*
*  BALANC  (input) CHARACTER*1
*          BALANC indicates how the input matrix should be diagonally
*          scaled and/or permuted to improve the conditioning of its
*          eigenvalues.
*             = 'N' Do not diagonally scale or permute.
*             = 'P' Perform permutations to make the matrix more nearly
*                   upper triangular. Do not diagonally scale.
*             = 'S' Diagonally scale the matrix, ie. replace A by
*                   D*A*D**(-1), where D is a diagonal matrix chosen
*                   to make the rows and columns of A more equal in
*                   norm. Do not permute.
*             = 'B' Both diagonally scale and permute A.
*          Computed reciprocal condition numbers will be for the matrix
*          after balancing and/or permuting. Permuting does not change
*          condition numbers (in exact arithmetic), but balancing does.
*
*  COMPVL  (input) CHARACTER*1
*          COMPVL specifies whether or not to compute the left
*          eigenvectors of A.
*             = 'N'  left eigenvectors are not computed.
*             = 'V'  left eigenvectors are computed.
*          If SENSE = 'E' or 'B', COMPVL must = 'V'.
*
*  COMPVR  (input) CHARACTER*1
*          COMPVR specifies whether or not to compute the right
*          eigenvectors of A.
*             = 'N'  right eigenvectors are not computed.
*             = 'V'  right eigenvectors are computed.
*          If SENSE = 'E' or 'B', COMPVR must = 'V'.
*
*  SENSE   (input) CHARACTER*1
*          SENSE determines which reciprocal condition numbers are
*          computed.
*             = 'N' None are computed.
*             = 'E' Computed for eigenvalues only.
*             = 'V' Computed for right eigenvectors only.
*             = 'B' Computed for eigenvalues and right eigenvectors.
*          If SENSE = 'E' or 'B', both left and right eigenvectors
*          must also be computed (COMPVL = 'V' and COMPVR = 'V').
*
*  N       (input) INTEGER
*          N specifies the number of rows and columns of the input
*          matrix A. N must be at least 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On input, A is the matrix whose eigenvalues and eigenvectors
*          are desired.
*          On output, A has been overwritten. If COMPVL = 'V' or
*          COMPVR = 'V', the upper Hessenberg part of A has been
*          overwritten with the Schur form of the scaled and balanced
*          version of A.
*
*  LDA     (input) INTEGER
*          LDA  specifies the leading dimension of A as
*          declared in the calling (sub) program. LDA must be at
*          least max(1,N).
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          On exit, WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues. The eigenvalues
*          will not be in any particular order, except that complex
*          conjugate pairs of eigenvalues will appear consecutively
*          with the eigenvalue having the positive imaginary part
*          first.
*
*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
*          Left eigenvectors of A are the same as the right
*          eigenvectors of transpose(A).
*          The left eigenvectors will be stored one after another in
*          the columns of VL, in the same order (but not necessarily
*          the same position) as their eigenvalues. An eigenvector
*          corresponding to a real eigenvalue will take up one column.
*          An eigenvector pair corresponding to a complex conjugate
*          pair of eigenvalues will take up two columns: the first
*          column will hold the real part, the second will hold the
*          imaginary part of the eigenvector corresponding to the
*          eigenvalue with positive imaginary part.
*
*          The eigenvectors will be normalized to have Euclidean
*          norm equal to 1 and largest component real.
*
*          If COMPVL = 'N',  VL will not be referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL must be at
*          least 1, and if COMPVL = 'V', at least N.
*
*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
*          The right eigenvectors will be stored one after another in
*          the columns of VR, in the same order (but not necessarily
*          the same position) as their eigenvalues. An eigenvector
*          corresponding to a real eigenvalue will take up one column.
*          An eigenvector pair corresponding to a complex conjugate
*          pair of eigenvalues will take up two columns: the first
*          column will hold the real part, the second will hold the
*          imaginary part of the eigenvector corresponding to the
*          eigenvalue with positive imaginary part.
*
*          The eigenvectors will be normalized to have Euclidean
*          norm equal to 1 and largest component real.
*
*          If COMPVR = 'N',  VR will not be referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR must be at
*          least 1, and if COMPVR = 'V', at least N.
*
*  ILO,IHI (output) INTEGER
*          On exit, ILO, IHI and SCALE describe how A was balanced.
*          The balanced A(i,j) is equal to zero if I is greater than
*          J and J = 1,...,ILO-1 or I = IHI+1,...,N.
*
*  SCALE   (output) DOUBLE PRECISION array, dimension (N)
*          On exit, SCALE contains information determining the
*          permutations and diagonal scaling factors used in balancing.
*          Suppose that the principal submatrix in rows ILO through
*          IHI has been balanced, that P(J) denotes the index inter-
*          changed with J during the permutation step, and that the
*          elements of the diagonal matrix used in diagonal scaling
*          are denoted by D(I,J).  Then
*          SCALE(J) = P(J),    for J = 1,...,ILO-1
*                   = D(J,J),      J = ILO,...,IHI
*                   = P(J)         J = IHI+1,...,N.
*          the order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  ABNRM   (output) DOUBLE PRECISION
*          On exit, the one-norm of the balanced matrix (the maximum
*          of the sum of absolute values of entries of any column)
*          is returned in ABNRM.
*
*  RCONDE  (output) DOUBLE PRECISION array, dimension (N)
*          RCONDE(i) is the reciprocal condition number of eigenvalue
*          WR(i) + sqrt(-1)*WI(i).
*
*  RCONDV  (output) DOUBLE PRECISION array, dimension (N)
*          RCONDV(i) is the reciprocal condition number of the i-th
*          right eigenvector in VR.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, WORK(1) contains the optimal workspace size LWORK
*          for high performance.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. It must be at least 1.
*          It must also be at least
*             2*N if COMPVL = 'N' and COMPVR = 'N', and
*             3*N otherwise
*          and at least N*N+6*N unless SENSE = 'N' or 'E'.
*          For good performance, LWORK must generally be larger.
*          The optimum value of LWORK for high performance is
*          returned in WORK(1).
*
*  IWORK   (workspace) INTEGER array, dimension (2*N-2)
*          If SENSE = 'N' or 'E', not referenced.
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*          >0 - DHSEQR failed to converge:
*            r*N+j  If the calculation of A, WR, WI, and VL or VR
*                   has failed. The eigenvalues in the WR and WI arrays
*                   should be correct for indices j+1,...,N.
*                   The value of "r" indicates the nature of the
*                   failure:
*               r=0   The block multishift method failed to find all
*                     eigenvalues in 30*N iterations.
*               r=1   The call to DLAHQR to find the shifts failed.
*               r=2   The call to DLAHQR to process a subblock of A
*                     failed.
*               r=3   DLAHQR was called because LWORK was too small,
*                     and DLAHQR failed to find all eigenvalues.
*            4*N    The 2 by 2 blocks could not be standardized
*
*           4*N+1  DTRSNA failed to standardize 2 by 2 blocks
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BALOK, SENSOK, WANTVL, WANTVR, WNTBLN, WNTSNB,
     $                   WNTSNE, WNTSNN, WNTSNV
      CHARACTER          JOB
      INTEGER            HSWORK, I, IERR, ISCL, ITAU, IWRK, K, MAXB,
     $                   MAXWRK, MINWRK, NOUT
      DOUBLE PRECISION   ANRM, BIGNUM, CS, EPS, R, SCL, SMLNUM, SN
*     ..
*     .. Local Arrays ..
      LOGICAL            SELECT( 1 )
      DOUBLE PRECISION   DUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEBAK, DGEBAL, DGEHRD, DHSEQR, DLACPY, DLARTG,
     $                   DLASCL, DORGHR, DROT, DSCAL, DTREVC, DTRSNA,
     $                   XERBLA
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX, ILAENV
      DOUBLE PRECISION   DLAMCH, DLANGE, DLAPY2, DNRM2
      EXTERNAL           LSAME, IDAMAX, ILAENV, DLAMCH, DLANGE, DLAPY2,
     $                   DNRM2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      WANTVL = LSAME( COMPVL, 'V' )
      WANTVR = LSAME( COMPVR, 'V' )
      WNTBLN = LSAME( BALANC, 'N' )
      BALOK = WNTBLN .OR. LSAME( BALANC, 'S' ) .OR.
     $        LSAME( BALANC, 'P' ) .OR. LSAME( BALANC, 'B' )
      WNTSNN = LSAME( SENSE, 'N' )
      WNTSNE = LSAME( SENSE, 'E' )
      WNTSNV = LSAME( SENSE, 'V' )
      WNTSNB = LSAME( SENSE, 'B' )
      SENSOK = WNTSNN .OR. WNTSNE .OR. WNTSNB .OR. WNTSNV
      IF( .NOT.BALOK ) THEN
         INFO = -1
      ELSE IF( ( .NOT.WANTVL ) .AND. ( .NOT.LSAME( COMPVL, 'N' ) ) )
     $          THEN
         INFO = -2
      ELSE IF( ( .NOT.WANTVR ) .AND. ( .NOT.LSAME( COMPVR, 'N' ) ) )
     $          THEN
         INFO = -3
      ELSE IF( ( .NOT.SENSOK ) .OR. ( ( WNTSNE .OR. WNTSNB ) .AND. .NOT.
     $         ( WANTVL .AND. WANTVR ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( WANTVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( WANTVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV.
*       HSWORK refers to the workspace preferred by DHSEQR, as
*       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
*       the worst case.)
*
      MINWRK = 1
      IF( INFO.EQ.0 .AND. LWORK.GE.1 ) THEN
         MAXWRK = N + N*ILAENV( 1, 'DGEHRD', ' ', N, 1, N, 0 )
         IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
            MINWRK = MAX( 1, 2*N )
            IF( .NOT.WNTSNN )
     $         MINWRK = MAX( MINWRK, N*N+6*N )
            MAXB = MAX( ILAENV( 3, 'DHSEQR', 'SN', N, 1, N, -1 ), 2 )
            IF( WNTSNN ) THEN
               K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'EN', N,
     $             1, N, -1 ) ) )
            ELSE
               K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'SN', N,
     $             1, N, -1 ) ) )
            END IF
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, 1, HSWORK )
            IF( .NOT.WNTSNN )
     $         MAXWRK = MAX( MAXWRK, N*N+6*N )
         ELSE
            MINWRK = MAX( 1, 3*N )
            IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) )
     $         MINWRK = MAX( MINWRK, N*N+6*N )
            MAXB = MAX( ILAENV( 3, 'DHSEQR', 'SN', N, 1, N, -1 ), 2 )
            K = MIN( MAXB, N, MAX( 2, ILAENV( 4, 'DHSEQR', 'EN', N, 1,
     $          N, -1 ) ) )
            HSWORK = MAX( K*( K+2 ), 2*N )
            MAXWRK = MAX( MAXWRK, 1, HSWORK )
            MAXWRK = MAX( MAXWRK, N+( N-1 )*
     $               ILAENV( 1, 'DORGHR', ' ', N, 1, N, -1 ) )
            IF( ( .NOT.WNTSNN ) .AND. ( .NOT.WNTSNE ) )
     $         MAXWRK = MAX( MAXWRK, N*N+6*N )
            MAXWRK = MAX( MAXWRK, 3*N, 1 )
         END IF
         WORK( 1 ) = MAXWRK
      END IF
      IF( LWORK.LT.MINWRK ) THEN
         INFO = -21
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEEVX', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) / EPS )
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max entry outside range [SMLNUM,BIGNUM]
*
      ANRM = DLANGE( 'M', N, N, A, LDA, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL DLASCL( 'G', 0, 0, ANRM, SMLNUM, N, N, A, LDA, INFO )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 2
         CALL DLASCL( 'G', 0, 0, ANRM, BIGNUM, N, N, A, LDA, INFO )
      END IF
*
*     Balance the matrix and compute ABNRM
*
      CALL DGEBAL( BALANC, N, A, LDA, ILO, IHI, SCALE, IERR )
      ABNRM = DLANGE( '1', N, N, A, LDA, DUM )
      DUM( 1 ) = ABNRM
      IF( ISCL.EQ.1 )
     $   CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, 1, 1, DUM, 1, INFO )
      IF( ISCL.EQ.2 )
     $   CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, 1, 1, DUM, 1, INFO )
      ABNRM = DUM( 1 )
      ITAU = 1
      IWRK = ITAU + N
*
*     Reduce to upper Hessenberg form
*     (Workspace: need 2*N, prefer N+N*NB)
*
      CALL DGEHRD( N, ILO, IHI, A, LDA, WORK( ITAU ), WORK( IWRK ),
     $             LWORK-IWRK+1, IERR )
*
*     Perform QR algorithm, compute eigenvectors if desired,
*     and unbalance
*
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
*
*        Compute eigenvalues only
*        If condition numbers desired, compute Schur form
*
         IF( WNTSNN ) THEN
            JOB = 'E'
         ELSE
            JOB = 'S'
         END IF
*
*        (Workspace: need 1, prefer HSWORK (see comments) )
         IWRK = ITAU
         CALL DHSEQR( JOB, 'N', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO, RSHIFT, CSHIFT )
*
*        Compute condition numbers if desired
*        (Workspace: need N*N+6*N unless SENSE = 'N')
*
         IF( .NOT.WNTSNN ) THEN
            CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR,
     $                   LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N,
     $                   IWORK, INFO )
            IF( INFO.EQ.1 )
     $         INFO = 4*N + 1
         END IF
*
      ELSE IF( ( .NOT.WANTVL ) .AND. WANTVR ) THEN
*
*        Compute right eigenvectors
*
*        Copy Householder vectors to VR
*
         CALL DLACPY( 'L', N, N, A, LDA, VR, LDVR )
*
*        Generate orthogonal matrix in VR
*        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VR, LDVR, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VR
*        (Workspace: need 1, prefer HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VR, LDVR,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
         IF( INFO.GT.0 )
     $      GO TO 80
*
*        Compute right eigenvectors
*        (Workspace: need 3*N)
*
         CALL DTREVC( 'R', 'O', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                N, NOUT, WORK( IWRK ), IERR )
*
*        Compute condition numbers if desired
*        (Workspace: need N*N+6*N unless SENSE = 'N')
*
         IF( .NOT.WNTSNN ) THEN
            CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR,
     $                   LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N,
     $                   IWORK, INFO )
            IF( INFO.EQ.1 )
     $         INFO = 4*N + 1
         END IF
*
*        Undo balancing
*
         CALL DGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR,
     $                IERR )
*
*        Normalize eigenvectors to have Euclidean norm 1
*        and largest component real
*
         DO 20 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ),
     $               DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 10 K = 1, N
                  WORK( K ) = VR( K, I )**2 + VR( K, I+1 )**2
   10          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
   20    CONTINUE
*
      ELSE IF( WANTVL .AND. ( .NOT.WANTVR ) ) THEN
*
*        Compute left eigenvectors
*
*        Copy Householder vectors to VL
*
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
*
*        Generate orthogonal matrix in VL
*        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (Workspace: need 1, prefer HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
         IF( INFO.GT.0 )
     $      GO TO 80
*
*        Compute left eigenvectors
*        (Workspace: need 3*N)
*
         CALL DTREVC( 'L', 'O', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                N, NOUT, WORK( IWRK ), IERR )
*
*        Compute condition numbers if desired
*        (Workspace: need N*N+6*N unless SENSE = 'N')
*
         IF( .NOT.WNTSNN ) THEN
            CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR,
     $                   LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N,
     $                   IWORK, INFO )
            IF( INFO.EQ.1 )
     $         INFO = 4*N + 1
         END IF
*
*        Undo balancing
*
         CALL DGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL,
     $                IERR )
*
*        Normalize eigenvectors to have Euclidean norm 1
*        and largest component real
*
         DO 40 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ),
     $               DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 30 K = 1, N
                  WORK( K ) = VL( K, I )**2 + VL( K, I+1 )**2
   30          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   40    CONTINUE
*
      ELSE IF( WANTVL .AND. WANTVR ) THEN
*
*        Compute right and left eigenvectors
*
*        Copy Householder vectors to VL
*
         CALL DLACPY( 'L', N, N, A, LDA, VL, LDVL )
*
*        Generate orthogonal matrix in VL
*        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
*
         CALL DORGHR( N, ILO, IHI, VL, LDVL, WORK( ITAU ), WORK( IWRK ),
     $                LWORK-IWRK+1, IERR )
*
*        Perform QR iteration, accumulating Schur vectors in VL
*        (Workspace: need 1, prefer HSWORK (see comments) )
*
         IWRK = ITAU
         CALL DHSEQR( 'S', 'V', N, ILO, IHI, A, LDA, WR, WI, VL, LDVL,
     $                WORK( IWRK ), LWORK-IWRK+1, INFO )
         IF( INFO.GT.0 )
     $      GO TO 80
*
*        Copy Schur vectors to VR
*
         CALL DLACPY( 'F', N, N, VL, LDVL, VR, LDVR )
*
*        Compute eigenvectors
*        (Workspace: need 3*N)
*
         CALL DTREVC( 'B', 'O', SELECT, N, A, LDA, VL, LDVL, VR, LDVR,
     $                N, NOUT, WORK( IWRK ), IERR )
*
*        Compute condition numbers if desired
*        (Workspace: need N*N+6*N unless SENSE = 'N' or 'E')
*
         IF( .NOT.WNTSNN ) THEN
            CALL DTRSNA( SENSE, 'A', SELECT, N, A, LDA, VL, LDVL, VR,
     $                   LDVR, RCONDE, RCONDV, N, NOUT, WORK( IWRK ), N,
     $                   IWORK, INFO )
            IF( INFO.EQ.1 )
     $         INFO = 4*N + 1
         END IF
*
*        Undo balancing
*
         CALL DGEBAK( BALANC, 'L', N, ILO, IHI, SCALE, N, VL, LDVL,
     $                IERR )
         CALL DGEBAK( BALANC, 'R', N, ILO, IHI, SCALE, N, VR, LDVR,
     $                IERR )
*
*        Normalize eigenvectors to have Euclidean norm 1
*        and largest component real
*
         DO 70 I = 1, N
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VR( 1, I ), 1 ),
     $               DNRM2( N, VR( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VR( 1, I ), 1 )
               CALL DSCAL( N, SCL, VR( 1, I+1 ), 1 )
               DO 50 K = 1, N
                  WORK( K ) = VR( K, I )**2 + VR( K, I+1 )**2
   50          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VR( K, I ), VR( K, I+1 ), CS, SN, R )
               CALL DROT( N, VR( 1, I ), 1, VR( 1, I+1 ), 1, CS, SN )
               VR( K, I+1 ) = ZERO
            END IF
            IF( WI( I ).EQ.ZERO ) THEN
               SCL = ONE / DNRM2( N, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
            ELSE IF( WI( I ).GT.ZERO ) THEN
               SCL = ONE / DLAPY2( DNRM2( N, VL( 1, I ), 1 ),
     $               DNRM2( N, VL( 1, I+1 ), 1 ) )
               CALL DSCAL( N, SCL, VL( 1, I ), 1 )
               CALL DSCAL( N, SCL, VL( 1, I+1 ), 1 )
               DO 60 K = 1, N
                  WORK( K ) = VL( K, I )**2 + VL( K, I+1 )**2
   60          CONTINUE
               K = IDAMAX( N, WORK, 1 )
               CALL DLARTG( VL( K, I ), VL( K, I+1 ), CS, SN, R )
               CALL DROT( N, VL( 1, I ), 1, VL( 1, I+1 ), 1, CS, SN )
               VL( K, I+1 ) = ZERO
            END IF
   70    CONTINUE
*
      END IF
*
*     Undo scaling if necessary
*
      IF( ISCL.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, N, 1, WR, N, INFO )
         CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, N, 1, WI, N, INFO )
         IF( ( WNTSNV .OR. WNTSNB ) .AND. INFO.NE.4*N+1 )
     $      CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, N, 1, RCONDV, N,
     $                   INFO )
      ELSE IF( ISCL.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, N, 1, WR, N, INFO )
         CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, N, 1, WI, N, INFO )
         IF( ( WNTSNV .OR. WNTSNB ) .AND. INFO.NE.4*N+1 )
     $      CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, N, 1, RCONDV, N,
     $                   INFO )
      END IF
*
   80 CONTINUE
      WORK( 1 ) = MAXWRK
*
      RETURN
*
*     End of DGEEVX
*
      END
