! This file is part of AGMG software package,
! Release 3.1.2 built on "Jan 17 2012"
!
!    AGMG is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AGMG is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AGMG.  If not, see <http://www.gnu.org/licenses/>.
!
! Up-to-date copies of the AGMG package can be obtained
! from the Web pages <http://homepages.ulb.ac.be/~ynotay/AGMG>.
!
! You can acknowledge, citing references [1] [2], and [3], the contribution
! of this package in any scientific publication dependent upon it use.
!
! [1] Y. Notay, An aggregation-based algebraic multigrid method,
!    Electronic Transactions on Numerical Analysis, vol. 37, pp. 123-146, 2010
!
! [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
!    convergence rate, to appear in SIAM J. Sci. Comput., 2012.
!
! [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
!    equations, Report GANMN 11-01, Universite Libre de Bruxelles, Brussels,
!    Belgium, 2011.
!
! See the accompanying userguide for more details on how to use the software,
! and the README file for installation instructions.
!
! AGMG Copyright (C) 2011 Yvan NOTAY
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file provides zagmg and dependencies; zagmg is a sequential
! implementation for complex matrices in double precision
! of the method presented in [1], where the used algorithms are described
! in detail. From realease 3.x, the coarsening has been modified according
! to the results in [2,3], see the release notes in the README file
! for more details.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! SUBROUTINE zagmg (MAIN DRIVER): see bottom of this file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! PARAMETERS DEFINITON -  may be tuned by expert users
!-----------------------------------------------------------------------
  MODULE zagmg_mem
    SAVE
!
! INTEGER
!
!  maxlev  is the maximal number of levels
!          (should be large enough - much larger than needed is armless).
!  complex_len is the length of 1 COMPLEX(kind(0.0d0)) in byte
!        (used only to display information on memory usage).
!  nsmooth  is the number of pre- and post-smoothing steps;
!  smoothtype indicates which smoother use:
!    if smoothtype==1, the smoothing scheme is
!        pre-smoothing: Forward Gauss-Seidel, then Backward GS, then Fwd GS,etc
!       post-smoothing: Backward Gauss-Seidel, then Forward GS, then Bwd GS,etc
!                       (nsmooth sweeps for both pre- and post-smoothing)
!  nstep   is the maximum number of coarsening steps
!          nstep<0 means that coarsening is stopped only according to
!          the matrix size, see parameter maxcoarsesize.
!  nlvcyc  is the number of coarse levels (from bottom) on which V-cycle
!          formulation is enforced (Rmk: K-cycle always allowed at first
!          coarse level).
!  npass   is the maximal number of pairwise aggregation passes for each
!          coarsening step, according to the algorithms in [2,3].
!  maxcoarsesize: when the size of the coarse grid matrix is less than or
!                 equal to maxcoarsesize*N^(1/3),  it is factorized exactly
!                 and coarsening is stopped;
!         maxcoarsesizeslow: in case of slow coarsening,
!                 exact factorization can be performed when the size of
!                 the coarse grid matrix is less than or equal to
!                 maxcoarsesizeslow*N^(1/3).
!         (N is the number of rows of the input matrix)
    INTEGER, PARAMETER :: maxlev=40, complex_len=16
    INTEGER, PARAMETER :: nsmooth=1, smoothtype=1, nstep=-1, nlvcyc=0
    INTEGER, PARAMETER :: npass=2,maxcoarsesize=40,maxcoarsesizeslow=400
!
! REAL
!
!  resi is the threshold t for the relative residual error in inner FCG & GCR
!       iterations, see Algorithm 3.2 in [1]
!  trspos is a threshold: if a row has a positive offdiagonal entry larger
!         than trspos times the diagonal entry, the corresponding node is
!         transferred unaggregated to the coarse grid
!  kaptg_ is the threshold used to accept or not a tentative aggregate
!         when applying the coarsening algorithms from [2,3];
!         kaptg_blocdia is used for control based on bloc diagonal smoother [2];
!         kaptg_dampJac is used for control based on Jacobi smoother [3].
!  checkdd is the threshold to keep outside aggregation nodes where
!         the matrix is strongly diagonally dominant (based on mean of row
!         and column);
!         In fact, AGMG uses the maximum of |checkdd| and of the default value
!            according to kaptg_ as indicated in [2,3]
!            (hence |checkdd| < 1 ensures that one uses this default value)
!         checkdd <0 : consider |checkdd|, but base the test on minus the
!                sum of offdiagonal elements, without taking the absolute value.
!  targetcoarsefac is the target coarsening factor (parameter tau in the main
!         coarsening algorithm in [2,3]): further pairwise aggregation passes
!         are omitted once the number of nonzero entries has been reduced by a
!         factor of at least targetcoarsefac.
!  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes,
!         where nl is the total number of nodes at that level, have
!         negative mean row and column sum, then the aggregation algorithm
!         of [2,3] is modified, exchanging all diagonal entries for the mean
!         row and column sum (that is, the algorithm is applied to a
!         modified matrix with mean row and colum sum enforced to be zero).
    REAL(kind(0.0d0)), PARAMETER :: resi=0.2, trspos=0.45
    REAL(kind(0.0d0)), PARAMETER :: kaptg_blocdia=8, kaptg_dampJac=10
    REAL(kind(0.0d0)), PARAMETER :: checkdd=0.5
    REAL(kind(0.0d0)), PARAMETER :: targetcoarsefac=2.0**npass
    REAL(kind(0.0d0)), PARAMETER :: fracnegrcsum=0.25
!!!!!!!!!!!!!!!!!!!! END of PARAMETERS DEFINITION -----------------
!!!!!!!!!!!!!!!!!!! Internal variables declaration
!
    TYPE InnerData
       COMPLEX(kind(0.0d0)), DIMENSION(:), POINTER :: a
       INTEGER, DIMENSION(:), POINTER :: ja
       INTEGER, DIMENSION(:), POINTER :: ia
       INTEGER, DIMENSION(:), POINTER :: il
       INTEGER, DIMENSION(:), POINTER :: iu
       COMPLEX(kind(0.0d0)), DIMENSION(:), POINTER :: p
       INTEGER, DIMENSION(:), POINTER :: idiag
       INTEGER, DIMENSION(:), POINTER :: ind
       INTEGER, DIMENSION(:), POINTER :: iext
       INTEGER, DIMENSION(:), POINTER :: ilstout
       INTEGER, DIMENSION(:), POINTER :: lstout
       INTEGER, DIMENSION(:), POINTER :: ilstin
       INTEGER, DIMENSION(:), POINTER :: lstin
       INTEGER, DIMENSION(:), POINTER :: iblockl
    END TYPE InnerData
!
    TYPE(InnerData) :: dt(maxlev)
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: scald(:)
    INTEGER :: nn(maxlev),kstat(2,maxlev)=0,innermax(maxlev)
    INTEGER :: nlev,nwrkcum,iout,nrst,nbblock
    INTEGER :: maxcoarset,maxcoarseslowt
    REAL(kind(0.0d0)) :: memi=0.0,memax=0.0,memr=0.0,mritr,rlenilen
    REAL(kind(0.0d0)) :: wcplex(4),fracnz(maxlev)
    LOGICAL :: spd,wfo,wff,allzero,trans,transint,zerors,gcrin
    REAL(kind(0.0d0)), PARAMETER :: cplxmax=3.0, xsi=0.6d0
    REAL(kind(0.0d0)), PARAMETER :: repsmach=SQRT(EPSILON(1.0d0))
    INTEGER :: nlc(2),nlcp(2),nlc1(2),icum,imult
    REAL(kind(0.0d0)) :: ngl(2),nglp(2),nlctot(2),ngl1(2),ngltot(2)
    INTEGER, DIMENSION(:), POINTER :: iextL1
    INTEGER, PARAMETER :: IRANK=-9999
    REAL(kind(0.0d0)), PARAMETER ::    &
               checkddJ=MAX(ABS(checkdd),kaptg_dampJac/(kaptg_dampJac-2))
    REAL(kind(0.0d0)), PARAMETER ::    &
               checkddB=MAX(ABS(checkdd),(kaptg_blocdia+1)/(kaptg_blocdia-1))
  END MODULE zagmg_mem
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of Internal variables declaration
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! TIMING
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_mestime(id,cputm,eltm)
    IMPLICIT NONE
    INTEGER, SAVE :: cpt_init(10)=-1,cpt_fin,cpt_max,freq,cpt
    REAL, SAVE :: t1(10), t2
    REAL(kind(0.0d0)) :: cputm,eltm
    INTEGER :: id
    IF (id>0) THEN
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t2)
       CALL SYSTEM_CLOCK(cpt_fin,freq,cpt_max)
       !
       cpt = cpt_fin - cpt_init(id)
       IF (cpt_fin < cpt_init(id)) cpt = cpt + cpt_max
       eltm = dble(cpt) / freq
       cputm = dble(t2 - t1(id))
       !
    ELSE
       !
       CALL SYSTEM_CLOCK(cpt_init(-id),freq,cpt_max)
       !Next line may be uncommented if FORTRAN 95 function
       !CPU_TIME is implemented
       !   CALL CPU_TIME(t1(-id))
       !
    END IF
    RETURN
  END SUBROUTINE zagmg_mestime
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of TIMING
!-----------------------------------------------------------------------
MODULE zagmg_ALLROUTINES
CONTAINS
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_relmem
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    !
    DEALLOCATE(scald)
    DO l=1,nlev-1
       IF(nn(l) > 0) DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%il,dt(l)%iu)
       IF (nn(l+1) > 0) DEALLOCATE(dt(l)%ind)
    END DO
    memi=0
    memr=0
    memax=0
    RETURN
  END SUBROUTINE zagmg_relmem
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_applyprec( N,f,X,a,ja,ia,flop)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER       :: N
    REAL(kind(0.0d0)) ::  flop
    COMPLEX(kind(0.0d0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: S(:)
    !
    mritr=nwrkcum+2*N
    ALLOCATE (S(nwrkcum+2*N))
    S(1:N)=f(1:n)*scald(1:N)
    CALL zagmg_prec_matv(1,S,X,S(N+1),flop,S(2*N+1),.FALSE.)
    X(1:N)=X(1:n)*scald(1:N)
    flop=flop+dble(2*N)
    kstat(2,1)=kstat(2,1)+1
    DEALLOCATE (S)
    !
    RETURN
  END SUBROUTINE zagmg_applyprec
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_partroword(n, a, ja, ia, idiag, w, iw, iext)
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER, TARGET :: iw(*)
    COMPLEX(kind(0.0d0)) :: a(*), p
    COMPLEX(kind(0.0d0)), TARGET :: w(*)
    INTEGER :: i, j, k, ipos, nzu, nze
    DO i = 1, n
       ipos = ia(i)
       nzu = 0
       nze = 0
       DO k = ia(i), ia(i+1) - 1
          j = ja(k)
          IF (j.GT.i) THEN
             nzu = nzu + 1
             w(nzu) = a(k)
             iw(nzu) = j
          ELSEIF (j.LT.i) THEN
             a(ipos) = a(k)
             ja(ipos) = j
             ipos = ipos + 1
          ELSE
             p = a(k)
          ENDIF
       ENDDO
       !
       idiag(i) = ipos
       a(ipos) = p
       ja(ipos) = i
       !
       DO k = 1, nzu
          ipos = ipos + 1
          a(ipos) = w(k)
          ja(ipos) = iw(k)
       ENDDO
       !
    ENDDO
    RETURN
  END SUBROUTINE zagmg_partroword
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_diagscal(n,a,ja,ia,idiag,as)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n)
    COMPLEX(kind(0.0d0)) :: a(*), val
    COMPLEX(kind(0.0d0)) :: as(*)
    INTEGER :: i,k,k1,idum(1),ier
    DO i=1,n
       val=a(idiag(i))
       scald(i)=sqrt(abs(val)/val)
    END DO
    DO i=1,n
       DO k=ia(i),ia(i+1)-1
          as(k)=a(k)*scald(i)*scald(ja(k))
       END DO
    END DO
    RETURN
  END SUBROUTINE zagmg_diagscal
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_csrdlu(n,a,ja,ia,idiag,ao,jao,il,iu,trans,iext,iexto)
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), idiag(n), jao(n+1:*), il(n+1), iu(n+1)
    INTEGER, OPTIONAL :: iext(n),iexto(n+1)
    COMPLEX(kind(0.0d0)) :: a(*),ao(*)
    INTEGER :: i,j,ili,iui,iei,ipos,ili0,iui0,iei0,nl,nu,k,next
    LOGICAL :: trans
    !
    !
    nl=0
    nu=0
    IF (.NOT.trans) THEN
       DO i=1,n
          ili=idiag(i)-ia(i)
          iui=ia(i+1)-idiag(i)-1
          iu(i)=iui
          il(i)=ili
          nl=nl+ili
          nu=nu+iui
       END DO
    ELSE
       il(2:n+1)=0
       iu(2:n+1)=0
       DO i=1,n
          iui=idiag(i)-ia(i)
          ili=ia(i+1)-idiag(i)-1
          nl=nl+ili
          nu=nu+iui
          DO k=ia(i),idiag(i)-1
             iu(ja(k)+1)=iu(ja(k)+1)+1
          END DO
          DO k=idiag(i)+1,ia(i+1)-1
             il(ja(k)+1)=il(ja(k)+1)+1
          END DO
       END DO
    END IF
    !
    !
    ipos=0
    ili=n+1
    iui=ili+nl
    iei=iui+nu
    IF (.NOT.trans) THEN
       DO i=1,n
          ili0=ili
          iui0=iui
          DO j=1,il(i)
             ipos=ipos+1
             ao(ili)=a(ipos)
             jao(ili)=ja(ipos)
             ili=ili+1
          END DO
          ipos=ipos+1
          ao(i)=a(ipos)
          DO j=1,iu(i)
             ipos=ipos+1
             ao(iui)=a(ipos)
             jao(iui)=ja(ipos)
             iui=iui+1
          END DO
          iu(i)=iui0
          il(i)=ili0
       END DO
       iu(n+1)=iui
       il(n+1)=ili
    ELSE
       il(1)=ili
       iu(1)=iui
       DO i=1,n
          il(i+1) = il(i) + il(i+1)
          iu(i+1) = iu(i) + iu(i+1)
       END DO
       DO i=1,n
          DO k=ia(i),idiag(i)-1
             j = ja(k)
             next = iu(j)
             ao(next) = a(k)
             jao(next) = i
             iu(j) = next+1
          END DO
          ao(i)=a(idiag(i))
          DO k=idiag(i)+1,ia(i+1)-1
             j = ja(k)
             next = il(j)
             ao(next) = a(k)
             jao(next) = i
             il(j) = next+1
          END DO
       END DO
       DO i=n,1,-1
          il(i+1) = il(i)
          iu(i+1) = iu(i)
       END DO
       il(1)=ili
       iu(1)=iui
    END IF
    RETURN
  END SUBROUTINE zagmg_csrdlu
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_csrmv(n, a, ja, ifja, ia, x, ifx, y, iad, flop, lstin)
!
!
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iad, ifx, i, k, j
    INTEGER, OPTIONAL :: lstin(0:*)
    COMPLEX(kind(0.0d0)) :: a(*), x(ifx:*), y(n), t
    REAL(kind(0.0d0)) :: flop
    !
       IF (iad .LT. -1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad .GT. 1) THEN
          DO i=1,n
             t=y(i)
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE IF (iad .EQ. -1) THEN
          DO i=1,n
             t=cmplx(0,0,kind(0.0d0))
             DO k=ia(i),ia(i+1)-1
                t = t - a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       ELSE
          DO i=1,n
             t=cmplx(0,0,kind(0.0d0))
             DO k=ia(i),ia(i+1)-1
                t = t + a(k)*x(ja(k))
             END DO
             y(i)=t
          END DO
       END IF
       !
    flop=flop+dble(2*(ia(n+1)-ia(1)))
    RETURN
  END SUBROUTINE zagmg_csrmv
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_csrlsolve(n, a, ja, ifja, ia, p, x, y, iunit, flop)
!
!
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    COMPLEX(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
!
    IF (iunit .LT. 0) THEN
       x(1) = y(1)
       DO i=2,n
          t = cmplx(0,0,kind(0.0d0))
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit .GT. 0) THEN
       x(1) = p(1)*y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(1) = y(1)
       DO i=2,n
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE zagmg_csrlsolve
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_csrusolve(n, a, ja, ifja, ia, p, x, y, iunit,flop)
!
!
!
!
!
    IMPLICIT NONE
    INTEGER :: n, ifja, ja(ifja:*), ia(n+1), iunit, i, k
    COMPLEX(kind(0.0d0)) :: a(*), p(n), x(n), y(n), t
    REAL(kind(0.0d0)) :: flop
!
    IF (iunit .LT. 0) THEN
       x(n) = y(n)
       DO i=n-1,1,-1
          t = cmplx(0,0,kind(0.0d0))
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = y(i) + p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE IF (iunit .GT. 0) THEN
       x(n) = p(n)*y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = p(i)*t
       END DO
       flop=flop+dble(n+2*(ia(n+1)-ia(1)))
    ELSE
       x(n) = y(n)
       DO i=n-1,1,-1
          t = y(i)
          DO k=ia(i),ia(i+1)-1
             t = t - a(k)*x(ja(k))
          END DO
          x(i) = t
       END DO
       flop=flop+dble(2*(ia(n+1)-ia(1)))
    END IF
    RETURN
  END SUBROUTINE zagmg_csrusolve
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_FlexCG(N,f,X,ITER,RESID,a,ja,ia,init,flop)
    !
    !
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, init
    REAL(kind(0.0d0)) :: flop
    COMPLEX(kind(0.0d0)) :: f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    INTEGER       :: MAXIT, ierr, kk
    REAL(kind(0.0d0)) :: TOL, BNORM, RESID, dum0, RESID0
    COMPLEX(kind(0.0d0)) :: ALPHA, BET0, RHO, dum(6), td
    COMPLEX(kind(0.0d0)) :: BET1
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: SD(:)
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: S(:),fsc(:)
    COMPLEX(kind(0.0d0)), external :: ZDOTC
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    !
    mritr=nwrkcum+4*N+ N
    ALLOCATE (S(2*N+1:nwrkcum+5*N),SD(N),fsc(N))
    !
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE(iout,940) IRANK
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (smoothtype == 1) THEN
          WRITE(iout,945) nsmooth
       ELSE
          WRITE(iout,946) nsmooth
       END IF
       WRITE(iout,947)
    END IF
    !
    TOL = RESID
    MAXIT = ITER
    RESID = 1.0d0
    ITER = 0
    dum(3) = DZNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    !
    !
    IF (init==1) THEN
       CALL zagmg_matv(N,x,SD,a,ja,ia,flop,trans )
       f(1:n)=f(1:n)-SD(1:N)
       dum(2) = DZNRM2(N, f, IONE)**2
       BNORM=SQRT(dble(dum(3)))
       RESID=SQRT(dble(dum(2)))
       RESID0=RESID
       IF (BNORM.EQ.0.0d0) THEN
          !
          !
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=cmplx(0,0,kind(0.0d0))
          RETURN
       END IF
       RESID=RESID/BNORM
       IF(wff.AND. (MAXIT <= 0 .OR. RESID <= TOL)) THEN
          WRITE(iout,900) 0, resid*bnorm, resid
       END IF
    END IF
    !
    !
    DO WHILE ( ITER < MAXIT .AND. RESID > TOL )
       ITER = ITER + 1
       !
       !
       fsc(1:N)=f(1:N)*scald(1:N)
       CALL zagmg_prec_matv(1,fsc,S(1+3*N),S(1+4*N),flop,S(1+5*N),.FALSE.)
       S(1+3*N:4*N)=S(1+3*N:4*N)*scald(1:N)
       flop=flop+dble(2*N)
       !
       !
       IF ( ITER > 1 ) THEN
          dum(1) = - ZDOTC(N, SD, IONE, S(1+3*N), IONE)
          BET0=dum(1)
          BET1=BET0/RHO
          CALL ZAXPY(N, BET1, S(1+2*N), IONE, S(1+3*N), IONE)
          flop=flop+dble(4*N)
       ENDIF
       CALL ZCOPY(N, S(1+3*N), IONE, S(1+2*N), IONE)
       CALL zagmg_matv(N,S(1+2*N),SD,a,ja,ia,flop,trans )
       !
       dum(1) =  ZDOTC(N, S(1+2*N), IONE, SD, IONE)
       dum(2) =  ZDOTC(N,S(1+2*N),IONE,f,IONE)
       flop=flop+dble(4*N)
       !
       IF (ITER==1) THEN
          IF (init == 0) THEN
             BNORM=SQRT(dble(dum(3)))
             RESID0=BNORM
             IF (BNORM.EQ.0.0d0) THEN
                !
                !
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=cmplx(0,0,kind(0.0d0))
                ITER=0
                RETURN
             END IF
          END IF
          IF(wff) THEN
             WRITE(iout,900) 0, resid*bnorm, resid
          END IF
       ELSE
       END IF
       !
       RHO=dum(1)
       ALPHA=dum(2)/RHO
       !
       IF (ITER == 1 .AND. init == 0) THEN
          CALL ZCOPY(N,S(1+2*N),IONE,X,IONE)
          CALL ZSCAL(N,ALPHA,X,IONE)
          flop=flop+dble(N)
       ELSE
          CALL ZAXPY(N, ALPHA, S(1+2*N), IONE, X, IONE)
          flop=flop+dble(2*N)
       END IF
       CALL ZAXPY(N, -ALPHA, SD, IONE, f, IONE)
       !
       dum0 = DZNRM2(N,f,IONE)**2
       RESID=dum0
       RESID=SQRT(RESID)
       RESID=RESID/BNORM
       IF (wff) THEN
          WRITE(iout,900) iter, resid*bnorm, resid
       END IF
       flop=flop+dble(4*N)
    END DO
    !
    IF( resid > tol ) THEN
       IF (IRANK <= 0) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*BNORM/RESID0
    !
    kstat(2,1)=ABS(iter)
    !
    DEALLOCATE(S,SD)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    RETURN
900 FORMAT('****  Iter=',i5,'        Resid=',e9.2,                 &
         '        Relat. res.=',e9.2)
940 FORMAT(i3,                                                     &
         '*SOLUTION: flexible conjugate gradient iterations (FCG(1))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
945 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre-and post-smoothing step(s)')
946 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' ILU(0) pre-and post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING!!!',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
    !
  END SUBROUTINE zagmg_FlexCG
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_GCR(N,f,X,ITER,RESID,a,ja,ia,init,flop)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER       :: N, ITER, init
    REAL(kind(0.0d0)) ::  flop
    COMPLEX(kind(0.0d0)) ::  f(N), X(N)
    INTEGER       :: ja(*), ia(N+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    INTEGER   :: MAXIT,i,itm,irst,ierr,j,k
    COMPLEX(kind(0.0d0)) ::  ALPHA, BET0
    REAL(kind(0.0d0)) ::  RESID, RESID2, BNORM2,  RHO, TRS, RESID0
    REAL(kind(0.0d0)) ::  TOL, TOL2BNORM2, TOLT, dum0
    COMPLEX(kind(0.0d0)) :: dum(6),xd,t
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: Sc(:),Siv(:),SiR(:)
    COMPLEX(kind(0.0d0)), ALLOCATABLE :: Su(:),R(:),fsc(:),Sw(:)
    COMPLEX(kind(0.0d0)), external :: ZDOTC
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    INTEGER  :: itm1, m, info
    !
    mritr=nwrkcum+N*2*nrst+((nrst+1)*nrst)/2+nrst+N
    ALLOCATE (Su(N*nrst),Sc(N*nrst)                  &
             ,SiR(((nrst+1)*nrst)/2),Siv(nrst)       &
             ,R(nwrkcum),fsc(N))
    !
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE(iout,938) IRANK,nrst
    END IF
    IF (wff) THEN
       IF (  trans) THEN
          WRITE(iout,941)
       END IF
       IF (smoothtype == 1) THEN
          WRITE(iout,945) nsmooth
       ELSE
          WRITE(iout,946) nsmooth
       END IF
       WRITE(iout,947)
    END IF
    !
    m=MIN(nrst,ITER)
    TRS=EPSILON(1.0d0)
    TRS=SQRT(SQRT(TRS))
    TOL = RESID
    TOL2BNORM2 = TOL
    MAXIT = ITER
    RESID2= 1.0d0
    ITER = 0
    dum(3) = DZNRM2(N, f, IONE)**2
    flop=flop+dble(2*N)
    !
    !
    IF (init==1) THEN
       CALL zagmg_matv(N,x,Sc,a,ja,ia,flop,trans )
       f(1:n)=f(1:n)-Sc(1:N)
       dum(2) = DZNRM2(N, f, IONE)**2
       BNORM2=dble(dum(3))
       RESID2=dble(dum(2))
       IF (BNORM2.EQ.0.0d0) THEN
          !
          !
          IF(wff) THEN
             WRITE(iout,998)
             WRITE(iout,999)
          END IF
          X(1:N)=cmplx(0,0,kind(0.0d0))
          RETURN
       END IF
       TOL2BNORM2 = TOL*TOL*BNORM2
       IF (wff.AND. (MAXIT <= 0 .OR. RESID2 <= TOL2BNORM2)) THEN
          WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
       END IF
       RESID0=SQRT(RESID2)
    END IF
    !
    itm  = -1
    irst = 0
    DO WHILE ( ITER < MAXIT .AND. RESID2 > TOL2BNORM2 )
       itm  = itm  + 1
       ITER = ITER + 1
       !
       IF (itm == m) THEN
          CALL ZTPTRS('U','N','U',m,IONE,SiR,Siv,m,info)
          IF (irst == 0 .AND. init == 0) THEN
             CALL ZGEMV('N',N,m,cmplx(1,0,kind(0.0d0)),Su,       &
                  N,Siv,IONE,cmplx(0,0,kind(0.0d0)),X,IONE)
             flop=flop+dble(2*m*N+m*(m+1))
          ELSE
             CALL ZGEMV('N',N,m,cmplx(1,0,kind(0.0d0)),Su,        &
                  N,Siv,IONE,cmplx(1,0,kind(0.0d0)),X,IONE)
             flop=flop+dble((2*m+1)*N+m*(m+1))
          END IF
          itm=0
          irst=irst+1
       END IF
       !
       !
       fsc(1:N)=f(1:N)*scald(1:N)
       CALL zagmg_prec_matv(1,fsc,Su(1+itm*N),Sc(1+itm*N)      &
            ,flop,R,.FALSE.)
       Su(1+itm*N:(1+itm)*N)=Su(1+itm*N:(1+itm)*N)*scald(1:N)
       flop=flop+dble(2*N)
       CALL zagmg_matv(N, Su(1+itm*N), Sc(1+itm*N)             &
            , a, ja, ia, flop,trans )
       !
       !
       IF (itm > 0) THEN
          DO i=0,itm-1
             dum(1)=ZDOTC(N,Sc(1+i*N),IONE,Sc(1+itm*N),IONE)
             bet0=dum(1)
             bet0=bet0/SiR(1+i+(i*(i+1))/2)
             SiR(1+i+(itm*(itm+1))/2)=bet0
             CALL ZAXPY(N,-bet0,Sc(1+i*N),IONE,Sc(1+itm*N),IONE)
             flop=flop+dble(4*N)
          END DO
       END IF
       !
       !
       dum(1)=DZNRM2(N,Sc(1+itm*N),IONE)**2
       dum(2)=ZDOTC(N, Sc(1+itm*N), IONE, f, IONE)
       IF (ITER == 1) THEN
          IF (init == 0) THEN
             BNORM2=dble(dum(3))
             RESID2=BNORM2
             RESID0=SQRT(BNORM2)
             IF (BNORM2.EQ.0.0d0) THEN
                !
                !
                IF(wff) THEN
                   WRITE(iout,998)
                   WRITE(iout,999)
                END IF
                X(1:N)=cmplx(0,0,kind(0.0d0))
                RETURN
             END IF
             TOL2BNORM2=TOL*TOL*BNORM2
          END IF
          IF (wff) THEN
             WRITE(iout,900) 0, 0,SQRT(resid2),SQRT(resid2/bnorm2)
          END IF
          TOLT = MAX(TOL2BNORM2,TRS*RESID2)
       ELSE
       END IF
       rho=dble(dum(1))
       alpha=dum(2)
       SiR(1+itm+(itm*(itm+1))/2)=rho
       bet0=alpha/rho
       Siv(1+itm)=bet0
       !
       CALL ZAXPY(N, -bet0, Sc(1+itm*N), IONE, f, IONE)
       flop=flop+dble(6*N)
       !
       RESID2 = RESID2 - conjg(alpha)*alpha/rho
       IF (RESID2 <= TOLT) THEN
          RESID2 = DZNRM2(N,f,IONE)**2
          flop=flop+dble(2*N)
          TOLT = MAX(TOL2BNORM2,TRS*RESID2)
       END IF
       IF (wff)THEN
          WRITE(iout,900) iter,irst,SQRT(ABS(resid2)),SQRT(ABS(resid2/bnorm2))
       END IF
       !
    END DO
    !
    IF (itm >= 0) THEN
       itm1=itm+1
       CALL ZTPTRS('U','N','U',itm1, IONE,SiR,Siv,m,info)
       IF (irst == 0 .AND. init == 0) THEN
          CALL ZGEMV('N',N,itm1,cmplx(1,0,kind(0.0d0)),Su,        &
               N,Siv,IONE,cmplx(0,0,kind(0.0d0)),X,IONE)
          flop=flop+dble(2*(itm+1)*N+(itm+1)*(itm+2))
       ELSE
          CALL ZGEMV('N',N,itm1,cmplx(1,0,kind(0.0d0)),Su,        &
               N,Siv,IONE,cmplx(1,0,kind(0.0d0)),X,IONE)
          flop=flop+dble((2*(itm+1)+1)*N+(itm+1)*(itm+2))
       END IF
    END IF
    !
    RESID=SQRT(RESID2/BNORM2)
    IF( resid > tol ) THEN
       IF (IRANK <= 0) THEN
          WRITE(iout,'()')
          WRITE(iout,950) iter
          WRITE(iout,951)
          WRITE(iout,'()')
       END IF
       iter=-iter
    ELSE
       IF (wff) THEN
          WRITE(iout,952) iter
          WRITE(iout,'()')
       END IF
    END IF
    RESID=RESID*SQRT(BNORM2)/RESID0
    !
    DEALLOCATE (Su,Sc,R,Siv,SiR)
    IF(ALLOCATED(fsc)) DEALLOCATE(fsc)
    IF(ALLOCATED(Sw)) DEALLOCATE(Sw)
    !
    kstat(2,1)=ABS(iter)
    !
    RETURN
900 FORMAT('****  Iter=',i5,' (',i2,' rest.)        Resid=',e9.2,    &
         '        Relat. res.=',e9.2)
938 FORMAT(i3,'*SOLUTION: GCR iterations (GCR(',i2,'))')
941 FORMAT('****     Rmk: solve system with the transpose of the input matrix')
945 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' Gauss-Seidel pre-and post-smoothing step(s)')
946 FORMAT(  '****     AMG preconditioner with',i2,             &
         ' ILU(0) pre-and post-smoothing step(s)')
947 FORMAT(  '****         at each level')
950 FORMAT('**** !!!   WARNING!!!',I5,' ITERATIONS WERE')
951 FORMAT('**** INSUFFICIENT TO ACHIEVE THE DESIRED ACCURACY')
952 FORMAT('****  - Convergence reached in',I5,' iterations -')
998 FORMAT('**** The norm of the right hand side is zero:')
999 FORMAT('****     set x equal to the zero vector and exit')
    !
  END SUBROUTINE zagmg_GCR
!--------------------------------------------------------------------
  SUBROUTINE zagmg_matv(n, x, y, a, ja, ia, flop, transpose,         &
       iext, lstout, ilstout, lstin, ilstin)
!
!
    IMPLICIT NONE
    INTEGER :: n, ja(*), ia(n+1), i, kk, k1, k2, ier, idum(1)
    INTEGER, OPTIONAL :: iext(n), lstout(*), ilstout(*), lstin(0:*), ilstin(*)
    COMPLEX(kind(0.0d0)) :: y(n), a(*), t, xx
    COMPLEX(kind(0.0d0)) :: x(n)
    REAL(kind(0.0d0)) :: flop
    LOGICAL :: transpose
    !
    !
    IF (.NOT.transpose) THEN
       DO i = 1, n
          k1 = ia(i)
          xx = x(ja(k1))
          t = a(k1) * xx
          k2 = ia(i+1)-1
          DO kk = k1+1, k2
             xx = x(ja(kk))
             t = t + a(kk)*xx
          ENDDO
          y(i) = t
       ENDDO
    ELSE
       y(1:n)=cmplx(0,0,kind(0.0d0))
       DO i = 1, n
          xx = x(i)
          DO kk = ia(i), ia(i+1)-1
             y(ja(kk)) = y(ja(kk)) + a(kk)*xx
          ENDDO
       ENDDO
    END IF
    !
    flop = flop + dble(2 *(ia(n+1)-1)-n)
    !
    RETURN
  END SUBROUTINE zagmg_matv
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE zagmg_prec_matv(l,B,X,AX,flop,R,matv)
!
!
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l
    REAL(kind(0.0d0)), OPTIONAL :: flop
    COMPLEX(kind(0.0d0)), OPTIONAL ::  B(*), X(*), AX(*), R(*)
    COMPLEX(kind(0.0d0)) ::  dum(1)
    LOGICAL, OPTIONAL :: matv
    LOGICAL :: update
    INTEGER :: is,n,nnext,XN,BN,AXN,RN,idum(1)
!
!
    INTEGER, PARAMETER :: nsmotot=nsmooth*smoothtype
    n=nn(l)
    nnext=nn(l+1)
    !
    IF (l == nlev) THEN
       !
       !
       X(1:N)=B(1:N)
       CALL zagmg_MUMPSseq(n,X,2,flop)
       !
       RETURN
       !
    END IF
    !
    !
    !
    IF (nsmotot == 1) THEN
       !
       CALL zagmg_fwGS(l,B,X,AX,flop,.FALSE.,R)
       !
    ELSE
       !
       update=.FALSE.
       DO is=2,nsmotot,2
          CALL zagmg_fwbwsolve1(l,dt(l)%a,B,X,AX,flop,update,R)
          update=.TRUE.
       END DO
       IF (mod(nsmotot,2) .EQ. 1) THEN
          !
          CALL zagmg_fwGS(l,B,X,AX,flop,.TRUE.,R)
          !
       END IF
    END IF
    !
    !
    IF (nnext > 0) THEN
       !
       !
       XN=1
       BN=XN+nnext
       IF (l+1 == nlev) BN=XN
       CALL zagmg_restag(N,nnext,AX,R(BN),dt(l)%ind,flop)
       !
       IF (l+1 == nlev) THEN
          !
          CALL zagmg_MUMPSseq(nnext,R,2,flop)
          !
       ELSE IF ( innermax(l+1) <= 1 ) THEN
          !
          AXN=BN+nnext
          RN=AXN+nnext
          CALL zagmg_prec_matv(l+1,R(BN),R(XN),R(AXN),flop,R(RN),.FALSE.)
          !
          kstat(1,l+1)=MAX(kstat(1,l+1),1)
          kstat(2,l+1)=kstat(2,l+1)+1
          !
       ELSE
          !
          CALL zagmg_inner_iter(nnext,R(XN),R(BN),l+1,flop)
          !
       END IF
       !
       CALL zagmg_prolag(N,nnext,X,R(XN),dt(l)%ind,flop)
       !
    END IF
    !
    !
    IF (nsmotot == 1) THEN
       !
       CALL zagmg_bwGS(l,B,X,AX,flop,matv)
       !
    ELSE
       IF (mod(nsmotot,2) .EQ. 1) THEN
          !
          CALL zagmg_bwGS(l,B,X,AX,flop,.FALSE.)
          !
       END IF
       !
       update=.FALSE.
       DO is=2,nsmotot,2
          IF (is .GE. nsmotot-1) update=matv
          CALL zagmg_fwbwsolve2(l,dt(l)%a,B,X,AX,flop,update,R)
       END DO
    END IF
    RETURN
  END SUBROUTINE zagmg_prec_matv
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_fwGS(l,B,X,AX,flop,update,R)
!
!
!
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) :: flop
    COMPLEX(kind(0.0d0)) ::  B(*), X(*), AX(*), R(*)
    LOGICAL :: update
    n=nn(l)
    !
    IF (update) THEN
       !
       !
       CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a   &
                            ,R,AX,1,flop)
       X(1:n)=X(1:n)+R(1:n)
       flop=flop+dble(n)
       !
       !
       !
       CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,R,1,AX,-1,flop)
       !
    ELSE
       !
       !
       CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,dt(l)%a,X,B,1,flop)
       !
       CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,-1,flop)
       !
    END IF
    RETURN
  END SUBROUTINE zagmg_fwGS
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_bwGS(l,B,X,AX,flop,matv)
!
!
!
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    REAL(kind(0.0d0)) :: flop
    COMPLEX(kind(0.0d0)) ::  B(*), X(*), AX(*)
    LOGICAL :: matv
    n=nn(l)
    !
    !
    AX(1:n)=B(1:n)
    CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
    !
    !
    CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,dt(l)%a,X,AX,1,flop)
    !
    IF (.NOT.matv) RETURN
    !
    !
    CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    RETURN
  END SUBROUTINE zagmg_bwGS
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_fwbwsolve1(l,p,B,X,AX,flop,update,R)
!
!
!
!
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    COMPLEX(kind(0.0d0)) ::  p(*), B(*), X(*), AX(*), R(*)
    REAL(kind(0.0d0)) :: flop
    LOGICAL :: update
    !
    n=nn(l)
    IF (.NOT.update) THEN
       !
       !
       CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,p    &
                            ,R,B,-1,flop)
       !
       !
       CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,p    &
                             ,X,R,1,flop)
       !
       !
       !
       AX(1:n)=B(1:n)-R(1:n)
       flop=flop+dble(n)
       !
       !
       CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,-2,flop)
       !
       IF (smoothtype == 2) AX(1:n)=AX(1:n)-dt(l)%a(1:n)*X(1:n)
    ELSE
       !
       !
       CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,p    &
                            ,R,AX,-1,flop)
       !
       !
       CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,p    &
                            ,R(n+1),R,1,flop)
       !
       !
       X(1:n)=X(1:n)+R(n+1:2*n)
       !
       !
       !
       AX(1:n)=AX(1:n)-R(1:n)
       flop=flop+dble(2*n)
       !
       !
       CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,R(n+1),1,AX  &
                         ,-2,flop)
       !
       IF (smoothtype == 2) AX(1:n)=AX(1:n)-dt(l)%a(1:n)*R(n+1:2*n)
       IF (smoothtype == 2) flop=flop+dble(2*n)
    END IF
    !
    RETURN
  END SUBROUTINE zagmg_fwbwsolve1
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_fwbwsolve2(l,p,B,X,AX,flop,matv,R)
!
!
!
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l, n, idum(1), ier
    COMPLEX(kind(0.0d0)) ::  p(*), B(*), X(*), AX(*), R(*)
    REAL(kind(0.0d0)) :: flop
    LOGICAL :: matv
    !
    n=nn(l)
    !
    !
    !
    CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,X,1,AX,0,flop)
    !
    R(1:n)=B(1:n)-AX(1:n)
    IF (smoothtype == 2) R(1:n)=R(1:n)-dt(l)%a(1:n)*X(1:n)
    IF (smoothtype == 2) flop=flop+dble(2*n)
    !
    !
    !
    CALL zagmg_csrlsolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,p    &
                         ,R(n+1),R,1,flop)
    !
    !
    R(1:n) = R(n+1:2*n) - X(1:n)
    CALL zagmg_csrusolve(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%iu,p    &
                         ,R,R,-1,flop)
    !
    !
    X(1:n)=X(1:n) + R(1:n)
    flop=flop+dble(3*n)
    !
    IF (.NOT.matv) RETURN
    !
    AX(1:n)=AX(1:n)+R(n+1:2*n)/P(1:n)
    flop=flop+dble(2*n)
    !
    !
    CALL  zagmg_csrmv(n,dt(l)%a,dt(l)%ja,n+1,dt(l)%il,X,1,AX,2,flop)
    IF (smoothtype == 2) AX(1:n)=AX(1:n)+dt(l)%a(1:n)*R(n+1:2*n)
    IF (smoothtype == 2) flop=flop+dble(2*n)
    RETURN
  END SUBROUTINE zagmg_fwbwsolve2
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE zagmg_inner_iter(n,X,R,l,flop)
    !
    !
    !
    !
    !
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER   :: N, ITER, l, ierr, k
    REAL(kind(0.0d0)) ::  RESID, flop, BNORM, det
    COMPLEX(kind(0.0d0)) :: X(N), R(N,*)
    COMPLEX(kind(0.0d0)) :: alpha1,alpha2,bet0,bet1,bet2,rho1,rho2,gamm0,gamm1,zeta
    COMPLEX(kind(0.0d0)) :: dum(10),y1,y2
    COMPLEX(kind(0.0d0)), external :: ZDOTC
    REAL(kind(0.0d0)), external :: DZNRM2
    INTEGER , parameter :: IONE=1
    !
    !
    dum(3)=DZNRM2(N, R, IONE)**2
    !
    !
    IF (dble(dum(3)) .EQ. 0.0d0) THEN
       X(1:N)=cmplx(0,0,kind(0.0d0))
       flop=flop+dble(2*N)
       RETURN
    END IF
    ITER = 1
    !
    !
    CALL zagmg_prec_matv(l,R,X,R(1,2),flop,R(1,3),.TRUE.)
    !
    !
    dum(1) = ZDOTC(N,X,IONE,R(1,2),IONE)
    dum(2) = ZDOTC(N,X,IONE,R,IONE)
    IF (ABS(dum(1)) .LE. repsmach*ABS(dum(2))) THEN
       flop=flop+dble(2*N)
       GOTO 100
    END IF
    BNORM = dble(dum(3))
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    !
    CALL ZAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    RESID = DZNRM2(N,R,IONE)**2
    IF (RESID <= resi*resi*BNORM) THEN
       CALL ZSCAL( N, bet0, X, IONE )
       flop=flop+dble(11*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    ITER = 2
    !
    !
    CALL zagmg_prec_matv(l,R,R(1,3),R(1,4),flop,R(1,5),.TRUE.)
    !
    !
    !
    dum(1) = ZDOTC(N,R(1,3),IONE,R(1,2),IONE)
    dum(2) = ZDOTC(N,R(1,3),IONE,R,IONE)
    dum(3) = ZDOTC(N,R(1,3),IONE,R(1,4),IONE)
    IF (.NOT.spd) dum(4) = ZDOTC(N,X,IONE,R(1,4),IONE)
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    IF (spd) THEN
       gamm1 = conjg(gamm0)
    ELSE
       gamm1 = dum(4)
    END IF
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    IF (ABS(bet1).LE.repsmach*ABS(alpha2) .OR.   &
        ABS(bet2)*repsmach.GE.1.0d0)            THEN
       flop=flop+dble(6*N)
       IF (.NOT.spd) flop=flop+dble(2*N)
      GOTO 200
    END IF
    CALL ZSCAL(N, bet2, X, IONE)
    CALL ZAXPY(N, zeta, R(1,3), IONE, X, IONE)
    !
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    IF (.NOT.spd) flop=flop+dble(2*N)
    !
    RETURN
    !
100 CONTINUE
    dum(1)=DZNRM2(N,R(1,2),IONE)**2
    dum(2)=ZDOTC(N, R(1,2), IONE, R, IONE )
    BNORM = dble(dum(3))
110 CONTINUE
    rho1=dum(1)
    alpha1=dum(2)
    !
    bet0=alpha1/rho1
    !
    RESID=BNORM-dble(conjg(alpha1)*bet0)
    IF (RESID <= resi*resi*BNORM) THEN
       CALL ZSCAL( N, bet0, X, IONE )
       flop=flop+dble(7*N)
       kstat(1,l)=MAX(kstat(1,l),iter)
       kstat(2,l)=kstat(2,l)+iter
       RETURN
    END IF
    !
    CALL ZAXPY(N, -bet0, R(1,2), IONE, R, IONE)
    !
    ITER = 2
    !
    !
    CALL zagmg_prec_matv(l,R,R(1,3),R(1,4),flop,R(1,5),.TRUE.)
    !
    !
    dum(1) = ZDOTC(N,R(1,4),IONE,R(1,2),IONE)
    dum(2) = ZDOTC(N,R(1,4),IONE,R,IONE)
    dum(3) = DZNRM2(N,R(1,4),IONE)**2
    gamm0 =dum(1)
    alpha2=dum(2)
    rho2  =dum(3)
    gamm1 = conjg(gamm0)
    !
    bet1=rho2-gamm0*gamm1/rho1
    bet2=(alpha1-alpha2*gamm1/bet1)/rho1
    zeta=alpha2/bet1
    !
    CALL ZSCAL(N, bet2, X, IONE)
    CALL ZAXPY(N, zeta, R(1,3), IONE, X, IONE)
    !
    !
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(19*N)
    !
    RETURN
    !
200 CONTINUE
    !
    !
    !
    !
    !
    !
    !
    !
    dum(1) = DZNRM2(N,R(1,2),IONE)**2
    dum(2) = ZDOTC(N,R(1,4),IONE,R(1,2),IONE)
    dum(3) = DZNRM2(N,R(1,4),IONE)**2
    dum(4) = ZDOTC(N,R(1,2),IONE,R,IONE)
    dum(5) = ZDOTC(N,R(1,4),IONE,R,IONE)
     !
     dum(4) = dum(4)+bet0*dum(1)
     dum(5) = dum(5)+bet0*dum(2)
     det = dble(dum(1)*dum(3)-dum(2)*conjg(dum(2)))
     y1  = (dum(3)*dum(4)-conjg(dum(2))*dum(5))/det
     y2  = (-dum(2)*dum(4)+dum(1)*dum(5))/det
     CALL ZSCAL( N, y1, X, IONE )
     CALL ZAXPY( N, y2, R(1,3), IONE, X, IONE )
    kstat(1,l)=MAX(kstat(1,l),iter)
    kstat(2,l)=kstat(2,l)+iter
    !
    flop=flop+dble(23*N)
    !
    RETURN
    !
  END SUBROUTINE zagmg_inner_iter
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_prolag(n, nc, V, B, ind, flop)
!
!
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    COMPLEX(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    !
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) THEN
          V(i) = V(i)+B(k)
       ENDIF
    ENDDO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE zagmg_prolag
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_restag(n, nc, V, B, ind, flop)
!
!
    IMPLICIT NONE
    INTEGER :: n, nc, ind(n), k, i
    COMPLEX(kind(0.0d0)) :: V(n), B(nc)
    REAL(kind(0.0d0)) :: flop
    !
    B(1:nc)=cmplx(0,0,kind(0.0d0))
    !
    DO i = 1, n
       k = ind(i)
       IF (k.GT.0) B(k) = B(k) + V(i)
    ENDDO
    flop = flop + dble(n)
    RETURN
  END SUBROUTINE zagmg_restag
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE zagmg_setupL1(n,a,ja,ia,listrank,ifl)
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n
    INTEGER :: ja(*),ia(n+1)
    COMPLEX(kind(0.0d0)) :: a(*)
    INTEGER :: ifl,listrank(ifl:*)
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)) :: ops,eta,dum(2)
    CHARACTER(len=13) :: prtint
    COMPLEX (kind(0.0d0)) :: fff(1)
    !
    nn(1)=n
    nlc(1)=n
    nlc(2)=ia(n+1)-ia(1)
    !
    ngl=nlc
       !
       wcplex=1.0d0
       nlctot=nlc
       ngltot=ngl
       ngl1=ngl
       nlc1=nlc
       icum=1
       fracnz(1)=1.0d0
       allzero=.FALSE.
       nglp=0.0d0
       maxcoarset=maxcoarsesize
       maxcoarseslowt=maxcoarsesizeslow
       maxcoarset=FLOOR(maxcoarset*(ngl(1)**(1.0d0/3)))
       maxcoarseslowt=FLOOR(maxcoarseslowt*(ngl(1)**(1.0d0/3)))
       !
       IF (wfo) THEN
          IF (n > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,'()')
          WRITE(iout,918) prtint(1:12)
          IF (nlc(2) > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,919) prtint(1:12),dble(nlc(2))/dble(n)
          WRITE(iout,'()')
       END IF
    IF ( 0 == nstep  .OR. 1 == maxlev .OR. ngl(1) <= maxcoarset ) nlev=1
    nlcp=nlc
    nglp=ngl
    !
    IF (1 /= nlev) THEN
       !
       !
       CALL zagmg_aggregation(1,n,a,ja,ia,nc)
       !
       !
       CALL zagmg_setup(2,nc)
       !
    ELSE
       !
       !
       DEALLOCATE(dt(1)%idiag)
       memi=memi-(n+1)
       CALL zagmg_MUMPSseq(n,fff,1,ops,a,ja,ia)
       IF (wfo) THEN
          WRITE(iout,911) IRANK,ops/dble(2*nlc1(2))
       END IF
       wcplex(4)=0.0d0
       IF (ngl(1) > maxcoarset)   wcplex(4)=-1.0d0
       IF (ngl(1) > maxcoarseslowt) wcplex(4)=ngl(1)/(1000*ngl1(1)**(1.0d0/3))
       !
       !
       RETURN
    END IF
    !
    !
       !
       nz=ia(n+1)-ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%il(n+1),dt(1)%iu(n+1))
          memi=memi+nz+n+2
          memr=memr+nz
       ELSE
          ALLOCATE(ap(nz),jap(nz-n),dt(1)%iu(n+1))
          dt(1)%il => dt(1)%idiag
          memi=memi+nz+1
          memr=memr+nz
       END IF
       memax=MAX(memax,memr+memi*rlenilen)
       !
       CALL zagmg_csrdlu(n,a,ja,ia,dt(1)%idiag              &
            ,ap,jap,dt(1)%il,dt(1)%iu,transint )
       IF (transint) THEN
          DEALLOCATE(dt(1)%idiag)
          memi=memi-(n+1)
       ELSE
          NULLIFY(dt(1)%idiag)
       END IF
       dt(1)%a  => ap
       dt(1)%ja => jap
       NULLIFY(ap,jap)
       !
       !
       innermax(nlev)=0
       innermax(1)=1
       eta=xsi/((1-xsi)*(cplxmax-1))
       icum=1
       DO i=2,nlev-1
          innermax(i)=min(2,floor(xsi**(i-1)/(eta*fracnz(i)*icum)))
          IF (nlev-i.LE.nlvcyc .AND. i.GT.2) innermax(i)=1
          icum=icum*innermax(i)
          wcplex(2)=wcplex(2)+icum*fracnz(i)
          wcplex(3)=wcplex(3)+(2**(i-1))*fracnz(i)
       END DO
       wcplex(2)=wcplex(2)+(2**(nlev-1))*fracnz(nlev)
       wcplex(3)=wcplex(3)+(2**(nlev-1))*fracnz(nlev)
       IF (nsmooth*smoothtype > 1)  THEN
          nwrkcum=2*nn(nlev-1)
       ELSE
          nwrkcum=nn(nlev)
       END IF
       nwrkcum=max(nwrkcum,1)
       DO i=nlev-2,1,-1
          nwrkcum=nwrkcum+3*nn(i+1)
          IF (innermax(i+1) > 1) nwrkcum=nwrkcum+2*nn(i+1)
          IF (nsmooth*smoothtype > 1)  nwrkcum=max(2*nn(i),nwrkcum)
       END DO
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,954) nlctot(1)/dble(nlc1(1))
          WRITE(iout,955) nlctot(2)/dble(nlc1(2))
       END IF
       IF (wff) THEN
          WRITE(iout,956) wcplex(3)
          WRITE(iout,957) wcplex(2)
          WRITE(iout,'()')
       END IF
    RETURN
911 FORMAT(i3,'*','        Exact factorization:',f12.3,' work units (*)')
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
954 FORMAT('****','                  Grid complexity:',f9.2)
955 FORMAT('****','              Operator complexity:',f9.2)
956 FORMAT('****','  Theoretical Weighted complexity:',f9.2, &
                ' (K-cycle at each level)' )
957 FORMAT('****','    Effective Weighted complexity:',f9.2, &
                ' (V-cycle enforced where needed)' )
  END SUBROUTINE zagmg_setupL1
!-----------------------------------------------------------------------
  RECURSIVE SUBROUTINE zagmg_setup(l,n,listrank)
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n
    INTEGER, OPTIONAL :: listrank(n+1:*)
    INTEGER :: nc,ierr,i,j,k,nz
    LOGICAL :: slcoarse
    INTEGER, POINTER, DIMENSION(:) :: jap
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    LOGICAL, SAVE :: slowcoarse
    REAL(kind(0.0d0)) :: ops,fw,eta,dum(2)
    CHARACTER(len=13) :: prtint
    COMPLEX (kind(0.0d0)) :: fff(1)
    !
    nn(l)=n
    nlc(1)=n
    IF (n > 0) THEN
       nlc(2)=dt(l)%ia(n+1)-dt(l)%ia(1)
    ELSE
       nlc(2)=0
    END IF
    !
    ngl=nlc
       !
       nlctot=nlctot+nlc
       ngltot=ngltot+ngl
       IF (wfo) THEN
          WRITE(iout,'()')
          WRITE(iout,914) IRANK,l
          IF (n > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(n)
          ELSE
             WRITE(prtint(1:12),'(i12)') n
          END IF
          WRITE(iout,920) prtint(1:12),dble(nlcp(1))/dble(n)
          IF (nlc(2) > 9.9e10) THEN
             WRITE(prtint(1:12),'(1pe12.5)') dble(nlc(2))
          ELSE
             WRITE(prtint(1:12),'(i12)') nlc(2)
          END IF
          WRITE(iout,921) prtint(1:12),dble(nlc(2))/dble(n),        &
               dble(nlcp(2))/dble(nlc(2))
       END IF
       fracnz(l)=ngl(2)/ngl1(2)
    IF (l==2) THEN
       wcplex(1)=ngl1(2)/ngl(2)
       slowcoarse=.FALSE.
    END IF
       !
    slcoarse = 2*nglp(1) < 3*ngl(1) .AND. 2*nglp(2) < 3*ngl(2)
    IF( l == nstep+1  .OR. l == maxlev                        &
         .OR. (  ngl(1) <= maxcoarset)                        &
         .OR. ( nglp(1) < 2*ngl(1) .AND. nglp(2) < 2*ngl(2)   &
                            .AND. ngl(1) <= maxcoarseslowt )  &
         .OR. ( slowcoarse .AND. slcoarse )                   &
         .OR. nglp(1) ==ngl(1) )                       THEN
       nlev=l
    END IF
    slowcoarse=slcoarse
    nlcp=nlc
    nglp=ngl
    IF (n == 0) THEN
       nc=0
       allzero=.TRUE.
       !
       !
       RETURN
    END IF
    IF (l /= nlev) THEN
       !
       !
       CALL zagmg_aggregation(l,n,dt(l)%a,dt(l)%ja,dt(l)%ia,nc)
       !
       !
       CALL zagmg_setup(l+1,nc)
       !
    ELSE
       !
       !
       DEALLOCATE(dt(l)%idiag)
       memi=memi-(n+1)
       CALL zagmg_MUMPSseq(n,fff,1,ops,dt(l)%a,dt(l)%ja,dt(l)%ia)
       !
       memi=memi-size(dt(l)%ja)-n-1
       memr=memr-size(dt(l)%a)
       DEALLOCATE(dt(l)%a,dt(l)%ja,dt(l)%ia)
       !
       IF (wfo) THEN
          WRITE(iout,911) IRANK,ops/dble(2*nlc1(2))
       END IF
       wcplex(4)=0.0d0
       IF (ngl(1) > maxcoarset)   wcplex(4)=-1.0d0
       IF (ngl(1) > maxcoarseslowt) wcplex(4)=ngl(1)/(1000*ngl1(1)**(1.0d0/3))
       !
       !
       RETURN
    END IF
    !
    !
       !
       nz=dt(l)%ia(n+1)-dt(l)%ia(1)
       IF (transint) THEN
          ALLOCATE(ap(nz),jap(nz-n),dt(l)%il(n+1),dt(l)%iu(n+1))
          memi=memi+nz+n+2
          memr=memr+nz
       ELSE
          ALLOCATE(ap(nz),jap(nz-n))
          dt(l)%iu => dt(l)%ia
          dt(l)%il => dt(l)%idiag
          memi=memi+nz-n
          memr=memr+nz
       END IF
       memax=MAX(memax,memr+memi*rlenilen)
       !
       CALL zagmg_csrdlu(n,dt(l)%a,dt(l)%ja,dt(l)%ia,dt(l)%idiag &
            ,ap,jap,dt(l)%il,dt(l)%iu,transint )
       memi=memi-size(dt(l)%ja)
       memr=memr-size(dt(l)%a)
       DEALLOCATE(dt(l)%a,dt(l)%ja)
       IF (transint) THEN
          DEALLOCATE(dt(l)%idiag,dt(l)%ia)
          memi=memi-2*(n+1)
       ELSE
          NULLIFY(dt(l)%idiag,dt(l)%ia)
       END IF
       dt(l)%a  => ap
       dt(l)%ja => jap
       NULLIFY(ap,jap)
    !
    RETURN
911 FORMAT(i3,'*','        Exact factorization:',f12.3,' work units (*)')
914 FORMAT(i3,'*','                      Level:',I12)
918 FORMAT('****','        Number of unknowns:', A12)
919 FORMAT('****','                 Nonzeros :', A12,                &
         ' (per row:',f7.2,')')
920 FORMAT('****','        Number of variables:',A12,              &
         '          (reduction ratio:',f5.2,')')
921 FORMAT('****','                   Nonzeros:',A12,              &
         ' (per row:',f4.1,    &
         '; red. ratio:',f5.2,')')
  END SUBROUTINE zagmg_setup
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_smoothsetup
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,i
    !
    !
    IF (smoothtype == 1) THEN
       DO l=1,nlev-1
          !
          DO i=1,nn(l)
             dt(l)%a(i)=1.0d0/dt(l)%a(i)
          END DO
       END DO
    END IF
    RETURN
  END SUBROUTINE zagmg_smoothsetup
!------------------------------------------------------------------
  SUBROUTINE zagmg_aggregation(l,n,a,ja,ia,nc,listrank)
!
!
!
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: l,n,nc
    INTEGER :: ja(*),ia(n+1)
    INTEGER, OPTIONAL :: listrank(*)
    COMPLEX (kind(0.0d0)) :: a(*), dum
    INTEGER :: ier,i,j,k,maxdg,np,kpass,nzc,m1,ndd,nzp,isize,nddp,npass1,nz,i0
    LOGICAL :: skipass
    !
    INTEGER, POINTER, DIMENSION(:) :: jan,ian,idiagn,iextn,ind2,lcg,lcgn
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: an
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sinn
    !
    INTEGER, POINTER, DIMENSION(:) :: jap,iap,idiagp,iextp,lcg1
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ap
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: sip
    !
    INTEGER, ALLOCATABLE, DIMENSION(:) :: ldd,iw,iperm,riperm
    REAL(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: si1,w
    COMPLEX(kind(0.0d0)), ALLOCATABLE, DIMENSION(:) :: wc
    !
    !
    IF (l .EQ. 1) THEN
       IF (wfo) THEN
          WRITE (iout,901) IRANK
       END IF
       IF (wff) THEN
          IF (spd) THEN
             WRITE (iout,906)
          ELSE IF (  transint) THEN
             WRITE (iout,908)
          END IF
          IF (.not.spd) then
             WRITE (iout,902) 'Jacobi',kaptg_dampJac,checkddJ
          ELSE
             WRITE (iout,902) 'BlockD',kaptg_blocdia,checkddB
          END IF
          IF (checkdd < 0) THEN
             WRITE (iout,904)
          END IF
          WRITE(iout,903) npass,targetcoarsefac
          WRITE (iout,905) trspos
       END IF
    END IF
    !
    IF (l .EQ. 1) THEN
       ALLOCATE(si1(n),ind2(n),iperm(n),riperm(n))
       memi=memi+4*n+1
       memr=memr+n
       CALL zagmg_setCMK(n,ja,ia,dt(l)%idiag,riperm,iperm)
    ELSE
       ALLOCATE(si1(n),ind2(n),iperm(n))
       memi=memi+2*n
       memr=memr+n
       iperm(1:n)=1
    END IF
    memax=MAX(memax,memr+memi*rlenilen)
    !
    call zagmg_prepareagg(n,a,ja,ia,dt(l)%idiag,ind2,iperm,si1,ndd,l )
    !
    IF (ndd .EQ. n) THEN
       nc=0
       nzc=0
       DEALLOCATE(si1,iperm,ind2)
       memi=memi-2*n
       memr=memr-n
       IF (l .EQ. 1) THEN
          DEALLOCATE(riperm)
          memi=memi-n
       END IF
       GOTO 999
    END IF
    !
    ALLOCATE(ldd(ndd),lcg(2*(n-ndd)))
    memi=memi+2*n-ndd
    memax=MAX(memax,memr+memi*rlenilen)
    !
    IF (dble(n) .GT. targetcoarsefac*(n-ndd)) THEN
       skipass=.TRUE.
       npass1=npass+1
    ELSE
       skipass=.FALSE.
       npass1=npass
    END IF
    !
    !
    IF (l > 1) THEN
       IF (spd) THEN
          CALL zagmg_findpairs_SI(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       ELSE
          CALL zagmg_findpairs_GI(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,iperm )
       END IF
       DEALLOCATE(iperm)
       memi=memi-n
    ELSE
       IF (spd) THEN
          CALL zagmg_findpairs_SI1(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       ELSE
          CALL zagmg_findpairs_GI1(n,a,ja,ia,dt(l)%idiag,si1  &
               ,ind2,lcg,nc,ndd,ldd,skipass,riperm,iperm )
       END IF
       DEALLOCATE(iperm,riperm)
       memi=memi-2*n
    END IF
10  CONTINUE
    nz=ia(n+1)-1
    !
    IF (npass1.GT.1) THEN
       isize=nc
    ELSE
       isize=1
    END IF
    ALLOCATE(an(nz-2*(n-nc)+ndd),jan(nz-2*(n-nc)+ndd)             &
         ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
    memi=memi+nz-2*(n-nc)+ndd+2*(nc+1)+2*nc+isize
    memr=memr+nz-2*(n-nc)+ndd+nc
    memax=MAX(memax,memr+memi*rlenilen)
    CALL zagmg_setcg(n,a,ja,ia,dt(l)%idiag,si1,ind2,lcg,nc,an,jan,ian      &
         ,idiagn,sinn,npass1.GT.1,maxdg,iw,wc )
    DEALLOCATE(wc,iw)
    memi=memi-2*nc
    memr=memr-nc
    nzc=ian(nc+1)-1
    IF (dble(nz).GT.targetcoarsefac*nzc .OR. npass1.LE.1) THEN
       DEALLOCATE(si1,sinn,lcg)
       IF(ALLOCATED(ldd)) DEALLOCATE(ldd)
       memi=memi-2*(n-ndd)
       memr=memr-n-isize
       dt(l)%ind  => ind2
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(ind2,an,jan,ian,idiagn)
       GOTO 999
    END IF
    !
    DEALLOCATE(ind2)
    memi=memi-n
    !
    lcg1 => lcg
    NULLIFY(lcg)
    m1=1
    !
    !
    !
    DO kpass=2,npass1
       m1=2*m1
       np  = nc
       nzp = nzc
       ap     => an
       jap    => jan
       iap    => ian
       idiagp => idiagn
       sip    => sinn
       NULLIFY(an,jan,ian,idiagn,sinn)
       ALLOCATE(lcg(2*np),ind2(np),w(maxdg),iw(maxdg))
       memi=memi+maxdg+3*np
       memr=memr+maxdg
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !
       ind2(1:np)=-1
       IF (spd) THEN
          CALL zagmg_findpairs_SF(np,ap,jap,iap,idiagp,sip  &
               ,ind2,lcg,nc                                  &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       ELSE
          CALL zagmg_findpairs_GF(np,ap,jap,iap,idiagp,sip     &
               ,ind2,lcg,nc                                     &
               ,m1,lcg1,a,ja,ia,dt(l)%idiag,si1,w,iw )
       END IF
       DEALLOCATE(w,iw)
       memi=memi-maxdg
       memr=memr-maxdg
       IF (kpass.NE.npass1) THEN
          isize=nc
       ELSE
          isize=1
          DEALLOCATE(si1)
          memr=memr-n
       END IF
       !
       !
       ALLOCATE(an(nzp-2*(np-nc)),jan(nzp-2*(np-nc))                 &
            ,ian(nc+1),idiagn(nc+1),sinn(isize),wc(nc),iw(2*nc))
       memi=memi+nzp-2*(np-nc)+2*(nc+1)+2*nc+isize
       memr=memr+nzp-2*(np-nc)+nc
       memax=MAX(memax,memr+memi*rlenilen)
       !
       !
       CALL zagmg_setcg(np,ap,jap,iap,idiagp,sip,ind2,lcg,nc,an     &
            ,jan,ian,idiagn,sinn,kpass.NE.npass1,maxdg,iw,wc )
       memi=memi-SIZE(jap)-2*(np+1)-np-2*nc
       memr=memr-SIZE(ap)-SIZE(sip)-nc
       DEALLOCATE(ap,jap,iap,idiagp,sip,ind2,wc,iw)
       !
       !
       ALLOCATE(lcgn(2*m1*nc))
       memi=memi+2*m1*nc
       memax=MAX(memax,memr+memi*rlenilen)
       CALL zagmg_lcgmix(nc,m1,lcg1,lcg,lcgn)
       memi=memi-SIZE(lcg)-SIZE(lcg1)
       DEALLOCATE(lcg,lcg1)
       lcg1 => lcgn
       NULLIFY(lcgn)
       nzc=ian(nc+1)-1
       IF ( kpass.NE.npass1 .AND. dble(nz).GT.targetcoarsefac*nzc ) THEN
          DEALLOCATE(si1)
          memr=memr-n
          EXIT
       END IF
    END DO
    !
    memr=memr-SIZE(sinn)
    DEALLOCATE(sinn)
    !
    ALLOCATE(dt(l)%ind(n))
    memi=memi+n
    memax=MAX(memax,memr+memi*rlenilen)
    CALL zagmg_setind(nc,ndd,ldd,lcg1,2*m1,dt(l)%ind)
    memi=memi-ndd-SIZE(lcg1)
    DEALLOCATE(lcg1,ldd)
    !
       dt(l+1)%a    => an
       dt(l+1)%ja   => jan
       dt(l+1)%ia   => ian
       dt(l+1)%idiag=> idiagn
       NULLIFY(an,jan,ian,idiagn)
999 CONTINUE
    !
    RETURN
901 FORMAT(i3,'*SETUP: Coarsening by multiple pairwise aggregation')
902 FORMAT('****       Quality threshold (',A6,'):',f6.2, &
         ' ;  Strong diag. dom. trs:',f5.2)
903 FORMAT('****         Maximal number of passes:',i3,     &
         '  ; Target coarsening factor:',f5.2)
904 FORMAT('****           Diag. dom. checked w.r.t. sum of offdiag', &
          ' (no absolute vaues)')
905 FORMAT('****',22x,'Threshold for rows with large pos. offdiag.:',f5.2)
906 FORMAT('****  Rmk: Setup performed assuming the matrix symmetric')
908 FORMAT('****  Rmk: Setup performed for the transpose of the input matrix')
  END SUBROUTINE zagmg_aggregation
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_setCMK(n,ja,ia,idiag,riperm,iperm,iext)
    !
    !
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),riperm(*),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    LOGICAL :: exc
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,mindg,kdim
    INTEGER :: ifirst,ilast,kb,i0
    !
    ifirst=1
    ilast=n
    i2=ifirst
    mindg=n+1
    DO i = ifirst,ilast
       dg=ia(i+1)-ia(i)
       IF (dg .GT. 1) THEN
          iperm(i)=-dg
          IF (dg.LT.mindg) THEN
             mindg=dg
             jj=i
          END IF
       ELSE
          riperm(i2)=i
          iperm(i)=i2
          i2=i2+1
       END IF
    ENDDO
    !
    ijs=ifirst-1
    i1=i2
15  CONTINUE
    !
    IF (i2 .LE. ilast) THEN
      riperm(i2)=jj
      iperm(jj)=i2
    END IF
    !
    DO WHILE (i1.LE.i2 .AND. i2.LT.ilast)
       !
       !
       i=riperm(i1)
       ijs1=i2+1
       !
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       DO kk = j1,jd-1
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       DO kk = jd+1,j2
          j=ja(kk)
          IF (iperm(j) .LT. 0) THEN
             i2=i2+1
             riperm(i2)=j
          END IF
       ENDDO
       !
       ijs2=i2
       exc=.TRUE. .AND. ijs2.GT.ijs1
       DO WHILE(exc)
          exc=.FALSE.
          DO kk=ijs1+1,ijs2
             IF( iperm(riperm(kk)) .GT. iperm(riperm(kk-1)) )THEN
                j=riperm(kk)
                riperm(kk)=riperm(kk-1)
                riperm(kk-1)=j
                exc=.TRUE.
             END IF
          END DO
       END DO
       DO kk=ijs1,ijs2
          iperm(riperm(kk))=kk
       END DO
       !
       i1=i1+1
    END DO
    IF (i2 .LT. ilast) THEN
       !
       jj=0
       DO WHILE (jj .EQ. 0)
          ijs=ijs+1
          IF (ijs .GT. ilast) THEN
             mindg=mindg+1
             ijs=ifirst
          END IF
          ijs1=ijs
          IF (iperm(ijs1).LT.0 .AND. ia(ijs1+1)-ia(ijs1).EQ.mindg) &
               jj=ijs1
       END DO
       i2=i2+1
       GOTO 15
    END IF
    !
 !
 !
    RETURN
  END SUBROUTINE zagmg_setCMK
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_prepareagg(n,a,ja,ia,idiag,ind2,iperm,si,ndd,l,iext)
    !
    !
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind2(n),iperm(n)
    INTEGER, OPTIONAL :: iext(*)
    INTEGER :: ndd, l
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) , TARGET :: si(n)
    REAL(kind(0.0d0)) :: checkddl,oda,odm,ods,vald
    INTEGER :: i,j,jj,jk,jd,k,kk,j1,j2,i1,i2,ijs,ijs1,ijs2,dg,kdim,nnegrcs
    REAL(kind(0.0d0)), POINTER, DIMENSION(:) :: odmax,odabs,osi
    !
    IF (.NOT.spd) THEN
       checkddl=checkddJ
    ELSE
       checkddl=checkddB
    END IF
    !
    IF (.NOT.spd) THEN
       !
       !
       IF (checkdd > 0) THEN
          ALLOCATE(odmax(n),odabs(n))
          memi=memi+2*n
          odabs(1:n)=0.0d0
       ELSE
          ALLOCATE(odmax(n))
          memi=memi+n
       END IF
       osi => si(1:n)
       si(1:n)=0.0d0
       odmax(1:n)=0.0d0
       memax=MAX(memax,memr+memi*rlenilen)
       DO i=n,1,-1
          j =ia(i)
          jd=idiag(i)
          jj=ia(i+1)-1
          DO k = j,jd-1
             jk=ja(k)
             osi(jk)=osi(jk)+dble(a(k))
             odmax(jk)=max(odmax(jk),dble(a(k)))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
          DO k = jd+1,jj
             jk=ja(k)
             osi(jk)=osi(jk)+dble(a(k))
             odmax(jk)=max(odmax(jk),dble(a(k)))
             IF (checkdd > 0) odabs(jk)=odabs(jk)+abs(a(k))
          ENDDO
       ENDDO
    END IF
    !
    ndd=0
    nnegrcs=0
    !
    DO i=1,n
       j1 =ia(i)
       jd=idiag(i)
       j2 = ia (i+1)-1
       vald = dble(a(jd))
       odm=0.0d0
       oda=0.0d0
       ods=0.0d0
       DO kk = j1,jd-1
          ods=ods+dble(a(kk))
          odm=max(odm,dble(a(kk)))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       DO kk = jd+1,j2
          ods=ods+dble(a(kk))
          odm=max(odm,dble(a(kk)))
          IF (checkdd > 0) oda=oda+abs(a(kk))
       ENDDO
       !
       IF (.NOT.spd) THEN
          ods=(osi(i)+ods)/2
          odm=max(odm,odmax(i))
          IF (checkdd > 0) oda=(oda+odabs(i))/2
       END IF
       !
       IF ((vald+ods) .LT. -repsmach*ABS(vald)) nnegrcs=nnegrcs+1
       !
       !
       si(i)=-ods
       IF ( (checkdd.GT.0 .AND. vald.GT.checkddl*oda)     &
            .OR. (checkdd.LT.0 .AND. vald.GT.checkddl*abs(ods)) ) THEN
          !
          ind2(i)=0
          ndd=ndd+1
       ELSE
          ind2(i)=-1
          IF (odm .GT. trspos*vald) iperm(i)=0
       ENDIF
    END DO
    !
    !
    zerors=.FALSE.
    IF (nnegrcs.GT.fracnegrcsum*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       zerors=.TRUE.
       ndd=0
       ind2(1:n)=-1
    END IF
    !
    IF (.NOT.spd) THEN
       IF (checkdd > 0) THEN
          DEALLOCATE(odmax,odabs)
          memi=memi-2*n
       ELSE
          DEALLOCATE(odmax)
          memi=memi-n
       END IF
    END IF
    RETURN
  END SUBROUTINE zagmg_prepareagg
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_findpairs_GF(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    COMPLEX(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(lcg1(2,j).EQ.0 .OR. ind(j).GE.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-dble(a(i))/2
          IF(kk .NE. 0) vals=vals-dble(a(kk))/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
             eta1=2*dble(a(idiag(isel)))
             eta2=2*dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       CALL zagmg_checktentagg_GF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_checktentagg_GF
!
!
!
!
!
    INTEGER, PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m,info, setdim1, setdim, l2, k2
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=dble(a1(idiag1(jj)))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=dble(a1(k))/2
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
       DO k=ia1(jj),idiag1(jj)-1
          DO l=1,j-1
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=dble(a1(k))/2
                W(j,l)=W(j,l)+alpha
                W(l,j)=W(j,l)
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
       v(j)=W(j,j)
   !
   !
   !
       W(j,j)=umdbndmum1*W(j,j)-abs(sig(j))
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
!
    beta=dbndmum1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
!
    acc=.FALSE.
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE zagmg_checktentagg_GF
  END SUBROUTINE zagmg_findpairs_GF
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_findpairs_SF(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,m1,lcg1,a1,ja1,ia1,idiag1,si1,rtent,jtent )
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: m1,ja1(*),ia1(*),idiag1(*),jtent(*),lcg1(m1,n)
    COMPLEX(kind(0.0d0)) :: a1(*)
    REAL(kind(0.0d0)) :: si1(*),rtent(*)
!
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (lcg1(2,isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       ntentleft=0
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(lcg1(2,j).EQ.0 .OR. ind(j).GE.0) CYCLE
          vals=-dble(a(i))
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          ntentleft=ntentleft+1
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          rtent(ntentleft)=tent
          jtent(ntentleft)=j
          CYCLE
9         CONTINUE
          rtent(ntentleft)=val
          jtent(ntentleft)=ipair
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       IF (ipair .EQ. 0) GOTO 25
20     CONTINUE
       CALL zagmg_checktentagg_SF
       IF (.NOT.acc) THEN
          ipair = 0
          IF (ntentleft .GT.0) THEN
             i=1
             j=1
             DO WHILE (i .LE. ntentleft)
                IF (jtent(j).GT.0) THEN
                   tent=rtent(j)
                   IF (ipair.EQ.0) GOTO 22
                   IF (16*(tent-val).LT.-1) GOTO 22
                   IF (16*(tent-val).LT.1 .AND. j.LT.ipair) GOTO 22
                   GOTO 23
22                 CONTINUE
                   val=tent
                   ipair=jtent(j)
                   ijtent=j
23                 CONTINUE
                   i=i+1
                END IF
                j=j+1
             END DO
             ntentleft=ntentleft-1
             jtent(ijtent)=0
             GOTO 20
          END IF
       END IF
       !
25     CONTINUE
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
    RETURN
  CONTAINS
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_checktentagg_SF
!
!
!
!
!
    INTEGER, PARAMETER :: mm=max(2**(npass+1),8)
    REAL(kind(0.0d0)) :: W(mm,mm), sig(mm), AGe(mm), v(mm)
    REAL(kind(0.0d0)) :: alpha, alp, tmp, beta, f1, f2
    INTEGER :: j,jj,k,l,m,info, setdim1, setdim, l2, k2
    INTEGER :: set(mm), l1, wdthT
    REAL(kind(0.0d0)) :: T
    LOGICAL :: exc
!
    IF (m1.eq.2) THEN
       IF (lcg1(2,isel) .LT. 0) THEN
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             setdim=2
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(1,ipair)
             set(3)=lcg1(2,ipair)
             setdim=3
          END IF
          l1=1
       ELSE
          IF (lcg1(2,ipair) .LT. 0) THEN
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             setdim=3
          ELSE
             set(1)=lcg1(1,isel)
             set(2)=lcg1(2,isel)
             set(3)=lcg1(1,ipair)
             set(4)=lcg1(2,ipair)
             setdim=4
          END IF
          l1=2
       END IF
    ELSE
       l1=m1
       IF (lcg1(m1,isel).LT.0) l1=-lcg1(m1,isel)
       set(1:l1)=lcg1(1:l1,isel)
       l2=m1
       IF (lcg1(m1,ipair).LT.0) l2=-lcg1(m1,ipair)
       set(l1+1:l1+l2)=lcg1(1:l2,ipair)
       setdim=l1+l2
    END IF
!
    exc=.TRUE.
    DO WHILE(exc)
       exc=.FALSE.
       DO l=2,SetDim
          IF( set(l)<set(l-1) )THEN
             jj=set(l)
             set(l)=set(l-1)
             set(l-1)=jj
             exc=.TRUE.
          END IF
       END DO
    END DO
!
    DO j=1,SetDim
       jj=Set(j)
       sig(j)=si1(jj)
       IF (zerors) THEN
          W(j,j)=sig(j)
          AGe(j)=0.0d0
       ELSE
          W(j,j)=dble(a1(idiag1(jj)))
          AGe(j)=W(j,j)-sig(j)
       END IF
       l2=j+1
       DO l=l2,SetDim
          W(j,l)=0.0d0
          W(l,j)=0.0d0
       END DO
       k2=ia1(jj+1)-1
       DO k=idiag1(jj)+1,k2
          DO l=l2,SetDim
             m=Set(l)
             IF(ja1(k)==m)THEN
                alpha=dble(a1(k))
                W(j,l)=alpha
                W(l,j)=alpha
                EXIT
             END IF
          END DO
       END DO
    END DO
!
    DO j=1,SetDim
       DO k=1,SetDim
          IF (j.ne.k) THEN
             sig(j)=sig(j)+W(j,k)
          END IF
       ENDDO
       IF (sig(j) < 0.0d0)  AGe(j)=AGe(j)+2*sig(j)
   !
   !
       W(j,j)=W(j,j)-abs(sig(j))
   !
   !
   !
   !
       tmp=2*abs(sig(j))
       W(j,j)=W(j,j)-bndmum1m1*tmp
   !
       v(j)=tmp+AGe(j)
       IF (j .eq. 1) THEN
          beta=v(j)
          alp=abs(AGe(j))
       ELSE
          beta=beta+v(j)
          alp=max(alp,abs(AGe(j)))
       END IF
    END DO
!
!
    beta=bndmum1m1/beta
    DO j=1,SetDim
       DO k=1,SetDim
          W(j,k)=W(j,k)+beta*v(j)*v(k)
       END DO
    END DO
!
!
    IF (alp.LT.repsmach*beta) THEN
       SetDim1=SetDim-1
    ELSE
       SetDim1=SetDim
    END IF
!
!
    acc=.FALSE.
!
    SELECT CASE (SetDim1)
    CASE (1)
       GOTO 11
    CASE (2)
       GOTO 12
    CASE (3)
       GOTO 13
    CASE (4)
       GOTO 14
    CASE (5)
       GOTO 15
    CASE (6)
       GOTO 16
    CASE (7)
       GOTO 17
    CASE (8)
       GOTO 18
    CASE DEFAULT
       CALL DPOTRF('U',SetDim1,W,mm,info)
       IF (info .NE. 0) RETURN
       GOTO 10
    END SELECT
18  CONTINUE
    IF (W(8,8) .LE. 0.0d0) RETURN
    W(7,7) = W(7,7) - (W(7,8)/W(8,8)) * W(7,8)
    T = W(6,8)/W(8,8)
    W(6,7) = W(6,7) - T * W(7,8)
    W(6,6) = W(6,6) - T * W(6,8)
    T = W(5,8)/W(8,8)
    W(5,7) = W(5,7) - T * W(7,8)
    W(5,6) = W(5,6) - T * W(6,8)
    W(5,5) = W(5,5) - T * W(5,8)
    T = W(4,8)/W(8,8)
    W(4,7) = W(4,7) - T * W(7,8)
    W(4,6) = W(4,6) - T * W(6,8)
    W(4,5) = W(4,5) - T * W(5,8)
    W(4,4) = W(4,4) - T * W(4,8)
    T = W(3,8)/W(8,8)
    W(3,7) = W(3,7) - T * W(7,8)
    W(3,6) = W(3,6) - T * W(6,8)
    W(3,5) = W(3,5) - T * W(5,8)
    W(3,4) = W(3,4) - T * W(4,8)
    W(3,3) = W(3,3) - T * W(3,8)
    T = W(2,8)/W(8,8)
    W(2,7) = W(2,7) - T * W(7,8)
    W(2,6) = W(2,6) - T * W(6,8)
    W(2,5) = W(2,5) - T * W(5,8)
    W(2,4) = W(2,4) - T * W(4,8)
    W(2,3) = W(2,3) - T * W(3,8)
    W(2,2) = W(2,2) - T * W(2,8)
    T = W(1,8)/W(8,8)
    W(1,7) = W(1,7) - T * W(7,8)
    W(1,6) = W(1,6) - T * W(6,8)
    W(1,5) = W(1,5) - T * W(5,8)
    W(1,4) = W(1,4) - T * W(4,8)
    W(1,3) = W(1,3) - T * W(3,8)
    W(1,2) = W(1,2) - T * W(2,8)
    W(1,1) = W(1,1) - T * W(1,8)
17  CONTINUE
    IF (W(7,7) .LE. 0.0d0) RETURN
    W(6,6) = W(6,6) - (W(6,7)/W(7,7)) * W(6,7)
    T = W(5,7)/W(7,7)
    W(5,6) = W(5,6) - T * W(6,7)
    W(5,5) = W(5,5) - T * W(5,7)
    T = W(4,7)/W(7,7)
    W(4,6) = W(4,6) - T * W(6,7)
    W(4,5) = W(4,5) - T * W(5,7)
    W(4,4) = W(4,4) - T * W(4,7)
    T = W(3,7)/W(7,7)
    W(3,6) = W(3,6) - T * W(6,7)
    W(3,5) = W(3,5) - T * W(5,7)
    W(3,4) = W(3,4) - T * W(4,7)
    W(3,3) = W(3,3) - T * W(3,7)
    T = W(2,7)/W(7,7)
    W(2,6) = W(2,6) - T * W(6,7)
    W(2,5) = W(2,5) - T * W(5,7)
    W(2,4) = W(2,4) - T * W(4,7)
    W(2,3) = W(2,3) - T * W(3,7)
    W(2,2) = W(2,2) - T * W(2,7)
    T = W(1,7)/W(7,7)
    W(1,6) = W(1,6) - T * W(6,7)
    W(1,5) = W(1,5) - T * W(5,7)
    W(1,4) = W(1,4) - T * W(4,7)
    W(1,3) = W(1,3) - T * W(3,7)
    W(1,2) = W(1,2) - T * W(2,7)
    W(1,1) = W(1,1) - T * W(1,7)
16  CONTINUE
    IF (W(6,6) .LE. 0.0d0) RETURN
    W(5,5) = W(5,5) - (W(5,6)/W(6,6)) * W(5,6)
    T = W(4,6)/W(6,6)
    W(4,5) = W(4,5) - T * W(5,6)
    W(4,4) = W(4,4) - T * W(4,6)
    T = W(3,6)/W(6,6)
    W(3,5) = W(3,5) - T * W(5,6)
    W(3,4) = W(3,4) - T * W(4,6)
    W(3,3) = W(3,3) - T * W(3,6)
    T = W(2,6)/W(6,6)
    W(2,5) = W(2,5) - T * W(5,6)
    W(2,4) = W(2,4) - T * W(4,6)
    W(2,3) = W(2,3) - T * W(3,6)
    W(2,2) = W(2,2) - T * W(2,6)
    T = W(1,6)/W(6,6)
    W(1,5) = W(1,5) - T * W(5,6)
    W(1,4) = W(1,4) - T * W(4,6)
    W(1,3) = W(1,3) - T * W(3,6)
    W(1,2) = W(1,2) - T * W(2,6)
    W(1,1) = W(1,1) - T * W(1,6)
15  CONTINUE
    IF (W(5,5) .LE. 0.0d0) RETURN
    W(4,4) = W(4,4) - (W(4,5)/W(5,5)) * W(4,5)
    T = W(3,5)/W(5,5)
    W(3,4) = W(3,4) - T * W(4,5)
    W(3,3) = W(3,3) - T * W(3,5)
    T = W(2,5)/W(5,5)
    W(2,4) = W(2,4) - T * W(4,5)
    W(2,3) = W(2,3) - T * W(3,5)
    W(2,2) = W(2,2) - T * W(2,5)
    T = W(1,5)/W(5,5)
    W(1,4) = W(1,4) - T * W(4,5)
    W(1,3) = W(1,3) - T * W(3,5)
    W(1,2) = W(1,2) - T * W(2,5)
    W(1,1) = W(1,1) - T * W(1,5)
14  CONTINUE
    IF (W(4,4) .LE. 0.0d0) RETURN
    W(3,3) = W(3,3) - (W(3,4)/W(4,4)) * W(3,4)
    T = W(2,4)/W(4,4)
    W(2,3) = W(2,3) - T * W(3,4)
    W(2,2) = W(2,2) - T * W(2,4)
    T = W(1,4)/W(4,4)
    W(1,3) = W(1,3) - T * W(3,4)
    W(1,2) = W(1,2) - T * W(2,4)
    W(1,1) = W(1,1) - T * W(1,4)
13  CONTINUE
    IF (W(3,3) .LE. 0.0d0) RETURN
    W(2,2) = W(2,2) - (W(2,3)/W(3,3)) * W(2,3)
    T = W(1,3)/W(3,3)
    W(1,2) = W(1,2) - T * W(2,3)
    W(1,1) = W(1,1) - T * W(1,3)
12  CONTINUE
    IF (W(2,2) .LE. 0.0d0) RETURN
    W(1,1) = W(1,1) - (W(1,2)/W(2,2)) * W(1,2)
11  CONTINUE
    IF (W(1,1) .LE. 0.0d0) RETURN
10  CONTINUE
!
    acc=.TRUE.
!
    RETURN
  END SUBROUTINE zagmg_checktentagg_SF
  END SUBROUTINE zagmg_findpairs_SF
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_findpairs_GI(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ipc(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-dble(a(i))/2
          IF(kk .NE. 0) vals=vals-dble(a(kk))/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
             eta1=2*dble(a(idiag(isel)))
             eta2=2*dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_GI
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_findpairs_SI(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,ipc )
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: ipc(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (ipc(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(ipc(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          vals=-dble(a(i))
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. j.LT.ipair)  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_SI
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_findpairs_GI1(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_dampJac
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(iperm(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          kk=0
          IF (i .LT. idiag(isel)) THEN
             j2=ia(j+1)-1
             DO jk=idiag(j)+1,j2
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ELSE
             DO jk=ia(j),idiag(j)-1
                IF (ja(jk) .EQ. isel) THEN
                   kk=jk
                   EXIT
                END IF
             END DO
          ENDIF
          vals=-dble(a(i))/2
          IF(kk .NE. 0) vals=vals-dble(a(kk))/2
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
             eta1=2*si(isel)
             eta2=2*si(j)
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
             eta1=2*dble(a(idiag(isel)))
             eta2=2*dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
          ELSE
             del1=rsi+2*sig1
          END IF
          IF (sig2 > 0.0d0) THEN
             del2=rsj
          ELSE
             del2=rsj+2*sig2
          END IF
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=((eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=(eta1*eta2)/(eta1+eta2)
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_dampJac
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_GI1
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_findpairs_SI1(n,a,ja,ia,idiag,si,ind,lcg,nc     &
    ,nddl,ldd,skipass                                      &
    ,riperm,iperm )
!
!
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),lcg(2,*),nc
    COMPLEX(kind(0.0d0)) :: a(*)
    REAL(kind(0.0d0)) :: si(n)
    INTEGER :: ldd(nddl),nddl
    LOGICAL :: skipass
    INTEGER :: iperm(n),riperm(n)
!
!
!
!
!
!
!
!
!
!
!
!----------------
!
    REAL(kind(0.0d0)) :: val,vals,valp,tent,rsi,rsj,epsr
    REAL(kind(0.0d0)) :: del1,del2,eta1,eta2,del12,sig1,sig2,rnd,vald
    INTEGER :: mindg,i,j,jj,k,kk,jk,isel,dg,ipair,nmark,idd
    INTEGER :: i1,i2,i3,ijs,ntentleft,npc,itrs,ijtent,j2
    LOGICAL :: acc,ifirst
    REAL(kind(0.0d0)) :: kaptg,bndmum1m1,dbndmum1,umdbndmum1
!
    ifirst=.TRUE.
  1 CONTINUE
!
    kaptg=imult*kaptg_blocdia
    bndmum1m1=1.0d0/(kaptg-1.0d0)
    dbndmum1=2*1.0d0/kaptg
    umdbndmum1=1.0d0-dbndmum1
    idd=0
    nmark=0
    nc=0
    ijs=1
    npc=0
    DO WHILE (nmark.LT.n)
       isel=ijs
       isel=riperm(ijs)
       ijs=ijs+1
       !
       !
       IF (ind(isel) .EQ. 0) THEN
          idd=idd+1
          nmark=nmark+1
          ldd(idd)=isel
          CYCLE
       END IF
       !
       IF (ind(isel) .GE. 0) CYCLE
       !
       nc = nc + 1
       lcg(1,nc) = isel
       nmark = nmark+1
       ind(isel) = nc
       ipair = 0
       !
       IF (iperm(isel) .EQ. 0) THEN
          lcg(2,nc) = 0
          npc=npc+1
          CYCLE
       END IF
       IF (skipass) THEN
          lcg(2,nc) = -1
          CYCLE
       END IF
       !
       !
       i2=ia(isel+1)-1
       DO i = ia(isel),i2
          IF (i .EQ. idiag(isel)) CYCLE
          j = ja (i)
          IF(iperm(j).EQ.0 .OR. ind(j).GE.0) CYCLE
          vals=-dble(a(i))
          IF (zerors) THEN
             rsi=0.0d0
             rsj=0.0d0
          ELSE
             rsi=-si(isel)+dble(a(idiag(isel)))
             rsj=-si(j)+dble(a(idiag(j)))
          END IF
          sig1=si(isel)-vals
          sig2=si(j)-vals
          !
          !
          IF (sig1 > 0.0d0) THEN
             del1=rsi
             eta1=rsi+2*sig1
          ELSE
             del1=rsi+2*sig1
             eta1=rsi
          END IF
          IF (eta1 < 0.0d0) CYCLE
          IF (sig2 > 0.0d0) THEN
             del2=rsj
             eta2=rsj+2*sig2
          ELSE
             del2=rsj+2*sig2
             eta2=rsj
          END IF
          IF (eta2 < 0.0d0) CYCLE
          IF (vals > 0.0d0) THEN
             epsr=repsmach*vals
             IF (ABS(del1) < epsr .AND. ABS(del2) < epsr) THEN
               valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del1) < epsr) THEN
                IF (del2 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE IF (ABS(del2) < epsr) THEN
                IF (del1 < -epsr) CYCLE
                valp=1.0d0+(eta1*eta2)/(vals*(eta1+eta2))
             ELSE
                del12=del1+del2
                IF (del12 < -epsr) CYCLE
                valp=vals+del1*del2/del12
                IF (valp < 0.0d0) CYCLE
                valp=(vals+(eta1*eta2)/(eta1+eta2))/valp
             END IF
          ELSE
             IF (del1 .LE. 0.0d0 .OR. del2 .LE. 0.0d0) CYCLE
             valp=vals+del1*del2/(del1+del2)
             IF (valp < 0.0d0) CYCLE
             vals=vals+(eta1*eta2)/(eta1+eta2)
             IF (vals < 0.0d0) CYCLE
             valp=vals/valp
          END IF
          IF (valp > kaptg) CYCLE
          !
          !
          tent=valp
          IF (ipair.EQ.0) GOTO 10
          IF (16*(tent-val).LT.-1) GOTO 9
          IF (16*(tent-val).LT.1 .AND. iperm(j).LT.iperm(ipair))  GOTO 9
          CYCLE
9         CONTINUE
10        CONTINUE
          ipair = j
          val = tent
       ENDDO
       !
       IF (ipair .EQ. 0) THEN
          lcg(2,nc) = -1
       ELSE
          lcg(2,nc) = ipair
          ind(ipair) = nc
          nmark = nmark+1
       END IF
    ENDDO
!
!
    IF (ifirst .AND. 3*nc.GT.2*n .AND. ngl(1).GT.maxcoarseslowt) THEN
       imult=2*imult
       IF (wfo) THEN
          WRITE (iout,901) IRANK,imult*kaptg_blocdia
       END IF
901    FORMAT(i3, &
       '*   Coarsening too slow: quality threshold increased to',f6.2)
       ind(1:n)=-1
       DO i=1,nddl
          ind(ldd(i))=0
       END DO
       ifirst=.FALSE.
       GOTO 1
    ENDIF
    RETURN
  END SUBROUTINE zagmg_findpairs_SI1
!------------------------------------------------------------------
  SUBROUTINE zagmg_lcgmix(nc,m,lcg1,lcg,lcgn)
    INTEGER :: nc,m,lcg1(m,*),lcg(2,*),lcgn(2*m,*),i,l,l1,l2
    IF (m.eq.2) THEN
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(4,i)=-1
          ELSE IF(lcg(2,i) .LT. 0) THEN
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=-1
                lcgn(4,i)=-1
             ELSE
                lcgn(1,i)=lcg1(1,lcg(1,i))
                lcgn(2,i)=lcg1(2,lcg(1,i))
                lcgn(4,i)=-2
             END IF
          ELSE
             IF (lcg1(2,lcg(1,i)) .LT. 0) THEN
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-2
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(1,lcg(2,i))
                   lcgn(3,i)=lcg1(2,lcg(2,i))
                   lcgn(4,i)=-3
                END IF
             ELSE
                IF (lcg1(2,lcg(2,i)) .LT. 0) THEN
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=-3
                ELSE
                   lcgn(1,i)=lcg1(1,lcg(1,i))
                   lcgn(2,i)=lcg1(2,lcg(1,i))
                   lcgn(3,i)=lcg1(1,lcg(2,i))
                   lcgn(4,i)=lcg1(2,lcg(2,i))
                END IF
             END IF
          END IF
       END DO
    ELSE
       DO i=1,nc
          IF(lcg(2,i) .EQ. 0) THEN
             lcgn(1,i)=lcg1(1,lcg(1,i))
             lcgn(2,i)=0
             lcgn(2*m,i)=-1
          ELSE
             lcgn(2,i)=-1
             l1=m
             IF (lcg1(m,lcg(1,i)).LT.0) l1=-lcg1(m,lcg(1,i))
             lcgn(1:l1,i)=lcg1(1:l1,lcg(1,i))
             IF(lcg(2,i) .LT. 0) THEN
                l=l1
             ELSE
                l2=m
                IF (lcg1(m,lcg(2,i)).LT.0) l2=-lcg1(m,lcg(2,i))
                lcgn(l1+1:l1+l2,i)=lcg1(1:l2,lcg(2,i))
                l=l1+l2
             END IF
             IF(l .LT. 2*m) lcgn(2*m,i)=-l
          END IF
       END DO
    END IF
    RETURN
  END SUBROUTINE zagmg_lcgmix
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_setind(nc,ndd,ldd,lcg,m,ind)
    INTEGER :: nc,m,lcg(m,*),nll,ldd(ndd),ind(*),i,k,l
    DO i=1,ndd
       ind(ldd(i))=0
    END DO
    DO i=1,nc
       l=m
       IF (lcg(m,i) .LT. 0) l=-lcg(m,i)
       DO k=1,l
          ind(lcg(k,i))=i
       END DO
    END DO
    RETURN
  END SUBROUTINE zagmg_setind
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_setcg(n,a,ja,ia,idiag,si,ind,lcg           &
       ,nc,a2,ja2,ia2,idiag2,si2,ysi,maxdg,iw,w,iext,iext2)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ja(*),ia(n+1),idiag(n),ind(n),nc,lcg(2,nc)
    INTEGER :: ja2(*),ia2(nc+1),idiag2(nc),maxdg
    INTEGER, TARGET :: iw(2*nc)
    INTEGER, OPTIONAL :: iext(*),iext2(*)
    COMPLEX(kind(0.0d0)) :: a(*),a2(*),w(nc),vald
    REAL(kind(0.0d0)) :: si(n),si2(*),sii
    LOGICAL :: ysi
    INTEGER :: nz,nzu,i,jj,jc,jcol,ki,kb,kf,jpos
    INTEGER, POINTER, DIMENSION(:) :: iw1, iw2
    !
    iw1 => iw(1:nc)
    iw2 => iw(nc+1:2*nc)
    !
    nz = 0
    iw1(1:nc)=0
    maxdg=0
    ia2(1)=1
    DO i = 1,nc
       sii=0.0d0
       vald=cmplx(0,0,kind(0.0d0))
       nzu=0
       DO ki= 1,2
          jj = lcg(ki,i)
          IF (ki.EQ.1 .OR. jj.GT.0) THEN
             IF (ysi) sii=sii+dble(si(jj))
             kf = ia(jj+1) - 1
             DO kb = ia(jj),kf
                jc = ja(kb)
                jcol = ind(jc)
                IF (jcol .GT. 0) THEN
                   IF (jcol .LT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nz = nz+1
                         ja2(nz) = jcol
                         iw1(jcol) = nz
                         a2(nz) = a(kb)
                      ELSE
                         a2(jpos) = a2(jpos) + a(kb)
                      ENDIF
                   ELSE IF (jcol .GT. i) THEN
                      jpos = iw1(jcol)
                      IF (jpos.EQ.0) THEN
                         nzu = nzu+1
                         iw2(nzu) = jcol
                         iw1(jcol) = nzu
                         w(nzu) = a(kb)
                      ELSE
                         w(jpos) = w(jpos) + a(kb)
                      ENDIF
                   ELSE
                      vald=vald+a(kb)
                      IF (ysi .AND. jc.NE.jj) sii=sii+dble(a(kb))
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
       ENDDO
       nz=nz+1
       a2(nz)=vald
       idiag2(i)=nz
       ja2(nz)=i
       a2(nz+1:nz+nzu)=w(1:nzu)
       ja2(nz+1:nz+nzu)=iw2(1:nzu)
       nz=nz+nzu
       maxdg=max(maxdg,nz-ia2(i))
       DO kb = ia2(i), nz
          iw1(ja2(kb))=0
       ENDDO
       IF (ysi) si2(i)=sii
       ia2(i+1)=nz+1
    ENDDO
    RETURN
  END SUBROUTINE zagmg_setcg
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_LAPACK(n,f,ijb,flop,a,ja,ia)
    USE zagmg_mem
    IMPLICIT NONE
    INTEGER :: n,ijb
    COMPLEX(kind(0.0d0)), TARGET :: f(n)
    REAL(kind(0.0d0)) :: flop
    COMPLEX(kind(0.0d0)), OPTIONAL, TARGET :: a(*)
    INTEGER, OPTIONAL :: ia(n+1)
    INTEGER, OPTIONAL, TARGET :: ja(*)
    !
    COMPLEX(kind(0.0d0)), ALLOCATABLE, SAVE :: ac(:,:)
    INTEGER , ALLOCATABLE, SAVE :: ipiv(:)
    REAL(kind(0.0d0)), SAVE :: iflop
    INTEGER :: i,kk,iext
    INTEGER  :: ierr
    INTEGER , parameter :: IONE=1
    !
    ierr=0
    IF (ijb == -2) THEN
       !
       DEALLOCATE (ac,ipiv)
       !
    ELSE IF (ijb == 1) THEN
       !
       ALLOCATE (ac(n,n),ipiv(n))
       memi=memi+n
       memr=memr+n*n
       memax=MAX(memax,memr+memi*rlenilen)
       ac=0.0d0
       DO i=1,n
          DO kk=ia(i),ia(i+1)-1
             ac(i,ja(kk))=a(kk)
          END DO
       END DO
       CALL ZGETRF(N,N,ac,N,ipiv,ierr)
       IF (ierr /= 0) THEN
          WRITE(iout, *) ' FATAL ERROR in GETRF: ierror=',ierr
          STOP
       END IF
       iflop=2*dble(n)**2-n
       flop=(2*1.0d0)/(3*1.0d0)*(dble(n)**3)
       !
    ELSE IF (ijb == 2) THEN
       !
       CALL ZGETRS('N',N,IONE,ac,N,ipiv,f,N,ierr)
       IF (ierr /= 0) THEN
          WRITE(iout, *) ' FATAL ERROR in GETRS: ierror=',ierr
          STOP
       END IF
       flop=flop+iflop
       !
    END IF
    !
    RETURN
  END SUBROUTINE zagmg_LAPACK
!-----------------------------------------------------------------------
  SUBROUTINE zagmg_MUMPSseq(n,f,ijob,flop,a,ja,ia)
    USE zagmg_mem
    IMPLICIT NONE
!!!!!!!
!
!
      TYPE ZMUMPS_ROOT_STRUC
        SEQUENCE
        INTEGER MBLOCK, NBLOCK, NPROW, NPCOL
        INTEGER MYROW, MYCOL
        INTEGER ROOT_SIZE, TOT_ROOT_SIZE
        INTEGER :: CNTXT_BLACS, truc
        INTEGER, DIMENSION(:), POINTER :: RG2L_ROW
        INTEGER, DIMENSION(:), POINTER :: RG2L_COL
        INTEGER , DIMENSION(:), POINTER :: IPIV
        INTEGER, DIMENSION( 9 ) :: DESCRIPTOR, DESCB
        LOGICAL yes, gridinit_done
        INTEGER LPIV, brol
!       Used to access Schur easily from root structure
        COMPLEX*16, DIMENSION(:), POINTER :: SCHUR_POINTER
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD, machin
!
!      Data for nullspace/QR
!
        COMPLEX*16, DIMENSION(:), POINTER :: QR_TAU
        DOUBLE PRECISION     QR_RCOND
!
!      Givens rotations
!
        INTEGER MAXG, GIND
        COMPLEX*16, DIMENSION(:),POINTER::GROW, GCOS, GSIN
!
!      RRRLU data
!
        INTEGER ELG_MAX,NULL_MAX
        INTEGER ELIND,EUIND,NLUPDATE,NUUPDATE
        INTEGER,DIMENSION(:),POINTER::PERM_ROW,PERM_COL
        INTEGER,DIMENSION(:),POINTER::ELROW, EUROW, PTREL, PTREU
        COMPLEX*16, DIMENSION(:), POINTER :: ELELG, EUELG, DL
!
      END TYPE ZMUMPS_ROOT_STRUC
      TYPE ZMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters
! for the interface to the user, plus internal
! information
!
! *****************
! INPUT PARAMETERS
! *****************
!    -----------------
!    MPI Communicator
!    -----------------
        INTEGER COMM
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite,
!        SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
        INTEGER SYM, PAR
        INTEGER JOB
!    --------------------
!    Order of Input matrix
!    --------------------
        INTEGER N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
        INTEGER NZ
        COMPLEX*16, DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN, JCN
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA, pad0
!
!       ------------------------------------
!       Case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        INTEGER NZ_loc, pad1
        INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
        COMPLEX*16, DIMENSION(:), POINTER :: A_loc, pad2
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
        INTEGER NELT, pad3
        INTEGER, DIMENSION(:), POINTER :: ELTPTR
        INTEGER, DIMENSION(:), POINTER :: ELTVAR
        COMPLEX*16, DIMENSION(:), POINTER :: A_ELT, pad4
!
!    ---------------------------------------------
!    Symmetric permutation :
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data
! ******************
!    --------------------------------------------------------
!    RHS / SOL_LOC
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        COMPLEX*16, DIMENSION(:), POINTER :: RHS, REDRHS
        COMPLEX*16, DIMENSION(:), POINTER :: RHS_SPARSE
        COMPLEX*16, DIMENSION(:), POINTER :: SOL_LOC
        INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
        INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
        INTEGER, DIMENSION(:), POINTER :: ISOL_LOC
        INTEGER LRHS, NRHS, NZ_RHS, LSOL_LOC, LREDRHS
        INTEGER pad5
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
        INTEGER ICNTL(40)
        INTEGER INFO(40)
        INTEGER INFOG(40)
        DOUBLE PRECISION COST_SUBTREES
        DOUBLE PRECISION CNTL(15)
        DOUBLE PRECISION RINFO(20)
        DOUBLE PRECISION RINFOG(20)
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!       SYM_PERM: Symmetric permutation
!       UNS_PERM: Column permutations (optionnal)
!    ---------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
!
!    -------------------------------------
!    Case of distributed matrix on entry:
!    ZMUMPS potentially provides mapping
!    -------------------------------------
        INTEGER, DIMENSION(:), POINTER :: MAPPING
!
!    -------------------------------
!    Deficiency and null space basis
!    -------------------------------
        COMPLEX*16, DIMENSION(:,:), POINTER :: NULL_SPACE
        INTEGER Deficiency, pad6
!    -----
!    Schur
!    -----
        INTEGER NPROW, NPCOL, MBLOCK, NBLOCK
        INTEGER SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER SIZE_SCHUR
        INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
        COMPLEX*16, DIMENSION(:), POINTER :: SCHUR
        COMPLEX*16, DIMENSION(:), POINTER :: SCHUR_CINTERFACE
!    --------------
!    Version number
!    --------------
        CHARACTER(LEN=16) VERSION_NUMBER
!    -----------
!    Out-of-core
!    -----------
        CHARACTER(LEN=256) :: OOC_TMPDIR
        CHARACTER(LEN=64) :: OOC_PREFIX
!    ------------------------------------------
!    To save the matrix in matrix market format
!    ------------------------------------------
        CHARACTER(LEN=256) WRITE_PROBLEM
!
!
! **********************
! INTERNAL Working data
! *********************
        INTEGER INST_Number
!       For MPI
        INTEGER COMM_NODES, MYID_NODES, COMM_LOAD
        INTEGER  MYID, NPROCS, NSLAVES
        INTEGER ASS_IRECV
        INTEGER, DIMENSION(:), POINTER :: POIDS
        INTEGER LBUFR
        INTEGER LBUFR_BYTES
        INTEGER, DIMENSION(:), POINTER ::  BUFR
!       For analysis/facto/solve phases
        INTEGER MAXIS1, pad7
        INTEGER KEEP(500)
        INTEGER*8 KEEP8(150)
!       IS is used for the factors + workspace for contrib. blocks
        INTEGER, DIMENSION(:), POINTER :: IS
!       is1 (maxis1) contains working arrays computed
!       and used only during analysis
        INTEGER, DIMENSION(:), POINTER :: IS1
!       The following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        INTEGER LNA
        INTEGER NBSA
        INTEGER,POINTER,DIMENSION(:)::STEP, NE_STEPS, ND_STEPS
!  Info for pruning tree
        INTEGER,POINTER,DIMENSION(:)::Step2node
!  ---------------------
        INTEGER,POINTER,DIMENSION(:)::FRERE_STEPS, DAD_STEPS
        INTEGER,POINTER,DIMENSION(:)::FILS, PTRAR, FRTPTR, FRTELT
        INTEGER,POINTER,DIMENSION(:)::NA, PROCNODE_STEPS
!       The two pointer arrays computed in facto and used by the solve
!          (except the factors) are PTLUST_S and PTRFAC.
        INTEGER, DIMENSION(:), POINTER :: PTLUST_S
        INTEGER(8), DIMENSION(:), POINTER :: PTRFAC
!       main real working arrays for factorization/solve phases
        COMPLEX*16, DIMENSION(:), POINTER :: S
!       Information on mapping
        INTEGER, DIMENSION(:), POINTER :: PROCNODE
!       Input matrix ready for numerical assembly
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
        INTEGER, DIMENSION(:), POINTER :: INTARR
        COMPLEX*16, DIMENSION(:), POINTER :: DBLARR
!       Element entry: internal data
        INTEGER NELT_LOC, LELTVAR, NA_ELT, pad8
        INTEGER, DIMENSION(:), POINTER :: ELTPROC
!       Candidates and node partitionning
        INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
        INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
        INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
        INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE
        LOGICAL, DIMENSION(:), POINTER :: I_AM_CAND
!       For heterogeneous architecture
        INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!       Compressed RHS
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP
        COMPLEX*16, DIMENSION(:), POINTER :: RHSCOMP
!       For C interface
!   Info on the subtrees to be used during factorization
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: MEM_SUBTREE
        INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
        INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
        INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: COST_TRAV
        INTEGER NBSA_LOCAL, zwave
        INTEGER(8) :: MAX_SURF_MASTER
        INTEGER :: LWK_USER, zozo
        COMPLEX*16, DIMENSION(:), POINTER :: WK_USER
!    For simulating parallel out-of-core stack.
        DOUBLE PRECISION, DIMENSION(:),POINTER ::CB_SON_SIZE
!   Instance number used/managed by the C/F77 interface
        INTEGER INSTANCE_NUMBER
!    OOC management data that must persist from factorization to solve.
        INTEGER OOC_MAX_NB_NODES_FOR_ZONE
        INTEGER, DIMENSION(:,:),   POINTER :: OOC_INODE_SEQUENCE
        INTEGER(8),DIMENSION(:,:), POINTER :: OOC_SIZE_OF_BLOCK
        INTEGER*8, DIMENSION(:,:),   POINTER :: OOC_VADDR
        INTEGER,DIMENSION(:), POINTER :: OOC_TOTAL_NB_NODES
        INTEGER,DIMENSION(:), POINTER :: OOC_NB_FILES
        CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES
        INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
!    Indices of nul pivots
        INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Internal control array
        DOUBLE PRECISION DKEEP(30)
!    Array needed to manage additionnal candidate processor
        INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC
!   ------------------------
!   Root structure(internal)
!   ------------------------
        TYPE (ZMUMPS_ROOT_STRUC) :: root
      END TYPE ZMUMPS_STRUC
!!!!!!!
    INTEGER :: n,ijob
    COMPLEX(kind(0.0d0)), TARGET :: f(n)
    REAL(kind(0.0d0)) :: flop
    COMPLEX(kind(0.0d0)), OPTIONAL, TARGET :: a(*)
    INTEGER, OPTIONAL :: ia(n+1)
    INTEGER, OPTIONAL, TARGET :: ja(*)
    !
    REAL(kind(0.0d0)) , SAVE :: iflop
    TYPE(ZMUMPS_STRUC), SAVE :: mumps_par
    INTEGER :: ierr, i, j, k
    !
    IF (ijob == -2) THEN
       !
       mumps_par%JOB = -2
       CALL ZAGMG_MUMPS(mumps_par)
       !
    ELSE IF (ijob == 1) THEN
       !
       mumps_par%COMM = 0
       mumps_par%JOB = -1
       mumps_par%SYM = 0
       mumps_par%PAR = 1
       CALL ZAGMG_MUMPS(mumps_par)
       mumps_par%ICNTL(2)=-1
       mumps_par%ICNTL(3)=-1
       mumps_par%ICNTL(4)=0
       mumps_par%ICNTL(14)=80
       mumps_par%ICNTL(18)=3
       mumps_par%ICNTL(24)=1
       mumps_par%NZ_loc=ia(n+1)-1
       IF (transint) THEN
          ALLOCATE( mumps_par%JCN_loc(mumps_par%NZ_loc) )
          do i=1,n
             do j=ia(i),ia(i+1)-1
                mumps_par%JCN_loc(j)=i
             end do
          end do
          mumps_par%IRN_loc => ja(1:mumps_par%NZ_loc)
       ELSE
          ALLOCATE( mumps_par%IRN_loc(mumps_par%NZ_loc) )
          do i=1,n
             do j=ia(i),ia(i+1)-1
                mumps_par%IRN_loc(j)=i
             end do
          end do
          mumps_par%JCN_loc => ja(1:mumps_par%NZ_loc)
       END IF
       memi=memi+mumps_par%NZ_loc
       memax=MAX(memax,memr+memi*rlenilen)
       mumps_par%A_loc => a(1:mumps_par%NZ_loc)
       !
       mumps_par%N=n
       mumps_par%JOB = 4
       CALL ZAGMG_MUMPS(mumps_par)
       flop=mumps_par%RINFO(3)
       i=mumps_par%INFO(27)
       IF (i .GT. 0) THEN
          iflop=4*dble(i)-n
       ELSE
          iflop=-4.0d6*dble(i)-n
       END IF
       IF (transint) THEN
          DEALLOCATE(mumps_par%JCN_loc)
          NULLIFY(mumps_par%IRN_loc,mumps_par%A_loc)
       ELSE
          DEALLOCATE(mumps_par%IRN_loc)
          NULLIFY(mumps_par%JCN_loc,mumps_par%A_loc)
       END IF
       memi=memi-mumps_par%NZ_loc
       !
    ELSE IF (ijob == 2) THEN
       !
       mumps_par%JOB = 3
       mumps_par%RHS => f(1:n)
       CALL ZAGMG_MUMPS(mumps_par)
       flop=flop+iflop
    END IF
    !
    RETURN
  END SUBROUTINE zagmg_MUMPSseq
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!
END MODULE zagmg_ALLROUTINES
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! MAIN DRIVER
!-----------------------------------------------------------------------
  SUBROUTINE zagmg( n,a,ja,ia,f,x,ijob,iprint,nrest,iter,tol )
    USE zagmg_mem
    USE zagmg_ALLROUTINES
    IMPLICIT NONE
    INTEGER    :: n,ia(n+1),ja(*),ijob,iprint,nrest,iter
    COMPLEX (kind(0.0d0)) :: a(*),f(n),x(n)
    REAL (kind(0.0d0)) :: tol
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Arguments
!  =========
!
!  N       (input) INTEGER.
!          The dimension of the matrix.
!
!  A       (input/output) COMPLEX (kind(0.0d0)). Numerical values of the matrix.
!  IA      (input/output) INTEGER. Pointers for every row.
!  JA      (input/output) INTEGER. Column indices.
!
!              AGMG ASSUMES THAT ALL DIAGONAL ENTRIES ARE POSITIVE
!
!          Detailed description of the matrix format
!
!              On input, IA(I), I=1,...,N, refers to the physical start
!              of row I. That is, the entries of row I are located
!              in A(K), where K=IA(I),...,IA(I+1)-1. JA(K) carries the
!              associated column indices. IA(N+1) must be defined in such
!              a way that the above rule also works for I=N (that is,
!              the last valid entry in arrays A,JA should correspond to
!              index K=IA(N+1)-1). According what is written
!              above, AGMG assumes that some of these JA(K) (for
!              IA(I)<= K < IA(I+1) ) is equal to I with corresponding
!              A(K) carrying the value of the diagonal element,
!              which must be positive.
!
!              A,IA,JA are "output" parameters because on exit the
!              entries of each row may occur in a different order (The
!              matrix is mathematically the same, but stored in
!              different way).
!
!  F       (input/output) COMPLEX (kind(0.0d0)).
!          On input, the right hand side vector f.
!          Overwritten on output.
!          Significant only if IJOB==0, 2, 3, 10, 12, 100, 102, 110, 112
!
!  X       (input/output) COMPLEX (kind(0.0d0)).
!          On input and if IJOB== 10, 12, 110, 112: initial guess
!             (for other values of IJOB, the default is used: the zero vector).
!          On output, the computed solution.
!
! IJOB     (input) INTEGER. Tells AGMG what has to be done.
!          0: performs setup + solve + memory release, no initial guess
!         10: performs setup + solve + memory release, initial guess in x(1:n)
!          1: performs setup only
!             (preprocessing: prepares all parameters for subsequent solves)
!          2: solves only (based on previous setup), no initial guess
!         12: solves only (based on previous setup), initial guess in x(1:n)
!          3: the vector returned in x(1:n) is not the solution of the linear
!                 system, but the result of the action of the multigrid
!                 preconditioner on the right hand side in f(1:n)
!         -1: erases the setup and releases internal memory
!
!   IJOB == 100,110,101,102,112: same as, respectively, IJOB==0,10,1,2,12
!       but, use the TRANSPOSE of the input matrix in A, JA, IA.
!
!   !!! IJOB==2,3,12,102,112 require that one has previously called AGMG
!       with IJOB==1 or IJOB==101
!
!   !!! (change with respect to versions 2.x) !!!
!       The preconditioner defined when calling AGMG
!         with IJOB==1 or IJOB==101 is entirely kept in internal memory.
!       Hence the arrays A, JA and IA are not accessed upon subsequent calls
!         with IJOB==3.
!       Upon subsequent calls with IJOB==2,12,102,112, a matrix needs to
!            be supplied in arrays A, JA, IA, but it will be used to
!            perform matrix vector product within the main iterative
!            solution process (and only for this).
!            Hence the system is solved with this matrix which
!            may differ from the matrix in A, JA, IA that was supplied
!            upon the previous call with IJOB==1 or IJOB==101;
!            then AGMG attempts to solve a linear system with the "new"
!            matrix (supplied when IJOB==2,12,102 or 112) using the
!            preconditioner set up for the "old" one (supplied when
!            IJOB==1 or 101).
!         The same remarks apply to IJOB >= 100 or not: the value IJOB==1
!            or 101 determines whether the preconditioner set up and stored
!            in internal memory is based on the matrix or its transpose;
!            the value IJOB==2,12 or 102,112 is used to determine whether
!            the linear system to be solved is with the matrix or its
!            transpose, independently of the set up.
!            Hence one may set up a preconditioner for a matrix and use it
!            for its transpose.
!       These functionalities (set up a preconditioner and use it for another
!            matrix) are provided for the sake of generality but should be
!            used with care; in general, set up is fast with AGMG and hence
!            it is recommended to rerun it even if the matrix changes only
!            slightly.
!
! IPRINT   (input) INTEGER.
!              Indicates the unit number where information is to be printed
!              (N.B.: 5 is converted to 6). If nonpositive, only error
!              messages are printed on standard output.
!
! NREST    (input) INTEGER.
!             Restart parameter for GCR (an implementation of GMRES)
!             Nonpositive values are converted to NREST=10 (default)
!
! !!  If NREST==1, Flexible CG is used instead of GCR (when IJOB==0,10,2,
!             12,100,110,102,112) and also (IJOB==0,1,100,101) performs some
!             simplification during the set up based on the assumption
!             that the matrix supplied in A, JA, IA is symmetric (there is
!             then no more difference between IJOB==1 and IJOB==101).
!
! !!!  NREST=1 Should be used if and only if the matrix is really SYMMETRIC
! !!!         (and positive definite).
!
!  ITER    (input/output) INTEGER.
!          On input, the maximum number of iterations. Should be positive.
!          On output, actual number of iterations.
!            If this number of iteration was insufficient to meet convergence
!            criterion, ITER will be returned negative and equal to the
!            opposite of the number of iterations performed.
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112
!
!  TOL     (input) REAL (kind(0.0d0)).
!          The tolerance on residual norm. Iterations are stopped whenever
!               || A*x-f || <= TOL* || f ||
!          Should be positive and less than 1.0
!          Significant only if IJOB==0, 2, 10, 12, 100, 102, 110, 112
!
!!!!! Remark !!!! Except insufficient number of iterations to achieve
!                 convergence (characterized by a negative value returned
!                 in ITER), all other detected errors are fatal and lead
!                 to a STOP statement.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
    INTEGER :: nza
    REAL(kind(0.0d0)), SAVE :: cputm=0.0d0,eltm=0.0d0,flop=0.0d0
    LOGICAL, SAVE :: preprocessed=.FALSE.,solve=.FALSE.
    INTEGER :: i,init,ijb,k
    REAL(kind(0.0d0)) :: cputmp,eltmp,memh
    REAL(kind(0.0d0)) :: resid,wupd
    INTEGER, POINTER, DIMENSION(:) :: iw
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: ascal
    COMPLEX(kind(0.0d0)), POINTER, DIMENSION(:) :: w
    COMPLEX(kind(0.0d0)) :: adum(1),fdum(1)
    INTEGER :: listrank(1),ifirstlistrank
!
    wfo=.TRUE.
    iout=iprint
    IF (iprint <= 0) THEN
       iout=6
       wfo=.FALSE.
    ELSE IF (iprint == 5) THEN
       iout=6
    END IF
    ijb=mod(ijob,100)
    trans=ijob.GE.100
    spd=nrest.EQ.1
    transint=trans.AND.(.NOT.spd)
!
    wff=wfo.AND.(IRANK<=0)
!
    IF (MOD(ijb,10) >= 2 .AND. .NOT.preprocessed) THEN
       WRITE (iout,1001) IRANK,ijob
       STOP
    END IF
    nza=ia(n+1)-ia(1)
!
    IF (ijb < 0 .AND. solve) GOTO 450
    IF (ijb < 0) GOTO 500
    IF (MOD(ijb,10) >= 2) GOTO 300
    IF (preprocessed) THEN
       CALL zagmg_relmem
       IF (.NOT.allzero)                        &
               CALL zagmg_MUMPSseq(nn(nlev),fdum,-2,flop)
       preprocessed=.FALSE.
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
    END IF
    CALL zagmg_mestime(-1,0.0d0,0.0d0)
!
!
    IF (HUGE(n) > 1.0e10) THEN
       rlenilen=dble(8)/dble(complex_len)
    ELSE
       rlenilen=dble(4)/dble(complex_len)
    END IF
!
!
!
    nlev=0
    imult=1
    gcrin=.FALSE.
    IF (wfo) THEN
       WRITE(iout,900) IRANK
    END IF
!
!
    ALLOCATE(dt(1)%idiag(n+1),w(n),iw(n))
    CALL zagmg_partroword(n,a,ja,ia,dt(1)%idiag,w,iw)
    DEALLOCATE(w,iw)
    !
    !
    !
    ALLOCATE(scald(n),ascal(ia(n+1)-ia(1)))
    memr=memr+n+ia(n+1)-ia(1)
    CALL zagmg_diagscal(n,a,ja,ia,dt(1)%idiag,ascal)
    CALL zagmg_setupL1(n,ascal,ja,ia,listrank,ifirstlistrank)
    DEALLOCATE(ascal)
    memr=memr-(ia(n+1)-ia(1))
    !
    !
    CALL zagmg_smoothsetup
    preprocessed=.TRUE.
    memh=memr+memi*rlenilen
!
    CALL zagmg_mestime(1,cputmp,eltmp)
    IF(wfo)THEN
       IF (nlev > 1) THEN
          WRITE(iout,960) IRANK, memax/nza, complex_len    &
                        , memax*complex_len/(2**20)
          IF(MOD(ijb,10) == 1) THEN
             WRITE(iout,961) IRANK, memh/nza, complex_len  &
                           , memh*complex_len/(2**20)
          END IF
       END IF
   !CPU_TIME: next line may be uncommented if implemented
       WRITE(iout,997) IRANK,eltmp
       WRITE(iout,'()')
    END IF
!
    IF (MOD(ijb,10) == 1) RETURN
!
!
300 CONTINUE
    CALL zagmg_mestime(-2,0.0d0,0.0d0)
    resid=tol
    nrst=nrest
    init=MAX(IJB/10,0)
    IF (nrst <= 0) nrst=10
   !
    IF (MOD(ijb,10) >= 3) THEN
       IF (wfo) THEN
          WRITE(iout,901) IRANK
       END IF
       CALL zagmg_applyprec(n,f,x,a,ja,ia,flop)
       CALL zagmg_mestime(2,cputmp,eltmp)
       cputm=cputm+cputmp
       eltm=eltm+eltmp
       solve=.TRUE.
       RETURN
   !
   !
    ELSE IF (nrst > 1) THEN
   !
       CALL zagmg_GCR(n,f,x,iter,resid,a,ja,ia,init,flop)
   !
   !
    ELSE
   !
       CALL zagmg_FlexCG(n,f,x,iter,resid,a,ja,ia,init,flop)
   !
   !
    END IF
!
    CALL zagmg_mestime(2,cputm,eltm)
    solve=.FALSE.
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          DO i=2,nlev-1
             WRITE(iout,955) i,kstat(2,i-1),kstat(2,i),         &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       IF (ITER .NE. 0) THEN
         flop=flop*log(1.0d0/10)/log(RESID)
         WRITE(iout,952) IRANK,flop/dble(2*nza)
       END IF
       WRITE(iout,962) IRANK,(memh+mritr)/nza, complex_len   &
                     , (memh+mritr)*complex_len/(2**20)
   !CPU_TIME: next line may be uncommented if implemented
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
    IF (MOD(ijb,10) > 0) THEN
       solve=.FALSE.
       eltm=0.0d0
       cputm=0.0d0
       flop=0.0d0
       kstat=0
       RETURN
    ELSE
       GOTO 500
    END IF
450 CONTINUE
    IF (wfo) THEN
       WRITE(iout,'()')
       WRITE(iout,990) IRANK
       WRITE(iout,'()')
    END IF
    IF (wfo) THEN
       IF (wff .AND. iter.NE.0) THEN
          DO i=2,nlev-1
             WRITE(iout,955) i,kstat(2,i-1),kstat(2,i),         &
                  dble(kstat(2,i))/dble(kstat(2,i-1)),kstat(1,i)
          END DO
       END IF
       WRITE(iout,'()')
       WRITE(iout,953) IRANK,flop/dble(2*nza)
       WRITE(iout,963) IRANK,(memh+mritr)/nza, complex_len   &
                     , (memh+mritr)*complex_len/(2**20)
   !CPU_TIME: next line may be uncommented if implemented
       WRITE(iout,999) IRANK,eltm
       WRITE(iout,'()')
    END IF
!
!
!
500 CONTINUE
    CALL zagmg_relmem
    IF (.NOT.allzero)                              &
            CALL zagmg_MUMPSseq(nn(nlev),fdum,-2,flop)
    preprocessed=.FALSE.
    solve=.FALSE.
    eltm=0.0d0
    cputm=0.0d0
    flop=0.0d0
    kstat=0
    IF (wfo) THEN
       WRITE (iout,902) IRANK
       WRITE (iout,903) IRANK
    END IF
!
    CONTINUE
!
    RETURN
900 FORMAT(i3,'*ENTERING AGMG **********************************',&
         '***************************')
901 FORMAT(i3,'*ONE APPLICATION OF AGMG PRECONDITIONER')
902 FORMAT(i3,            &
  ' (*) 1 work unit represents the cost of 1 (fine grid) resiudal evaluation ')
903 FORMAT(i3,'*LEAVING AGMG * (MEMORY RELEASED) ***************',&
         '***************************')
952 FORMAT(i3,'*','       Number of work units:',f9.2,             &
              ' per digit of accuracy (*)')
953 FORMAT(i3,'*',' Total number of work units:',f9.2, '(*)')
955 FORMAT('****     level',i2,'   #call=',i6,'   #cycle=',i6,    &
         '   mean=',f7.2,'    max=',i3)
960 FORMAT(  i3,'*','         memory used (peak):',f9.2,        &
         ' cplx(',i2,') words per nnz (',f8.2,' Mb)')
961 FORMAT(  i3,'*','     memory still allocated:',f9.2,        &
         ' cplx(',i2,') words per nnz (',f8.2,' Mb)')
962 FORMAT(  i3,'*','   memory used for solution:',f9.2,        &
         ' cplx(',i2,') words per nnz (',f8.2,' Mb)')
963 FORMAT(  i3,'*','                memory used:',f9.2,        &
         ' cplx(',i2,') words per nnz (',f8.2,' Mb)')
990 FORMAT(i3,'*GLOBAL STATISTICS for preconditioner application:')
996 FORMAT(i3,'*','           Setup time (CPU):   ',1PE10.2,     &
         ' seconds')
997 FORMAT(i3,'*','       Setup time (Elapsed):   ',1PE10.2,     &
         ' seconds')
998 FORMAT(i3,'*','        Solution time (CPU):   ',1PE10.2,     &
         ' seconds')
999 FORMAT(i3,'*','    Solution time (Elapsed):   ',1PE10.2,     &
         ' seconds')
1001 FORMAT(i3,'*',' FATAL ERROR: setup not done: ijob=',i3, &
         ' is not allowed')
1002 FORMAT(i3,'*',' FATAL ERROR: ijob=',i3, &
         ' (i.e. >= 100: work with transpose) is not allowed in the || case')
  END SUBROUTINE zagmg
!-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!! END of MAIN DRIVER
!-----------------------------------------------------------------------
!------------------ END of source file ---------------------------------
