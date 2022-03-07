

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
        program example_par
!
!  Solves the discrete Laplacian on the unit square by simple call to agmg.
!  The right-hand-side is such that the exact solution is the vector of all 1.
!  Uses a strip partitioning of the domain, with internal boundaries parallel
!       to the x direction.
!
        implicit none
        include 'mpif.h'
        real (kind(0d0)),allocatable :: a(:),f(:),x(:)
        integer,allocatable :: ja(:),ia(:),listrank(:)
        integer :: n,iter,iprint,nhinv,NPROC,IRANK,mx,my,ifirstlistrank,ierr
        real (kind(0d0)) :: tol
        character*10 filename
!
!       set inverse of the mesh size (feel free to change)
        nhinv=1000
!
!       maximal number of iterations
        iter=50
!
!       tolerance on relative residual norm
        tol=1.e-6
!
!  Initialize MPI
!
       call MPI_INIT(ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,ierr)
!
!       unit number for output messages (alternative: iprint=10+IRANK)
        iprint=10
        filename(1:8)='res.out_'
        write (filename(9:10),'(i2.2)') IRANK        ! processor dependent
        open(iprint,file=filename,form='formatted')
!
!    calculate local grid size
!
       mx=nhinv-1
       my=(nhinv-1)/NPROC
       if (IRANK < mod(nhinv-1,NPROC)) my=my+1
!
!       generate the matrix in required format (CSR)
!
!         first allocate the vectors with correct size
            N=mx*my
            allocate (a(5*N),ja(5*N),ia(N+1),f(N),x(n),listrank(2*mx))
!           external nodes connected with local ones on top and bottom
!           internal boundaries will receive numbers [N+1,...,N+2*mx]
            ifirstlistrank=N+1
!         next call subroutine to set entries
!           before, initialize listrank to zero so that entries
!           that do not correspond to a nonlocal variable present
!           in ja are anyway properly defined
         listrank(1:2*mx)=0
         call uni2dstrip(mx,my,f,a,ja,ia,IRANK,NPROC,listrank,ifirstlistrank)
!
!       call agmg
!         argument 5 (ijob)  is 0 because we want a complete solve
!         argument 7 (nrest) is 1 because we want to use flexible CG
!                            (the matrix is symmetric positive definite)
!
         call dagmgpar(N,a,ja,ia,f,x,0,iprint,1,iter,tol,         &
                       MPI_COMM_WORLD,listrank,ifirstlistrank)
!
!      uncomment the following to write solution on disk for checking
!
!        filename(1:8)='sol.out_'
!        write (filename(9:10),'(i2.2)') IRANK         ! processor dependent
!        open(11,file=filename,form='formatted')
!        write(11,'(e22.15)') x(1:n)
!        close(11)

      END program example_par
!----------------------------------------------------------------------
    subroutine uni2dstrip(mx,my,f,a,ja,ia,IRANK,NPROC,listrank,ifirstlistrank)
!
! Fill a matrix in CSR format corresponding to a constant coefficient
! five-point stencil on a rectangular grid
! Bottom boundary is an internal boundary if IRANK > 0, and
!    top boundary is an internal boundary if IRANK < NPROC-1
!
      implicit none
      real (kind(0d0)) :: f(*),a(*)
      integer :: mx,my,ia(*),ja(*),ifirstlistrank,listrank(ifirstlistrank:*)
      integer :: IRANK,NPROC,k,l,i,j
      real (kind(0d0)), parameter :: zero=0.0d0,cx=-1.0d0,cy=-1.0d0, cd=4.0d0
!
      k=0
      l=0
      ia(1)=1
      do i=1,my
        do j=1,mx
          k=k+1
          l=l+1
          a(l)=cd
          ja(l)=k
          f(k)=zero
          if(j < mx) then
             l=l+1
             a(l)=cx
             ja(l)=k+1
            else
             f(k)=f(k)-cx
          end if
          if(i < my) then
             l=l+1
             a(l)=cy
             ja(l)=k+mx
            else if (IRANK == NPROC-1) then
             f(k)=f(k)-cy             !real boundary
            else
             l=l+1                  !internal boundary (top)
             a(l)=cy                ! these external nodes are given the
             ja(l)=k+mx             !  numbers [mx*my+1,...,mx*(my+1)]
             listrank(k+mx)=IRANK+1 !Thus listrank(mx*my+1:mx*(my+1))=IRANK+1
          end if
          if(j > 1) then
             l=l+1
             a(l)=cx
             ja(l)=k-1
            else
             f(k)=f(k)-cx
          end if
          if(i >  1) then
             l=l+1
             a(l)=cy
             ja(l)=k-mx
            else if (IRANK == 0) then
             f(k)=f(k)-cy             !real boundary
            else
             l=l+1                    !internal boundary (bottom)
             a(l)=cy                  ! these external nodes are given the
             ja(l)=k+mx*(my+1)        ! numbers [mx*(my+1)+1,...,mx*(my+2)]
             listrank(k+mx*(my+1))=IRANK-1
                              !Thus listrank(mx*(my+1)+1:mx*(my+2))=IRANK-1
          end if
          ia(k+1)=l+1
        end do
      end do

      return
    end subroutine uni2dstrip
