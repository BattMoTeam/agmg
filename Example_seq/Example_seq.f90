

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
      program example_seq
!
!  Solves the discrete Laplacian on the unit square by simple call to agmg.
!  The right-hand-side is such that the exact solution is the vector of all 1.
!
        implicit none
        real (kind(0d0)),allocatable :: a(:),f(:),x(:)
        integer,allocatable :: ja(:),ia(:)
        integer :: n,iter,iprint,nhinv
        real (kind(0d0)) :: tol
!
!       set inverse of the mesh size (feel free to change)
        nhinv=500
!
!       maximal number of iterations
        iter=50
!
!       tolerance on relative residual norm
        tol=1.e-6
!
!       unit number for output messages: 6 => standard output
        iprint=6
!
!       generate the matrix in required format (CSR)
!
!         first allocate the vectors with correct size
            N=(nhinv-1)**2
            allocate (a(5*N),ja(5*N),ia(N+1),f(N),x(N))
!         next call subroutine to set entries
            call uni2d(nhinv-1,f,a,ja,ia)
!
!       call agmg
!         argument 5 (ijob)  is 0 because we want a complete solve
!         argument 7 (nrest) is 1 because we want to use flexible CG
!                            (the matrix is symmetric positive definite)
!
         call dagmg(N,a,ja,ia,f,x,0,iprint,1,iter,tol)
!
!      uncomment the following lines to write solution on disk for checking
!
!       open(10,file='sol.out',form='formatted')
!       write(10,'(e22.15)') f(1:n)
!       close(10)

      END program example_seq
!----------------------------------------------------------------------
      subroutine uni2d(m,f,a,ja,ia)
!
! Fill a matrix in CSR format corresponding to a constant coefficient
! five-point stencil on a square grid
!
      implicit none
      real (kind(0d0)) :: f(*),a(*)
      integer :: m,ia(*),ja(*)
      integer :: k,l,i,j
      real (kind(0d0)), parameter :: zero=0.0d0,cx=-1.0d0,cy=-1.0d0, cd=4.0d0
!
      k=0
      l=0
      ia(1)=1
      do i=1,m
        do j=1,m
          k=k+1
          l=l+1
          a(l)=cd
          ja(l)=k
          f(k)=zero
          if(j < m) then
             l=l+1
             a(l)=cx
             ja(l)=k+1
            else
             f(k)=f(k)-cx
          end if
          if(i < m) then
             l=l+1
             a(l)=cy
             ja(l)=k+m
            else
             f(k)=f(k)-cy
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
             ja(l)=k-m
            else
             f(k)=f(k)-cy
          end if
          ia(k+1)=l+1
        end do
      end do
!
      return
      end subroutine uni2D
