/*

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
*/
#include "mex.h"
#include "matrix.h"

/* computational subroutines */
extern void dag2l_twolev__(int *, int *, double *, int *, int *, int *,
			   int *, int *, double *, double *, double *,
			   double *, double *, int *);


/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    mwIndex *irs,*jcs;
    double  *as, *kap, *checkd, *targetcf, *fracnz, *trspos;
    double  *npas1, *sym1, *l1, *x , *iprint1;
    mwSize  nn, nsiz, ncol;
    int     *ja, *ia, *ind, n, nz, i, l, sym, npas, iprint;
  
    /* get pointers to input*/
    nn      = mxGetN(prhs[0]);
    as      = mxGetPr(prhs[0]);
    irs     = mxGetIr(prhs[0]);
    jcs     = mxGetJc(prhs[0]);
    l1      = mxGetPr(prhs[1]);
    sym1    = mxGetPr(prhs[2]);
    npas1   = mxGetPr(prhs[3]);
    kap     = mxGetPr(prhs[4]);
    checkd  = mxGetPr(prhs[5]);
    targetcf= mxGetPr(prhs[6]);
    fracnz  = mxGetPr(prhs[7]);
    trspos  = mxGetPr(prhs[8]);
    iprint1 = mxGetPr(prhs[9]);
    n       = nn;
    
    /* transform to integer */
    l=(int)*l1;
    sym=(int)*sym1;
    npas=(int)*npas1;
    iprint=(int)*iprint1;
  
    /* create workspace, get pointers to it*/
    nz = jcs[n];
    nsiz = nz*sizeof(int);
    ja  = mxMalloc(nsiz);
    nsiz = (n+1)*sizeof(int);
    ia  = mxMalloc(nsiz);
        for (i=0 ; i<=n ; i++) {
            ia[i] = jcs[i]+1;}
        for (i=0 ; i<jcs[n] ; i++) {
            ja[i] = irs[i]+1;}
    nsiz = n*sizeof(int);
    ind    = mxMalloc(nsiz);
    /* call the aggregation procedure */
    dag2l_twolev__(&l,&n,as,ja,ia,ind,&sym,&npas,kap,checkd,targetcf,fracnz,trspos,&iprint);

    /* create a new array and set the output pointer to it */
    ncol=1;
    plhs[0] = mxCreateDoubleMatrix(nn, ncol, mxREAL);
    x = mxGetPr(plhs[0]);
        for (i=0 ; i<n ; i++) {
            x[i] = ind[i];
        }
}
