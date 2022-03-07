function [P ind] = agtwolev(A,l,sym,npass,kappa,checkdd,targetcf,fracnz,trspos,verbose)
%
% Compute aggregation and associated prolongation matrix according to the
% algorithms in [2] (symmetric matrices) or [3] (general matrices) 
%
% [2] A. Napov and Y. Notay, An algebraic multigrid method with guaranteed
%    convergence rate, Report GANMN 10-03, Universite Libre de Bruxelles,
%    Brussels, Belgium, 2010.
%
% [3] Y. Notay, Aggregation-based algebraic multigrid for convection-diffusion
%    equations, Report GANMN 11-01, Universite Libre de Bruxelles, Brussels,
%    Belgium, 2011.
%
% USAGE:
% [P Ind] = agtwolev(A)
% [P Ind] = agtwolev(A,l,sym,npass,kappa,checkdd,targetcoarsefac
%                   ,fracnegrcsum,trspos,verbose)
%
% INPUT:
%
%  A: sparse matrix to deal with.
%
%  l:  l==1: (default) top level aggregation 
%            (priority rule according to CMK permutation)
%      l==2: further level aggregation (priority according to input ordering).
%
%  sym: sym==0: general matrix (default); sym==1: symmetric matrix
%
%  npass   is the maximal number of pairwise aggregation passes (default:2).
%
%  kappa is the threshold used to accept or not a tentative aggregate
%        (default: 10 (sym==0) or 8 (sym==1)).
%
%  checkdd is the threshold to keep outside aggregation nodes where
%         the matrix is strongly diagonally dominant (based on mean of row
%         and column);
%         In fact, uses the maximum of |checkdd| and of the default value 
%            according to kappa as indicated in [1,2]
%            (hence |checkdd| < 1 ensures that one uses this default value)
%         checkdd <0 : consider |checkdd|, but base the test on minus the
%               sum of offdiagonal elements, without taking the absolute value
%         (default: 0.5).
%
%  targetcoarsefac is the target coarsening factor (parameter tau in the main
%         coarsening algorithm in [1,2]): further pairwise aggregation passes
%         are omitted once the number of nonzero entries has been reduced by a
%         factor of at least targetcoarsefac (default: 2^npass).
%
%  fracnegrcsum: if, at some level, more than fracnegrcsum*nl nodes, 
%         where nl is the total number of nodes at that level, have
%         negative mean row and column sum, then the aggregation algorithm
%         of [2,3] is modified, exchanging all diagonal entries for the mean
%         row and column sum (that is, the algorithm is applied to a
%         modified matrix with mean row and colum sum enforced to be zero);
%         (default: 0.25).
%
%  trspos is a threshold: if a row has a positive offidiagonal entry larger
%         than trspos times the diagonal entry, the corresponding node is
%         transferred unaggregated to the coarse grid (default: 0.45).
%
% OUTPUT:
%  ind: vector of length N, where N is the numbers of rows & columns in A; 
%       ind(i) is the index of the aggregates to which i
%       belongs; ind(i)=0 iff i has been kept outside aggregation (see checkdd). 
%    P: associated prologation matrix; sparse N x Nc matrix, where Nc is the
%       number of aggregates;        
%       P(i,j)=1 iff ind(i)=j and P(i,j)=0 otherwise (i=1,N ; j=1,Nc).
%
% AGMG Copyright (C) 2011 Yvan NOTAY
% This function is part of AGMG software package distributed under the terms
%      of the GNU General Public License <http://www.gnu.org/licenses/>.
% Enter agtwolev('v') for detailed condition of use and version number.
% Check the web site <http://homepages.ulb.ac.be/~ynotay/agmg> for
%       up-to-date copies and documentation.
%
if (size(A,1) == 1 && size(A,2) ==1)
if (A=='v' || A=='V' || A=='version' || A=='Version')
      type agmglicense,return
end		            
end
if (nargin<1 || isempty(A) || any(class(A)~='double') || ~issparse(A) || ~isreal(A))
     error('MATLAB:agmg:NonSparseMatrix', ...
    'Input Matrix A must be a sparse nonempty array of double real');
end

if nargin < 2,    l   =  1; end
if nargin < 3,    sym   =0; end; 
if sym       ,    sym   =1; else, sym=0; end  
if nargin < 4,    npass =2; end
if nargin < 5,    if (sym), kappa=8; else, kappa=10; end; end
if nargin < 6,    checkdd = 0.5;  end
if nargin < 7,    targetcf= 2^npass;    end
if nargin < 8,    fracnz  = 0.25; end
if nargin < 9,    trspos  = 0.45; end
if nargin < 10,    verbose = 6;    end
if verbose   ,    verbose = 6;    end

ind = dmtlagtwolev(A,l,sym,npass,kappa,checkdd,targetcf,fracnz,trspos,verbose);

%%%%%%%%%%%%%%%%%%%%%%%
ind2 = find(ind);
P = sparse(ind2, ind(ind2), sign(ind(ind2)), size(A, 1), max(ind));
%%%%%%%%%%%%%%%%%%%%%%%

return;
