%%  InducedSchattenNormModHerm
%   This function is derived from the function InducedSchattenNorm in
%   QETLAB version 0.9 (http://qetlab.com).
%
%   This function is nearly the same as InducedSchattenNormMod (and so all of
%   the comments in the header of that function apply here). The only
%   difference with this function is that it searches over Hermitian
%   inputs, and as such uses spectral decompositions instead of singular
%   value decompositions. Note that the use of the spectral decomposition 
%   when searching over Hermitian inputs in this algorithm is only valid
%   when Phi is Hermiticity preserving (which the maps for which this
%   algorithm was written to investigate are). It is not checked whether
%   this is the case, so this function will not behave properly if
%   non-Hermiticity preserving maps are passed to it.

function [nrm,X] = InducedSchattenNormModHerm(Phi,PhiA,dim_in,dim_out)

% fix the optional inputs for my purposes
dim = [dim_in, dim_out];
tol = 10^(-5);

% In all other cases, we iterate to compute the induced Schatten norm.

% generate a random starting matrix X0, if appropriate
X = randn(dim(1),dim(1));
% make the starting matrix hermitian
X = X+ X';

X = X/SchattenNorm(X,1); % normalize the starting matrix

% Preparation is done; now do the actual iteration.
it_err = 2*tol+1;
Y = Phi(X);
nrm = SchattenNorm(Y,1);

while it_err > tol
    % Get the eigenvalue decomposition of Y
    [V,D] = eig(Y);
    
    % construct the U such that <U,Y> = ||Y||_1
    U=V*sign(D)*V';
    
    % Next, find the best right matrix X, keeping the left marix Y fixed.
    X = PhiA(U);
    
    [V,D] = eig(X);
    
    D = diag(D);
    [v,i] = max(abs(D));
    
    X = V(:,i)*V(:,i)';
    
    % Check to see if we made any progress; if so, keep iterating.
    Y = Phi(X);
    new_nrm = SchattenNorm(Y,1);
    it_err = abs(new_nrm - nrm);
    nrm = new_nrm;
end