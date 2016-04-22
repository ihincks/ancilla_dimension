%%  InducedSchattenNormMod
%   This function is derived from the function InducedSchattenNorm in
%   QETLAB version 0.9 (http://qetlab.com).
%   The arguments are:
%       - Phi - a linear map taking square matrices to square matrices,
%               passed as a function
%       - PhiA - The adjoint of Phi, also passed as a function
%       - dim_in - the input dimension of the linear map Phi 
%                   (e.g. dim_in = n for nxn matrices)
%       - dim_out - the output dimension of Phi
%
%   This function computes a lower bound on the induced 1->1 norm of the map
%   Phi. Aside from taking the linear maps as functions (as opposed to
%   the representations used in QETLAB), this function is just a pared down
%   version of InducedSchattenNorm. See the documentation for
%   InducedSchattenNorm for an explanation of the algorithm used. The
%   output is the computed lower bound, nrm, along with the operator X
%   which achieves this value.

function [nrm,X] = InducedSchattenNormMod(Phi,PhiA,dim_in,dim_out)

% fix the optional inputs for my purposes
dim = [dim_in, dim_out];
tol = 10^(-5);

X = randn(dim(1),dim(1));
X = X/SchattenNorm(X,1); % normalize the starting matrix

% Preparation is done; now do the actual iteration.
it_err = 2*tol+1;
Y = Phi(X);
nrm = SchattenNorm(Y,1);

while it_err > tol
    % First, find the best left matrix Y, keeping the right matrix X fixed.
    [U,S,V] = svd(Y);
    S = diag(S); % only want the diagonal part of S
    
    S = S/max(S); % pre-process in this way first for numerical reasons
    S = S.^0; % this is the equality condition from the Schatten Holder inequality
    S = S/norm(S,Inf);
    
    Y = U*diag(S)*V'; % reconstruct the optimal Y from the new SVD that we just computed
    
    % Next, find the best right matrix X, keeping the left marix Y fixed.
    X = PhiA(Y);
    [U,S,V] = svd(X);
    S = diag(S); % only want the diagonal part of S
    
    [~,ind] = max(S);
    S = zeros(dim(1),1);
    S(ind) = 1;
    
    X = U*diag(S)*V'; % reconstruct the optimal X from the new SVD that we just computed
    
    % Check to see if we made any progress; if so, keep iterating.
    Y = Phi(X);
    new_nrm = SchattenNorm(Y,1);
    it_err = abs(new_nrm - nrm);
    nrm = new_nrm;
end