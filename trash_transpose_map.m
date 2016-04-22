function sig = trash_channel( rho, n, m)
%TRASH_TRANSPOSE_MAP
%   For a given n and m, acts the linear map Psi_{n,2} with an identity
%   ancilla of dimension m on the matrix rho (i.e. rho is a (n x n x m) x
%   (n x n x m) matrix). 
%   This function depends on the PartialTranspose and PartialTrace functions in
%   QETLAB version 0.9 (http://qetlab.com)
    
    sig = (kron([1, 0; 0, 0], PartialTranspose(PartialTrace(rho,2,[n,n,m]), 1, [n,m])) +kron([0,0;0,1],  PartialTranspose(PartialTrace(rho,1,[n,n,m]),1,[n,m])));

end

