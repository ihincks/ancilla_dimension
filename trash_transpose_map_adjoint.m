function sig = trash_transpose_map_adjoint( Y, n, m)
%TRASH_TRANSPOSE_MAP_ADJOINT 
%   This is the adjoint of the map performed in trash_transpose_map.
%   This function depends on the PartialTranspose and PartialTrace functions in
%   QETLAB version 0.9 (http://qetlab.com)
    
    YT = PartialTranspose(Y, 2, [2,n,m]);
    sig = (Swap(kron(eye(n),  PartialTrace(kron([1,0;0,0],eye(n*m))*YT,1,[2,n,m])), [1,2], [n,n,m]) ...
                        + kron(eye(n),  PartialTrace(kron([0,0;0,1],eye(n*m))*YT,1,[2,n,m])));

end

