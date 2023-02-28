#
# DSS
#
function DSS_matrix(SD::NSD_1D, QT::Exact, Me::AbstractArray, conn, nelem, npoin, N, T)
    
    M    = zeros(npoin, npoin)
    Minv = zeros(npoin, npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            for j=1:N+1
                J = conn[j,iel]
                M[I,J] = M[I,J] + Me[i,j,iel]                
            end
        end
    end
    Minv = inv(M)
    
    return M , Minv
end


function DSS_vector(SD::NSD_1D, QT::Inexact, Ae::AbstractArray, conn, nelem, npoin, N, T)

    V    = zeros(npoin)
    Vinv = zeros(npoin)
    
    for iel=1:nelem
        for i=1:N+1
            I = conn[i,iel]
            V[I] = V[I] + Ve[i,iel]
        end
    end
    Vinv = 1.0./V
    
    return V, Vinv
end
