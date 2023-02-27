#
# Element Mass and Differentiation matrices in 1D
#
function mt_build_element_matrices!(SD::NSD_1D, QT::Exact, ψ, dψdξ, ω, order, nelem, Δx, T::Float64)

    N  = order + 1    
    Nq = order + 2
    M = zeros(T, N, N, nelem)
    D = zeros(T, N, N, nelem)
    
    for iel=1:nelem
        Jac = Δx[iel]/2
        
        for iq=1:Nq
            for i=1:N
                for j=1:N
                    M[i,j,iel] = M[i,j,iel] + Jac*ω[iq]*ψ[i,iq]*ψ[j,iq]
                    D[i,j,iel] = D[i,j,iel] +     ω[iq]*ψ[i,iq]*dψdξ[j,iq]
                end
            end
        end
    end
    #show(stdout, "text/plain", M)
    
    return M, D
end

function mt_build_element_matrices!(SD::NSD_1D, QT::Inexact ψ, dψdξ, ω, order, nelem, Δx, T::Float64)
    
    N  = order + 1    
    Nq = order + 1 
    M = zeros(T, N,    nelem)
    D = zeros(T, N, N, nelem)

    for iel=1:nelem
        Jac = Δx[iel]/2
        
        for iq=1:Q+1
            for i=1:order+1
                M[i,iel] = M[i,iel] + Jac*ω[iq] #Store only the diagonal elements
                for j=1:order+1
                    D[i,j,iel] = D[i,j,iel] + ω[iq]*ψ[i,iq]*dψdξ[j,iq] #Sparse
                end
            end
        end
    end
    #show(stdout, "text/plain", D)
    
    return M, D
    
end
