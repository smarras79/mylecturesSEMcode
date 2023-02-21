function mt_LagrangeBasis!()
    
    p = 2
    Q = 100
    
    (L, dLdξ) = mt_buildLagrange(; order=p, Nq=Q)
      
end

function mt_buildLagrange(; order=1, Nq=1)

    #------------------------------------------------------------------
    # A note on notation:
    #
    # L[i,l] is the digital equivalent of Lᵢ(xk) in the book
    # where 'i' indicates the iᵗʰ Lagrange basis that is evaluated
    # at the kᵗʰ quadrature point xk
    #------------------------------------------------------------------
    
    #
    # Nq: number of quadrature points
    #
    ξ = Array{Float64}(undef, 1)
    ξ = range(-1, 1,  order + 1) # NOTICE: use LegendreGaussNodesAndWeights() for LGL points instead of equally spaced.

    x = Array{Float64}(undef, 1)
    x = range(-1, 1, Nq)
    
    L = dLdx = Array{Float64}(undef, 2)
    L = dLdx = zeros(order+1, Nq)
    
    for l=1:Nq
        xl = ξq[l]

        for i=1:order+1
            
            xi        = ξ[i]
            L[i,l]    = 1.0
            dLdx[i,l] = 0.0
            for j=1:order+1
                xj = ξ[j]

                #L
                if (j != i)
                    L[i,l] = L[i,l]*(xl - xj)/(xi - xj)
                end
                
                ddL=1
                if (j != i)
                    for k=1:order+1
                        xk = ξ[k]
                        
                        #dL/dx
                        if (k !=i && k !=j)
                            ddL = ddL*(xl - xk)/(xi - xk)
                        end
                    end
                    dLdx[i, l] = dLdx[i, l] + ddL/(xi - xj)
                end
            end
        end
    end

    return (L, dLdx)
end
