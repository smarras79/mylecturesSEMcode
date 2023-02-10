using Plots

function mt_LagrangeBasis!()
    
    p = 2
    Q = 100
    
    (L, dLdξ) = mt_buildLagrange(; order=p, Nq=Q)
      
end

function mt_buildLagrange(; order=1, Nq=1)

    #
    # Nq: number of quadrature points
    #
    ξ = Array{Float64}(undef, 1)
    ξ = range(-1, 1,  order + 1)

    x = Array{Float64}(undef, 1)
    x = range(-1, 1, Nq)
    
    L = dLdx = Array{Float64}(undef, 2)
    L = dLdx = zeros(order+1, Nq)
    
    for l=1:Nq
        xl = x[l]

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
                
                #=ddL=1
                if (j != i)
                    for k=1:order+1
                        xk = ξ[k]
                        
                        #dL/dx
                        if (k !=i && k !=j)
                            ddL = ddL*(xl - xk)/(xi - xk)
                        end
                    end
                    dLdx[i, l] = dLdx[i, l] + ddL/(xi - xj)
                end=#
            end
        end
    end

    
    #Plot to PDF
    plt = plot()
    for p=1:order+1
        L_legend = string("L(x) p=", p)

        plt = plot!(x, L[p,:], label=[L_legend], lw=4)
        plt = scatter!(ξ, ξ.*0)
        display(plt)
    end
    #fname = string("tm_Lagrange.pdf")
    #savefig(fname)
    #savefig(plt, fname)
    
    return (L, dLdx)
end

mt_LagrangeBasis!()
