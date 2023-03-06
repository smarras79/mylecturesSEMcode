using("mt_Chebyshev.jl")

function ChebyshevGaussNodesAndWeights!(cg::St_lg, nop::Int64)
       for j=0:nop
        cg.ξ[j+1]=-cospi((2*j+1)/(2*nop+2))
        cg.ω[j+1]=π/(nop+1)
    end
end

function ChebyshevGaussLobattoNodesAndWeights!(cgl::St_cgl,nop::Int64)
    for j=0:nop
        cgl.ξ[j+1]=-cospi(j/nop)
        cgl.ω[j+1]=π/nop
    end
    cgl.ω[1]=cgl.ω[1]/2
    cgl.ω[nop+1]=cgl.ω[nop+1]/2
end
