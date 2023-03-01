using("mt_Legendre.jl")

#
# LG points and weights
#
function mt_LegendreGaussNodesAndWeights!(Legendre::St_Legendre, lg::St_lg, p::Int64)
    
    Δ::Float64 = 0.0
    
    NITER = 100
    TOL = 4*eps(Float64)
    if (p == 0)
        ξ[1] = 0.0
        ω[1] = 2.0
    elseif (p == 1)
        ξ[1] = -sqrt(1/3)
        ω[1] =  1.0
        ξ[2] = -ξ[1]
        ω[2] =  ω[1]
    else
        for j=0:floor(Int64,(p+1)/2)-1
            ξ[j+1] = - cospi((2*j+1)/(2*p+2))
            for k=0:NITER
                (ϕ, dϕdx, _, _) = mt_LegendreAndDerivativeAndQ!(p+1, ξ[j+1])
                Δ = -ϕ/dϕdx
                ξ[j+1] = ξ[j+1]+Δ
                if (abs(Δ)≤TOL)
                    break
                end
            end
            (ϕ, dϕdx, _, _) = mt_LegendreAndDerivativeAndQ!(p+1, ξ[j+1])
            ξ[p-j+1] = -ξ[j+1]
            ω[j+1]   = 2.0/(1-ξ[j+1]^2)/(dϕdx)^2
            ω[p-j+1] = ω[j+1]
        end
    end
    if (mod(p,2) == 0)
        (ϕ, dϕdx, _, _) = mt_LegendreAndDerivativeAndQ!(p+1, 0.0)
        ξ[Int64(p/2)] = 0
        ω[Int64(p/2)] = 2/dϕdx^2
    end

    return ξ, ω
    
end


#
# LLG points and weights
#
function mt_LegendreGaussLobattoNodesAndWeights!(Legendre::St_Legendre, lgl::St_llg, nop::Int64)
    
    NITER = 100
    TOL = 4*eps()
    
    ξ0 ::Float64=0.0
    ξ1 ::Float64=1.0
    ξP ::Float64=1.0
    ξj ::Float64=0.0
    ξj2::Float64=0.0
    L2 ::Float64=0.0
    
    ω0 ::Float64=0.0
    ω1 ::Float64=0.0
    ωP ::Float64=0.0
    
    Δ  ::Float64=0.0    
    for j=1:p+1
	ξ[j] = 0.0;
	ω[j] = 1.0;
    end
    
    if (p == 1)
	ξllg0      = -1.0
	ω0         =  1.0
	ξ1         =  1.0
	ω1         =   ω0
        
	ξ[1]   = ξ0
	ξ[p+1] = ξ1
	ω[1]   = ω0
	ω[p+1] = ω1
    else 
	ξ0         = -1.0
	ω0         =  Float64(2.0/(p*(p + 1)))
	ξP         =  1.0
	ωP         =  ω0
        
	ξ[1]   =  ξ0
	ξ[p+1] =  ξP
	ω[1]   =  ω0
	ω[p+1] =  ωP
	
        for jj = 2:floor(Int,(p + 1)/2) 
	    j = jj - 1
	    ξj = -cos((j + 0.25)*π/p - 3.0/(8.0*p*π*(j + 0.25)))
	    ξ[jj] = ξj;
            
            for k = 0:NITER
	        (ϕ, dϕdx, q, dqdx) = LegendreAndDerivativeAndQ!(p, ξj)
	        Δ = -q/dqdx
	        ξj    =  ξj + Δ
	        
	        if (abs(Δ) <= TOL*abs(ξj))
                    break
                end
	    end
            (ϕ, dϕdx, _, _) = LegendreAndDerivativeAndQ!(p, ξj)
	    ξ[jj]      =  ξj
	    ξ[p+1-j] = -ξj
	    xj2            =  ξj*ξj
	    L2             = ϕ^2
	    ω[jj]      = 2.0/(p*(p + 1.0)*L2)
	    ω[p+1-j] =  ω[jj]
            
        end
    end
    
    if (mod(p,2) == 0)
	(ϕ, dϕdx, _, _) = LegendreAndDerivativeAndQ!(p, 0.0);
	ξ[Int64(p/2)+1] = 0.0;
	
	L2           = ϕ^2
	ω[Int64(p/2)+1] = 2.0/(p*(p + 1.0)*L2);
    end

    return ξ, ω
    
end
