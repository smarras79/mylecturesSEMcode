function mt_LegendreAndDerivativeAndQ(p::Int64, x::Float64)
    
    """
             Evaluate by recursion, the Legendre polynomial of order p
             its Derivatives at coordinate x, and q:
     
              L_{p}  --> legendre of order p
              L'_{p} --> dlegendr of order p
              q  = L_{p+1}  -  L_{p-1}
              q' = L'_{p+1} -  L'_{p-1} 
     
              Algorithm 24 from Kopriva, D. Implementing spectral methods for PDEs, Springer
                   
              Note that this algorithm looses precision at high enough polynomial order p
              if we intendt to use particularly large nop we should examine Yakimiw 1996
         """

    ϕ   ::Float64=0.0
    dϕ  ::Float64=0.0
    ϕp1 ::Float64=0.0
    dϕp1::Float64=0.0
    ϕm1 ::Float64=0.0
    dϕm1::Float64=0.0
    ϕm2 ::Float64=0.0
    dϕm2::Float64=0.0
    
    #st_legendre Legendre;  
    if (p == 0) #Order 0 case
        ϕ    = 1.0
        dϕ   = 0.0
        q    = x
        dq   = 1.0
    elseif (p == 1)
        ϕ    = x
        dϕ   = 1.0
        q    = 0.5*(3*x^2 - 2) - 1
        dq   = 3*x        
    else
        ϕm2  = 1.0
	ϕm1  = x
	dϕm2 = 0.0
	dϕm1 = 1.0

        #Construct Nth Order Legendre Polynomial
	for k=2:p
             
            ϕ    = x*ϕm1*(2.0*k - 1.0)/k - ϕm2*(k - 1.0)/k
	    dϕ   = dϕm2 + (2.0*k - 1.0)*ϕm1
	    
	    ϕp1  = x*ϕ*(2.0*k + 1.0)/(k+1) - ϕm1*k/(k+1)
	    dϕp1 = dϕm1 + (2.0*k - 1.0)*ϕ
            
            q  =  ϕp1 -  ϕm1
            dq = dϕp1 - dϕm1
            
	    ϕm2 = ϕm1
	    ϕm1 = ϕ

	    dϕm2 = dϕm1
	    dϕm1 = dϕ
	end
    end
    
    return ϕ, dϕ, q, dq
    
end
