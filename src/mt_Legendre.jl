function LegendreAndDerivativeAndQ(p, x)
    
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
    
    ϕ   = dϕdx   = 0.0
    ϕp1 = dϕp1dx = 0.0
    ϕm1 = dϕm1dx = 0.0
    ϕm2 = dϕm2dx = 0.0
    if (p == 0) #Order 0 case
	ϕ    = 1.0
        dϕdx = 0.0
        q    = x
        dqdx = 1.0
    elseif (p == 1)
	ϕ    = x
        dϕdx = 1.0
        q    = 0.5*(3*x^2 - 2) - 1
        dqdx = 3*x
    else
	ϕm2    = 1.0
	ϕm1    = x
	dϕm2dx = 0.0
	dϕm1dx = 1.0
	
	#Construct pᵗʰ Order Legendre Polynomial
	for k=2:p

	    ϕ    = x*ϕm1*(2.0*k - 1.0)/k - ϕm2*(k - 1.0)/k
	    dϕdx = dϕm2dx + (2.0*k - 1.0)*ϕm1
            
	    ϕp1  = x*ϕ*(2.0*k + 1.0)/(k+1) - ϕm1*(k+1 - 1.0)/(k+1)
	    dϕp1 = dϕm1dx + (2.0*k - 1.0)*ϕ
	    
	    ϕm2 = ϕm1
	    ϕm1 = ϕ

	    dϕm2dx = dϕm1dx
	    dϕm1dx = dϕdx
	end
        q    = ϕp1 - ϕm2
        dqdx = dϕp1dx - dϕm2dx
    end

    return ϕ, dϕdx, q, dqdx
    
end
