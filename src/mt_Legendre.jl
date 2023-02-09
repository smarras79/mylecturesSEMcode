function LegendreAndDerivativeAndQ!(Legendre::St_Legendre, nop::TInt, x::TFloat)
    
    """
             Evaluate by recursion, the Legendre polynomial of order p
             its Derivatives at coordinate x, and q:
     
              L_{p}  --> legendre of order p
              L'_{p} --> dlegendr of order p
              q  = L_{p+1}  -  L_{p-1}
              q' = L'_{p+1} -  L'_{p-1} 
     
              Algorithm 24 of Kopriva's book
     
              Simone Marras, October 2021
              
              Note that this algorithm looses precision at high enough nop
              if we intendt to use particularly large nop we should examine Yakimiw 1996 (Yassine) 
         """
    TFloat=Float64
    
    a   ::TFloat=0.0
    b   ::TFloat=0.0
    L   ::TFloat=0.0
    dL  ::TFloat=0.0
    Lp1 ::TFloat=0.0
    dLp1::TFloat=0.0
    Lm1 ::TFloat=0.0
    dLm1::TFloat=0.0
    Lm2 ::TFloat=0.0
    dLm2::TFloat=0.0
    
    
    #st_legendre Legendre;
    
    if (nop == 0) #Order 0 case
	Legendre.legendre  = 1.0
	Legendre.dlegendre = 0.0
        Legendre.q = x
        Legendre.dq = 1.0
    elseif (nop == 1)
	Legendre.legendre  = x
	Legendre.dlegendre = 1.0
        Legendre.q = 0.5*(3*x^2-2)-1
        Legendre.dq = 3*x
    else
	Lm2  = 1.0
	Lm1  = x
	dLm2 = 0.0
	dLm1 = 1.0
	
	#Construct Nth Order Legendre Polynomial
	for k=2:nop
            
	    a = TFloat((2.0*k - 1.0)/k)
	    b = TFloat((k - 1.0)/k)
            
	    L  = a*x*Lm1 - b*Lm2
	    dL = dLm2 + (2.0*k - 1.0)*Lm1
	    
	    a = TFloat((2.0*(k+1) - 1.0)/(k+1))
	    b = TFloat((k+1 - 1.0)/(k+1))
	    Lp1  = a*x*L - b*Lm1
	    dLp1 = dLm1 + (2.0*k - 1.0)*L
	    
	    Lm2 = Lm1
	    Lm1 = L

	    dLm2 = dLm1
	    dLm1 = dL
	end
        Legendre.legendre = L
        Legendre.dlegendre = dL
        Legendre.q = Lp1 - Lm2
        Legendre.dq = dLp1 - dLm2
    end
end
