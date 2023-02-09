function ChebyshevPolynomial!(Chebyshev::St_Chebyshev,nop::TInt,x::TFloat,Ks::TInt)
    """
          Evaluate by recursion the Chebyshev Polynomial
          and switch to direct evaluation when nop > Ks
          uses Algorithm 21 of Kopriva's book, 
          NOTE Ks depends should never exceed 70 regardless of machine architecture
          Yassine Tissaoui March 2022
       """
    T2::TFloat=0.0
    T1::TFloat=0.0
    T::TFloat=0.0
    if (nop == 0) #case for order 0
        Chebyshev.chebyshev = 1
    elseif (nop ==1) #case for order 1
        Chebyshev.chebyshev = x
    elseif (nop ≤ Ks)
        T2=1
        T1=x
        for j=2:nop
            T=2*x*T1 - T2
            T2=T1
            T1=T
        end
        Chebyshev.chebyshev = T
    else
        T=cos(nop*acos(x))
        Chebyshev.chebyshev = T
    end
end
