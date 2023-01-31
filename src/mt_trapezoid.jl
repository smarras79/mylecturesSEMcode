function quad_trapezoid(f; a, b, Nq)
    
    Δx = (b - a)/Nq
    ∫ = Δx*(f(a) + f(b))/2
    for k = 1:Nq - 1
        xk = (b - a)*k/Nq + a
        ∫ = ∫ + Δx*f(xk)
    end
    return ∫
end
