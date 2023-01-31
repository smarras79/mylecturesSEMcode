function quad_trapezoid(f; a, b, Nq)   
    Δx = (b - a)/Nq
    If = Δx * (f(a) + f(b))/2
    for k = 1:Nq - 1
        xk = (b - a) * k/Nq + a
        If = If + Δx*f(xk)
    end
    return If
end
