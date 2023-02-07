using("mt_quadTrapezoid.jl")

function mt_FourierExpansion(N, Nq, Nx)

    #
    # N wave modes (1:∞)
    # Nq number of quadrature points
    # Nx number of coordinate points along the period domain
    #
    F = zeros(Float64, Nx)
    
    #       1
    # âo = ----∫_{a,b} u(x)cos(x)dx
    #       2π
    u(x) = 3.0/(5.0 - 4.0*cos(x))
    a0 = mt_quad_Trapezoid(u; a=0.0, b=2π, Nq)/(2π)
    
    for i = 1:Nx
        xi = 2π*i/Nx
        
        sumcos = sumsin = 0.0
        an     = bn     = 0.0
        for n = 1:N
            #       1
            # ân = ----∫_{a,b} u(x)cos(x)dx
            #      cn⋅π
            uxcos(x) = u(x)*cos(n*x)
            an = mt_quadTrapezoid(uxcos; a=0.0, b=2π, Nq)/π
            
            #       1
            # b̂n = ---∫_{a,b} u(x)sin(x)dx
            #       π
            uxsin(x) = u(x)*sin(n*x)
            bn = mt_quadTrapezoid(uxsin; a=0.0, b=2π, Nq)/π
            
            sumcos = sumcos + an*cos(n*xi)
            sumsin = sumsin + bn*sin(n*xi)
            
        end
        F[i] = a0 + sumcos + sumsin
    end
end

