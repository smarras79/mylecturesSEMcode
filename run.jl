using Plots
using Printf

include("./src/mt_trapezoid.jl")

function mt_fourier(N, Nq, Nx, expansion_type)
    
    pu     = zeros(Float64, Nx)
    uexact = zeros(Float64, Nx)
    xcoord = zeros(Float64, Nx)
    
    if (expansion_type == "complex")
        #=  pu = zeros(Complex, Nx)
        for j=1:Nx
        u = 0.0
        for i=1:N
        un = 1.0/2.0^(abs(n[i]))
        #un = (2.0/π)/(1.0 - 4*n[i]*n[i])
        u = u + un*exp(n[i]*x[j]im)
        end
        pu[j] = u
        end=#
        nothing
    elseif (expansion_type == "trig")
        
        #       1
        # âo = ----∫_{a,b} u(x)cos(x)dx
        #       2π
        u(x) = 3.0/(5.0 - 4.0*cos(x))
        a0 = mt_quadTrapezoid(u; a=0.0, b=2π, Nq)/(2π)
        
        for i = 1:Nx
            #x-oordinates
            xi = 2π*i/Nx
            xcoord[i] = xi
            
            ∑ancosx = ∑bnsinx = 0.0
            an = bn = 0.0
            for n = 1:N
                #n-modes
                
                #       1
                # an = ----∫_{a,b} u(x)cos(x)dx
                #      cn⋅π
                uxcos(x) = u(x)*cos(n*x)
                an = mt_quadTrapezoid(uxcos; a=0.0, b=2π, Nq)/π
                
                #       1
                # bn = ---∫_{a,b} u(x)sin(x)dx
                #       π
                uxsin(x) = u(x)*sin(n*x)
                bn = mt_quadTrapezoid(uxsin; a=0.0, b=2π, Nq)/π
                
                ∑ancosx += an*cos(n*xi)
                ∑bnsinx += bn*sin(n*xi)
                
            end
            
            pu[i]     = a0 + ∑ancosx + ∑bnsinx
            uexact[i] = 3.0/(5.0 - 4.0*cos(xi))
        end
    end
    x = range(0,2π, Nx)

    #Plot to PDF
    fourier_legend = string("F(u) n=", N)
    if (N == 1)
        plt = plot()
        plt = plot!(x/2π, [pu uexact], label=[fourier_legend "exact"], lw=[4 2])
    else
        plt = plot!(x/2π, pu, label=fourier_legend, lw=[4 2])
    end
    fname = string("tm_Fourier.pdf")
    savefig(fname)
    savefig(plt, fname)
end


N  = 2  #number of wave modes in the Fourier series
Nq = 1000   #quadrature point to find ân and b̂n
Nx = 1000 #point along the x-axis

expansion_type = "trig"
for N in [1 2 4 8 16]
    mt_fourier(N, Nq, Nx, expansion_type)
end
