using Plots

N = 16
n = range(-N/2, N/2, N)
@info n
Nx = 100
x  = range(0, 2π, Nx)
pu = zeros(Complex, Nx)
for j=1:100
    u = 0.0
    for i=1:N
        un = 1.0/2.0^(abs(n[i]))
        #un = (2.0/π)/(1.0 - 4*n[i]*n[i])
        u = u + un*exp(n[i]*x[j]im)
    end
    pu[j] = u
end


#uex = sin.(x/2)
uex =  3.0./(5.0.-4.0*cos.(x))
plot(x, [real(pu) uex], label=["P(u)" "exact"], lw=[2 1])
