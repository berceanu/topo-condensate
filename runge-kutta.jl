using ODE


μ(l::Int,N::Int) = div(l-1,N) - div(N-1,2)
ν(l::Int,N::Int) = div(N-1,2) - rem(l-1,N)


cl(m::Int,n::Int,m₀::Int,n₀::Int) = (m - m₀)^2 + (n - n₀)^2

function F(t::Float64,a::Array{Complex{Float64}, 1}, N::Int,m₀::Int,n₀::Int,σ::Float64)
    d = N^2
    adot = Array(Complex{Float64}, d) 

    for l = 1:d
        m = μ(l::Int,N::Int)
        n = ν(l::Int,N::Int)

        p₃ = -(4σ + 1/σ*cl(m,n,m₀,n₀))

        adot[l] = p₃*a[l]
        if l > N
            adot[l] += σ*a[l-N] + σ*a[l-1]
        elseif l > 1
            adot[l] += σ*a[l-1]
        end 
        if l < d+1-N
            adot[l] += σ*a[l+1] + σ*a[l+N]
        elseif l < d
            adot[l] += σ*a[l+1]
        end
    end
    adot
end 



@time tout, aout = ode45((t,a)->F(t,a, 20, 0,0, 25.),
ones(Complex{Float64}, N^2), [0., 1.], points=:specified)


### --->
#plotting
using PyPlot

# matplotlib parameters
matplotlib["rcParams"][:update](["axes.labelsize" => 22,
                                 "axes.titlesize" => 20,
                                 "font.size" => 18,
                                 "legend.fontsize" => 14,
                                 "axes.linewidth" => 1.5,
                                 "font.family" => "serif",
                                 "font.serif" => "Computer Modern Roman",
                                 "xtick.labelsize" => 20,
                                 "xtick.major.size" => 5.5,
                                 "xtick.major.width" => 1.5,
                                 "ytick.labelsize" => 20,
                                 "ytick.major.size" => 5.5,
                                 "ytick.major.width" => 1.5,
                                 "text.usetex" => true,
                                 "figure.autolayout" => true])

edge = div(N-1,2)

fig, ax = plt.subplots(figsize=(5, 5))

ax[:imshow](abs2(reshape(aout[end], N,N)),
            origin="upper", ColorMap("gist_heat_r"), interpolation="none",
            extent=[-edge, edge, -edge, edge])

ax[:set_xlabel](L"$m$")
ax[:set_ylabel](L"$n$")

fig[:savefig]("caca.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(fig)
### <---


### <--- profiling
## Profile.clear()
## @profile (for i=1:100;  F(1.,ones(Complex{Float64}, 512^2), 512,0,0,1.) ; end)
## Profile.print()
