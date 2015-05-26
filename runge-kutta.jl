using ODE


μ(l::Int,N::Int) = div(l-1,N) - div(N-1,2)
ν(l::Int,N::Int) = div(N-1,2) - rem(l-1,N)


cl(m::Int,n::Int,m₀::Int,n₀::Int) = (m - m₀)^2 + (n - n₀)^2

function F(t::Float64,a::Array{Float64, 1}, N::Int,m₀::Int,n₀::Int,σ::Float64)
    d = N^2
    adot = Array(Float64, d) 

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

const N = 20



@time tout, aout = ode45((t,a)->F(t,a, N, 0,0, 25.),
ones(Float64, N^2)/N^2, linspace(0,15,16), points=:specified)



function kinetic(N::Int, a::Array{Float64, 1}, σ::Float64)
    d = N^2
    s = 0.
    for l = 1:d #<-- TODO: check this part
        if (l>1) && (l < d+1-N)
            s += abs2(a[l+N] - a[l]) + abs2(a[l-1] - a[l])
        elseif l>1
            s += abs2( - a[l]) + abs2(a[l-1] - a[l])
        elseif l < d+1-N
            s += abs2(a[l+N] - a[l]) + abs2( - a[l])
        else 
            s += abs2(- a[l]) + abs2(- a[l])
        end
    end ##
    return σ*s/sum(abs2(a))
end 

# TODO: constant norm
# TODO: plot energies vs time (kin, pot, tot)
# TODO: try different lattice sizes


function potential(N::Int, a::Array{Float64, 1}, σ::Float64, m₀::Int, n₀::Int)
    s = 0.
    d = N^2 
    for l = 1:d
        m = μ(l::Int,N::Int)
        n = ν(l::Int,N::Int)

        s += abs2(a[l])*cl(m,n,m₀,n₀)
    end
    return 1/σ*s/sum(abs2(a))
end 

kinetic(N, aout[end], 25.)

potential(N, aout[end], 25., 0, 0)



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
## @profile (for i=1:100;  F(1.,ones(Float64, 512^2), 512,0,0,1.) ; end)
## Profile.print()
