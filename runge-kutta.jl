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



const N = 20
const m₀ = 0
const n₀ = 0

h = 0.1
t₀ = 0.
q = 100

kin = Array(Float64, q+1)
pot = Array(Float64, q+1)
tout = Array(Float64, q+1)

aout = Array[ones(Float64, N^2), ones(Float64, N^2)]
for k in 0:q
    tk = t₀ + k*h
    t, aout = ode45((t,a)->F(t,a, N, m₀,n₀, 25.), aout[end]/norm(aout[end]), [tk, tk+h], points=:specified)
    kin[k+1] = kinetic(N, aout[end]/norm(aout[end]), 25.)
    pot[k+1] = potential(N, aout[end]/norm(aout[end]), 25., m₀, n₀)
    tout[k+1] = tk+h
end

aout = aout[end]/norm(aout[end])



# TODO: try different lattice sizes



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

ax[:imshow](abs2(reshape(aout, N,N)),
            origin="upper", ColorMap("gist_heat_r"), interpolation="none",
            extent=[-edge, edge, -edge, edge])

ax[:set_xlabel](L"$m$")
ax[:set_ylabel](L"$n$")

fig[:savefig]("psisq.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(10, 4))
ax[:plot](tout, kin, label="kinetic energy")
ax[:plot](tout, pot, label="potential energy")


ax[:legend](loc="lower right")



ax[:set_xlim](tout[1], tout[end])
#ax[:set_ylim](0, 70)

ax[:set_xlabel](L"$t\,[\frac{2}{\omega_0}]$")
ax[:set_ylabel](L"$E\,[\frac{\hbar\omega_0}{2}]$ ")

fig[:savefig]("energies.pdf", transparent=true, pad_inches=0.0, bbox_inches="tight")

plt.close(fig)
### <---


### <--- profiling
## Profile.clear()
## @profile (for i=1:100;  F(1.,ones(Float64, 512^2), 512,0,0,1.) ; end)
## Profile.print()
