import Pkg; Pkg.add("ODE"); Pkg.add("PyPlot")
using ODE


m(l::Int,N::Int) = div(l-1,N) - div(N-1,2)
n(l::Int,N::Int) = div(N-1,2) - rem(l-1,N)


cl(N::Int,l::Int,m₀::Int,n₀::Int) = (m(l,N) - m₀)^2 + (n(l,N) - n₀)^2

## # symmetric gauge
## function Jl(N::Int,l::Int,p::Int,q::Int,J::Float64)
##     α = p/q
##     (-J * exp(-π*α*n(l,N)*im), -J * exp(π*α*m(l,N)*im))
## end

# Landau gauge
function Jl(N::Int,l::Int,p::Int,q::Int,J::Float64)
    α = p/q
    (-J,-J * exp(2π*α*m(l,N)*im))
end


function F(t::Float64,a::Array{Complex{Float64}, 1}, N::Int,p::Int,q::Int,m₀::Int,n₀::Int,J::Float64,σ::Float64,Ω::Float64,κ::Float64,f::Array{Float64,1})
    d = N^2
    adot = Array{Complex{Float64}}(undef, d)

    for l = 1:d
        (J1, J2) = Jl(N,l,p,q,J)
        #p3 and p4 are always present
        denominator = -σ/Ω*f[l] + im
        p₃ = (κ*cl(N,l,m₀,n₀) + σ*(f[l] - 1)*im) / denominator
        p₄ = (1 - σ*im) / denominator
        adot[l] = p₃*a[l] + p₄*abs2(a[l])*a[l]
        if l > N
            p₁ = J1 / denominator
            p₂ = conj(J2) / denominator
            adot[l] += p₁*a[l-N] + p₂*a[l-1]
        elseif l > 1
            p₂ = conj(J2) / denominator
            adot[l] += p₂*a[l-1]
        end
        if l < d+1-N
            p₅ = J2 / denominator
            p₆ = conj(J1) / denominator
            adot[l] += p₅*a[l+1] + p₆*a[l+N]
        elseif l < d
            p₅ = J2 / denominator
            adot[l] += p₅*a[l+1]
        end
    end
    adot
end

## N = 3

## γ0 = 1.
## f0 = ones(N^2)
## f = f0./γ0


## Γ0 = 1.
## a0 = ones(Complex{Float64}, N^2)
## a₀ = sqrt(Γ0/γ0)*a0

## U0 = 1.
## σ = Γ0/U0

## J0 = 1.
## J = Γ0*J0/(U0*γ0)

## κ0 = 1.
## κ = Γ0*κ0/(U0*γ0)

## ℏ = 1.
## Ω0 = 1./(2ℏ)
## Ω = 2ℏ*Γ0*Ω0/(U0*γ0)

## m₀ = 0
## n₀ = 0
## p = 1
## q = 11

# N = 65
N = 3

f = ones(N^2)
Ω = 1.


a₀ = ones(Complex{Float64}, N^2)

σ = 0.
J = 1.
κ = 10.

m₀ = 0
n₀ = 0

p = 0
q = 1


## tout, aout = odeXX((t,a)->F(t,a), a0, tspan; keywords...)
## F(t,a, N,p,q,m₀,n₀,J,σ,Ω,κ,f)



tout, aout = ode45((t,a)->F(t,a, N,p,q,m₀,n₀,J,σ,Ω,κ,f), a₀, [0., 10])


#a = sqrt(γ0/Γ0)*aout
#t = 2ℏ*Γ0*tout/(U0*γ0)

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



f, ax = plt.subplots(figsize=(5, 5))
ax[:imshow](abs2.(reshape(aout[end], N,N)),
            origin="upper", ColorMap("gist_heat_r"), interpolation="none",
            extent=[-edge, edge, -edge, edge])

ax[:set_xlabel](L"$m$")
ax[:set_ylabel](L"$n$")
ax[:set_xticks]([-edge,0, edge])
ax[:set_yticks]([-edge,0, edge])

savefig("squares.png")

#testing
ttst, atst = ode45((t,a)->F(t,a, 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)), ones(Complex{Float64}, 3^2), [0., 1])

using Test

const A = [-0.8090528248566107 - 0.3496936423122082im,
           -0.1587464671688188 + 0.8412535328311812im,
           -0.6587464671688188 + 0.34125353283118115im,
           0.5 + 1.5im,
           1.0 + 2.0im,
           0.5 + 1.5im,
           -0.6587464671688188 + 0.3412535328311812im,
           -0.15874646716881885 + 0.8412535328311812im,
           -0.8090528248566107 - 0.3496936423122081im];
const B = [0.5590362122824426 + 0.2612536645663328im,
           0.5843645134936586 + 0.725410327353141im,
           0.5304330500743882 + 0.6028800777195348im,
           0.5862026223731139 + 0.9334549437696245im,
           0.41646252857869204 + 1.2971224463101734im,
           0.586202622373114 + 0.9334549437696245im,
           0.5304330500743882 + 0.6028800777195347im,
           0.5843645134936588 + 0.725410327353141im,
           0.5590362122824425 + 0.26125366456633276im];
const C = [0.0,
           0.028768110347909062,
           0.12566178426740987,
           0.2228075252877766,
           0.32359991817436035,
           0.4261488432928628,
           0.5303748765226354,
           0.6369499640459351,
           0.7466894493187844,
           0.8601933256409379,
           0.977728108995228,
           1.0];

@test isapprox(F(1.,ones(Complex{Float64}, 3^2), 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)), A; atol=1e-16)
@test isapprox(atst[end], B; atol=1e-6)
@test isapprox(ttst, C; atol=1e-11)

using Profile


# @time for i=1:10^6; F(1.,ones(Complex{Float64}, 3^2), 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)) ; end
# 2.6s

Profile.clear()
@profile (for i = 1:10^6; F(1.,ones(Complex{Float64}, 3^2), 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)); end)

@time for i=1:100; F(1.,ones(Complex{Float64}, 512^2), 512,1,11,0,0,1.,1.,1.,1.,ones(512^2)) ; end
# 10.9s

## Profile.clear()
## @profile (for i=1:100; F(1.,ones(Complex{Float64}, 512^2), 512,1,11,0,0,1.,1.,1.,1.,ones(512^2)) ; end)
## Profile.print()
