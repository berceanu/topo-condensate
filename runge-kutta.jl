using ODE


μ(l::Int,N::Int) = div(l-1,N) - div(N-1,2)
ν(l::Int,N::Int) = div(N-1,2) - rem(l-1,N)


cl(m::Int,n::Int,m₀::Int,n₀::Int) = (m - m₀)^2 + (n - n₀)^2

function F(t::Float64,a::Array{Complex{Float64}, 1}, N::Int,p::Int,q::Int,m₀::Int,n₀::Int,J::Float64,σ::Float64,Ω::Float64,κ::Float64,f::Array{Float64,1})
    d = N^2
    adot = Array(Complex{Float64}, d) 

    ϖ = im*J 
    for l = 1:d
        m = μ(l::Int,N::Int)
        n = ν(l::Int,N::Int)

        p₃ = -im*κ*cl(m,n,m₀,n₀)

        adot[l] = p₃*a[l] - im*abs2(a[l])*a[l]
        if l > N
            adot[l] += ϖ*a[l-N] + ϖ*a[l-1]
        elseif l > 1
            adot[l] += ϖ*a[l-1]
        end 
        if l < d+1-N
            adot[l] += ϖ*a[l+1] + ϖ*a[l+N]
        elseif l < d
            adot[l] += ϖ*a[l+1]
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



### --->
## N = 3
# N = 65

## f = ones(N^2)
## Ω = 1.


## a₀ = ones(Complex{Float64}, N^2)

## σ = 0.
## J = 1.
## κ = 10.

## m₀ = 0
## n₀ = 0

## p = 0
## q = 1
### <---

## tout, aout = odeXX((t,a)->F(t,a), a0, tspan; keywords...)
## F(t,a, N,p,q,m₀,n₀,J,σ,Ω,κ,f)



### tout, aout = ode45((t,a)->F(t,a, N,p,q,m₀,n₀,J,σ,Ω,κ,f), a₀, [0., 10], points=:specified)


#a = sqrt(γ0/Γ0)*aout
#t = 2ℏ*Γ0*tout/(U0*γ0)

### --->
#plotting
## using PyPlot

## # matplotlib parameters
## matplotlib["rcParams"][:update](["axes.labelsize" => 22,
##                                  "axes.titlesize" => 20,
##                                  "font.size" => 18,
##                                  "legend.fontsize" => 14,
##                                  "axes.linewidth" => 1.5,
##                                  "font.family" => "serif",
##                                  "font.serif" => "Computer Modern Roman",
##                                  "xtick.labelsize" => 20,
##                                  "xtick.major.size" => 5.5,
##                                  "xtick.major.width" => 1.5,
##                                  "ytick.labelsize" => 20,
##                                  "ytick.major.size" => 5.5,
##                                  "ytick.major.width" => 1.5,
##                                  "text.usetex" => true,
##                                  "figure.autolayout" => true])

## edge = div(N-1,2)


## f, ax = plt.subplots(figsize=(5, 5))
## ax[:imshow](abs2(reshape(aout[end], N,N)),
##             origin="upper", ColorMap("gist_heat_r"), interpolation="none",
##             extent=[-edge, edge, -edge, edge])

## ax[:set_xlabel](L"$m$")
## ax[:set_ylabel](L"$n$")
## ax[:set_xticks]([-edge,0, edge])
## ax[:set_yticks]([-edge,0, edge])
### <---


#testing
ttst, atst = ode45((t,a)->F(t,a, 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)), ones(Complex{Float64}, 3^2), [0., 1])

using Base.Test

const A = [0.0-1.0im
           0.0+1.0im
           0.0+0.0im
           0.0+2.0im
           0.0+3.0im
           0.0+2.0im
           0.0+0.0im
           0.0+1.0im
           0.0-1.0im]

const B =[-0.00714133+0.388327im
          0.533606+0.842441im
          0.455033+0.816523im
          0.721952+0.912371im
          1.34469+0.668148im
          0.721952+0.912371im
          0.455033+0.816523im
          0.533606+0.842441im
          -0.00714133+0.388327im]


const C = [0.0      
           0.0278208
           0.095038 
           0.162282 
           0.229763 
           0.298015 
           0.367175 
           0.436445 
           0.504675 
           0.571481 
           0.636676 
           0.699467 
           0.759766 
           0.81879  
           0.878097 
           0.93911  
           1.0]


@test_approx_eq_eps(F(1.,ones(Complex{Float64}, 3^2), 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)), A, 1e-16)

@test_approx_eq_eps(atst[end], B, 1e-5)

@test_approx_eq_eps(ttst, C, 1e-6)

# 2.6s; 0.8s
#@time for i=1:10^6; F(1.,ones(Complex{Float64}, 3^2), 3,1,11,0,0,1.,1.,1.,1.,ones(3^2)) ; end

# 10.9s; 1.8s
#@time for i=1:100; F(1.,ones(Complex{Float64}, 512^2), 512,1,11,0,0,1.,1.,1.,1.,ones(512^2)) ; end


## Profile.clear()
## @profile (for i=1:100; F(1.,ones(Complex{Float64}, 512^2), 512,1,11,0,0,1.,1.,1.,1.,ones(512^2)) ; end)
## Profile.print()
