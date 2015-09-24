# This is second basic test

include("cauchy.jl")
include("short_step.jl")

σ = 0.01

n = 25
#Λ = [σ*ones(n);(σ+1.0)/2;ones(n)]
#Λ = [σ;uniform(0.8*σ+0.2,0.2*σ+0.8,n-2);1.0]
Λ = [σ*ones(n);1.0*ones(n)]
n = 2*n

x₀ = 1./sqrt(Λ)

x, iter = cauchy_method(Λ, x₀, max_iter=100)

x, iter = short_step_method(Λ, x₀, max_iter=100)
