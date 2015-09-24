# This is second basic test

include("cauchy.jl")
include("short_step.jl")

σ = 0.01

n = 50
Λ = exp_normal(σ, 1.0, n; weigh_left=false)

x₀ = 1./sqrt(Λ)

x, iter = cauchy_method(Λ, x₀, max_iter=100)

x, iter = short_step_method(Λ, x₀, max_iter=100)
