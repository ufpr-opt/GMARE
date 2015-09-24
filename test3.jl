# This is second basic test

σ = 0.01

n = 50
Λ = linspace(σ, 1.0, n)

include("cauchy.jl")

x₀ = 1./sqrt(Λ)
x, iter = cauchy_method(Λ, x₀, max_iter=100)

include("short_step.jl")
x₀ = 1./sqrt(Λ)
x, iter = short_step_method(Λ, x₀, max_iter=100)
