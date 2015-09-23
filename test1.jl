# This is the basic test

σ = 0.01

n = 50
Λ = linspace(σ, 1.0, n)

x₀ = 1./sqrt(Λ)

include("plots.jl")

plot_g(Λ, x₀, filename="iter-0")
d = -Λ.*x₀
λ = dot(d,d)/dot(d,Λ.*d)
x₁ = x₀ + λ*d
plot_g(Λ, x₀, x₁, filename="iter-0-1")
d = -Λ.*x₁
λ = dot(d,d)/dot(d,Λ.*d)
x₂ = x₁ + λ*d
plot_g(Λ, x₀, x₁, x₂, 1/λ, filename="iter-0-1-2")
