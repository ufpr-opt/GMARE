include("../src/universal_gradients.jl")

f(x) = x[1]^2 + 4*x[2]^2
∇f(x) = [2*x[1]; 8*x[2]]
x₀ = [1.0; 0.5]

println(fgm(f, ∇f, x₀))

f(x) = x[1]^4 + 4*x[2]^4
∇f(x) = [4*x[1]^3; 16*x[2]^3]
 
println(fgm(f, ∇f, x₀))
