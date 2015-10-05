include("methods.jl")
include("aux.jl")

methods = [cauchy, barzilai_borwein, short_step]

σ = 0.001
Λ = [σ;1.0]

function f(θ,mtd)
  x₀ = [cos(θ)/σ;sin(θ)]
  x₀ = x₀/norm(x₀)
  x, iter, nMV = mtd(diagm(Λ), [0.0;0.0], x₀)
  #println("x₀, it, nmv = $x₀, $iter, $nMV")
  return iter
end

N = 1000
#θ = fmax(x->f(x,cauchy), 0, pi/2)
fmax = 0
xmax = -1
for x = linspace(0,pi/2,N)
  fx = f(x,cauchy)
  if fx > fmax
    fmax = fx
    xmax = x
  end
end
θ = xmax
println("θ = $θ")
x₀ = [cos(θ)/σ;sin(θ)]
println("x = $(x₀/sin(θ))")
x₀ = x₀/norm(x₀)
println("x = $x₀")
println("fmax = $fmax")
