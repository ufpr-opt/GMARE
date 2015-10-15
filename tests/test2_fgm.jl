include("../src/universal_gradients.jl")

a = 1
b = 1
function fmax(x)
  f1 = x[1]^2 + 4*x[2]^2
  f2 = 9*(x[1]-a)^2 + 5*(x[2]-b)^2
  return max(f1,f2)
end
function ∇fmax(x)
  f1 = x[1]^2 + 4*x[2]^2
  f2 = 9*(x[1]-a)^2 + 5*(x[2]-b)^2
  if f1 > f2
    return [2*x[1]; 8*x[2]]
  else
    return [18*(x[1]-a); 10*(x[2]-b)]
  end
end
x₀ = [-2.0; -1.0]

println(fgm(fmax, ∇fmax, x₀))
