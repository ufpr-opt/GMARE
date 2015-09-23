# Cauchy method for quadratic functions

include("plots.jl")

function cauchy_method(G::Matrix, g::Vector, x₀::Vector)
  x = copy(x₀)
  d = -(G*x + g)
  iter = 0
  while norm(d) > 1e-5
    x = x + dot(d,d)*d/dot(d,G*d)
    d = -(G*x + g)
    iter += 1
    if iter > 10000
      break
    end
  end
  return x, iter
end

function cauchy_method(Λ::Vector, x₀::Vector)
  x = copy(x₀)
  iter = 0
  while norm(Λ.*x) > 1e-5
    d = -Λ.*x
    λ = dot(d,d)/dot(d,Λ.*d)
    x = x + λ*d
    iter += 1
    if iter > 20
      break
    end
    if iter >= 2
      plot_g(Λ, x₀, x₁, x, 1/λ, filename="cauchy-$(iter-2)-$(iter-1)-$(iter)")
      copy!(x₀, x₁)
      copy!(x₁, x)
    else
      x₁ = copy(x)
    end
  end
  return x, iter
end
