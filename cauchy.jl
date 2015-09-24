# Cauchy method for quadratic functions

include("plots.jl")
include("aux.jl")

#=
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
=#

function cauchy_method(Λ::Vector, x₀::Vector; max_iter = 20)
  x = copy(x₀)
  xpp = copy(x)
  iter = 0
  while norm(Λ.*x) > 1e-5
    d = -Λ.*x
    λ = dot(d,d)/dot(d,Λ.*d)
    x = x + λ*d
    iter += 1
    if iter > max_iter
      break
    end
    if iter >= 2
      fname = "cauchy-$(getname(iter, max_iter))"
      plot_g(Λ, xpp, xp, x, 1/λ, filename=fname)
      copy!(xpp, xp)
      copy!(xp, x)
    else
      xp = copy(x)
    end
  end
  return x, iter
end
