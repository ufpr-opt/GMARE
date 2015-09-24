# Short step method for quadratic functions

include("plots.jl")
include("aux.jl")

function short_step_method(Λ::Vector, x₀::Vector;
    max_iter=20, Ki = 10, Ks = 2, Kc = 8, S = 1e4)
  x = copy(x₀)
  xpp = copy(x)
  iter = 0
  first_sstep = false
  while norm(Λ.*x) > 1e-5
    d = -Λ.*x
    if iter < Ki || (iter-Ki)%(Ks+Kc) >= Ks
      # Cauchy step
      λ = dot(d,d)/dot(d,Λ.*d)
      first_sstep = true
    else
      # Short step
      if first_sstep
        x⁺ = x - S*Λ.*x
        g⁺ = Λ.*x⁺
        λ = dot(g⁺,g⁺)/dot(g⁺,Λ.*g⁺)
        first_sstep = false
      end
    end
    x = x + λ*d

    iter += 1
    if iter > max_iter
      break
    end
    if iter >= 2
      fname = "sstep-$(getname(iter, max_iter))"
      plot_g(Λ, xpp, xp, x, 1/λ, filename=fname)
      copy!(xpp, xp)
      copy!(xp, x)
    else
      xp = copy(x)
    end
  end
  return x, iter
end
