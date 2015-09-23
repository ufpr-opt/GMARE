# Short step method for quadratic functions

include("plots.jl")

function short_step_method(Λ::Vector, x₀::Vector)
  x = copy(x₀)
  Ks = 2
  Kc = 8
  iter = 0
  first_sstep = false
  while norm(Λ.*x) > 1e-5
    d = -Λ.*x
    if iter%(Ks+Kc) < Kc
      # Cauchy step
      λ = dot(d,d)/dot(d,Λ.*d)
      first_sstep = true
    else
      # Short step
      if first_sstep
        S = 1e4
        x⁺ = x - S*Λ.*x
        g⁺ = Λ.*x⁺
        λ = dot(g⁺,g⁺)/dot(g⁺,Λ.*g⁺)
        first_sstep = false
      end
    end
    x = x + λ*d

    iter += 1
    if iter > 20
      break
    end
    if iter >= 2
      plot_g(Λ, x₀, x₁, x, 1/λ, filename="sstep-$(iter-2)-$(iter-1)-$(iter)")
      copy!(x₀, x₁)
      copy!(x₁, x)
    else
      x₁ = copy(x)
    end
  end
  return x, iter
end
