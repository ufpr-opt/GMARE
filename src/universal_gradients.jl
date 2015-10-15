# ξ(x,y) = 0.5‖x-y‖²
function fgm(f::Function, ∇f::Function, x₀::Vector;
    tol = 1e-6, max_iter = 1000, history=true)
  A = 0
  x = copy(x₀)
  if history
    X = zeros(length(x),max_iter+1)
    X[:,1] = x₀
  end
  y = copy(x₀)
  v = copy(x₀)
  L = 1.0
  iter = 0
  ef = 0
  a = 0.0
  Δf = Inf
  while iter < max_iter
    while true
      a = (1 + sqrt(1 + 4*L*A))/(2*L)
      A⁺ = A + a
      τ = a/A⁺
      x⁺ = τ*v + (1-τ)*y
      xhat = v - a*∇f(x⁺)
      y⁺ = τ*xhat + (1-τ)*y
      if f(y⁺) <= f(x⁺) + dot(∇f(x⁺),y⁺-x⁺) + L*norm(y⁺-x⁺)^2/2 + tol*τ/2
        Δf = f(y⁺) - f(y)
        (x, y, A, L) = (x⁺, y⁺, A, L/2)

        break
      else
        L = L*2
        if L == Inf
          error("INF")
        end
      end
    end
    iter += 1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      ef = 1
      break
    end
    v = v - a*∇f(x)
  end
  if history
    return x, f(x), ∇f(x), iter, ef, X[:,1:iter+1]
  else
    return x, f(x), ∇f(x), iter, ef
  end
end
