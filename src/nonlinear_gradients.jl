function cauchy(f::Function, ∇f::Function, x₀::Vector;
    line_search::Function = armijo,
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    A = zeros(max_iter)
    X[:,1] = x₀
  end
  fx = f(x)
  nf = 1
  d = -∇f(x)
  ng = 1
  iter = 0
  while norm(d) > tol
    x, λ, fx, nfp = line_search(f, x, fx, d)
    nf += nfp
    d = -∇f(x)
    iter += 1
    if history
      X[:,iter+1] = x
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nf, ng, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nf, ng
  end
end

function barzilai_borwein(f::Function, ∇f::Function, x₀::Vector;
    line_search::Function = armijo,
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    A = zeros(max_iter)
    X[:,1] = x₀
  end
  fx = f(x)
  nf = 1
  d = -∇f(x)
  ng = 1
  iter = 0
  dold = d
  Δx = zeros(length(x))
  if norm(d) > tol
    x, λ, fx, nfp = line_search(f, x, fx, d)
    nf += nfp
    Δx = λ*d
    d = -∇f(x)
    ng += 1
    iter += 1
    if history
      X[:,iter+1] = x
      A[iter] = λ
    end
  end
  while norm(d) > tol
    Δg = dold-d
    λ = dot(Δx, Δg)/dot(Δg,Δg)
    x = x + λ*d
    Δx = λ*d
    dold = d
    d = -∇f(x)
    ng += 1
    iter += 1
    if history
      X[:,iter+1] = x
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nf, ng, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nf, ng
  end
end
