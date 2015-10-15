function coord_cauchy(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 1000, history = false)
  x = copy(x₀)
  n = length(x)
  if history
    X = zeros(length(x₀),max_iter+2)
    A = zeros(max_iter+1)
    X[:,1] = x₀
  end
  iter = 0
  nMV = 0
  d = fill(Inf,n)
  while norm(d) > tol
    i = iter%n + 1
    d[i] = -rowdot(G,i,x) - v[i]
    λ = 1/G[i,i]
    x[i] = x[i] + λ*d[i]
    iter += 1
    if history && iter%n == 0
      nMV += 1
      X[:,nMV+1] = x
      A[nMV] = λ
    end
    if iter >= max_iter*n
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:nMV+1], A[1:nMV]
  else
    return x, iter, nMV
  end
end

function coord_cauchy_rand(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 1000, history = false)
  x = copy(x₀)
  n = length(x)
  if history
    X = zeros(length(x₀),max_iter+2)
    A = zeros(max_iter+1)
    X[:,1] = x₀
  end
  L = maximum([G[i,i] for i = 1:n])
  iter = 0
  nMV = 0
  d = fill(Inf,n)
  while norm(d) > tol
    i = rand(1:n)
    d[i] = -rowdot(G,i,x) - v[i]
    λ = 1/L
    x[i] = x[i] + λ*d[i]
    iter += 1
    if history && iter%n == 0
      nMV += 1
      X[:,nMV+1] = x
      A[nMV] = λ
    end
    if iter >= max_iter*n
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:nMV+1], A[1:nMV]
  else
    return x, iter, nMV
  end
end
