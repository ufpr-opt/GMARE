
function cauchy(G::Matrix, g::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    if max_iter > 100
      error("Don't use such a large max_iter")
    end
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
  end
  Gx = G*x
  nMV = 1
  d = -(Gx + g)
  iter = 0
  while norm(d) > tol
    λ = dot(d,d)/dot(d,G*d)
    nMV += 1
    x = x + λ*d
    Gx = G*x
    nMV += 1
    d = -(Gx + g)
    iter += 1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter]
  else
    return x, iter, nMV
  end
end

function random_decrease(G::Matrix, g::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000)
  β = 1e-4
  rnd() = rand()*(2-2*β)+β
  x = copy(x₀)
  Gx = G*x
  nMV = 1
  d = -(Gx + g)
  iter = 0
  while norm(d) > tol
    Gd = G*d
    nMV += 1
    λ = rnd()*dot(d,d)/dot(d,Gd)
    x = x + λ*d
    Gx = G*x
    nMV += 1
    d = -(Gx + g)
    iter += 1
    if iter > max_iter
      break
    end
  end
  return x, iter, nMV
end

function barzilai_borwein(G::Matrix, g::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000)
  x = copy(x₀)
  Gx = G*x
  nMV = 1
  d = -(Gx + g)
  iter = 0
  while norm(d) > tol
    if iter%2 == 0
      λ = dot(d,d)/dot(d,G*d)
      nMV += 1
    end
    x = x + λ*d
    Gx = G*x
    nMV += 1
    d = -(Gx + g)
    iter += 1
    if iter > max_iter
      break
    end
  end
  return x, iter, nMV
end


function short_step(G::Matrix, g::Vector, x₀::Vector;
    tol = 1e-6, max_iter=10000, Ki = 10, Ks = 2, Kc = 8, S = 1e4)
  x = copy(x₀)
  Gx = G*x
  nMV = 1
  iter = 0
  first_sstep = false
  d = -(Gx + g)
  while norm(d) > tol
    if iter < Ki || (iter-Ki)%(Ks+Kc) >= Ks
      # Cauchy step
      λ = dot(d,d)/dot(d,G*d)
      nMV += 1
      first_sstep = true
    else
      # Short step
      if first_sstep
        x⁺ = x + S*d
        g⁺ = G*x⁺+g
        λ = dot(g⁺,g⁺)/dot(g⁺,G*g⁺)
        nMV += 2
        first_sstep = false
      end
    end
    x = x + λ*d
    Gx = G*x
    nMV += 1
    d = -(Gx + g)

    iter += 1
    if iter > max_iter
      break
    end
  end
  return x, iter, nMV
end
