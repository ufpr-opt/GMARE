
function cauchy(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
  end
  d = -(G*x + v)
  nMV = 1
  iter = 0
  while norm(d) > tol
    Gd = G*d
    nMV += 1
    λ = dot(d,d)/dot(d,Gd)
    x = x + λ*d
    d = d - λ*Gd
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

function random_decrease(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  β = 1e-4
  rnd() = rand()*(2-2*β)+β
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
  end
  d = -(G*x + v)
  nMV = 1
  iter = 0
  while norm(d) > tol
    Gd = G*d
    nMV += 1
    λ = rnd()*dot(d,d)/dot(d,Gd)
    x = x + λ*d
    d = d - λ*Gd
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

function barzilai_borwein(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
  end
  iter = 0
  d = -(G*x + v)
  nMV = 1
  λ = dot(d,d)/dot(d,G*d)
  λ₊ = λ
  while norm(d) > tol
    Gd = G*d
    nMV += 1
    x = x + λ*d
    d = d - λ*Gd
    λ = λ₊
    λ₊ = dot(d,d)/dot(d,G*d)
    iter +=1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X
  else
    return x, iter, nMV
  end
end

function alternate_cauchy(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false, interval = 2)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
  end
  d = -(G*x + v)
  nMV = 1
  iter = 0
  while norm(d) > tol
    Gd = G*d
    nMV += 1
    if iter%interval == 0
      λ = dot(d,d)/dot(d,Gd)
    end
    x = x + λ*d
    d = d - λ*Gd
    iter += 1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X
  else
    return x, iter, nMV
  end
end


function short_step(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter=10000, Ki = 10, Ks = 4, Kc = 8, S = 1e4,
    history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
  end
  iter = 0
  first_sstep = false
  d = -(G*x + v)
  nMV = 1
  while norm(d) > tol
    λ = 0.0
    if iter < Ki || (iter-Ki)%(1+Kc) >= 1
      # Cauchy step
      Gd = G*d
      nMV += 1
      λ = dot(d,d)/dot(d,Gd)
      first_sstep = true
      x = x + λ*d
      d = d - λ*Gd
    else
      # Short step
      if first_sstep
        Gd = G*d
        nMV += 1
        v⁺ = v + S*Gd
        λ = dot(v⁺,v⁺)/dot(v⁺,G*v⁺)
        nMV += 1
        first_sstep = false
      end
      for i = 1:Ks
        x = x + λ*d
        d = d - λ*G*d
        nMV += 1
      end
    end

    iter += 1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X
  else
    return x, iter, nMV
  end
end
