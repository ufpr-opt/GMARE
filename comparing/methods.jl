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
    return x, iter, nMV, X[:,1:iter+1]
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
    return x, iter, nMV, X[:,1:iter+1]
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
    return x, iter, nMV, X[:,1:iter+1]
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
    return x, iter, nMV, X[:,1:iter+1]
  else
    return x, iter, nMV
  end
end


function short_step(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter=10000, Ki = 10, Ks = 2, Kc = 8, S = 1e4,
    history = false, hist_nmv = true)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),2*max_iter+1)
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
      if history && hist_nmv
        X[:,nMV] = x
      end
      d = d - λ*Gd
    else
      # Short step
      if first_sstep
        Gd = G*d
        nMV += 1
        d⁺ = d - S*Gd
        λ = dot(d⁺,d⁺)/dot(d⁺,G*d⁺)
        nMV += 1
        first_sstep = false
        x = x + λ*d
        if history && hist_nmv
          X[:,nMV-1] = x
          X[:,nMV] = x
        end
        d = d - λ*Gd
      end
      for i = 2:Ks
        x = x + λ*d
        d = d - λ*G*d
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
        end
      end
    end

    iter += 1
    if history && !hist_nmv
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    if hist_nmv
      return x, iter, nMV, X[:,1:nMV]
    else
      return x, iter, nMV, X[:,1:iter+1]
    end
  else
    return x, iter, nMV
  end
end

function alternate_short_step(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter=10000, Ki = 10, Ks = 2, Kc = 8, S = 1e4,
    history = false, hist_nmv = true)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),10*max_iter+1)
    X[:,1] = x₀
  end
  iter = 0
  first_sstep = false
  d = -(G*x + v)
  nMV = 1
  λs = 0.0
  while norm(d) > tol
    if iter < Ki || (iter-Ki)%(1+Kc) >= 1
      # Cauchy step
      Gd = G*d
      nMV += 1
      λ = dot(d,d)/dot(d,Gd)
      first_sstep = true
      x = x + λ*d
      if history && hist_nmv
        X[:,nMV] = x
      end
      d = d - λ*Gd
      if iter >= Ki
        x = x + λs*d
        d = d - λs*G*d
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
        end
      end
    else
      # Short step
      if first_sstep
        Gd = G*d
        nMV += 1
        d⁺ = d - S*Gd
        λs = dot(d⁺,d⁺)/dot(d⁺,G*d⁺)
        nMV += 1
        first_sstep = false
        x = x + λs*d
        if history && hist_nmv
          X[:,nMV-1] = x
          X[:,nMV] = x
        end
        d = d - λs*Gd
      end
      for i = 2:Ks
        x = x + λs*d
        d = d - λs*G*d
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
        end
      end
    end

    iter += 1
    if history && !hist_nmv
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    if hist_nmv
      return x, iter, nMV, X[:,1:nMV]
    else
      return x, iter, nMV, X[:,1:iter+1]
    end
  else
    return x, iter, nMV
  end
end

function dai_yuan(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  d = -(G*x + v)
  Gd = G*d
  dot_dd = dot(d,d)
  λ = dot_dd/dot(d,Gd)
  λp = 0.0
  x = x + λ*d
  d = d - λ*Gd
  nMV = 2
  iter = 1
  if history
    X = zeros(length(x₀),max_iter+1)
    X[:,1] = x₀
    X[:,2] = x
  end
  while norm(d) > tol
    Gd = G*d
    old_dot_dd = dot_dd
    dot_dd = dot(d,d)
    nMV += 1
    λp = λ
    λ = dot_dd/dot(d,Gd)
    λnow = 0.0
    if iter%4 == 1 || iter%4 == 2
      λnow = 2/(sqrt((1/λp-1/λ)^2 + 4*dot_dd/(λp^2*old_dot_dd)) + 1/λp + 1/λ)
    else
      λnow = λ
    end
    x = x + λnow*d
    d = d - λnow*Gd
    iter += 1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1]
  else
    return x, iter, nMV
  end
end
