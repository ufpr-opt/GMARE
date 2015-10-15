function cauchy(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    A = zeros(max_iter)
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
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1], A[1:iter]
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
    A = zeros(max_iter)
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
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nMV
  end
end

function barzilai_borwein(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    A = zeros(max_iter)
    X[:,1] = x₀
  end
  iter = 0
  d = -(G*x + v)
  Gd = G*d
  λ = dot(d,d)/dot(d,Gd)
  nMV = 2
  λ₊ = λ
  while norm(d) > tol
    #println("($λ,$(dot(d,d)/dot(d,G*d)))")
    x = x + λ*d
    d = d - λ*Gd
    if history
      A[iter+1] = λ
    end
    λ = λ₊
    Gd = G*d
    nMV += 1
    λ₊ = dot(d,d)/dot(d,Gd)
    iter +=1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nMV
  end
end

function alternate_cauchy(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false, interval = 2)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),max_iter+1)
    A = zeros(max_iter)
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
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nMV
  end
end


function short_step(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter=10000, Ki = 10, Ks = 2, Kc = 6,
    history = false, hist_nmv = true)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),2*max_iter+1)
    A = zeros(2*max_iter)
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
      λ *= 0.8 + rand()*0.4
      first_sstep = true
      x = x + λ*d
      if history && hist_nmv
        X[:,nMV] = x
        A[nMV-1] = λ
      end
      d = d - λ*Gd
    else
      # Short step
      if first_sstep
        Gd = G*d
        nMV += 1
        # S is infinite
        d⁺ = -Gd
        λ = dot(d⁺,d⁺)/dot(d⁺,G*d⁺)
        nMV += 1
        first_sstep = false
        x = x + λ*d
        if history && hist_nmv
          X[:,nMV-1] = x
          X[:,nMV] = x
          A[nMV-1] = λ
        end
        d = d - λ*Gd
      end
      for i = 2:Ks
        x = x + λ*d
        d = d - λ*G*d
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
          A[nMV-1] = λ
        end
      end
    end

    iter += 1
    if history && !hist_nmv
      X[:,iter+1] = x
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    if hist_nmv
      return x, iter, nMV, X[:,1:nMV], A[1:nMV-1]
    else
      return x, iter, nMV, X[:,1:iter+1], A[1:iter]
    end
  else
    return x, iter, nMV
  end
end

function alternate_short_step(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter=10000, Ki = 10, Ks = 2, Kc = 6,
    history = false, hist_nmv = true)
  x = copy(x₀)
  if history
    X = zeros(length(x₀),10*max_iter+1)
    A = zeros(10*max_iter)
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
      λ *= 0.8 + rand()*0.4
      first_sstep = true
      x = x + λ*d
      if history && hist_nmv
        X[:,nMV] = x
        A[nMV-1] = λ
      end
      d = d - λ*Gd
      if iter >= Ki
        x = x + λs*d
        d = d - λs*G*d
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
          A[nMV-1] = λs
        end
      end
    else
      # Short step
      Gd = G*d
      nMV += 1
      # S is infinite
      d⁺ = -Gd;
      λs = dot(d⁺,d⁺)/dot(d⁺,G*d⁺)
      nMV += 1
      x = x + λs*d
      if history && hist_nmv
        X[:,nMV-1] = x
        X[:,nMV] = x
        A[nMV-1] = λs
      end
      d = d - λs*Gd
      for i = 2:Ks
        x = x + λs*d
        d = d - λs*G*d
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
          A[nMV-1] = λs
        end
      end
    end

    iter += 1
    if history && !hist_nmv
      X[:,iter+1] = x
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    if hist_nmv
      return x, iter, nMV, X[:,1:nMV], A[1:nMV-1]
    else
      return x, iter, nMV, X[:,1:iter+1], A[1:iter]
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
    A = zeros(max_iter)
    X[:,1] = x₀
    X[:,2] = x
    A[1] = λ
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
      A[iter] = λnow
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nMV
  end
end

function alternate_dai_yuan(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false, hist_nmv = true,
    Ki = 8, Kc = 8, Kdy = 6)
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
    X = zeros(length(x₀),10*max_iter+1)
    A = zeros(10*max_iter)
    X[:,1] = x₀
    X[:,2] = x
    A[1] = λ
  end
  λdy = 0.0
  while norm(d) > tol
    if iter < Ki || (iter-Ki)%(1+Kc) >= 1
      # Cauchy step
      Gd = G*d
      nMV += 1
      old_dot_dd = dot_dd
      dot_dd = dot(d,d)
      λp = λ
      λ = dot_dd/dot(d,Gd)
      first_sstep = true
      x = x + λ*d
      if history && hist_nmv
        X[:,nMV] = x
        A[nMV-1] = λ
      end
      d = d - λ*Gd
    else
      λdy = 2/(sqrt((1/λp-1/λ)^2 + 4*dot_dd/(λp^2*old_dot_dd)) + 1/λp + 1/λ)
      for i = 1:Kdy
        Gd = G*d
        x = x + λdy*d
        d = d - λdy*Gd
        nMV += 1
        if history && hist_nmv
          X[:,nMV] = x
          A[nMV-1] = λdy
        end
      end
    end
    iter += 1
    if history && !hist_nmv
      X[:,iter+1] = x
      A[iter] = λ
    end
    if iter >= max_iter
      break
    end
  end
  if history
    if hist_nmv
      return x, iter, nMV, X[:,1:nMV], A[1:nMV-1]
    else
      return x, iter, nMV, X[:,1:iter+1], A[1:iter]
    end
  else
    return x, iter, nMV
  end
end

function conjugate_gradient(G::Matrix, v::Vector, x₀::Vector;
    tol = 1e-6, max_iter = 10000, history = false)
  r = G*x₀ + v
  x = copy(x₀)
  if history
    X = zeros(length(x₀), max_iter+1)
    A = zeros(max_iter)
    X[:,1] = x₀
  end
  p = -r
  nMV = 1
  iter = 0
  dot_rr = dot(r,r)
  while dot_rr > tol*tol
    Gp = G*p
    nMV += 1
    α = dot_rr/dot(p,Gp)
    x = x + α*p
    r = r + α*Gp
    dot_rrp = dot(r,r)
    β = dot_rrp/dot_rr
    dot_rr = dot_rrp
    p = β*p - r
    iter += 1
    if history
      X[:,iter+1] = x
    end
    if iter >= max_iter
      break
    end
  end
  if history
    return x, iter, nMV, X[:,1:iter+1], A[1:iter]
  else
    return x, iter, nMV
  end
end
