module QuadGen

"""
Generation of quadratic functions in the form
    f(x) = ¹/₂xᵀQx + vᵀx + c
  In general, these rules are valid:

  - every function accepts 3 values:
      λmin - the smallest eigenvalue
      λmax - the largest eigenvalue
      n - the dimension

  - every function outputs either
      Q - Matrix
      v - Vector
      c - Real
    or
      f  - Objective function
      ∇f - Gradient of the objective function
      H  - Hessian of the objective function
    The first is the default.
    To choose, use the keyword `output` which can have the values
    `:matrix` or `:function`.

  - every problem can accept an eigenvalue distribution in the form of
    the keyword argument `distribution`. This paramers must be
    a function (a,b,n)->f(a,b,n), such as linspace(a,b,n), that
    returns an array of size n with values in [a,b].
    This distribution will be scaled to fit into [λmin, λmax], but it won't be
    sorted.
"""

macro sanitize_input()
  quote
    outputs = [:matrix, :function]
    if λmax < λmin
      error("λmax ($λmax) can't be smaller than λmin ($λmin)")
    elseif !(output in outputs)
      error("output `$($output)` not accepted. Use one of $outputs")
    end
  end
end

function fit_on_interval(x, λmin, λmax)
  @sanitize_input
  m = minimum(x); M = maximum(x)
  if M != m
    x = (x-m)*(λmax-λmin)/(M-m) + λmin
  elseif λmin != λmax
    error("Can't scale [$m,$M] to [$λmin,$λmax]")
  else
    x = fill(λmin, length(x))
  end
  return x
end

function diagonal(λmin, λmax, n;
    distribution=(a,b,n)->linspace(a,b,n),
    output=:matrix)
  @sanitize_input
  Λ = fit_on_interval(distribution(λmin, λmax, n), λmin, λmax)
  if output == :matrix
    return diagm(Λ), zeros(n), 0.0
  elseif output == :function
    return x->0.5*dot(x,Λ.*x), x->Λ.*x, x->diagm(Λ)
  end
end
diagonal_rand(λmin, λmax, n; output=:matrix) = diagonal(λmin, λmax, n,
    distribution=(a,b,n)->rand(n)*(b-a)+a, output=output)

function general(λmin, λmax, n;
    distribution=(a,b,n)->linspace(a,b,n),
    output=:matrix)
  @sanitize_input
  Λ = fit_on_interval(distribution(λmin, λmax, n), λmin, λmax)
  (G,R) = qr(rand(n,n))
  Q = G'*diagm(Λ)*G
  sol = distribution(-1.0, 1.0, n)
  if output == :matrix
    v = -Q*sol
    c = 0.5*dot(sol,Q*sol)
    return Q, v, c
  elseif output == :function
    return x->0.5*dot(x-sol,Q*(x-sol)), x->Q*(x-sol), x->Q
  end
end

end
