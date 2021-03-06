# The two most important line searches
# All these assume d = -∇f(x)

function armijo(f::Function, x::Vector, fx::Real, d::Vector;
    α::Real = 0.5, η::Real = 0.7)
  dot_dd = dot(d,d)
  t = 1.0
  fx⁺ = f(x + t*d)
  nf = 1
  while fx⁺ > fx - α*t*dot_dd
    t *= η
    fx⁺ = f(x + t*d)
    nf += 1
  end
  return x+t*d, t, fx⁺, nf
end

function golden_search(f::Function, x::Vector, fx::Real, v::Vector;
    tol = 1e-3)
  ϕ(t) = f(x + t*v)
  g = 1/golden

  (a, b) = (0, 1)
  ϕa, ϕb = (fx, ϕ(b))
  if ϕb >= ϕa
    c, d = b - g*(b-a), a + g*(b-a)
    ϕc = ϕ(c); ϕd = ϕ(d)
    nf = 2
  else
    d, b = b, golden
    ϕd, ϕb = ϕb, ϕ(d)
    nf = 1
    while ϕb < ϕd
      a, d, b = d, b, golden*b
      ϕa, ϕd, ϕb = ϕd, ϕb, ϕ(b)
      nf += 1
      if b > 1e6
        error("Golden search failed")
      end
    end
    c = b - g*(b-a)
    ϕc = ϕ(c)
    nf += 1
  end
  while abs(c-d) > tol
    if ϕc < ϕd
      b, d = d, c
      ϕb, ϕd = ϕd, ϕc
      c = b - g*(b-a)
      ϕc = ϕ(c)
    else
      a, c = c, d
      ϕa, ϕc = ϕc, ϕd
      d = a + g*(b-a)
      ϕd = ϕ(d)
    end
    nf += 1
  end
  return x+c*v, c, ϕc, nf
end
