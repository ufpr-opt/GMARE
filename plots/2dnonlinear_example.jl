include("../src/methods.jl")
include("../src/aux.jl")

using Gadfly

function example2d_nonlinear()
  methods = [cauchy]
  line_searches = [armijo, golden_search]

  path = createpath("2d-nonlinear-examples")
  σ = 10

  f(x) = 1-exp(-x[1]^2-σ*x[2]^2)
  ∇f(x) = -exp(-x[1]^2-σ*x[2]^2)*[-2*x[1];-2*σ*x[2]]
  x₀ = ones(2)

  for ls in line_searches
    for mtd in methods
      x, iter, nMV, X = mtd(f, ∇f, x₀, history=true,
          line_search=ls, max_iter=100, tol=0.02)
      m = size(X)[2]
      F = [0.5*dot(X[:,i],Λ.*X[:,i]) for i=1:m]

      N = 200
      t = linspace(-2.0,2.0,N);
      p = Gadfly.plot(layer(z=(x,y)->f([x;y]), x=t, y=t,
          Geom.contour(levels=F)),
          layer(x=X[1,:], y=X[2,:], Geom.line))
      draw(PDF("$path/$mtd-$ls.pdf", 4inch, 3inch), p)
    end
  end
end

example2d()
