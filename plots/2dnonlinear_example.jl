include("../src/methods.jl")
include("../src/aux.jl")

using Gadfly

function example2d_nonlinear()
  methods = [cauchy, barzilai_borwein]
  line_searches = [armijo, golden_search]

  path = createpath("2d-nonlinear-examples")
  σ = 10

  #f(x) = 1-exp(-x[1]^2-σ*x[2]^2)
  #∇f(x) = -exp(-x[1]^2-σ*x[2]^2)*[-2*x[1];-2*σ*x[2]]
  f(x) = x[1]^4 + σ*x[2]^4
  ∇f(x) = [4*x[1]^3; 4σ*x[2]^3]
  x₀ = ones(2)./[1.0;σ]

  for ls in line_searches
    for mtd in methods
      x, iter, nf, ng, X = mtd(f, ∇f, x₀, history=true,
          line_search=ls, max_iter=100, tol=1e-4)
      m = size(X)[2]
      F = [f(X[:,i]) for i=1:m]

      N = 200

      # Axis
      xm = min(minimum(X[1,:]),-2.0); xM = max(maximum(X[1,:]),2.0)
      ym = min(minimum(X[2,:]),-2.0); yM = max(maximum(X[2,:]),2.0)
      w = 0.01
      xm, xM = xm - w*(xM-xm), xM + w*(xM-xm)
      ym, yM = ym - w*(yM-ym), yM + w*(yM-ym)

      tx = linspace(xm,xM,N);
      ty = linspace(ym,yM,N);
      p = Gadfly.plot(layer(z=(x,y)->f([x;y]), x=tx, y=ty,
          Geom.contour(levels=F)),
          layer(x=X[1,:], y=X[2,:], Geom.line))
      println("$mtd $ls")
      draw(PDF("$path/$mtd-$ls.pdf", 4inch, 3inch), p)
    end
  end
end

example2d_nonlinear()
