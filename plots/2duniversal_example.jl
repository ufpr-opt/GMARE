include("../src/methods.jl")
include("../src/aux.jl")

using Gadfly

function example2d_universal()
  methods = [fgm]

  path = createpath("2d-universal-examples")
  σ = 10

  #f(x) = 1-exp(-x[1]^2-σ*x[2]^2)
  #∇f(x) = -exp(-x[1]^2-σ*x[2]^2)*[-2*x[1];-2*σ*x[2]]
  #f(x) = x[1]^4 + σ*x[2]^4
  #∇f(x) = [4*x[1]^3; 4σ*x[2]^3]
  #x₀ = ones(2)./[1.0;σ]

  function f(x)
    return max(x[1]^2+4*x[2]^2, 9*(x[1]-1)^2+5*(x[2]-1)^2)
  end

  function ∇f(x)
    f1 = x[1]^2 + 4*x[2]^2
    f2 = 9*(x[1]-1)^2 + 5*(x[2]-1)^2
    if f1 > f2
      return [2*x[1]; 8*x[2]]
    else
      return [18*(x[1]-1); 10*(x[2]-1)]
    end
  end
  
  x₀ = [1.5; 0.3]

  for mtd in methods
    x, fx, ∇fx, iter, ef, X = mtd(f, ∇f, x₀, history=true,
        max_iter=1000, tol=1e-1)
    println(iter)
    println(x)
    #X = X[:,10:20]
    m = size(X, 2)
    F = [f(X[:,i]) for i=5:100:m]

    N = 200

    # Axis
    xm = minimum(X[1,:]); xM = maximum(X[1,:])
    ym = minimum(X[2,:]); yM = maximum(X[2,:])
    w = 1
    xm, xM = xm - w*(xM-xm), xM + w*(xM-xm)
    ym, yM = ym - w*(yM-ym), yM + w*(yM-ym)

    tx = linspace(xm,xM,N);
    ty = linspace(ym,yM,N);
    p = Gadfly.plot(layer(z=(x,y)->f([x;y]), x=tx, y=ty,
        Geom.contour(levels=F)),
        layer(x=X[1,:], y=X[2,:], Geom.line))
    draw(PDF("$path/$mtd.pdf", 4inch, 3inch), p)
  end
end

example2d_universal()
