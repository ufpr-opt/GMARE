include("../src/methods.jl")
include("../src/aux.jl")

using Gadfly

function example2d()
  methods = [cauchy, barzilai_borwein, alternate_cauchy, dai_yuan,
      conjugate_gradient]

  path = createpath("2d-examples")
  σ = 0.1
  Λ = [σ;1.0]

  x₀ = ones(Λ)./sqrt(Λ)
  x₀ = 2*x₀/norm(x₀)
  println("x₀ = $x₀")

  for mtd in methods
    x, iter, nMV, X = mtd(diagm(Λ), [0.0;0.0], x₀, history=true,
    max_iter=100, tol=0.02)

    m = size(X)[2]
    F = [0.5*dot(X[:,i],Λ.*X[:,i]) for i=1:m]

    N = 200
    t = linspace(-2.0,2.0,N);
    p = Gadfly.plot(layer(z=(x,y)->0.5*x^2*Λ[1]+0.5*y^2*Λ[2], x=t, y=t,
    Geom.contour(levels=F)),
    layer(x=X[1,:], y=X[2,:], Geom.line))
    draw(PDF("$path/$mtd.pdf", 4inch, 3inch), p)
  end
end

example2d()
