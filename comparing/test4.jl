include("methods.jl")
include("aux.jl")

using Gadfly

#methods = [cauchy, barzilai_borwein, short_step]
methods = [cauchy]

σ = 0.1
Λ = [σ;1.0]

x₀ = ones(Λ)./sqrt(Λ)
x₀ = 2*x₀/norm(x₀)
println("x₀ = $x₀")
x, iter, nMV, X = cauchy(diagm(Λ), [0.0;0.0], x₀, history=true,
    max_iter=100, tol=0.02)

m = size(X)[2]
F = [0.5*dot(X[:,i],Λ.*X[:,i]) for i=1:m]

N = 200
t = linspace(-2.0,2.0,N);
#p = plot(z=(x,y)->0.5*x^2*Λ[1]+0.5*y^2*Λ[2], x=t, y=t,
    #Geom.contour(levels=linspace(0.0, 2.0, 10)))
p = plot(layer(z=(x,y)->0.5*x^2*Λ[1]+0.5*y^2*Λ[2], x=t, y=t,
      Geom.contour(levels=F)),
#      Geom.contour(levels=linspace(0.0, 2.0, 10))),
    layer(x=X[1,:], y=X[2,:], Geom.line))
draw(PDF("test.pdf", 4inch, 3inch), p)

#for i = 1
