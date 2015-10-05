include("methods.jl")
include("aux.jl")

methods = [cauchy, random_decrease, barzilai_borwein, short_step]
dims = [2 3 4 5 10 50 100 150 200]
iters = zeros(Int, length(dims),length(methods))
nMVs = zeros(Int, length(dims),length(methods))
x0funs = [Λ->(ones(Λ)./Λ)/norm(ones(Λ)./Λ), Λ->ones(Λ)/norm(ones(Λ)),
          Λ->(ones(Λ)./sqrt(Λ))/norm(ones(Λ)./sqrt(Λ))]
#x0funs = [Λ->ones(Λ)./Λ, Λ->ones(Λ), Λ->ones(Λ)./sqrt(Λ)]
σ = 0.001

for x0fun = x0funs
  println("### New x₀")
  for (idim, n) = enumerate(dims)
    Λ = linspace(σ,1.0,n)
    x₀ = x0fun(Λ)
    for (imethod, method) in enumerate(methods)
      x, iter, nMV = method(diagm(Λ), zeros(n), x₀)
      iters[idim,imethod] = iter
      nMVs[idim,imethod] = nMV
    end
  end

  print_table("iters", "n", dims, "methods", methods, iters)
  println("")
  print_table("nMVs", "n", dims, "methods", methods, nMVs)
end

