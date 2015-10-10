include("methods.jl")
include("aux.jl")

using Winston
using Formatting

methods = [cauchy, barzilai_borwein, alternate_cauchy, short_step,
    random_decrease]
colors = ["k", "r", "b", "k--", "r--", "b--"]

σ = 0.001
n_values = [5 10 50 100 200 400]

for n in n_values
  Λ = linspace(σ, 1.0, n)

  x₀ = ones(Λ)./sqrt(Λ)

  M = 0
  hold(false)
  for (i_mtd, mtd) in enumerate(methods)
    x, iter, nMV, X = mtd(diagm(Λ), zeros(n), x₀, history=true, max_iter = 10*n)
    if i_mtd != 1 && iter > M
      M = iter
    end

    m = size(X, 2)
    F = zeros(m)
    for i = 1:m
      F[i] = 0.5*dot(X[:,i],Λ.*X[:,i])
    end

    t = 1:m
    if i_mtd == 1
      semilogy(t, F, colors[i_mtd])
    else
      oplot(t, F, colors[i_mtd], ylog=true)
    end
  end
  xlim(1, M)
  # Ignore cauchy iterations
  legend([replace(string(x),"_"," ") for x in methods], (0.6,0.95))
  # Filename generation
  L = length(string(maximum(n_values)))
  number = format("{1:>0$(L)d}", n)
  savefig("function-decrease-$number.png")

end
