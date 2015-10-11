include("methods.jl")
include("aux.jl")

using Winston
using Formatting

methods = [cauchy, barzilai_borwein, alternate_cauchy, short_step,
    alternate_short_step, dai_yuan]
#methods = [cauchy, barzilai_borwein, alternate_cauchy, short_step,
#    random_decrease, alternate_short_step]
colors = ["k", "r", "b", "k--", "r--", "b--"]

rnd_seed = 1

σ = 0.0001
n_values = [5 10 50 100 200 400 800]
#n_values = [100]

for hist_nmv in [true, false]
  for n in n_values
    Λ = linspace(σ, 1.0, n)

    x₀ = ones(n)./sqrt(Λ)

    M = 0
    hold(false)
    for (i_mtd, mtd) in enumerate(methods)
      srand(rnd_seed)
      if mtd == short_step || mtd == alternate_short_step
        x, iter, nMV, X = mtd(diagm(Λ), zeros(n), x₀, history=true,
            max_iter=10*n, hist_nmv=hist_nmv, Ks = 2, Kc = 8,
            Ki = 10)
      else
        x, iter, nMV, X = mtd(diagm(Λ), zeros(n), x₀, history=true, max_iter = 10*n)
      end
      if i_mtd != 1
        if hist_nmv && nMV > M
          M = nMV
        elseif !hist_nmv && iter+1> M
          M = iter + 1
        end
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
    filename = "function-decrease"
    if hist_nmv
      filename = "$(filename)-nmv"
      title("function value by number of Matrix-Vector mult")
    else
      filename = "$(filename)-iter"
      title("function value by number of iterations")
    end
    savefig("$filename-$number.png")
  end
end
