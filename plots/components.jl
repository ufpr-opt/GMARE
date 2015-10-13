include("../src/methods.jl")
include("../src/aux.jl")

using Winston

function plot_components()
  #methods = [barzilai_borwein, short_step, alternate_short_step,
      #dai_yuan, alternate_dai_yuan, conjugate_gradient]
  methods = [cauchy, barzilai_borwein, short_step, alternate_cauchy,
      alternate_short_step, dai_yuan, alternate_dai_yuan, conjugate_gradient]
  #methods = [cauchy, barzilai_borwein, alternate_cauchy, short_step,
  #    random_decrease, alternate_short_step]
  colors = ["black", "red", "blue"]
  linekinds = ["solid", "dashed", "dotted"]

  path="output/components"
  if !ispath("output")
    mkdir("output")
  elseif ispath(path)
    rm(path, recursive=true)
  end
  mkdir(path)
  rnd_seed = 1

  tol = 1e-4
  σ = 0.01
  #n_values = [5 10 50 100 200 400 800]
  n_values = [50]
  L = maximum(n_values)

  for hist_nmv in [true, false]
    for n in n_values
      Λ = linspace(σ, 1.0, n)

      x₀ = ones(n)./sqrt(Λ)
      max_iter = 100

      for (i_mtd, mtd) in enumerate(methods)
        # Each method has a plot
        srand(rnd_seed)
        if mtd == short_step || mtd == alternate_short_step ||
            mtd == alternate_dai_yuan
          x, iter, nMV, X, A = mtd(diagm(Λ), zeros(n), x₀, history=true,
              max_iter=max_iter, hist_nmv=hist_nmv, tol=tol)
        else
          x, iter, nMV, X, A = mtd(diagm(Λ), zeros(n), x₀, history=true,
              max_iter=max_iter, tol=tol)
        end

        m = size(X, 2)

        for t = 3:m
          hold(false)
          semilogy(Λ, abs(Λ.*X[:,t-2]), ":k")
          hold(true)
          μ = 1/A[t-2]
          semilogy([μ; μ], [1.0; 1e-10], "--b")
          semilogy(Λ, abs(Λ.*X[:,t-1]), "--b")
          μ = 1/A[t-1]
          semilogy([μ; μ], [1.0; 1e-10], "-r")
          semilogy(Λ, abs(Λ.*X[:,t]), "o-r")
          xlim(0.0, maximum(Λ))
          ylim(1e-10,10)
          number = getname(t, L)
          if hist_nmv
            xtype = "nmv"
          else
            xtype = "iter"
          end
          savefig("$path/component-$mtd-$xtype-$n-$number.png")
        end
      end
    end
  end
end

plot_components()
