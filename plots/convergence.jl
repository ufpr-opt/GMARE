include("../src/methods.jl")
include("../src/aux.jl")

using Winston
using Formatting

function plots()
  methods = [barzilai_borwein, short_step, alternate_short_step,
      dai_yuan, alternate_dai_yuan, conjugate_gradient]
  #methods = [cauchy, barzilai_borwein, alternate_cauchy, alternate_short_step,
      #dai_yuan, alternate_dai_yuan, conjugate_gradient]
  #methods = [barzilai_borwein, alternate_short_step, short_step]
  colors = ["black", "red", "blue"]
  linekinds = ["solid", "dashed", "dotted"]

  rnd_seed = 1
  path = createpath("convergence")

  σ = 0.001
  #n_values = [10 100 1000]
  n_values = [1000]
  tol = 5e-5

  #for hist_nmv in [true, false]
  for hist_nmv in true
    for n in n_values
      Λ = linspace(σ, 1.0, n)

      x₀ = ones(n)./sqrt(Λ)

      M = 300
      hold(false)
      f_plt = FramedPlot()
      g_plt = FramedPlot()
      for (i_mtd, mtd) in enumerate(methods)
        srand(rnd_seed)
        if mtd == short_step || mtd == alternate_short_step ||
            mtd == alternate_dai_yuan
          x, iter, nMV, X = mtd(diagm(Λ), zeros(n), x₀, history=true,
          max_iter=10*n, hist_nmv=hist_nmv, tol=tol)
        else
          x, iter, nMV, X = mtd(diagm(Λ), zeros(n), x₀, history=true, max_iter =
          10*n, tol=tol)
        end
        if mtd != cauchy
          if hist_nmv && nMV > M
            M = nMV
          elseif !hist_nmv && iter+1> M
            M = iter + 1
          end
        end

        m = size(X, 2)
        F = zeros(m)
        G = zeros(m)
        for i = 1:m
          F[i] = 0.5*dot(X[:,i],Λ.*X[:,i])
          G[i] = norm(Λ.*X[:,i])
        end

        t = 1:m

        i = (i_mtd-1)%length(colors) + 1
        j = div(i_mtd-1, length(linekinds)) + 1

        c = Curve(t, F, "color", colors[i], "linekind", linekinds[j])
        setattr(c, "label", replace(string(mtd), "_", " "))
        add(f_plt, c)

        c = Curve(t, G, "color", colors[i], "linekind", linekinds[j])
        setattr(c, "label", replace(string(mtd), "_", " "))
        add(g_plt, c)
      end
      setattr(f_plt, "xrange", (1, M))
      setattr(g_plt, "xrange", (1, M))
      setattr(f_plt, "ylog", true)
      setattr(g_plt, "ylog", true)
      # Ignore cauchy iterations
      add(f_plt, Legend(0.6, 0.95, f_plt.content1.components))
      add(g_plt, Legend(0.6, 0.95, g_plt.content1.components))
      # Filename generation
      L = length(string(maximum(n_values)))
      number = format("{1:>0$(L)d}", n)
      if hist_nmv
        filename = "nmv"
        #title("function value by number of Matrix-Vector mult")
      else
        filename = "iter"
        #title("function value by number of iterations")
      end
      savefig(f_plt, "$path/function-decrease-$filename-$number.png", "width",
      800, "height", 600)
      savefig(g_plt, "$path/gradient-$filename-$number.png", "width", 800,
      "height", 600)
    end
  end
end
plots()
