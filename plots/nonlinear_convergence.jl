include("../src/methods.jl")
include("../src/aux.jl")

using Winston
using Formatting

function nonlinear_plots()
  methods = [cauchy, barzilai_borwein]
  line_searches = [armijo, golden_search]
  colors = ["black", "red", "blue"]
  linekinds = ["solid", "dashed", "dotted"]
  width = 800
  height = 600

  rnd_seed = 1
  path = createpath("nonlinear-convergence")
  n = 1000
  σ = 0.001
  Λ = linspace(σ, 1.0, n)
  p = 2

  #f(x) = 1-exp(-x[1]^2 - σ*x[2]^2)
  #∇f(x) = -exp(-x[1]^2 - σ*x[2]^2)*[-2*x[1];-2*σ*x[2]]
  #∇²f(x) = -exp(-x[1]^2 - σ*x[2]^2) * (
    #[4*x[1]^2 4*σ*x[1]*x[2]; 4*σ*x[1]*x[2] 4*σ^2*x[2]^2]
    #+ [-2 0;0 -2*σ])
  #f(x) = x[1]^4 + σ*x[2]^4
  #∇f(x) = [4*x[1]^3; 4σ*x[2]^3]
  f(x) = sum(Λ.*x.^p)/p
  ∇f(x) = Λ.*x.^(p-1)
  x₀ = ones(n)./sqrt(Λ)
  tol = 1e-6

  M = 0
  hold(false)
  f_plt = FramedPlot()
  g_plt = FramedPlot()

  count = 0
  for line_search in line_searches
    for mtd in methods
      srand(rnd_seed)
      x, iter, nf, ng, X, A = mtd(f, ∇f, x₀, history=true, max_iter = 10000,
          tol=tol, line_search=line_search)
      if mtd != cauchy && iter + 1 > M
        M = iter + 1
      end

      m = size(X, 2)
      F = zeros(m)
      G = zeros(m)
      for i = 1:m
        F[i] = f(X[:,i])
        G[i] = norm(∇f(X[:,i]))
      end

      t = 1:m

      i = count%length(colors) + 1
      j = div(count, length(linekinds)) + 1
      count += 1

      c = Curve(t, F, "color", colors[i], "linekind", linekinds[j])
      name = replace(string("$mtd-$line_search"), "_", " ")
      setattr(c, "label", name)
      add(f_plt, c)

      c = Curve(t, G, "color", colors[i], "linekind", linekinds[j])
      setattr(c, "label", name)
      add(g_plt, c)
    end
  end
  setattr(f_plt, "xrange", (1, M))
  setattr(g_plt, "xrange", (1, M))
  setattr(f_plt, "ylog", true)
  setattr(g_plt, "ylog", true)
  # Ignore cauchy iterations
  add(f_plt, Legend(0.6, 0.95, f_plt.content1.components))
  add(g_plt, Legend(0.6, 0.95, g_plt.content1.components))
  # Filename generation
  xtype = "iter"
  savefig(f_plt, "$path/function-decrease.png", "width",
      width, "height", height)
  savefig(g_plt, "$path/gradient.png", "width", width,
      "height", height)
end

nonlinear_plots()
