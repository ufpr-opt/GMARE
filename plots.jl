using Winston

# Function to plot |gᵢ|
#=
function plot_g(Λ::Vector, x::Vector; filename = "test1")
  n = length(x)
  σ = minimum(Λ)
  hold(false)
  semilogy(linspace(σ,1,n), abs(Λ.*x), "o-r")
  ylim(1e-10,10)
  savefig("$filename.png")
end

function plot_g(Λ::Vector, x::Vector, y::Vector; filename = "test1")
  n = length(x)
  σ = minimum(Λ)
  hold(false)
  semilogy(linspace(σ,1,n), abs(Λ.*x), "--b")
  hold(true)
  semilogy(linspace(σ,1,n), abs(Λ.*y), "o-r")
  ylim(1e-10,10)
  savefig("$filename.png")
end
=#

function plot_g(Λ::Vector, x::Vector, y::Vector, z::Vector, μ::Real;
    filename = "tmp")
  n = length(x)
  hold(false)
  semilogy(Λ, abs(Λ.*x), ":k")
  hold(true)
  μ = Λ[1]*x[1]/(x[1] - y[1])
  semilogy([μ; μ], [1.0; 1e-10], "--b")
  semilogy(Λ, abs(Λ.*y), "--b")
  μ = Λ[1]*y[1]/(y[1] - z[1])
  semilogy([μ; μ], [1.0; 1e-10], "-r")
  semilogy(Λ, abs(Λ.*z), "o-r")
  xlim(minimum(Λ), maximum(Λ))
  ylim(1e-10,10)
  savefig("$filename.png")
end
