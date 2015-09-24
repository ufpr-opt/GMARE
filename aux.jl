# Auxiliary functions

using Formatting

function getname(i::Int, N::Int)
  s = length(string(N)) # or floor(Int, log10(N))
  return format("{1:0$(s)d}-{2:0$(s)d}-{3:0$(s)d}", i-2, i-1, i)
end

function uniform(a::Real, b::Real, n::Int)
  v = rand(n)
  v = a + (b-a)*(v-minimum(v))/(maximum(v) - minimum(v))
  return sort(v)
end

function normal(a::Real, b::Real, n::Int)
  v = randn(n)
  v = a + (b-a)*(v-minimum(v))/(maximum(v) - minimum(v))
  return sort(v)
end

function exp_normal(a::Real, b::Real, n::Int; weigh_left::Bool = true)
  if weigh_left
    v = exp(randn(n))
  else
    v = -exp(randn(n))
  end
  v = a + (b-a)*(v-minimum(v))/(maximum(v) - minimum(v))
  return sort(v)
end
