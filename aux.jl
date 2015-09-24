# Auxiliary functions

using Formatting

function getname(i::Int, N::Int)
  s = length(string(N)) # or floor(Int, log10(N))
  return format("{1:0$(s)d}-{2:0$(s)d}-{3:0$(s)d}", i-2, i-1, i)
end
