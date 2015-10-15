using Formatting

function print_table(title, rowlabel, rowvalues, collabel, colvalues, table)
  (m,n) = size(table)
  max_len_col(col) = maximum(map(x->length(string(x)),col))
  str_table = fill("", m+2,n+2)
  # Column 1
  c = length(rowlabel)
  for i = 1:m+2
    str_table[i,1] = format("{1:>$c}", "")
  end
  str_table[2+floor((m+1)/2),1] = rowlabel
  # Column 2
  c = max_len_col(rowvalues)
  for i = 1:2
    str_table[i,2] = format("{1:>$c}", "")
  end
  for i = 1:m
    str_table[i+2,2] = format("{1:>$c}", rowvalues[i])
  end
  # Other columns
  for j = 1:n
    c = max(max_len_col(table[:,j]), length(string(colvalues[j])))
    str_table[1,j+2] = format("{1:>$c}", "")
    str_table[2,j+2] = format("{1:>$c}", colvalues[j])
    for i = 1:m
      str_table[i+2,j+2] = format("{1:>$c}",string(table[i,j]))
    end
  end
  str_table[1,2+floor((n+1)/2)] = collabel

  # Printing
  println(title)
  bar = [join(fill("-",length(str_table[2,j]))) for j = 1:n+2]
  println(join(bar,"-|-"))
  print(join(str_table[1,1:2],"   "))
  print(" | ")
  println(join(str_table[1,3:n+2],"   "))

  print(join(str_table[1,1:2],"   "))
  print(" | ")
  println(join(str_table[2,3:n+2]," | "))
  bar = [join(fill("-",length(str_table[2,j]))) for j = 1:n+2]
  println(join(bar,"-|-"))
  println(join([join(str_table[i,:]," | ") for i=3:m+2],'\n'))
end

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

function createpath(name)
  dir = "output"
  path = "$dir/$name"
  if !ispath("output")
    mkdir("output")
  elseif ispath(path)
    rm(path, recursive=true)
  end
  mkdir(path)
  return path
end

function rowdot(A::Matrix, row::Int, x::Vector)
  (m,n) = size(A)
  if n != length(x)
    error("size(A,2) != length(x)")
  end
  if row < 0 || row > m
    error("1 ≦ row ≦ size(A,1)")
  end
  s = 0.0
  for j = 1:n
    @inbounds s += A[row,j]*x[j]
  end
  return s
end
