using Formatting

function print_table(title, rowlabel::String, rowvalues, collabel::String, colvalues, table)
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
