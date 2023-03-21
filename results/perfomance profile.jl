using Plots

function perfprof(A)
  (m, n) = size(A)
  a = minimum(A, dims=n)# mostra os minimos de cada problema
  R = A ./ a #deixar a coluna de quem for melhor com entrada 1 
  sum(R .== 1, dims=1) / m #ver a porcentagem em que o método é melhor
  sum(R .<= 1.5, dims=1) / m #verifica até 50% do tempo extra comparado ao minimo pra resolver o problema
  ρ(t, k) = sum(R[:, k] .<= t) / m
  plot(leg=:bottomright)
  plot!(t -> ρ(t, 1), 1.0, 100.0, lab="CGA")
  plot!(t -> ρ(t, 2), 1.0, 100.0, lab="LM")
  ylims!(0, 1)


end