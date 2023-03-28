using Plots

function perfprof(A)
  (m, n) = size(A)
  a = minimum(A, dims=2)# mostra os minimos de cada problema
  R = A ./ a #deixar a coluna de quem for melhor com entrada 1 
  sum(R .== 1, dims=1) / m #ver a porcentagem em que o método é melhor
  sum(R .<= 1.5, dims=1) / m #verifica até 50% do tempo extra comparado ao minimo pra resolver o problema
  if maximum(R) < Inf
    r = maximum(R)
  else
    r = unique(sort(R[:]))[end-1]
  end
  ρ(t, k) = sum(R[:, k] .<= t) / m
  plot(leg=:bottomright, xaxis=:log)
  plot!(t -> ρ(t, 1), 1.0, 1.1r, lab="LMCLASSICO")
  plot!(t -> ρ(t, 2), 1.0, 1.1r, lab="LMSeq")
  plot!(t -> ρ(t, 3), 1.0, 1.1r, lab="CGA")

  ylims!(0, 1.1)

  savefig("3alg3d_40.png")
end

function ss(A,x)
  (m,n) = size(A)
  E = 1.0e-4
  y = [x[1], x[2]]
  h = zeros(m)
  for i = 1:m
    h[i] = norm(A[i,:] - y)^2 - x[3]^2
    if abs(h[i]) < E
      println(h[i], i)
    end
  end
  return h
end
