function geradorcircle()
  n = 1000  # tamanho da lista desejada
  a = -120  # limite inferior do intervalo
  b = 120  # limite superior do intervalo
  x = -50
  y = 50
  c1 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
  c2 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
  c3 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
  v1 = round.(float(rand(n) .* (y - x) .+ x), digits=1)
  v2 = round.(float(rand(n) .* (y - x) .+ x), digits=1)
  v3 = round.(float(rand(n) .* (y - x) .+ x), digits=1)
  r = float(rand(5:170, 1000))
  npts = float(rand(10:2000, 401))
  for i = 1:400
    build_problem("circle3d", [1.0, 1.0], [c1[i], c2[i], c3[i], r[i], v1[i], v2[i], v3[i], -v2[i], v1[i], 0.0, npts[i], 0])
  end
end

function geradoraut(h)
  n = 120  # tamanho da lista desejada
  a = -120  # limite inferior do intervalo
  b = 120  # limite superior do intervalo
  c1 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
  c2 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
  c3 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
  r = float(rand(5:170, 1000))
  npts = float(rand(8:50, 120))
  nout = float([floor(Int, h * x) for x in npts])
  for i = 1:100
    build_problem("sphere3D", [1.0, 1.0], [c1[i], c2[i], c3[i], r[i], npts[i], nout[i]])
  end
end