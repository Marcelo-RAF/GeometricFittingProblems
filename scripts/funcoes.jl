function fline(x, data)
  (m, n) = size(data)
  f = zeros(m)
  for i = 1:m
    f[i] = x[1] * data[i, 1] + x[2] - data[i, 2]
  end
  return f
end

function jline(x, data)
  (m, n) = size(data)
  tam = length(x)
  J = zeros(m, tam)
  J[:, 1] = data[:, 1]
  J[:, 2] = ones(m)
  return J
end

function fcubic(x, data)
  (m, n) = size(data)
  f = zeros(m)
  for i = 1:m
    f[i] = x[1] * data[i, 1]^3 + x[2] * data[i, 1]^2 + x[3] * data[i, 1] + x[4] - data[i, 2]
  end
  return f
end

function jcubic(x, data)
  (m, n) = size(data)
  tam = length(x)
  J = zeros(m, tam)
  J[:, 1] = data[:, 1] .^ 3
  J[:, 2] = data[:, 1] .^ 2
  J[:, 3] = data[:, 1]
  J[:, 4] = ones(m)
  return J
end

function fsphere(xinit, data)
  h = data
  (m, n) = size(data)
  r = zeros(m)
  for i = 1:m
    for j = 1:n
      r[i] = (h[i, j] - xinit[j])^2 + r[i]
    end
    r[i] = r[i] - xinit[end]^2
  end
  return r
end

function jsphere(xinit, data)
  h = data
  (m, n) = size(data)
  J = zeros(m, n + 1)
  for i = 1:m
    for j = 1:n
      J[i, j] = -2 * (h[i, j] - xinit[j])
    end
    J[i, end] = -2 * xinit[end]
  end
  return J
end

function fplane(xi, data)
  h = data
  (m, n) = size(data)
  r = zeros(m)
  for i = 1:m
    for j = 1:n
      r[i] = r[i] + xi[j] * h[i, j]
    end
    r[i] = r[i] - xi[end]
  end
  return r
end

function jplane(xi, data)
  h = data
  (m, n) = size(data)
  J = zeros(m, n + 1)
  for i = 1:m
    for j = 1:n
      J[i, j] = h[i, j]
    end
    J[i, end] = -1.0
  end
  return J
end

function fcircle(x, P)
  (m, n) = size(P)
  r = zeros(m)
  a = zeros(m)
  for i = 1:m
    a[i] = (dot(P[i, :] - x[4:6], x[1:3]))^2
    for j = 1:n
      r[i] = r[i] + (P[i, j] - x[3+j])^2
    end
    r[i] = (r[i] - x[7]^2)^2 + a[i]
  end
  return r
end


function jcircle(x, P)
  (m, n) = size(P)
  J = zeros(m, 7)
  a = zeros(m)
  h = zeros(m)
  for i = 1:m
    for j = 1:n
      h[i] = h[i] + (P[i, j] - x[3+j])^2
    end
    h[i] = h[i] - x[end]^2
  end
  for i = 1:m
    a[i] = (dot(P[i, :] - x[4:6], x[1:3]))
    for j = 1:n
      J[i, j] = 2 * (P[i, j] - x[j+3]) * a[i]
      J[i, j+3] = -4 * (P[i, j] - x[j+3]) * h[i] - 2 * x[j] * a[i]
    end
    J[i, end] = -4 * x[end] * h[i]
  end
  return J #sum(J, dims=1)[:, :]#J
end


function fexponencial(x, data)
  (m, n) = size(data)
  f = zeros(m)
  # x[1] * exp(- x[2] * t)
  for i = 1:m
    f[i] = x[1] * exp(-x[2] * data[i, 1]) - data[i, 2]
  end
  return f
end

function jexponencial(x, data)
  (m, n) = size(data)
  tam = length(x)
  J = zeros(m, tam)
  for i = 1:m
    J[i, 1] = exp(-x[2] * data[i, 1])
    J[i, 2] = -data[i, 1] * x[2] * exp(-x[2] * data[i, 1])
  end
  return J
end

function sort_funcion_res(x, model, data, nout)
  P = data
  (n, m) = size(data)
  v = zeros(n)
  for i = 1:n
    v[i] = (prob.model(x, prob.data[i, :]))^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux
        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end

function sort_exponencial(P, x, nout)
  n = length(P[:, 1])
  m = length(P[1, :])
  v = zeros(n)
  for i = 1:n
    v[i] = (x[1] * exp(-x[2] * P[i, 1]) - P[i, 2])^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux
        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end


function sort_line(P, x, nout)
  n = length(P[:, 1])
  m = length(P[1, :])
  v = zeros(n)
  for i = 1:n
    v[i] = (x[1] * P[i, 1] + x[2] - P[i, 2])^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux
        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end

function sort_cubic(P, x, nout)
  n = length(P[:, 1])
  m = length(P[1, :])
  v = zeros(n)
  for i = 1:n
    v[i] = (x[1] * P[i, 1]^3 + x[2] * P[i, 1]^2 + x[3] * P[i, 1] + x[4] - P[i, 2])^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux
        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end

function sort_sphere_res(P, x, nout)
  n = length(P[:, 1])
  m = length(P[1, :])
  v = zeros(n)
  for i = 1:n
    for j = 1:m
      v[i] = v[i] + (P[i, j] - x[j])^2
    end
    v[i] = (v[i] - x[end]^2)^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux

        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end

function sort_plane_res(P, x, nout)
  n = length(P[:, 1])
  m = length(P[1, :])
  v = zeros(n)
  for i = 1:n
    for j = 1:m
      v[i] = v[i] + P[i, j] * x[j]
    end
    v[i] = (v[i] - x[end])^2
  end
  indtrust = [1:n;]
  for i = 1:n-nout+1
    for j = i+1:n
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux

        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  #    println(indtrust[n-nout+1:n])
  return P[indtrust[1:n-nout], :], sum(v[1:n-nout])
end

function sort_circle_res(P, x, nout)
  N = length(P[:, 1])
  M = length(P[1, :])
  v = zeros(N)
  a = zeros(N)
  #for i = 1:N
  #       v[i] = (norm(P[i, :] - x[4:6])^2-x[7]^2)^2 + (dot(P[i,:]-x[4:6],x[1:3]))^2  
  #end
  for i = 1:N
    a[i] = abs(dot(P[i, :] - x[4:6], x[1:3]))
    for j = 1:M
      v[i] = v[i] + (P[i, j] - x[3+j])^2 #corrigir aqui
    end
    v[i] = abs(v[i] - x[7]^2) + a[i]
  end
  indtrust = [1:N;]
  for i = 1:N-nout+1    #1:N-nout+1
    for j = i+1:N    #i+1:N
      if v[i] > v[j]
        aux = v[j]
        v[j] = v[i]
        v[i] = aux

        aux2 = indtrust[j]
        indtrust[j] = indtrust[i]
        indtrust[i] = aux2
      end
    end
  end
  return P[indtrust[1:N-nout], :], sum(v[1:N-nout])
end






