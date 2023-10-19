using PyPlot

function gradtangential(x, data)
    (N, n) = size(data)
    v = [0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    h = [norm(data[i, :], 2)^2 for i = 1:N]
    pp = zeros(N)
    for i = 1:N
        pp[i] = dot(data[i, :], x[1:end-2]) - x[end-1] - x[end] * v[i]
    end
    J = zeros(N, n + 2)
    for i = 1:N
        for j = 1:n
            J[i, j] = 2 * pp[i] * data[i, j]
        end
        J[i, end-1] = -2 * pp[i]
        J[i, end] = (pp[i])
    end
    return sum(J, dims=1)[:, :]
end



function gradgeometric(x, h)
    (m, n) = size(h)
    J = zeros(m, n + 1)
    p = zeros(m)
    r = zeros(m)
    for i = 1:m
        for j = 1:n
            p[i] = p[i] + (h[i, j] - x[j])^2
        end
        r[i] = sqrt(p[i])
        p[i] = sqrt(p[i]) - x[end]
    end
    J = zeros(m, n + 1)
    for i = 1:m
        for j = 1:n
            J[i, j] = -4 * p[i] * (h[i, j] - x[j]) * (1 / r[i])
        end
        J[i, end] = -2 * p[i]
    end
    return norm(sum(J, dims=1)[:, :])
end

function gradalgebric(x, h)
    (m, n) = size(h)
    r = zeros(m)
    J = zeros(m, n + 1)

    for i = 1:m
        for j = 1:n
            r[i] = r[i] + (h[i, j] - x[j])^2
        end
        r[i] = r[i] - x[end]^2
    end

    for i = 1:m
        for j = 1:n
            J[i, j] = -4 * r[i] * (h[i, j] - x[j])
        end
        J[i, end] = -4 * x[end] * r[i]
    end
    return norm(sum(J, dims=1)[:, :])
end




function residtangential(x, data)
    (N, n) = size(data)
    v = [0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    pp = zeros(N)
    for i = 1:N
        pp[i] = sqrt(2 * abs(dot(data[i, :], x[1:end-2]) - x[end-1] - x[end] * v[i]))
    end
    return sum(pp)
end

function residalgebric(x, h)
    (m, n) = size(h)
    r = zeros(m)
    for i = 1:m
        for j = 1:n
            r[i] = r[i] + (h[i, j] - x[j])^2
        end
        r[i] = (r[i] - x[end]^2)^2
    end
    sum(r[1:m])
end

function residgeometric(x, h)
    (m, n) = size(h)
    r = zeros(m)
    for i = 1:m
        for j = 1:n
            r[i] = r[i] + (h[i, j] - x[j])^2
        end
        r[i] = abs(sqrt(r[i]) - x[end])
    end
    return sum(r[1:m])
end

function residcirc(s1, s2, data)
    (N, n) = size(data)
    D = [data'; -ones(1, N)]
    v = [-0.5 * norm(D[1:n, i], 2)^2 for i = 1:N]
    D = [D; v']
    aux = s1[end-1]
    s1[end] = s1[end-1]
    s1[end-1] = aux
    aux2 = s2[end-1]
    s2[end] = s2[end-1]
    s2[end-1] = aux2
    h = zeros(N)
    for i = 1:N
        h[i] = norm((dot(D'[i, :], s1) * s2 - dot(D'[i, :], s2) * s1))
    end
    return h
end
#pyplot()
function visusu(prob, a, b)
    #pyplot()#beckendpyplot
    plt = plot()
    if prob.name == "circle3d" || prob.name == "\tcircle3d"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], prob.data[:, 3], line=:scatter, aspect_ratio=:equal, lab="pontos do problema")
        n = 20
        u = range(0, stop=2 * pi, length=n)
        v = range(0, stop=pi, length=n)
        h1 = zeros(n)
        h2 = zeros(n)
        h3 = zeros(n)
        h4 = zeros(n)
        h5 = zeros(n)
        h6 = zeros(n)
        for i = 1:n
            h1[i] = b[1]
            h2[i] = b[2]
            h3[i] = b[3]
            h4[i] = a[1]
            h5[i] = a[2]
            h6[i] = a[3]
        end
        x = h1 .+ b[4] * cos.(u) * sin.(v)'
        y = h2 .+ b[4] * sin.(u) * sin.(v)'
        z = h3 .+ b[4] * cos.(v)'
        xs = h4 .+ a[4] * cos.(u) * sin.(v)'
        ys = h5 .+ a[4] * sin.(u) * sin.(v)'
        zs = h6 .+ a[4] * cos.(v)'

        plot!(plt, xs, ys, zs, aspect_ratio=:equal, color=:blue, label="solução do algoritmo")
        plot!(plt, x, y, z, aspect_ratio=:equal, color=:red, label="solução perfeita")
        #wireframe!(xs, ys, zs, aspect_ratio=:equal, color=:blue, label="solução do algoritmo")
        #wireframe!(x, y, z, aspect_ratio=:equal, color=:red, label="solução perfeita")
        return plt
    end
end