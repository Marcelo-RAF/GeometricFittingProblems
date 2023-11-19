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
    return sum(J, dims=1)[:, :]
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
    return sum(J, dims=1)[:, :]
end

function residexp(x, P)
    (m, n) = size(P)
    r = zeros(m)
    for i = 1:m
        r[i] = ((x[1] * exp(-x[2] * P[i, 1]) - P[i, 2])^2)^2
    end
    return sum(r)
end

function residcubic(x, P)
    (m, n) = size(P)
    r = zeros(m)
    for i = 1:m
        r[i] = (x[1] * P[i, 1]^3 + x[2] * P[i, 1]^2 + x[3] * P[i, 1] + x[4] - P[i, 2])^2
    end
    return sum(r)
end

function residline(x, P)
    (m, n) = size(P)
    r = zeros(m)
    for i = 1:m
        r[i] = (x[1] * P[i, 1] + x[2] - P[i, 2])^2
    end
    return sum(r)
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
    return sum(r)
end

function residplane(x, P)
    n = length(P[:, 1])
    m = length(P[1, :])
    v = zeros(n)
    for i = 1:n
        for j = 1:m
            v[i] = v[i] + P[i, j] * x[j]
        end
        v[i] = (v[i] - x[end])^2
    end
    return sum(v)
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

function residconformal(x, P)
    (m, n) = size(P)
    h = zeros(m)
    D = [P'; ones(1, m)]
    v = [0.5 * norm(D[1:n, i], 2)^2 for i = 1:m]
    D = [D; v']'
    for i = 1:m
        for j = 1:n
            h[i] = h[i] + D[i, j] * x[j]
        end
        h[i] = (h[i] - x[end] - x[end-1] * v[i])^2
    end
    return sum(h)
end


function residcircle(x, P)
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
    return sum(r)
end