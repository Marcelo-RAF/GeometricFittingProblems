function gradtangential(x, data)
    (N, n) = size(data)
    v = [0.5 * norm(data[i, :],2)^2 for i = 1:N]
    h = [norm(data[i, :],2)^2 for i = 1:N]
    pp = zeros(N)
    for i = 1:N
        pp[i] = dot(data[i, :], x[1:end-2]) - x[end-1] - x[end] * v[i]
    end
    J = zeros(N, n + 2)
    for i = 1:N
        for j = 1:n
            J[i, j] = 2*pp[i]*data[i,j] 
        end
        J[i, end-1] = -2*pp[i]
        J[i, end] = -pp[i]*h[i]
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
            J[i, j] = -4 * p[i] * (h[i, j] - x[j])*(1/r[i])
        end
        J[i, end] = -2 * p[i]
    end
    return sum(J, dims=1)[:, :]
end

function gradalgebric(x, h)
    (m, n) = size(h)
    r = zeros(m)
    J = zeros(m, n + 1)
    
    for i=1:m
        for j=1:n
            r[i] = r[i] + (h[i, j] - x[j])^2 
        end
        r[i] = r[i]-x[end]^2
    end

    for i = 1:m
        for j = 1:n
            J[i, j] = -4 * r[i] * (h[i, j] - x[j])
        end
        J[i, end] = -4 * x[end] * r[i]
    end
    return sum(J, dims=1)[:, :]
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
