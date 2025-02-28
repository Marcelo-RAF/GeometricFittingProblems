module GeometricFittingProblems

using DelimitedFiles, LinearAlgebra, Plots

export load_problem, solve, build_problem, visualize, FitProbType, FitOutputType, sort_funcion_res, AACGA, AGCGA, ICGA, LOVOCGA

import Base.show


"""
    FitProbType

It is an immutable type used by main functions of this package

"""
struct FitProbType
    name::String
    data::Array{Float64,2}
    npts::Int
    nout::Int
    model::Function
    dim::Int
    cluster::Bool
    noise::Bool
    solution::Array{Float64,1}
    description::String
    sandbox::Vector{Vector{Any}}
    function FitProbType(_name, _data, _npts, _nout, _model, _dim, _cluster, _noise, _solution, _description)
        new(_name, _data, _npts, _nout, _model, _dim, _cluster, _noise, _solution, _description, [["NONE"]])
    end
end

struct FitOutputType
    status::Bool
    solution::Vector{Float64}
    niter::Int
    minimum::Float64
    feval::Int
end



"""
    load_problem(filename::String)

This function is used to load a problem from a csv file and convert to FitProbType. It is an important function because FitProbType is the unique supported format in this package. 

# Examples
```
julia-repl
julia> load_problem("toy.csv")

returns a FitProbType
```
"""
function load_problem(filename::String)
    prob_matrix = readdlm(filename, ':')
    (m, n) = size(prob_matrix)
    if m == 10
        return FitProbType(prob_matrix[1, 2], eval(Meta.parse(prob_matrix[2, 2])), prob_matrix[3, 2], prob_matrix[4, 2], eval(Meta.parse(prob_matrix[5, 2])), prob_matrix[6, 2], prob_matrix[7, 2], prob_matrix[8, 2], eval(Meta.parse(prob_matrix[9, 2])), prob_matrix[10, 2])
    elseif m == 11
        return FitProbType(prob_matrix[1, 2], eval(Meta.parse(prob_matrix[2, 2])), prob_matrix[3, 2], prob_matrix[4, 2], eval(Meta.parse(prob_matrix[5, 2])), prob_matrix[6, 2], prob_matrix[7, 2], prob_matrix[8, 2], eval(Meta.parse(prob_matrix[9, 2])), prob_matrix[10, 2], eval(Meta.parse(prob_matrix[11, 2])))
    else
        error("No type identified!!")
    end
end

"""
    AGCGA(data::Matrix, Method::String)

Essa função ajusta pontos do espaço R^n, dentro do conjunto {hiperplano, hiperesfera, hipercírculo, hiperreta} por um modelo geométrico aproximado de um método de Álgebra Geométrica no Epaço Conformal que minimiza as distâncias dos pontos dados a um dos objetos supracitados. Os dados de entrada são: uma matriz em que cada linha é um vetor do espaço euclidiano e o objeto a ser encontrado. 

# Examples
```
julia-repl
julia> prob = load_problem("toy.csv")
julia> AGCGA(prob.data)

returns a vector
```
"""


function AGCGA(data, object::String, ε=1.0e-5) #algoritmo dorst esferas
    (N, n) = size(data)
    v = [0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    D = [data'; v'; ones(1, N)]
    J = copy(D')
    H = -copy(J[:, n+1])
    J[:, n+1] = -J[:, n+2]
    J[:, n+2] = H
    DDt = D * D'
    aux = -copy(DDt[:, n+1])
    DDt[:, n+1] = -DDt[:, n+2]
    DDt[:, n+2] = aux
    p = (1.0 / N)
    P = p .* (DDt)
    F = eigen(P)
    if object == "sphere" || object == "plane"
        return F.vectors[:, 2], F.values[2], F.values[3]
    end
    if object == "line" || object == "circle"
        return F.vectors[:, 2], F.vectors[:, 3]
    end#hhh push!(center, √(norm(center, 2)^2 - 2.0 * xnorm[end]))
end

function simetrica(D)
    (m, n) = size(D)
    H = zeros(m, m)
    for i = 1:m
        for j = i:m
            H[i, j] = dot(D[i, :], D[j, :])
        end
    end
    for j = 1:m-1
        for i = j+1:m
            H[i, j] = H[j, i]
        end
    end
    return H
end

"""
    AGCGA(data::Matrix)

Essa função ajusta pontos do espaço R^n, dentro do conjunto {hiperplano, hiperesfera, hipercírculo, hiperreta} por um modelo algébrico aproximado de um método de Álgebra Geométrica no Epaço Conformal que minimiza as distâncias dos pontos dados a um dos objetos supracitados. Os dados de entrada são: uma matriz em que cada linha é um vetor do espaço euclidiano e o objeto a ser encontrado. 

# Examples
```
julia-repl
julia> prob = load_problem("toy.csv")
julia> AACGA(prob.data)

returns a vector
```
"""


function AACGA(data, object::String, ε=1.0e-5)
    (N, n) = size(data)
    v = [-0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    D = [data'; -ones(1, N); v']
    return D
    Dd = simetrica(D)
    F = eigen(Dd)
    indmin = 1
    valmin = F.values[1]
    for i = 2:n
        if abs(valmin) > abs(F.values[i])
            if F.values[i] > -ε
                indmin = i
                valmin = F.values[i]
            end
        end
    end
    if valmin < -ε
        error("P does not have postive eigen value!")
    end
    if nullspace(Dd) == zeros(n + 2, 0)
        if object == "sphere" || object == "plane"
            return F.vectors[:, indmin]
        end
        if object == "line" || object == "circle"
            return F.vectors[:, indmin], F.vectors[:, indmin+1]
        end
    else
        return nullspace(Dd)
    end
end

"""
    ICGA(data::Matrix, object::String)

Essa função ajusta pontos do espaço R^n, dentro do conjunto {hiperplano, hiperesfera, hipercírculo, hiperreta} que minimiza as distâncias tangenciais para esferas e círculos e as distâncias ortogonais para planos e retas. Os dados de entrada são: conjunto de pontos do espaço euclidiano e o objeto que deve ser encontrado. 

# Examples
```
julia-repl
julia> prob = load_problem("toy.csv")
julia> ICGA(prob.data, "sphere")

returns a vector
```
"""

function ICGA(data, object::String)
    (N, n) = size(data)
    v = [-0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    D = [data'; -ones(1, N); v']
    Dd = 1 / N * simetrica(D)
    if object == "sphere"
        Dd[end, :] .= 0.0
        np = nullspace(Dd)
        #p = np / np[end]
        #centernp = np[1:end-2]
        #s = push!(centernp, √(norm(centernp, 2)^2 - 2.0 * np[end-1]))
        return np
    end
    if object == "plane"
        B = Dd[1:end-2, 1:end-2]
        u = Dd[1:end-2, end-1]
        a = Dd[end-1, end-1]
        a2 = Dd[end-1, end]
        H = B - (u * u') / a
        F = eigen(H)
        vn = F.vectors[:, 1]
        λ1 = F.values[1]
        λ2 = F.values[2]
        d = -(u' * vn) / a
        π = [vn; d; 0]
        #coef = (-w'*vn)/a2

        return π, λ1, λ2#, coef
    end
    if object == "reta"
        B = Dd[1:end-2, 1:end-2]
        u = Dd[1:end-2, end-1]
        a = Dd[end-1, end-1]
        a2 = Dd[end-1, end]
        H = B - (u * u') / a
        F = eigen(H)
        display(F)
        vn = F.vectors[:, 1]
        vn2 = F.vectors[:, 2]
        λ1 = F.values[1]
        λ2 = F.values[2]
        d = -(u' * vn) / a
        d2 = -(u' * vn2) / a
        π = [vn; d; 0]
        pi2 = [vn2; d2; 0]
        return hcat(π, pi2), λ1, λ2
    end
    if object == "circle"
        B = Dd[1:end-2, 1:end-2]
        u = Dd[1:end-2, end-1]
        a = Dd[end-1, end-1]
        H = B - (u * u') / a
        F = eigen(H)
        vn = F.vectors[:, 1]
        d = -(u' * vn) / a
        π = [vn; d; 0]
        Px = copy(Dd)
        Px[end, :] .= 0.0
        np = nullspace(Px)
        return hcat(np, π)
    end
end

function conformalsort(P, x, nout)
    (m, n) = size(P)
    #display((m,n))
    h = zeros(m)
    D = [P'; ones(1, m)]
    v = [0.5 * norm(D[1:n, i], 2)^2 for i = 1:m]
    D = [D; v']'
    for i = 1:m
        for j = 1:n
            h[i] = h[i] + D[i, j] * x[j]
        end
        h[i] = (h[i] - x[end-1] - x[end] * v[i])^2
    end
    indtrust = [1:m;]
    for i = 1:m-nout+1
        for j = i+1:m
            if h[i] > h[j]
                aux = h[j]
                h[j] = h[i]
                h[i] = aux

                aux2 = indtrust[j]
                indtrust[j] = indtrust[i]
                indtrust[i] = aux2
            end
        end
    end
    #    println(indtrust[n-nout+1:n])

    return P[indtrust[1:m-nout], :], sum(h[1:m-nout]), P[indtrust[(m-nout)+1:m], :]
end

function conformalsort2(P, x, y, nout)
    (m, n) = size(P)
    #display((m,n))
    h1 = zeros(m)
    h2 = zeros(m)
    h = zeros(m)
    D = [P'; ones(1, m)]
    v = [0.5 * norm(D[1:n, i], 2)^2 for i = 1:m]
    D = [D; v']'
    for i = 1:m
        for j = 1:n
            h1[i] = h1[i] + D[i, j] * x[j]
            h2[i] = h2[i] + D[i, j] * y[j]
        end
        h[i] = ((h1[i] - x[end-1] - x[end] * v[i])^2 + (h2[i] - y[end-1] - y[end] * v[i])^2) / 2
    end
    indtrust = [1:m;]
    for i = 1:m-nout+1
        for j = i+1:m
            if h[i] > h[j]
                aux = h[j]
                h[j] = h[i]
                h[i] = aux

                aux2 = indtrust[j]
                indtrust[j] = indtrust[i]
                indtrust[i] = aux2
            end
        end
    end
    #    println(indtrust[n-nout+1:n])

    return P[indtrust[1:m-nout], :], sum(h[1:m-nout]), P[indtrust[(m-nout)+1:m], :]
end



"""
LOVOCGA(data, nout, θ, nome, ε=1.0e-5)

Essa função recebe de parâmetros de entrada pontos de um espaço euclidiano, uma quantidade de pontos não confiáveis, um chute inicial e aproximação desejada, se é algébrica ou geométrica. A função funciona eliminando os outliers por meio do número de pontos não confiáveis do problema e retornando a melhor aproximação possível dentro do conjunto {hiperplanos, hipercírculos, hiperesferas, hiperetas.} 

# Examples
```
julia-repl
julia> prob = load_problem("toy.csv")
julia> LOVOCGA(prob.data, prob.nout, x_0, "algebraic")

returns a vector
```
"""


function LOVOCGA(data, nout, θ, name, object, ε=1.0e-4)
    ordres = conformalsort(data, θ, nout)
    k = 1
    antres = 0.0
    eig1 = 0.0
    eig2 = 0.0
    while abs(ordres[2] - antres) > ε && k < 100
        antres = ordres[2]
        if name == "AACGA"
            θ = AACGA(ordres[1], object)
            ordres = conformalsort(data, θ, nout)
        end
        if name == "AGCGA"
            θ, eig1, eig2 = AGCGA(ordres[1], object)
            ordres = conformalsort(data, θ, nout)
        end
        if name == "ICGA"
            θ = ICGA(ordres[1], object)
            if size(θ, 2) == 2
                ordres = conformalsort2(data, θ[:, 1], θ[:, 2], nout)
            else
                ordres = conformalsort(data, θ, nout)
            end
        end
        k = k + 1
    end
    return θ, eig1, eig2, k, ordres[2], ordres[1], ordres[3]
end

function LOVOAGCGA(data, nout, θ, object, ε=1.0e-4)
    ordres = conformalsort(data, θ, nout)
    k = 1
    antres = 0.0
    eig1 = 0.0
    eig2 = 0.0
    while abs(ordres[2] - antres) > ε && k < 100
        antres = ordres[2]
        θ, eig1, eig2 = AGCGA(ordres[1], object)
        ordres = conformalsort(data, θ, nout)
        k = k + 1
    end
    return θ, eig1, eig2, ordres[2]
end

function LOVOICGA(data, nout, θ, object, ε=1.0e-4)
    ordres = conformalsort(data, θ, nout)
    k = 1
    antres = 0.0
    λ1 = 0.0
    λ2 = 0.0
    while abs(ordres[2] - antres) > ε && k < 100
        antres = ordres[2]
        θ, λ1, λ2 = ICGA(ordres[1], object)
        if size(θ, 2) == 2
            ordres = conformalsort(data, θ[:, 1], nout)
        else
            ordres = conformalsort2(data, θ[:, 1], θ[:,2], nout)
        end
    end
    return θ, λ1, λ2#, ordres[2]
end


function transphere(np)
    npnorm = np / np[end]
    centernp = npnorm[1:end-2]
    hj = push!(centernp, √(norm(centernp, 2)^2 - 2.0 * npnorm[end-1]))
    return hj
end

"""
build_problem(probtype::String, limit::Vector{Float64}, params::Vector{Float64})

Essa função recebe de parâmetros de entrada pontos de um espaço euclidiano, uma quantidade de pontos não confiáveis, um chute inicial e aproximação desejada, se é algébrica ou geométrica. A função funciona eliminando os outliers por meio do número de pontos não confiáveis do problema e retornando a melhor aproximação possível dentro do conjunto {hiperplanos, hipercírculos, hiperesferas, hiperetas.} 

# Examples
```
julia-repl
julia> prob = load_problem("toy.csv")
julia> LOVOConformal(prob.data, prob.nout, x_0, "algebraic")

returns a vector
```
"""

function build_problem(probtype::String, params::Vector{Float64}, noise::Bool)
    if probtype == "line2D"
        println("params need to be setup as [vector, npts, nout]")
        p = [params[1], params[2]]
        npts = Int(params[3])
        nout = Int(params[4])
        t = range(-15.0, stop=15.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        sgn = sign(randn())
        ruid = randn(1, npts)
        if noise == true
            for i = 1:npts
                x[i] = t[i]
                y[i] = p[1] * x[i] + p[2] + ruid[1, i] #+ (1.0 + 2 * rand()) * 7.0 * sgn
            end
        else
            for i = 1:npts
                x[i] = t[i]
                y[i] = p[1] * x[i] + p[2]
            end
        end
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        pt = 50.0
        for k = 1:nout
            y[iout[k]] = p[1] * x[iout[k]] + p[2] + rand([-pt:0.1:pt;]) #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end

        FileMatrix = ["name :" "line2d"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> t[1]*x + t[2]"; "dim :" 2; "cluster :" "false"; "noise :" noise; "solution :" [push!(p)]; "description :" "type2: line model"]

        open("line2d_$(p[1])_$(p[2])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end
    if probtype == "line3D"
        println("params need to be setup as [point,direction,npts,nout]")
        p0 = [params[1], params[2], params[3]]
        u = [params[4], params[5], params[6]]
        npts = Int(params[7])
        pp = range(-50.0, stop=50.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        ruid = randn(3, npts)
        if noise == true
            for i = 1:npts
                λ = rand(pp)
                x[i] = p0[1] + λ * u[1] + ruid[1, i]
                y[i] = p0[2] + λ * u[2] + ruid[2, i]
                z[i] = p0[3] + λ * u[3] + ruid[3, i]
            end
        else
            for i = 1:npts
                λ = rand(pp)
                x[i] = p0[1] + λ * u[1]
                y[i] = p0[2] + λ * u[2]
                z[i] = p0[3] + λ * u[3]
            end
        end
        nout = Int(params[8])
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        pt = 50.0
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-pt:0.1:pt;])
            y[iout[k]] = y[iout[k]] + rand([-pt:0.1:pt;])
            z[iout[k]] = z[iout[k]] + rand([-pt:0.1:pt;])
        end
        #FileMatrix = ["name :" "line3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :"[push!(u)]; "description :" [[p0]]]

        FileMatrix = ["name :" "line3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - t[4]^2"; "dim :" 3; "cluster :" "false"; "noise :" noise; "solution :" [push!(u, pt)]; "description :" [[p0]]]

        open("line3d_$(u[1])_$(u[2])_$(u[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "plane"
        println("params need to be setup as [point,directions,npts,nout]")
        p0 = [params[1], params[2], params[3]]
        u = [params[4], params[5], params[6]]
        v = [params[7], params[8], params[9]]
        npts = Int(params[10])
        pp = range(-50.0, stop=50.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        w = zeros(npts)
        vn = cross(u, v)
        vn = vn / norm(vn)
        d = dot(vn, p0)
        vn = push!(vn, d)
        ruid = randn(3, npts)
        sgn = sign(randn())
        if noise == true
            for i = 1:npts
                λ = rand(pp)
                μ = rand(pp)
                x[i] = p0[1] + λ * u[1] + μ * v[1] + ruid[1, i]
                y[i] = p0[2] + λ * u[2] + μ * v[2] + ruid[2, i]
                z[i] = p0[3] + λ * u[3] + μ * v[3] + ruid[3, i]
            end
        else
            for i = 1:npts
                λ = rand(pp)
                μ = rand(pp)
                x[i] = p0[1] + λ * u[1] + μ * v[1]
                y[i] = p0[2] + λ * u[2] + μ * v[2]
                z[i] = p0[3] + λ * u[3] + μ * v[3]
            end
        end
        nout = Int(params[11])
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        pt = 500
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-pt:0.1:pt;])
            y[iout[k]] = y[iout[k]] + rand([-pt:0.1:pt;])
            z[iout[k]] = z[iout[k]] + rand([-pt:0.1:pt;])
        end
        FileMatrix = ["name :" "plane"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> x[1]*t[1] + x[2]*t[2] + x[3]*t[3] + t[4]"; "dim :" 4; "cluster :" "false"; "noise :" noise; "solution :" [push!(vn)]; "description :" [[p0, p0]]]

        open("plane_$(vn[1])_$(vn[2])_$(vn[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "circle"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1], params[2], params[3]]
        r = params[4]
        u = [params[5], params[6], params[7]]
        v = [params[8], params[9], params[10]]
        npts = Int(params[11])
        u = u / norm(u)
        h = v - (dot(v, u) / norm(u)^2) * u
        v = h / norm(h)
        u = round.(u, digits=1)
        v = round.(v, digits=1)
        vn = cross(u, v) / norm(cross(u, v))
        vnc = vcat(vn, c)
        λ = [0:4/npts:1;]
        w = zeros(Int(3.0), npts)
        #h = zeros(Int(3.0), npts)
        nn = zeros(npts)
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        for i = 1:Int(round(npts / 4))
            w[:, i] = c + r * ((λ[i] * u + (1 - λ[i]) * v) / (norm(λ[i] * u + (1 - λ[i]) * v)))
        end
        for i = (Int(round(npts / 4))+1):Int(round(npts / 2))
            w[:, i] = c + r * ((λ[i-(Int(round(npts / 4)))] * (-u) + (1 - λ[i-(Int(round(npts / 4)))]) * v) / (norm(λ[i-(Int(round(npts / 4)))] * (-u) + (1 - λ[i-(Int(round(npts / 4)))]) * v)))
        end
        for i = (Int(round(npts / 2))+1):Int(round(3 * npts / 4))
            w[:, i] = c + r * ((λ[i-(Int(round(npts / 2)))] * u + (1 - λ[i-(Int(round(npts / 2)))]) * (-v)) / (norm(λ[i-(Int(round(npts / 2)))] * u + (1 - λ[i-(Int(round(npts / 2)))]) * (-v))))
        end
        for i = (Int(round(3 * npts / 4))+1):npts
            w[:, i] = c + r * ((λ[i-(Int(round(3 * npts / 4)))] * (-u) + (1 - λ[i-(Int(round(3 * npts / 4)))]) * (-v)) / (norm((λ[i-(Int(round(3 * npts / 4)))] * (-u) + (1 - λ[i-(Int(round(3 * npts / 4)))]) * (-v)))))
        end
        nout = Int(params[12])
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        pt = 50.0
        for k = 1:nout
            w[:, iout[k]] = w[:, iout[k]] + [rand([-pt:0.1:pt;]), rand([-pt:0.1:pt;]), rand([-pt:0.1:pt;])]
        end
        G = randn(3, npts)
        if noise == true
            for i = 1:npts
                x[i] = w[1, i] + G[1, i]
                y[i] = w[2, i] + G[2, i]
                z[i] = w[3, i] + G[3, i]
            end
        else
            for i = 1:npts
                x[i] = w[1, i]
                y[i] = w[2, i]
                z[i] = w[3, i]
            end
        end
        FileMatrix = ["name :" "circle3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> ( (x[1] - t[4])*t[1] +(x[2]-t[5])*t[2] +(x[3]-t[6])*t[3])^2 + ((x[1]-t[4])^2 + (x[2]-t[5])^2 + (x[3]-t[6])^2 - t[7]^2)^2"; "dim :" 7; "cluster :" "false"; "noise :" noise; "solution :" [push!(vnc, r)]; "description :" [[u, v]]]

        open("circle3D_$(c[1])_$(c[2])_$(c[3])_$(r)_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "sphere2D"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1], params[2]]
        r = params[3]
        npts = Int(params[4])
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        ruid = randn(2, npts)
        θ = range(0, stop=2π, length=npts) #Int(ceil(npts/2)))
        #θ2 = range(5*π/4, stop=7*π/4, length= 2*npts)#Int(ceil(npts/2)))
        if noise == false
            for k = 1:npts
                x[k] = c[1] + r * cos(θ[k]) #+ ruid[1, k]
                y[k] = c[2] + r * sin(θ[k]) #+ ruid[2, k]
            end
        else
            for k = 1:npts
                x[k] = c[1] + r * cos(θ[k]) + ruid[1, k]
                y[k] = c[2] + r * sin(θ[k]) + ruid[2, k]
            end
        end
        nout = Int(params[5])
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        #dx = rand() * 0.1 * r # deslocamento aleatório em x
        #dy = rand() * 0.1 * r # deslocamento aleatório em y
        pt = 50.0
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-pt:0.1:pt;])
            y[iout[k]] = y[iout[k]] + rand([-pt:0.1:pt;])  #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end
        FileMatrix = ["name :" "sphere2D"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 3; "cluster :" "false"; "noise :" noise; "solution :" [push!(c, r)]; "description :" "type3: test sphere2d with noise and outliers"]

        open("sphere2D_$(c[1])_$(c[2])_$(c[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end
    if probtype == "sphere3D2"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1], params[2], params[3]]
        r = params[4]
        npts = Int(params[5])
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)

        for k = 1:npts
            θ = 2π * rand()  # Ângulo azimutal
            φ = acos(2 * rand() - 1)  # Ângulo polar

            x[k] = c[1] + r * sin(φ) * cos(θ)
            y[k] = c[2] + r * sin(φ) * sin(θ)
            z[k] = c[3] + r * cos(φ)
        end
        if noise == true
            rd = randn(3, npts)
            for k = 1:npts
                θ = 2π * rand()  # Ângulo azimutal
                φ = acos(2 * rand() - 1)  # Ângulo polar#forma de espiral - ao criar outro forma, se obtem metade dos circulos máximos
                x[k] = c[1] + r * cos(θ) * sin(φ) + rd[1, k]
                y[k] = c[2] + r * sin(θ) * sin(φ) + rd[2, k]
                z[k] = c[3] + r * cos(φ) + rd[3, k]
            end
        end
        nout = Int(params[6])
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        pt = 50
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-pt:0.1:pt;])
            y[iout[k]] = y[iout[k]] + rand([-pt:0.1:pt;])
            z[iout[k]] = z[iout[k]] + rand([-pt:0.1:pt;])
        end
        FileMatrix = ["name :" "sphere3D"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - t[4]^2"; "dim :" 4; "cluster :" "false"; "noise :" noise; "solution :" [push!(c, r)]; "description :" [[c, c]]]

        open("sphere3D_$(c[1])_$(c[2])_$(c[3])_$(c[4])_$(nout).csv", "w") do io #o que essa linha faz exatamente?
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "sphere3D"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1], params[2], params[3]]
        r = params[4]
        npts = Int(params[5])
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        w = zeros(npts)
        θ = range(0, stop=2π, length=npts)
        φ = range(0, stop=π, length=npts)
        rd = randn(3, npts)
        if noise == true
            for k = 1:npts #forma de espiral - ao criar outro forma, se obtem metade dos circulos máximos
                x[k] = c[1] + r * cos(θ[k]) * sin(φ[k]) + rd[1, k]
                y[k] = c[2] + r * sin(θ[k]) * sin(φ[k]) + rd[2, k]
                z[k] = c[3] + r * cos(φ[k]) + rd[3, k]
            end
        else
            for k = 1:npts #forma de espiral - ao criar outro forma, se obtem metade dos circulos máximos
                x[k] = c[1] + r * cos(θ[k]) * sin(φ[k])
                y[k] = c[2] + r * sin(θ[k]) * sin(φ[k])
                z[k] = c[3] + r * cos(φ[k])
            end
        end
        nout = Int(params[6])
        k = 1
        iout = []
        while k <= nout
            i = rand([1:npts;])
            if i ∉ iout
                push!(iout, i)
                k = k + 1
            end
        end
        pt = 50
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-pt:0.1:pt;])
            y[iout[k]] = y[iout[k]] + rand([-pt:0.1:pt;])
            z[iout[k]] = z[iout[k]] + rand([-pt:0.1:pt;])
        end
        FileMatrix = ["name :" "sphere3D"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - t[4]^2"; "dim :" 4; "cluster :" "false"; "noise :" noise; "solution :" [push!(c, r)]; "description :" [[c, c]]]

        open("sphere3D_$(c[1])_$(c[2])_$(c[3])_$(c[4])_$(nout).csv", "w") do io #o que essa linha faz exatamente?
            writedlm(io, FileMatrix)
        end
    end
end


function LOVOCLASSIFICATION(data, nout, θ, ε1, ε2, ε=1.0e-4)
    A = LOVOAGCGA(data, nout, θ, "sphere")
    B, λ1, λ2 = LOVOICGA(prob.data, prob.nout, ones(5), "reta")
    sphere = transphere(A[1])
    display(λ1)
    display(λ2)
    if 1 / sphere[end] < ε2
        if λ2 - λ1 < ε1
            return B[:, 1], B[:, 2], "line"
        else
            return B[:, 1], 0.0, "plane"
        end
    end
    if abs(A[2] - λ1) < ε1
        return A[1], B[1], "circle"
    else
        return A[1], 0.0, "sphere"
    end
end


function fittingclass(data, ε1, ε2)
    (N, n) = size(data)
    v = [0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    D = [data'; v'; ones(1, N)]
    J = copy(D')
    H = -copy(J[:, n+1])
    J[:, n+1] = -J[:, n+2]
    J[:, n+2] = H
    DDt = D * D'
    p = (1.0 / N)
    aux = -copy(DDt[:, n+1])
    DDt[:, n+1] = -DDt[:, n+2]
    DDt[:, n+2] = aux
    P = p .* (DDt)
    IM = copy(P)
    IM[[end - 1, end], :] = -P[[end, end - 1], :]
    F = eigen(P)
    λ1 = F.values[2]
    λ2 = F.values[3]
    #display(λ1)
    #display(λ2)
    v1 = real(F.vectors[:, 2])
    v2 = real(F.vectors[:, 3])
    obj = ""
    s1 = transphere(v1)

    if λ2 - λ1 > ε1
        if 1 / s1[end] > ε2
            obj = "sphere"
            return v1, obj
        else
            B = IM[1:end-2, 1:end-2]
            u = IM[1:end-2, end-1]
            w = IM[1:end-2, end]
            a = IM[end-1, end-1]
            Hm = B - (u * u') / a
            Fv = eigen(Hm)
            vn = Fv.vectors[:, 1]
            d = -(u' * vn) / a
            π = [vn; d; 0]
            #println("hyperplane")
            obj = "plane"
            return π, obj
        end
    else
        s2 = transphere(v2)
        #display(s1)
        if 1 / s1[end] < ε2
            B = IM[1:end-2, 1:end-2]
            u = IM[1:end-2, end-1]
            w = IM[1:end-2, end]
            a = IM[end-1, end-1]
            Hm = B - (u * u') / a
            Fv = eigen(Hm)
            vn = Fv.vectors[:, 1]
            vn2 = Fv.vectors[:, 2]
            d = -(u' * vn) / a
            d2 = -(u' * vn2) / a
            π = [vn; d; 0]
            π2 = [vn2; d2; 0]
            #println("line")
            obj = "line"
            return hcat(π, π2), obj
        else
            B = IM[1:end-2, 1:end-2]
            u = IM[1:end-2, end-1]
            w = IM[1:end-2, end]
            a = IM[end-1, end-1]
            Hm = B - (u * u') / a
            Fv = eigen(Hm)
            vn = Fv.vectors[:, 1]
            d = -(u' * vn) / a
            π = [vn; d; 0]
            #println("hypercircle")
            obj = "circle"
            return hcat(v1, π), obj
        end
    end
end

function show(io::IO, fout::FitOutputType)

    print(io, "  ▶ Output ◀ \n")
    if Bool(fout.status) == true
        print(io, "  ↳ Status (.status) = Convergent \n")
    else
        print(io, "  ↳ Status (.status) = Divergent \n")
    end
    print(io, "  ↳ Solution (.solution) = $(fout.solution) \n")
    print(io, "  ↳ Number of iterations (.niter) = $(fout.niter) \n")
    print(io, "  ↳ Minimum (.minimum) = $(fout.minimum) \n")
    print(io, "  ↳ Number of function calls (.feval) = $(fout.feval) \n")
end

#comentando 

end # module

