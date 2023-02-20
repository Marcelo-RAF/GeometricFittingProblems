module GeometricFittingProblems

using DelimitedFiles, LinearAlgebra, Plots
pyplot()

export load_problem, solve, build_problem, inverse_power_method, solve2, visualize

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
    description::Vector{Vector{Any}}
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
    return FitProbType(prob_matrix[1, 2], eval(Meta.parse(prob_matrix[2, 2])), prob_matrix[3, 2], prob_matrix[4, 2], eval(Meta.parse(prob_matrix[5, 2])), prob_matrix[6, 2], prob_matrix[7, 2], prob_matrix[8, 2], eval(Meta.parse(prob_matrix[9, 2])), eval(Meta.parse(prob_matrix[10, 2])))
end


function CGAHypersphere(data; ε=1.0e-4) #algoritmo dorst esferas
    (N, n) = size(data)
    D = [data'; ones(1, N)]
    v = [0.5 * norm(D[1:n, i], 2)^2 for i = 1:N]
    D = [D; v']
    #println("First D")
    #display(D)
    J = copy(D')
    H = -copy(J[:, n+1])
    J[:, n+1] = -J[:, n+2]
    J[:, n+2] = H
    #display(J)

    DDt = D * D'
    #println("Second D")
    #M = zeros(n+2,n+2)
    #for i=1:n
    # M[i,i] = 1.0
    #end
    #M[n+1,n+2] = -1.0
    #M[n+2,n+1] = -1.0
    aux = -copy(DDt[:, n+1])
    DDt[:, n+1] = -DDt[:, n+2]
    DDt[:, n+2] = aux
    #display(DDt)
    p = (1.0 / N)
    P = p .* (DDt)
    F = eigen(P)
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
    xnorm = (1.0 / (F.vectors[:, indmin][end-1])) * F.vectors[:, indmin]
    center = xnorm[1:end-2]

    #u = 0.5*norm(xnorm[1:end-2])^2
    #y =[xnorm[1] xnorm[2] xnorm[3] u]
    #display(J*y')

    return push!(center, √(norm(center, 2)^2 - 2.0 * xnorm[end]))
end

function sort_sphere_res(P, x, nout)
    n = length(P[:, 1])
    m = length(P[1, :])
    v = zeros(n)
    for i = 1:n
        for j = 1:m
            v[i] = v[i] + (P[i, j] - x[j])^2
        end
        v[i] = abs(v[i] - x[end]^2)
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

function LOVOCGAHypersphere(data, nout, θ, ε=1.0e-8)

    ordres = sort_sphere_res(data, θ, nout)
    k = 1
    antres = 0.0
    while abs(ordres[2] - antres) > ε
        display(ordres[2])
        antres = ordres[2]
        θ = CGAHypersphere(ordres[1])
        println(θ)
        ordres = sort_sphere_res(data, θ, nout)
        k = k + 1
    end
    #   display(k)
    return θ

end




"""
    solve :: Function

This functions is able to solve a fitting problem previous loaded.

# Example
```
julia-repl
julia> prob =  load_problem("sphere2D_50.0_50.0_8.0_10.csv")

julia> solve(prob,"CGA-Hypersphere")

```
return the smallest positive eigen and eigen vector associated which give to us the desired parameters.

# Example
```
julia-repl
julia> prob =  load_problem("sphere2D_50.0_50.0_8.0_10.csv")

julia> solve(prob,"LOVO-CGA-Hypersphere")

```
return the desired parameters.

another way to run is using a particular initial guess.

# Example
```
julia-repl
julia> prob =  load_problem("sphere2D_50.0_50.0_8.0_10.csv")

julia> solve(prob,"LOVO-CGA-Hypersphere",rand(3))

```

"""
function solve(prob::FitProbType, method::String, initθ=CGAHypersphere(prob.data))
    if method == "CGA-Hypersphere"
        return CGAHypersphere(prob.data)
    end
    if method == "LOVO-CGA-Hypersphere"
        return LOVOCGAHypersphere(prob.data, prob.nout, initθ)
    end
end

function solve2(prob::FitProbType, method::String, initθ=CGAHypercircle(prob.data))
    if method == "CGA-Hypercircle"
        return CGAHypercircle(prob.data)
    end
    if method == "LOVO-CGA-Hypercircle"
        return LOVOCGAHypercircle(prob.data, prob.nout, initθ)
    end
end

function build_problem(probtype::String, limit::Vector{Float64}, params::Vector{Float64})
    if probtype == "circle3d"
        println("params need to be setup as [center,radious,npts,nout]")
        c = [params[1], params[2], params[3]]
        r = params[4]
        u = [params[5], params[6], params[7]]
        v = [params[8], params[9], params[10]]
        npts = Int(params[11])
        u = u / norm(u)
        h = v - (dot(v, u) / norm(u)^2) * u
        v = h / norm(h)
        display(u)
        display(v)
        λ = [0:4/npts:1;]
        println(λ)
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
        for k = 1:nout
            w[:, iout[k]] = w[:, iout[k]] + [rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;]), rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;]), rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])]
        end
        G = randn(3, npts)
        for i = 1:npts
            x[i] = w[1, i] + G[1, i]
            y[i] = w[2, i] + G[2, i]
            z[i] = w[3, i] + G[3, i]
        end
        FileMatrix = ["name :" "circle3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c, r)]; "description :" [[u, v]]]

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
        #xr = randn(npts)
        #yr = randn(npts)
        ruid = randn(2, npts)
        θ = [0.0:2*π/(npts-1):2*π;]
        for k = 1:npts
            x[k] = c[1] + r * cos(θ[k]) + ruid[1, k]#xr[k]
            y[k] = c[2] + r * sin(θ[k]) + ruid[2, k] #yr[k]
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
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])   #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end
        FileMatrix = ["name :" "sphere2D"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c, r)]; "description :" [[c, c]]]

        open("sphere2D_$(c[1])_$(c[2])_$(c[3])_$(nout).csv", "w") do io
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
        θ = [0.0:2*π/(npts-1):2*π;]
        φ = [0.0:π/(npts-1):π;]
        rd = randn(3, npts)
        for k = 1:npts #forma de espiral - ao criar outro forma, se obtem metade dos circulos máximos
            x[k] = c[1] + r * cos(θ[k]) * sin(φ[k]) + rd[1, k]
            y[k] = c[2] + r * sin(θ[k]) * sin(φ[k]) + rd[2, k]
            z[k] = c[3] + r * cos(φ[k]) + rd[3, k]
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
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            z[iout[k]] = z[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
        end
        FileMatrix = ["name :" "sphere3D"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - t[4]^2"; "dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c, r)]; "description :" [[c, c]]]

        open("sphere3D_$(c[1])_$(c[2])_$(c[3])_$(c[4])_$(nout).csv", "w") do io #o que essa linha faz exatamente?
            writedlm(io, FileMatrix)
        end
    end

end

"""
    inverse_power_method :: function

This functions implements the inverse power method to find the smallest eigen value associated to an array A.

# Examples
```
julia-repl

julia> A = [1.0 2.0 0.0; 2.0 -5.0 3.0; 0.0 3.0 4.0]

julia> inverse_power_method(A,[1.0,1.0,1.0])

returns ???
```
"""
function CGAHypercircle(data; ε=1.0e-4)
    (N, n) = size(data)
    p = (1.0 / N)
    H1 = zeros(3, 3)
    H2 = zeros(3, 3)
    H3 = 0.0
    H4 = zeros(3, 3)
    H5 = zeros(3)'
    H6 = 0.0
    H7 = zeros(3)'
    Id = [1 0 0; 0 1 0; 0 0 1]
    SI = zeros(3, 3)
    P = zeros(10, 10)
    Nf = convert(AbstractFloat, N)

    for i = 1:N
        H1 = H1 + [0.0 -data[i, 3] data[i, 2]; data[i, 3] 0.0 -data[i, 1]; -data[i, 2] data[i, 1] 0.0]
        H2 = H2 + ([0.0 -data[i, 3] data[i, 2]; data[i, 3] 0.0 -data[i, 1]; -data[i, 2] data[i, 1] 0.0])^2
        H3 = H3 + norm(data[i, :])^2
        H4 = H4 + (norm(data[i, :])^2) * [0.0 -data[i, 3] data[i, 2]; data[i, 3] 0.0 -data[i, 1]; -data[i, 2] data[i, 1] 0.0]
        H5 = H5 + (norm(data[i, :])^2) * (data[i, :])'
        H6 = H6 + norm(data[i, :])^4
        H7 = H7 + data[i, :]'
        SI = SI + [1 0 0; 0 1 0; 0 0 1]
    end
    P[1:3, 1:3] = -H2
    P[4:6, 1:3] = -H1
    P[7:9, 1:3] = (1 / 2) * H4 #antes era só h1
    P[1:3, 4:6] = -(1 / 2) * H4
    P[4:6, 4:6] = H2 + (1 / 2) * H3 * Id
    P[7:9, 4:6] = (1 / 4) * H6 * Id
    P[10, 4:6] = (1 / 2) * H5
    P[1:3, 7:9] = H1
    P[4:6, 7:9] = SI
    P[7:9, 7:9] = H2 + (1 / 2) * H3 * Id
    P[10, 7:9] = H7
    P[4:6, 10] = -H7'
    P[7:9, 10] = -(1 / 2) * H5'
    P[10, 10] = -H3
    P = p .* (P)
    #display(P)
    F = eigen(P)
    indmin = 1
    valmin = F.values[1]
    for i = 2:10
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
    A = F.vectors[:, indmin]
    n1 = -A[4:6]
    d1, d2, d3 = n1[1], n1[2], n1[3]
    C = [d1 d2 d3; 0 d3 -d2; -d3 0 d1; d2 -d1 0]
    α = norm(A[4:6])
    center = [A[10]; A[1:3]]' / -C'    #talvez tenha que alterar muito
    n1 = n1 / α
    radius = (center * center' - 2 * n1' * A[7:9] / α - 2 * (n1' * center')^2)


    #α = norm(n1)
    #daqui pra baixo ta diferente
    #H = A/α
    #n = -H[4:6]
    #B0, B1, B2, B3 = -H[10], H[1], H[2], H[3]
    #c = [B0 -B3 B2; B3 B0 -B1; -B2 B1 B0]*n
    #vinf = H[7:9]
    #r = sqrt(norm(c)^2 -2*(n1'*vinf)/α - 2*((B0)^2))


    #n = n1/α
    #B0, B1, B2, B3 = -A[10], A[1], A[2], A[3]
    #c = [B0 -B3 B2; B3 B0 -B1; -B2 B1 B0] * (n1/α^2)
    #vinf = A[7:9]
    #r = sqrt(norm(c)^2 -2*(n1'*vinf/α^2) - 2*((B0)^2)/α^2)

    # push!(center,√(norm(center,2)^2 -2.0*xnorm[end]))
    #return push!(n1, center, sqrt(radius))
    u = [n1[1], n1[2], n1[3], center[1], center[2], center[3], sqrt(radius)]
    return u
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

function LOVOCGAHypercircle(data, nout, θ, ε=1.0e-12)
    ordres = sort_circle_res(data, θ, nout)
    k = 1
    antres = 0.0
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        θ = CGAHypercircle(ordres[1])
        ordres = sort_circle_res(data, θ, nout)
        k = k + 1
    end
    display(k)
    return θ
end


function inverse_power_method(A::Array{Float64}; q0=ones(size(A)[1]), ε=10.0^(-4), limit=100)
    stop_criteria = 1000.0
    F = lu(A)
    B = inv(A)
    k = 1
    s = 0.0
    q = zeros(length(q0))
    while stop_criteria > ε && k < limit
        s = norm(q0, Inf)
        q = B * (q0 / s)
        stop_criteria = norm(abs.(q) - abs.(q0), Inf)
        q0 = copy(q)
        k = k + 1
    end
    if k == limit
        error("iteration limit of inverse power method was reached")
    else
        return q, 1.0 / s
    end
end

function teste(prob)
    plt = plot()
    if prob.name == "sphere3D" || prob.name == "\tsphere3D"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], prob.data[:, 3], line=:scatter, aspect_ratio=:equal, lab="pontos do problema")
        n = 100
        h1 = zeros(n)
        h2 = zeros(n)
        h3 = zeros(n)
        h4 = zeros(n)
        for i = 1:n
            h1[i] = prob.solution[1]
            h2[i] = prob.solution[2]
            h3[i] = prob.solution[3]
            h4[i] = prob.solution[4]
        end
        u = range(-π, π; length=n)
        v = range(0, π; length=n)
        x = h1 .+ prob.solution[4] * cos.(u) * sin.(v)'
        y = h2 .+ prob.solution[4] * sin.(u) * sin.(v)'
        z = h3 .+ prob.solution[4] * cos.(v)'
        xs = h1 .+ cos.(u) * sin.(v)'
        ys = h1 .+ sin.(u) * sin.(v)'
        zs = h1 .+ ones(n) * cos.(v)'
        plot!(plt, x, y, z, xs, ys, zs, st=:wireframe, camera=(-50, 50))
        display(plt)
    end
end



function visualize(prob, a)
    plt = plot()
    if prob.name == "sphere2D" || prob.name == "\tsphere2D"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], line=:scatter, aspect_ratio=:equal, lab="pontos do problema")
        θ = [0.0:2*π/360:2*π;]
        xs = a[1] .+ a[3] * cos.(θ)
        ys = a[2] .+ a[3] * sin.(θ)
        x = prob.solution[1] .+ prob.solution[3] * cos.(θ)
        y = prob.solution[2] .+ prob.solution[3] * sin.(θ)
        plot!(plt, xs, ys, color=:red, lab="solução do algoritmo")
        plot!(plt, x, y, color=:green, lab="solução perfeita")
        display(plt)
    end
    if prob.name == "sphere3D" || prob.name == "\tsphere3D"
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
            h1[i] = prob.solution[1]
            h2[i] = prob.solution[2]
            h3[i] = prob.solution[3]
            h4[i] = c[1]
            h5[i] = c[2]
            h6[i] = c[3]
        end
        x = h1 .+ prob.solution[4] * cos.(u) * sin.(v)'
        y = h2 .+ prob.solution[4] * sin.(u) * sin.(v)'
        z = h3 .+ prob.solution[4] * cos.(v)'
        xs = h4 .+ c[4] * cos.(u) * sin.(v)'
        ys = h5 .+ c[4] * sin.(u) * sin.(v)'
        zs = h6 .+ c[4] * cos.(v)'
        wireframe!(xs, ys, zs, aspect_ratio=:equal, lab="solução do algoritmo")
        wireframe!(x, y, z, aspect_ratio=:equal, color=:red, lab="solução perfeita")
        display(plt)
    end
    if prob.name == "circle3d" || prob.name == "\tcircle3d"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], prob.data[:, 3], line=:scatter, aspect_ratio=:equal, lab="pontos do problema")
        vn = [a[1], a[2], a[3]]
        u = [-a[2], a[1], 0.0]
        v = [0.0, a[3], -a[2]]
        u = u / norm(u)
        h = v - (dot(v, u) / norm(u)^2) * u
        v = h / norm(h)
        uprob = prob.description[1]
        vprob = prob.description[2]
        sol = prob.solution
        θ = [0.0:2*π/360:2*π;]
        xprob = sol[1] .+ sol[4] * (cos.(θ)) * uprob[1] .+ sol[4] * (sin.(θ)) * vprob[1]
        yprob = sol[2] .+ sol[4] * (cos.(θ)) * uprob[2] .+ sol[4] * (sin.(θ)) * vprob[2]
        zprob = sol[3] .+ sol[4] * (cos.(θ)) * uprob[3] .+ sol[4] * (sin.(θ)) * vprob[3]
        x = a[4] .+ a[7] * (cos.(θ)) * u[1] .+ a[7] * (sin.(θ)) * v[1]
        y = a[5] .+ a[7] * (cos.(θ)) * u[2] .+ a[7] * (sin.(θ)) * v[2]
        z = a[6] .+ a[7] * (cos.(θ)) * u[3] .+ a[7] * (sin.(θ)) * v[3]
        plot!(plt, xprob, yprob, zprob, camera=(20, 50), color=:green, lab="solução perfeita")
        plot!(plt, x, y, z, camera=(20, 50), color=:red, lab="solução do algoritmo") #usar camera=(100,40) pro nout=20
        display(plt)
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


end # module

