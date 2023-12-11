module GeometricFittingProblems

using DelimitedFiles, LinearAlgebra, Plots

export load_problem, solve, build_problem, solve2, visualize, plot_plane, Levenberg, LovoLM, LMPersistent

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


function CGAHypersphere(data, method::String, ε=1.0e-5) #algoritmo dorst esferas
    (N, n) = size(data)
    D = [data'; ones(1, N)]
    v = [0.5 * norm(D[1:n, i], 2)^2 for i = 1:N]
    D = [D; v']
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
    if method == "eigenvector"
        F = eigen(P)
        indmin = 1
        #println(F.values)
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
        #xnorm = (1.0 / (F.vectors[:, indmin][end-1])) * F.vectors[:, indmin]
        #center = xnorm[1:end-2]
        #display(F)
        return F.vectors[:, indmin] #push!(center, √(norm(center, 2)^2 - 2.0 * xnorm[end]))
    elseif method == "nullspace"

        P[end, :] = zeros(n + 2)
        np = nullspace(P)

        #npnorm = np / np[end-1]
        #centernp = npnorm[1:end-2]
        #push!(centernp, √(norm(centernp, 2)^2 - 2.0 * npnorm[end]))
        return np
    end
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


function hildebran(data, method::String, ε=1.0e-5)
    (N, n) = size(data)
    v = [-0.5 * norm(data[i, :], 2)^2 for i = 1:N]
    D = [data'; v'; -ones(1, N)]
    Dd = simetrica(D)
    Dd = 1.0 / N * Dd #apagar
    if method == "eigenvector"
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
        xnorm1 = (1.0 / (F.vectors[:, indmin][end-1])) * F.vectors[:, indmin]
        center1 = xnorm1[1:end-2]

        #xnorm = (1.0 / (F.vectors[:, indmin][end-1])) * F.vectors[:, indmin]
        #center = xnorm[1:end-2]

        #display(F)
        return F.vectors[:, indmin] #push!(center, √(norm(center, 2)^2 - 2.0 * xnorm[end]))
    elseif method == "nullspace"
        Dd[end-1, :] = zeros(n + 2)
        #Dd[2, :] = zeros(n + 2)
        np = nullspace(Dd)

        #npnew = [np[1], np[2], np[3], np[5], np[4]]
        #npnorm = np / np[3]
        #centernp = npnorm[1:end-2]
        #display(Dd)
        return np #push!(centernp, √(norm(centernp, 2)^2 - 2.0 * npnorm[end]))
    end
end


function conformalsort(P, x, nout)
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
    return P[indtrust[1:m-nout], :], sum(h[1:m-nout])
end



function LOVOConformal(data, nout, θ, nome, ε=1.0e-5)
    ordres = conformalsort(data, θ, nout)
    k = 1
    antres = 0.0
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        if nome == "algebraic"
            θ = hildebran(ordres[1], "nullspace")
        end
        if nome == "geometric"
            θ = CGAHypersphere(ordres[1], "nullspace")
        end
        ordres = conformalsort(data, θ, nout)
        k = k + 1
    end
    return θ, k
end


function solve(prob::FitProbType, method::String)
    if method == "CGA-Hypersphere"
        return CGAHypersphere(prob.data)
    end
    if method == "LOVO-CGA-Geometric"
        initθ = CGAHypersphere(prob.data, "nullspace")
        return LOVOConformal(prob.data, prob.nout, initθ, "geometric")
    end
    if method == "LOVO-CGA-Algebraic"
        initθ = hildebran(prob.data, "nullspace")
        return LOVOConformal(prob.data, prob.nout, initθ, "algebraic")
    end
end



function build_problem(probtype::String, limit::Vector{Float64}, params::Vector{Float64})
    if probtype == "logistic"
        println("params need to be setup as [vector, npts, nout]")
        p = [params[1], params[2], params[3], params[4]]
        npts = Int(params[5])
        nout = Int(params[6])
        t = range(-50.0, stop=50.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        sgn = sign(randn())
        for i = 1:npts
            x[i] = t[i]
            y[i] = p[1] + p[2] / (1 + exp(-p[3] * x[i] + p[4]))
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
        for k = 1:nout
            y[iout[k]] = p[1] + p[2] / (1 + exp(-p[3] * x[iout[k]] + p[4])) + randn() * 200 #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end
        FileMatrix = ["name :" "logistic"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> x[1] + x[2]/(1 + exp(-x[3]*t + x[4]))"; "dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [push!(p)]; "description :" "type2: cubic model"]

        open("logistic_$(p[1])_$(p[2])_$(p[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end



    end
    if probtype == "exponencial"
        println("params need to be setup as [vector, npts, nout]")
        p = [params[1], params[2]]
        npts = Int(params[3])
        nout = Int(params[4])
        t = range(0.0, stop=15.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        ruid = randn(npts)
        for i = 1:npts
            x[i] = t[i]
            y[i] = p[1] * exp(-p[2] * x[i]) + ruid[i]
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
        r = 3
        for k = 1:nout
            y[iout[k]] = y[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])#rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
        end

        FileMatrix = ["name :" "exponencial"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> x[1]*exp(-x[2]*t) "; "dim :" 2; "cluster :" "false"; "noise :" "false"; "solution :" [push!(p)]; "description :" "type: exponencial function"]

        open("exponencial_$(p[1])_$(p[2])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "cubic"
        println("params need to be setup as [vector, npts, nout]")
        p = [params[1], params[2], params[3], params[4]]
        npts = Int(params[5])
        nout = Int(params[6])
        t = range(-5.0, stop=5.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        sgn = sign(randn())
        for i = 1:npts
            x[i] = t[i]
            y[i] = p[1] * x[i]^3 + p[2] * x[i]^2 + p[3] * x[i] + p[4] #+ (1.0+2*rand()) * 7.0*sgn
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
        for k = 1:nout
            y[iout[k]] = p[1] * x[iout[k]]^3 + p[2] * x[iout[k]]^2 + p[3] * x[iout[k]] + p[4] + randn() * 200 #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end

        FileMatrix = ["name :" "cubic"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> x[1]*t^3 + x[2]*t^2 + x[3]*t + x[4]"; "dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [push!(p)]; "description :" "type2: cubic model"]

        open("cubic_$(p[1])_$(p[2])_$(p[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
    if probtype == "line2d"
        println("params need to be setup as [vector, npts, nout]")
        p = [params[1], params[2]]
        npts = Int(params[3])
        nout = Int(params[4])
        t = range(-15.0, stop=15.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        sgn = sign(randn())
        for i = 1:npts
            x[i] = t[i]
            y[i] = p[1] * x[i] + p[2] + (1.0 + 2 * rand()) * 7.0 * sgn
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
        for k = 1:nout
            y[iout[k]] = p[1] * x[iout[k]] + p[2] + randn() * 200 #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end

        FileMatrix = ["name :" "line2d"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> x[1]*t + x[2]"; "dim :" 2; "cluster :" "false"; "noise :" "true"; "solution :" [push!(p)]; "description :" "type2: line model"]

        open("line2d_$(p[1])_$(p[2])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end

    end
    if probtype == "line3d"
        println("params need to be setup as [point,direction,npts,nout]")
        p0 = [params[1], params[2], params[3]]
        u = [params[4], params[5], params[6]]
        npts = Int(params[7])
        pp = range(-50.0, stop=50.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        for i = 1:npts
            for j = 1:npts
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
        r = 3.0
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            z[iout[k]] = z[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
        end
        FileMatrix = ["name :" "line3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> p0 + λ*u"; "dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" "description :" [[p0, p0]]]

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
        # λ = range(0, stop = 66, length=npts)
        # μ = range(-50, stop = 5, length=npts)
        pp = range(-10.0, stop=10.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        vn = cross(u, v)
        vn = vn / norm(vn)
        d = dot(vn, p0)
        vn = push!(vn, d)
        ruid = randn(3, npts)
        for i = 1:npts
            for j = 1:npts
                λ = rand(pp)
                μ = rand(pp)
                x[i] = p0[1] + λ * u[1] + μ * v[1] #+ ruid[1,i]
                y[i] = p0[2] + λ * u[2] + μ * v[2] #+ ruid[2,i]
                z[i] = p0[3] + λ * u[3] + μ * v[3] #+ ruid[3,i]
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
        r = 5.0
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            z[iout[k]] = z[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
        end
        FileMatrix = ["name :" "plane"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> p0 + λ*u + μ*v"; "dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [push!(vn)]; "description :" [[p0, p0]]]

        open("plane_$(vn[1])_$(vn[2])_$(vn[3])_$(nout).csv", "w") do io
            writedlm(io, FileMatrix)
        end
    end
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
        for k = 1:nout
            w[:, iout[k]] = w[:, iout[k]] + [rand([-0.5*r:0.1:0.5*r;]), rand([-0.5*r:0.1:0.5*r;]), rand([-0.5*r:0.1:0.5*r;])]
        end
        #G = randn(3, npts)
        for i = 1:npts
            x[i] = w[1, i] #+ G[1, i]
            y[i] = w[2, i] #+ G[2, i]
            z[i] = w[3, i] #+ G[3, i]
        end
        FileMatrix = ["name :" "circle3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 7; "cluster :" "false"; "noise :" "false"; "solution :" [push!(vnc, r)]; "description :" [[u, v]]]

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
        ruid = randn(2, npts)
        θ = range(0, stop=2π, length=npts) #Int(ceil(npts/2)))
        #θ2 = range(5*π/4, stop=7*π/4, length= 2*npts)#Int(ceil(npts/2)))
        for k = 1:npts
            x[k] = c[1] + r * cos(θ[k]) + ruid[1, k]
            y[k] = c[2] + r * sin(θ[k]) + ruid[2, k]
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
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-(0.5)*r:0.1:(0.5)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(0.5)*r:0.1:(0.5)*r;])   #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end
        FileMatrix = ["name :" "sphere2D"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - x[3]^2"; "dim :" 3; "cluster :" "false"; "noise :" "true"; "solution :" [push!(c, r)]; "description :" "type3: test sphere2d with noise and outliers"]

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
        θ = range(0, stop=2π, length=npts)
        φ = range(0, stop=π, length=npts)
        #φ2 = range(5π/6, stop=π, length=npts)
        rd = randn(3, npts)
        #if iseven(npts)==false
        #    l = Int(round(npts/2))
        #   h = Int(ceil(npts/2))
        #else
        #   l = Int(npts/2)
        #  h = Int(npts/2) + 1
        # end
        for k = 1:npts #forma de espiral - ao criar outro forma, se obtem metade dos circulos máximos
            x[k] = c[1] + r * cos(θ[k]) * sin(φ[k]) #+ rd[1, k]
            y[k] = c[2] + r * sin(θ[k]) * sin(φ[k]) #+ rd[2, k]
            z[k] = c[3] + r * cos(φ[k]) #+ rd[3, k]
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
            dx = rand() * 0.2 * r
            dy = rand() * 0.2 * r
            dz = rand() * 0.2 * r
            x[iout[k]] = x[iout[k]] + dx #rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;]) 
            y[iout[k]] = y[iout[k]] + dy #rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
            z[iout[k]] = z[iout[k]] + dz #rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
        end
        FileMatrix = ["name :" "sphere3D"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - x[4]^2"; "dim :" 4; "cluster :" "false"; "noise :" "true"; "solution :" [push!(c, r)]; "description :" [[c, c]]]

        open("sphere3D_$(c[1])_$(c[2])_$(c[3])_$(c[4])_$(nout).csv", "w") do io #o que essa linha faz exatamente?
            writedlm(io, FileMatrix)
        end
    end
end

function Levenberg(Function, Jacobian, x, data, ε=10e-4, λ_min=1e-2)
    k = 0
    F = Function(x, data)
    J = Jacobian(x, data)
    (m, n) = size(J)
    xn = zeros(length(x))
    Id = Matrix{Float64}(I, n, n)
    λ = 1.0#norm((J') * F, 2) / (norm(F, 2)^2)
    k1 = 2.5
    k2 = 3.5
    while norm((J') * F) > ε && k < 50
        d = (J' * J + λ * Id) \ ((-J') * F)
        xn = x + d
        if 0.5 * norm(Function(xn, data), 2)^2 < 0.5 * norm(Function(x, data), 2)^2
            x = xn
            if λ < λ_min
                λ = λ_min
            else
                λ = λ / k1
            end
            F = Function(x, data)
            J = Jacobian(x, data)
        else
            λ = λ * k2
        end
        k = k + 1
    end
    x[1:3] = x[1:3] / norm(x[1:3])
    return x, k
end

function LovoLM(Function, Jacobian, Ord, xk, data, nout, ε=1.0e-4, MAXIT=50)
    newdata = Ord(data, xk, nout)
    R = Function(xk, newdata[1])
    J = Jacobian(xk, newdata[1])
    (m, n) = size(J)
    Id = Matrix{Float64}(I, n, n)
    k = 0
    λ_up = 2.0
    λ_down = 2.0
    λ = 1.0
    μ = 0.7
    dk = 0.0
    while norm(J' * R, 2) >= ε && k < MAXIT
        dk = (J' * J + λ * Id) \ ((-J') * R)
        md = 0.5 * (norm((R + J * dk), 2))^2 + λ * norm(dk, 2)^2
        Rd = Function(xk + dk, newdata[1])
        ρk = (0.5 * norm(R, 2)^2 - 0.5 * norm(Rd, 2)^2) / (0.5 * norm(R, 2)^2 - md)
        if ρk < μ
            λ = λ * λ_up
        else
            λ = λ / λ_down
            xk = xk + dk
            newdata = Ord(data, xk, nout)
            R = Function(xk, newdata[1])
            J = Jacobian(xk, newdata[1])
            k = k + 1
        end
    end
    #xk = xk/norm(xk[1:3])
    return xk, k, newdata[1]
end

#Levenberg(Function, Jacobian, x, data, ε=10e-5, λ_min=1e-4)

function LMPersistent(Function, Jacobian, Ord, xk, data, nout, ε=1.0e-4)
    ordres = Ord(data, xk[1], nout)
    antres = 0.0
    k = 1
    kk = 0
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        xk = Levenberg(Function, Jacobian, xk[1], ordres[1])
        kk = kk + xk[2]
        ordres = Ord(data, xk[1], nout)
        k = k + 1
        #display(abs(ordres[2] - antres))
    end
    #x = xk[1]
    #x[1:3] = x[1:3]/norm(x[1:3])
    return xk[1], kk, k, ordres[1]
end


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

    display(radius)

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

function circleag(s1, s2)
    s = externo(s1, s2)
    A = [s[5], -s[2], s[1], -s[3], -s[6], -s[8], s[4], s[7], s[9], s[10]]
    n1 = -A[4:6]
    d1, d2, d3 = n1[1], n1[2], n1[3]
    C = [d1 d2 d3; 0 d3 -d2; -d3 0 d1; d2 -d1 0]
    α = norm(A[4:6])
    center = [A[10]; A[1:3]]' / -C'    #talvez tenha que alterar muito
    n1 = n1 / α
    radius = abs((center * center' - 2 * n1' * A[7:9] / α - 2 * (n1' * center')^2))
    u = [n1[1], n1[2], n1[3], center[1], center[2], center[3], sqrt(radius)]
    return u
end


function LOVOCGAHypercircle(data, nout, θ, ε=1.0e-6)
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


function externo(u, v)
    m = length(u)
    p = Int((m * (m - 1)) / 2)
    w = zeros(p)
    k = 1
    for i = 1:m
        for j = (i+1):m
            w[k] = u[i] * v[j] - u[j]v[i]
            k = k + 1
        end
    end
    return w
end

function visualize(prob, a)
    pyplot()#beckendpyplot
    plt = plot()
    if prob.name == "sphere2D" || prob.name == "\tsphere2D"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], line=:scatter, aspect_ratio=:equal, lab="pontos do problema")
        θ = [0.0:2*π/360:2*π;]
        xs = a[1] .+ a[3] * cos.(θ)
        ys = a[2] .+ a[3] * sin.(θ)
        x = prob.solution[1] .+ prob.solution[3] * cos.(θ)
        y = prob.solution[2] .+ prob.solution[3] * sin.(θ)
        plot!(plt, xs, ys, color=:red, lab="solução do algoritmo")
        #plot!(plt, x, y, color=:green, lab="solução perfeita")
        display(plt)
    end
    if prob.name == "sphere3D" || prob.name == "\tsphere3D"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], prob.data[:, 3], line=:scatter, aspect_ratio=:equal)#, lab="pontos do problema")
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
            h4[i] = a[1]
            h5[i] = a[2]
            h6[i] = a[3]
        end
        #x = h1 .+ prob.solution[4] * cos.(u) * sin.(v)'
        #y = h2 .+ prob.solution[4] * sin.(u) * sin.(v)'
        #z = h3 .+ prob.solution[4] * cos.(v)'
        xs = h4 .+ a[4] * cos.(u) * sin.(v)'
        ys = h5 .+ a[4] * sin.(u) * sin.(v)'
        zs = h6 .+ a[4] * cos.(v)'
        wireframe!(xs, ys, zs, aspect_ratio=:equal, color=:red)#, label="CGA")
        #wireframe!(x, y, z, aspect_ratio=:equal, color=:red, label="LOVO-CGA")
        return plt
    end
    if prob.name == "circle3d" || prob.name == "\tcircle3d"
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

function comparsol(prob, a1, a2, a3, a4)
    pyplot()#beckendpyplot
    plt = plot()
    if prob.name == "sphere2D" || prob.name == "\tsphere2D"
        plot!(plt, prob.data[:, 1], prob.data[:, 2], line=:scatter, aspect_ratio=:equal, lab="pontos do problema")
        θ = [0.0:2*π/360:2*π;]
        xs = a1[1] .+ a1[3] * cos.(θ)
        ys = a1[2] .+ a1[3] * sin.(θ)
        xss = a2[1] .+ a2[3] * cos.(θ)
        yss = a2[2] .+ a2[3] * sin.(θ)
        xk = a3[1] .+ a3[3] * cos.(θ)
        yk = a3[2] .+ a3[3] * sin.(θ)
        xkk = a4[1] .+ a4[3] * cos.(θ)
        ykk = a4[2] .+ a4[3] * sin.(θ)


        plot!(plt, xs, ys, color=:red, lab="solução L1")# legend=:outerbottomright)
        plot!(plt, xss, yss, color=:green, lab="solução L2")# legend=:outerbottomright)
        #plot!(plt, xk, yk, color=:blue, lab="solução L3")# legend=:outerbottomright)
        # plot!(plt, xkk, ykk, color=:black, lab="solução L4") # legend=:outerbottomright)
        plot!(plt, legend=:outerright)


        #plot!(plt, x, y, color=:green, lab="solução perfeita")
        display(plt)
    end
end



function plot_plane(prob, m::Vector{Float64}, d::Float64, m2::Vector{Float64}, d2::Float64)
    pyplot()#beckendpyplot
    plt = plot()
    plot!(plt, prob.data[:, 1], prob.data[:, 2], prob.data[:, 3], line=:scatter, aspect_ratio=:equal)
    u1 = [m[3], 0.0, m[1]]
    u2 = [m[1], m[3], 0.0]
    v1 = [-m2[3], 0.0, m2[1]]
    v2 = [m2[2], -m2[1], 0.0]
    v1 = v1 / norm(v1)
    v2 = v2 - (dot(v2, v1) / norm(v1)^2) * v1
    v2 = v2 / norm(v2)
    p0 = [d, d, d]
    p1 = [d2, d2, d2]
    #u = u / norm(u)
    #h = v - (dot(v, u) / norm(u)^2) * u
    #v = h / norm(h)
    # Cria um meshgrid para o plano
    x = range(-100, stop=100, length=10)
    y = range(-100, stop=100, length=10)


    xx, yy = [xi for xi in x, yi in y], [yi for xi in x, yi in y]





    # Calcula a equação do plano
    #z = (-m[1] .* xx .- m[2] .* yy .- d) / m[3]
    xs = p0[1] .+ xx .* u1[1] + yy .* u2[1]
    ys = p0[2] .+ xx .* u1[2] + yy .* u2[2]
    zs = p0[3] .+ xx .* u1[3] + yy .* u2[3]
    xn = p1[1] .+ xx .* v1[1] + yy .* v2[1]
    yn = p1[2] .+ xx .* v1[2] + yy .* v2[2]
    zn = p1[3] .+ xx .* v1[3] + yy .* v2[3]


    # Plota o plano
    wireframe!(xs, ys, zs, aspect_ratio=:equal, color=:blue, label="CGA")
    wireframe!(xn, yn, zn, aspect_ratio=:equal, color=:red, label="CGdA")
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

