module GeometricFittingProblems

using DelimitedFiles, LinearAlgebra, Plots

export load_problem, solve, build_problem, inverse_power_method, solve2, visualize, LMsphere, LMcircle, LMClass, LMClassCirc, geradoraut

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
    #println("First D")
    #display(D)
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
        xnorm = (1.0 / (F.vectors[:, indmin][end-1])) * F.vectors[:, indmin]
        center = xnorm[1:end-2]
        return push!(center, √(norm(center, 2)^2 - 2.0 * xnorm[end]))
    elseif method == "nullspace"
        P[end, :] = zeros(n + 2)
        np = nullspace(P)
        npnorm = np / np[end-1]
        centernp = npnorm[1:end-2]
        push!(centernp, √(norm(centernp, 2)^2 - 2.0 * npnorm[end]))
    end
end

function hildebran(data, method::String, ε=1.0e-5)
    (N, n) = size(data)
    D = [data'; -ones(1, N)]
    v = [-0.5 * norm(D[1:n, i], 2)^2 for i = 1:N]
    D = [D; v']
    Dd = D * D'
    #println(F.values)
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
        xnorm = (1.0 / (F.vectors[:, indmin][end])) * F.vectors[:, indmin]
        center = xnorm[1:end-2]
        return push!(center, √(norm(center, 2)^2 - 2.0 * xnorm[end-1]))
    elseif method == "nullspace"
        Dd[end, :] = zeros(n + 2)
        np = nullspace(Dd)
        npnorm = np / np[end]
        centernp = npnorm[1:end-2]
        return push!(centernp, √(norm(centernp, 2)^2 - 2.0 * npnorm[end-1]))
    end
end

function sort_plane_res(P, x, nout)
    n = length(P[:, 1])
    m = length(P[1, :])
    v = zeros(n)
    println(x)
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

function LMCircsq(data, nout, θ, ε=1.0e-6)
    ordres = sort_circle_res(data, θ, nout)
    antres = 0.0
    k = 1
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        θ = LMcircle(ordres[1], θ)
        ordres = sort_circle_res(data, θ, nout)
        k = k + 1
    end
    return θ, k
end


function LMSORT(data, nout, θ, ε=1.0e-8)
    ordres = sort_sphere_res(data, θ, nout)
    antres = 0.0
    k = 1
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        θ = LMsphere(ordres[1], θ)
        ordres = sort_sphere_res(data, θ, nout)
        k = k + 1
    end
    return θ, k
end



function LOVOCGAHypersphere(data, nout, θ, ε=1.0e-6)
    ordres = sort_sphere_res(data, θ, nout)
    k = 1
    antres = 0.0
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        θ = CGAHypersphere(ordres[1])
        ordres = sort_sphere_res(data, θ, nout)
        k = k + 1
    end
    return θ, k
end

function LOVOCGAHyperplane(data, nout, θ, ε=1.0e-6)
    ordres = sort_plane_res(data, θ, nout)
    k = 1
    antres = 0.0
    while abs(ordres[2] - antres) > ε
        antres = ordres[2]
        θ = CGAHypersphere(ordres[1])
        ordres = sort_plane_res(data, θ, nout)
        k = k + 1
    end
    return θ#, k
end



function solve(prob::FitProbType, method::String, initθ=CGAHypersphere(prob.data))
    if method == "LMsphere"
        return LMsphere(prob.data, initθ)
    end
    if method == "LMSORT"
        return LMSORT(prob.data, prob.nout, initθ)
    end
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
    if method == "LMCircsq"
        return LMCircsq(prob.data, prob.nout, initθ)
    end
end

function build_problem(probtype::String, limit::Vector{Float64}, params::Vector{Float64})
    if probtype == "line"
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
        FileMatrix = ["name :" "plane"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> p0 + λu + μv"; "dim :" 3; "cluster :" "false"; "noise :" "false"; "solution :" [push!(u)]; "description :" [[p0, p0]]]

        open("line_$(u[1])_$(u[2])_$(u[3])_$(nout).csv", "w") do io
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
        pp = range(-50.0, stop=50.0, length=npts)
        x = zeros(npts)
        y = zeros(npts)
        z = zeros(npts)
        vn = cross(u, v)
        vn = vn / norm(vn)
        for i = 1:npts
            for j = 1:npts
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
        r = 3.0
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
            z[iout[k]] = z[iout[k]] + rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])
        end
        FileMatrix = ["name :" "plane"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> p0 + λu + μv"; "dim :" 4; "cluster :" "false"; "noise :" "false"; "solution :" [push!(vn)]; "description :" [[p0, p0]]]

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
        vn = cross(u, v) / norm(cross(u, v))
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
            w[:, iout[k]] = w[:, iout[k]] + [rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;]), rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;]), rand([-(1 + 0.25)*r:0.1:(1+0.25)*r;])]
        end
        G = randn(3, npts)
        for i = 1:npts
            x[i] = w[1, i] + G[1, i]
            y[i] = w[2, i] + G[2, i]
            z[i] = w[3, i] + G[3, i]
        end
        FileMatrix = ["name :" "circle3d"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 7; "cluster :" "false"; "noise :" "false"; "solution :" [push!(c, r)]; "description :" [[u, v]]]

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
            x[iout[k]] = x[iout[k]] + rand([-(0.2)*r:0.1:(0.2)*r;])
            y[iout[k]] = y[iout[k]] + rand([-(0.2)*r:0.1:(0.2)*r;])   #rand([0.25*r:0.1*(r); (1 + 0.25) * r])
        end
        FileMatrix = ["name :" "sphere2D"; "data :" [[x y]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 - t[3]^2"; "dim :" 3; "cluster :" "false"; "noise :" "true"; "solution :" [push!(c, r)]; "description :" "type3: test sphere2d with noise and outliers"]

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
        dx = randn() * 0.1 * r
        dy = randn() * 0.1 * r
        dz = randn() * 0.1 * r
        for k = 1:nout
            x[iout[k]] = x[iout[k]] + dx #rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;]) 
            y[iout[k]] = y[iout[k]] + dy #rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
            z[iout[k]] = z[iout[k]] + dz #rand([-(1 + 0.15)*r:0.1:(1+0.15)*r;])
        end
        FileMatrix = ["name :" "sphere3D"; "data :" [[x y z]]; "npts :" npts; "nout :" nout; "model :" "(x,t) -> (x[1]-t[1])^2 + (x[2]-t[2])^2 +(x[3]-t[3])^2 - t[4]^2"; "dim :" 4; "cluster :" "false"; "noise :" "true"; "solution :" [push!(c, r)]; "description :" [[c, c]]]

        open("sphere3D_$(c[1])_$(c[2])_$(c[3])_$(c[4])_$(nout).csv", "w") do io #o que essa linha faz exatamente?
            writedlm(io, FileMatrix)
        end
    end

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
        J[i, end] = 1.0
    end
    return J
end


function LMClass(prob, xk, ε=1.0e-5, MAXIT=100)
    newdata = sort_sphere_res(prob.data, xk, prob.nout)
    R = fsphere(xk, newdata[1])
    J = jsphere(xk, newdata[1])
    (m, n) = size(J)
    Id = Matrix{Float64}(I, n, n)
    k = 1
    λ_up = 2.0
    λ_down = 2.0
    λ = 1.0
    μ = 0.7
    dk = 0.0
    while norm(J' * R, 2) > ε && k < MAXIT
        dk = (J' * J + λ * Id) \ ((-J') * R)
        md = 0.5 * (norm((R + J * dk), 2))^2 + λ * norm(dk, 2)^2
        Rd = fsphere(xk + dk, newdata[1])
        ρk = (0.5 * norm(R, 2)^2 - 0.5 * norm(Rd, 2)^2) / (0.5 * norm(R, 2)^2 - md)
        if ρk < μ
            λ = λ * λ_up
        else
            λ = λ / λ_down
            xk = xk + dk
            newdata = sort_sphere_res(prob.data, xk, prob.nout)
            R = fsphere(xk, newdata[1])
            J = jsphere(xk, newdata[1])
            k = k + 1
        end
    end
    return xk, k
end


function LMClassCirc(prob, xk, ε=1.0e-5, MAXIT=100)
    newdata = sort_circle_res(prob.data, xk, prob.nout)
    R = fcircle(xk, newdata[1])
    J = jcircle(xk, newdata[1])
    (m, n) = size(J)
    Id = Matrix{Float64}(I, n, n)
    k = 1
    λ_up = 2.0
    λ_down = 2.0
    λ = 1.0
    μ = 0.7
    dk = 0.0
    while norm(J' * R, 2) > ε && k < MAXIT
        dk = (J' * J + λ * Id) \ ((-J') * R)
        md = 0.5 * (norm((R + J * dk), 2))^2 + λ * norm(dk, 2)^2
        Rd = fcircle(xk + dk, newdata[1])
        ρk = (0.5 * norm(R, 2)^2 - 0.5 * norm(Rd, 2)^2) / (0.5 * norm(R, 2)^2 - md)
        if ρk < μ
            λ = λ * λ_up
        else
            λ = λ / λ_down
            xk = xk + dk
            newdata = sort_circle_res(prob.data, xk, prob.nout)
            R = fcircle(xk, newdata[1])
            J = jcircle(xk, newdata[1])
            k = k + 1
        end
    end
    return xk, k
end


function LMsphere(data, x0, ε=1.0e-5, λ_min=1e-4)
    k = 1
    x = x0
    R = fsphere(x, data)
    J = jsphere(x, data)
    (m, n) = size(J)
    xn = zeros(length(x))
    Id = Matrix{Float64}(I, n, n)
    λ = norm((J') * R, 2) / (norm(R, 2)^2)
    k1 = 2
    k2 = 1.5
    while norm((J') * R) > ε && k < 500
        d = (J' * J + λ * Id) \ ((-J') * R)
        xn = x + d
        if 0.5 * norm(fsphere(xn, data), 2)^2 < 0.5 * norm(fsphere(x, data), 2)^2
            x = xn
            if λ < λ_min
                λ = λ_min
            else
                λ = λ / k1
            end
            R = fsphere(x, data)
            J = jsphere(x, data)
        else
            λ = λ * k2
        end
        k = k + 1
    end
    return x, k
end

function LMplane(data, x0, ε=1.0e-5, λ_min=1e-4)
    k = 1
    x = x0
    R = fplane(x, data)
    J = jplane(x, data)
    (m, n) = size(J)
    xn = zeros(length(x))
    Id = Matrix{Float64}(I, n, n)
    λ = norm((J') * R, 2) / (norm(R, 2)^2)
    k1 = 2
    k2 = 1.5
    while norm((J') * R) > ε && k < 500
        d = (J' * J + λ * Id) \ ((-J') * R)
        xn = x + d
        if 0.5 * norm(fplane(xn, data), 2)^2 < 0.5 * norm(fplane(x, data), 2)^2
            x = xn
            if λ < λ_min
                λ = λ_min
            else
                λ = λ / k1
            end
            R = fplane(x, data)
            J = jplane(x, data)
        else
            λ = λ * k2
        end
        k = k + 1
    end
    x = x / norm(x)
    return x, k
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
    return sum(r)
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
    return sum(J, dims=1)[:, :]#J
end



function LMcircle(data, x0, ε=1.0e-6, λ_min=1e-4)
    k = 1
    x = x0
    R = fcircle(x, data)
    J = jcircle(x, data)
    (m, n) = size(J)
    xn = zeros(length(x))
    Id = Matrix{Float64}(I, n, n)
    λ = norm((J') * R, 2) / (norm(R, 2)^2)
    k1 = 2
    k2 = 1.5
    while norm((J') * R) > ε && k < 1000
        d = (J' * J + λ * Id) \ ((-J') * R)
        xn = x + d
        if 0.5 * norm(fcircle(xn, data), 2)^2 < 0.5 * norm(fcircle(x, data), 2)^2
            x = xn
            if λ < λ_min
                λ = λ_min
            else
                λ = λ / k1
            end
            R = fcircle(x, data)
            J = jcircle(x, data)
        else
            λ = λ * k2
        end
        k = k + 1
    end
    return x
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

function findcirc(s1, s2)
    vn = s1[1:end-1] - s2[1:end-1]
    vn = vn / norm(vn)
    r = norm(s1[1:end-1] - s2[1:end-1])
    r1 = s1[end]
    r2 = s2[end]
    t = (r2^2 - r1^2 - r^2) / (2 * r^2)
    rc = sqrt(r1^2 - (r1^2 - r2^2 + r^2)^2 / (4 * r^2))
    #rc = r^2*t + t*(r1^2 -r2^2 + r^2) + 
    d = 0.5 * (r1^2 - r2^2 + r^2) + dot(s1[1:end-1], s2[1:end-1] - s1[1:end-1])
    cc = s1[1:3] - t * (s2[1:3] - s1[1:3])
    sol = [vn; d; cc; rc]
    return sol
end

function geradoraut(h)
    n = 1000  # tamanho da lista desejada
    a = -120  # limite inferior do intervalo
    b = 120  # limite superior do intervalo
    c1 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
    c2 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
    c3 = round.(float(rand(n) .* (b - a) .+ a), digits=1)
    r = float(rand(5:170, 1000))
    npts = float(rand(10:2000, 401))
    nout = float([floor(Int, h * x) for x in npts])
    for i = 1:400
        build_problem("sphere2D", [1.0, 1.0], [c1[i], c2[i], r[i], npts[i], nout[i]])
    end
end

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


function visualize(prob, a)
    pyplot() #beckendpyplot
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
        plt = plot()
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
        x = h1 .+ prob.solution[4] * cos.(u) * sin.(v)'
        y = h2 .+ prob.solution[4] * sin.(u) * sin.(v)'
        z = h3 .+ prob.solution[4] * cos.(v)'
        xs = h4 .+ a[4] * cos.(u) * sin.(v)'
        ys = h5 .+ a[4] * sin.(u) * sin.(v)'
        zs = h6 .+ a[4] * cos.(v)'
        wireframe!(xs, ys, zs, aspect_ratio=:equal, color=:green, label="solução do algoritmo")
        wireframe!(x, y, z, aspect_ratio=:equal, color=:red, label="solução perfeita")
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

