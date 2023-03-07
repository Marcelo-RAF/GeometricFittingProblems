
using GeometricFittingProblems, DelimitedFiles

function run_sphere_2D(file::String,method::String,pinit=[0,0,1.0])
    set_problem = String.(readdlm(file))
    for probname ∈ set_problem
        prob = load_problem(probname)
        s = solve(prob,method,pinit)
        try solve(prob,method,pinit)
            solved = true
            #escreve probs_results.csv
        catch
            solved = false
        end
        if solved == true 
            @bench
            #escreve a resposta do benchmark em arquivo... bench_results.csv 
        else
            #escreva no arquivo de log
        end
    end
end
