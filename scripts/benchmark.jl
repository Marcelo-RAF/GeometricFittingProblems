using BenchmarkTools, CSV, DataFrames




function testes(x) #esse arquivo está na pasta resultsLM\Cubic\semruido
  prob = load_problem("cubic_0.1_-7.7_-4.0_10.csv")
  (m, n) = size(prob.data)
  s = zeros(m)
  cl2(x) = invokelatest(prob.model(x, t))
  for i = 1:m
    s[i] = invokelatest(prob.model, x, prob.data[i, :])
  end
  k = LMPers(([1.0, 1.0, 1.0, 1.0], 0), prob.model, prob.data, prob.dim, sort_funcion_res, prob.nout)
  return k
end


function benchmarklovo(namecsv::String, file::String, method::String, Ord, xk)
  if method == "LMPersistent"
    set_problem = String.(readdlm(file))
    csv_file = open(namecsv, "w")
    #csv_file_benchmark = open("benchlm3_40.csv", "w")
    df = DataFrame()
    k = 0
    #benchmark_df = DataFrame()
    for probname ∈ set_problem
      log_file = open("logpersistentlovo.txt", "w")
      prob = load_problem(probname)
      solved = false
      try
        #x = Levenberg(fcubic, jcubic, [1.0, 1.0, 1.0, 1.0], prob.data)
        s = LMPers(xk, prob.model, prob.data, prob.dim, Ord, prob.nout)
        a = @benchmark LMPers($xk, $prob.model, $prob.data, $prob.dim, $Ord, $prob.nout) samples = 100 seconds = 70
        #ndif = norm(prob.solution - s[1])
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], s[4], median(a.times) / 1e9)])
        df = vcat(df, row)
        #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)])
        #benchmark_df = vcat(benchmark_df, benchmark_row)
        #df = DataFrame(solution_LOVOCGA = [s], prob_solution = [prob.solution])
        CSV.write(csv_file, df)
        #CSV.write(csv_file_benchmark, benchmark_df)
      catch e
        println("erro: ", e)
        solved = false
        write(log_file, "$probname\n")
      end
      close(log_file)
    end
    close(csv_file)
    #close(csv_file_benchmark)
  end
  if method == "LMClass"
    set_problem = String.(readdlm(file))
    csv_file = open(namecsv, "w")
    #csv_file_benchmark = open("benchlm3_40.csv", "w")
    df = DataFrame()
    k = 0
    #benchmark_df = DataFrame()
    for probname ∈ set_problem
      log_file = open("loglmclass.txt", "w")
      prob = load_problem(probname)
      solved = false
      try
        #x = Levenberg(fcubic, jcubic, [1.0, 1.0, 1.0, 1.0], prob.data)
        s = LMLOVO(xk, prob.model, prob.data, prob.dim, Ord, prob.nout)
        a = @benchmark LMLOVO($xk, $prob.model, $prob.data, $prob.dim, $Ord, $prob.nout) samples = 100 seconds = 70
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], median(a.times) / 1e9)])
        df = vcat(df, row)
        #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)])
        #benchmark_df = vcat(benchmark_df, benchmark_row)
        #df = DataFrame(solution_LOVOCGA = [s], prob_solution = [prob.solution])
        CSV.write(csv_file, df)
        #CSV.write(csv_file_benchmark, benchmark_df)
      catch e
        println("erro: ", e)
        solved = false
        write(log_file, "$probname\n")
      end
      close(log_file)
    end
    close(csv_file)
    #close(csv_file_benchmark)
  end
end


function benchmk(namecsv::String, file::String, method::String, funcao, jacobiana, Ord, xk)
  if method == "LMPersistent"
    set_problem = String.(readdlm(file))
    csv_file = open(namecsv, "w")
    #csv_file_benchmark = open("benchlm3_40.csv", "w")
    df = DataFrame()
    k = 0
    #benchmark_df = DataFrame()
    for probname ∈ set_problem
      log_file = open("logpersistentlovo.txt", "w")
      prob = load_problem(probname)
      solved = false
      try
        #x = Levenberg(fcubic, jcubic, [1.0, 1.0, 1.0, 1.0], prob.data)
        s = LMPersistent(funcao, jacobiana, Ord, xk, prob.data, prob.nout)
        a = @benchmark LMPersistent($funcao, $jacobiana, $Ord, $xk, $prob.data, $prob.nout) samples = 100 seconds = 70
        #ndif = norm(prob.solution - s[1])
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], s[4], median(a.times) / 1e9)])
        df = vcat(df, row)
        #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)])
        #benchmark_df = vcat(benchmark_df, benchmark_row)
        #df = DataFrame(solution_LOVOCGA = [s], prob_solution = [prob.solution])
        CSV.write(csv_file, df)
        #CSV.write(csv_file_benchmark, benchmark_df)
      catch e
        println("erro: ", e)
        solved = false
        write(log_file, "$probname\n")
      end
      close(log_file)
    end
    close(csv_file)
    #close(csv_file_benchmark)
  end
  if method == "LMClass"
    set_problem = String.(readdlm(file))
    csv_file = open(namecsv, "w")
    #csv_file_benchmark = open("benchlm3_40.csv", "w")
    df = DataFrame()
    k = 0
    #benchmark_df = DataFrame()
    for probname ∈ set_problem
      log_file = open("loglmclass.txt", "w")
      prob = load_problem(probname)
      solved = false
      try
        #x = Levenberg(fcubic, jcubic, [1.0, 1.0, 1.0, 1.0], prob.data)
        s = LovoLM(funcao, jacobiana, Ord, xk, prob.data, prob.nout)
        a = @benchmark LovoLM($funcao, $jacobiana, $Ord, $xk, $prob.data, $prob.nout) samples = 100 seconds = 70
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], median(a.times) / 1e9)])
        df = vcat(df, row)
        #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)])
        #benchmark_df = vcat(benchmark_df, benchmark_row)
        #df = DataFrame(solution_LOVOCGA = [s], prob_solution = [prob.solution])
        CSV.write(csv_file, df)
        #CSV.write(csv_file_benchmark, benchmark_df)
      catch e
        println("erro: ", e)
        solved = false
        write(log_file, "$probname\n")
      end
      close(log_file)
    end
    close(csv_file)
    #close(csv_file_benchmark)
  end
end






