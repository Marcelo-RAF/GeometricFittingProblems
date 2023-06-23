using BenchmarkTools, CSV, DataFrames

function benchLM(file::String, method::String, pinit=[0.0, 0.0, 1.0])
  set_problem = String.(readdlm(file))
  csv_file = open("LMSORT_20.csv", "w")
  csv_file_benchmark = open("benchLMSORT_20.csv", "w")
  df = DataFrame()
  k = 0
  benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("logLMSORT_20.txt", "w")
    prob = load_problem(probname)
    #pinit = CGAHypersphere(prob.data)
    x0 = LMsphere(prob.data, pinit)
    solved = false
    try
      s = solve(prob, method, x0)
      a = @benchmark solve($prob, $method, $x0) samples = 5000 #usa
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2])])
      df = vcat(df, row)
      benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)]) #usa
      benchmark_df = vcat(benchmark_df, benchmark_row) #usa
      CSV.write(csv_file, df)
      CSV.write(csv_file_benchmark, benchmark_df)
    catch e
      println("erro: ", e)
      solved = false
      write(log_file, "$probname\n")
    end
    close(log_file)
  end
  close(csv_file)
  close(csv_file_benchmark)
end

function benchCGA(file::String, method::String)
  set_problem = String.(readdlm(file))
  csv_file = open("NULLSPACE.csv", "w")
  #csv_file_benchmark = open("benchCGA_20.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("logcga.txt", "w")
    prob = load_problem(probname)
    #pinit = CGAHypersphere(prob.data)
    solved = false
    try
      s = solve(prob, method)#, pinit)
      #a = @benchmark solve($prob, $method, $pinit) samples = 5000 #usa
      k = k + 1
      println(k)
      ch = fsphere(s, prob.data)
      row = DataFrame([(probname, prob.npts, prob.solution, ch, s)]) #depois do prob.npts voltar o prob.nout
      df = vcat(df, row)
      #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)]) #usa
      #benchmark_df = vcat(benchmark_df, benchmark_row) #usa
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

function testeclass(file::String, pinit=[0.0, 0.0, 1.0])
  set_problem = String.(readdlm(file))
  csv_file = open("LMCLASS20.csv", "w")
  csv_file_benchmark = open("benchLMCLASS20.csv", "w")
  df = DataFrame()
  benchmark_df = DataFrame()
  k = 0
  for probname ∈ set_problem
    log_file = open("logLMCLASS_20.txt", "w")
    prob = load_problem(probname)
    solved = false
    try
      s = LMClass(prob, pinit)
      a = @benchmark LMClass($prob, $pinit) samples = 5000
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2])])
      df = vcat(df, row)
      benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)]) #usa
      benchmark_df = vcat(benchmark_df, benchmark_row) #usa
      CSV.write(csv_file, df)
      CSV.write(csv_file_benchmark, benchmark_df)
    catch e
      println("erro: ", e)
      solved = false
      write(log_file, "$probname\n")
    end
    close(log_file)
  end
  close(csv_file)
  close(csv_file_benchmark)
end

function testeLM(file::String, method::String, pinit=[[0, 0, 0.0, 1.0], 0.0])
  set_problem = String.(readdlm(file))
  csv_file = open("sollm3_40.csv", "w")
  csv_file_benchmark = open("benchlm3_40.csv", "w")
  df = DataFrame()
  k = 0
  benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("loglm4.txt", "w")
    prob = load_problem(probname)
    solved = false
    try
      s = solve(prob, method, pinit)
      a = @benchmark solve($prob, $method, $pinit) samples = 5000
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3])])
      df = vcat(df, row)
      benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)])
      benchmark_df = vcat(benchmark_df, benchmark_row)
      #df = DataFrame(solution_LOVOCGA = [s], prob_solution = [prob.solution])
      CSV.write(csv_file, df)
      CSV.write(csv_file_benchmark, benchmark_df)
    catch e
      println("erro: ", e)
      solved = false
      write(log_file, "$probname\n")
    end
    close(log_file)
  end
  close(csv_file)
  close(csv_file_benchmark)
end

function calcgradresid(file::String)
  set_problem = String.(readdlm(file))
  csv_file = open("hildebranautovalor.csv", "w")
  #csv_file_benchmark = open("benchLMSORT_20.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("log.txt", "w")
    prob = load_problem(probname)
    s = hildebran(prob.data)
    solved = false
    try
      ch1 = residalgebric(s, prob.data)
      ch2 = residgeometric(s, prob.data)
      ch3 = gradalgebric(s, prob.data)
      ch4 = gradgeometric(s, prob.data)
      #a = @benchmark solve($prob, $method, $x0) samples = 5000 #usa
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, s, ch1, ch2, ch3, ch4)])
      df = vcat(df, row)
      #benchmark_row = DataFrame([(probname, prob.npts, prob.nout, minimum(a.times) / 1e9, median(a.times) / 1e9, maximum(a.times) / 1e9)]) #usa
      #benchmark_df = vcat(benchmark_df, benchmark_row) #usa
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