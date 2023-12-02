using BenchmarkTools, CSV, DataFrames

function difer(u, v)
  h = 0.0
  for j = 1:length(u)-1
    h = h + sqrt((u[j] - v[j])^2)
  end
  h = h + sqrt((u[end] - abs(v[end]))^2)
  return h
end

function solqual(u, v)
  h = 0.0
  h = norm(u - v)
  return h
end

function benchLOVO(namecsv::String, file::String, method::String, Function, Jacobian, Ord, qsol, xk)
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
        s = LMPersistent(Function, Jacobian, Ord, xk, prob.data, prob.nout)
        a = @benchmark LMPersistent($Function, $Jacobian, $Ord, $xk, $prob.data, $prob.nout) samples = 100 seconds = 70
        ndif = qsol(s[1], s[4])
        #ndif = norm(prob.solution - s[1])
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], ndif, median(a.times) / 1e9)])
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
        s = LovoLM(Function, Jacobian, Ord, xk, prob.data, prob.nout)
        a = @benchmark LovoLM($Function, $Jacobian, $Ord, $xk, $prob.data, $prob.nout) samples = 100 seconds = 70
        ndif = qsol(s[1], s[3])
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], ndif, median(a.times) / 1e9)])
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

function benchLOVOwnoise(namecsv::String, file::String, method::String, Function, Jacobian, Ord, qsol, xk)
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
        s = LMPersistent(Function, Jacobian, Ord, xk, prob.data, prob.nout)
        a = @benchmark LMPersistent($Function, $Jacobian, $Ord, $xk, $prob.data, $prob.nout) samples = 100 seconds = 70
        ndif = qsol(prob.solution, s[1])
        #ndif = norm(prob.solution - s[1])
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], ndif, median(a.times) / 1e9)])
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
        s = LovoLM(Function, Jacobian, Ord, xk, prob.data, prob.nout)
        a = @benchmark LovoLM($Function, $Jacobian, $Ord, $xk, $prob.data, $prob.nout) samples = 100 seconds = 70
        ndif = qsol(prob.solution, s[1])
        k = k + 1
        println(k)
        row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], ndif, median(a.times) / 1e9)])
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


function chamadasbench() #comecar 2d com ruido em Levtest
  benchLOVO("LMPersnew.csv", "testenomes.txt", "LMPersistent", fsphere, jsphere, sort_sphere_res, residalgebric, [[1.0, 1.0, 1.0], 0])
  cd("..")
  cd("semruido")
  benchLOVOwnoise("LMPersnew.csv", "testenomes.txt", "LMPersistent", fsphere, jsphere, sort_sphere_res, difer, [[1.0, 1.0, 1.0], 0])
  cd("..")
  cd("..")
  cd("sphere3d\\semruido\\20")
  benchLOVOwnoise("LMPersnew.csv", "testenomes.txt", "LMPersistent", fsphere, jsphere, sort_sphere_res, difer, [[1.0, 1.0, 1.0, 1.0], 0])
  cd("..")
  cd("30")
  benchLOVOwnoise("LMPersnew.csv", "testenomes.txt", "LMPersistent", fsphere, jsphere, sort_sphere_res, difer, [[1.0, 1.0, 1.0, 1.0], 0])
  cd("..")
  cd("..")
  cd("ruido\\20")
  benchLOVO("LMPersnew.csv", "testenomes.txt", "LMPersistent", fsphere, jsphere, sort_sphere_res, residalgebric, [[1.0, 1.0, 1.0, 1.0], 0])
  cd("..")
  cd("30")
  benchLOVO("LMPersnew.csv", "testenomes.txt", "LMPersistent", fsphere, jsphere, sort_sphere_res, residalgebric, [[1.0, 1.0, 1.0, 1.0], 0])
end

function spherebench()
  benchLOVOwnoise("LMClass.csv", "testenomes.txt", "LMClass", fline, jline, sort_line, solqual, [1.0, 1.0])
  benchLOVOwnoise("LMPers.csv", "testenomes.txt", "LMPersistent", fline, jline, sort_line, solqual, [[1.0, 1.0], 0])
  cd("..")
  cd("ruido")
  benchLOVO("LMClass.csv", "testenomes.txt", "LMClass", fline, jline, sort_line, residline, [1.0, 1.0])
  benchLOVO("LMPers.csv", "testenomes.txt", "LMPersistent", fline, jline, sort_line, residline, [[1.0, 1.0], 0])
end


function testeLM(file::String, method::String)
  set_problem = String.(readdlm(file))
  csv_file = open("persistentLM2.csv", "w")
  #csv_file_benchmark = open("benchlm3_40.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("loglm4.txt", "w")
    prob = load_problem(probname)
    solved = false
    try
      s = solve(prob, method)
      a = @benchmark solve($prob, $method) samples = 100 seconds = 60
      #ndif = norm(prob.solution - s[1])
      ndif = difer(prob.solution, s[1])
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], s[3], ndif, median(a.times) / 1e9)])
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

function testeLMClass(file::String)
  set_problem = String.(readdlm(file))
  csv_file = open("LMClass2.csv", "w")
  #csv_file_benchmark = open("benchlm3_40.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("loglm4.txt", "w")
    prob = load_problem(probname)
    solved = false
    try
      s = LMClass(prob, [0.0, 0.0, 0.0, 1.0])
      a = @benchmark LMClass($prob, $[0.0, 0.0, 0.0, 1.0]) samples = 100 seconds = 60
      ndif = difer(prob.solution, s[1])
      #ndif = norm(prob.solution - s[1])
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, prob.nout, prob.solution, s[1], s[2], ndif, median(a.times) / 1e9)])
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

function benchcircle(file::String)
  set_problem = String.(readdlm(file))
  csv_file = open("hildnullspace.csv", "w")
  #csv_file_benchmark = open("benchLMSORT_20.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("log.txt", "w")
    prob = load_problem(probname)
    s = hildebrancircle(prob.data, "nullspace")
    solved = false
    try
      ch1 = residcircle(s, prob.data)
      a = @benchmark hildebrancircle($prob.data, $"nullspace") samples = 5000 seconds = 30#usa
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, s, ch1, median(a.times) / 1e9)])
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

function withoutres(file::String)
  set_problem = String.(readdlm(file))
  csv_file = open("hildebran.csv", "w")
  #csv_file_benchmark = open("benchLMSORT_20.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("log.txt", "w")
    prob = load_problem(probname)
    s = solve(prob, "LOVO-HildSphere")
    h = norm(s[1] - prob.solution)
    solved = false
    try
      a = @benchmark solve($prob, $"LOVO-HildSphere")# samples = 5000 seconds = 20 #usa
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, s[1], s[2], h, median(a.times) / 1e9)])
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

function withres(file::String)
  set_problem = String.(readdlm(file))
  csv_file = open("hildebran.csv", "w")
  #csv_file_benchmark = open("benchLMSORT_20.csv", "w")
  df = DataFrame()
  k = 0
  #benchmark_df = DataFrame()
  for probname ∈ set_problem
    log_file = open("log.txt", "w")
    prob = load_problem(probname)
    s = solve(prob, "LOVO-HildSphere")
    solved = false
    try
      ch1 = residalgebric(s[1], prob.data)
      ch2 = residgeometric(s[1], prob.data)
      a = @benchmark solve($prob, $"LOVO-HildSphere")# samples = 5000 seconds = 20 #usa
      k = k + 1
      println(k)
      row = DataFrame([(probname, prob.npts, s[1], s[2], ch1, ch2, median(a.times) / 1e9)])
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
