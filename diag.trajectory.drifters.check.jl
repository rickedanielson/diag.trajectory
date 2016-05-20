#=
 = Identify all files like 2010-06-17-0000.00092861 in the current dir and check
 = that they tracked for longer than a day (or that their corresponding simulated
 = trajectory files also exist) - RD August 2015
 =#

if (argc = length(ARGS)) != 0
  print("\nUsage: jj $(basename(@__FILE__))\n\n")
  exit(1)
end

files = filter(x -> ismatch(r"^\d{4}-\d{2}-\d{2}-\d{4}\.\d{8}$", x), readdir("."))
# (count,) = size(files) ; print("found $count files\n")

for (a, file) in enumerate(files)
  filf = file * ".gif"
# filf = file * ".traj." * file[1:15]
  if file[1:3] == "201" && !isfile(filf)  print("$file missing $filf\n")  end
# if filesize(file) < 1000  print("mv $file ../limbo\n")  end
# if filesize(file) < 1000  print("wc $file\n")  end
end
