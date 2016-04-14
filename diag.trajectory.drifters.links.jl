#=
 = Create local links to all files in a directory and any subdirectories
 = (note that a link is created only to the first file found, where two
 =  files from different subdirs have the same name) - RD September 2015
 =#

if size(ARGS) != (1,)
  print("\nUsage: jjj $(basename(@__FILE__)) /home/cercache/project/oceanheatflux/data/references/ifremerflux\n\n")
  exit(1)
end

function walkdir(dir)
  absdir = abspath(dir)
  files  = readdir(absdir)
  for file in files
    full = joinpath(absdir, file)
    if isfile(full) && endswith(full, ".nc")
      if !islink(file)
        symlink(full, file)
      else
        print("$file link already exists, so not replacing with a link to $full\n")
      end
    elseif isdir(full)
      walkdir(full)
    end
  end
end

walkdir(ARGS[1])
