#=
 = Split a set of observation locations into locations of good and poor temporal coverage
 = (e.g., as given by the minimum ucur coverage among all analyses) throughout a core period
 = (i.e., most of 2012-2014) and save the two location subsets as separate files.  Call these
 = locations of good coverage "the collocations".  Then create one-sided NFFT power spectral
 = coefficients for the collocations.  Note that all data gaps in any one analysis are mapped
 = to the other analyses and that a nonequidistant FFT algorithm (NFFT; www.nfft.org) is thus
 = employed - RD May 2016.
 =#

using My, NFFT
const UCUR             = 10                             # indecies of all data variables
const VCUR             = 11

const TIMTOT           = 3408                           # length of the full timeseries
const TIMSTA           =  180                           #  start of the mid-October 2012-2014 period of interest
const TIMLEN           = 2920                           # length of the mid-October 2012-2014 period of interest (and number of spectral coefficients)
const CUTOFF           = 2000                           # minimum number from all analyses during period of interest
const MISS             = -9999.0                        # generic missing value
const CALIB            = 0                              # flag determining whether to calibrate (1 = yes)

if (argc = length(ARGS)) != 2
  print("\nUsage: jjj $(basename(@__FILE__)) buoydata_1993_2014_drogON.asc.nonmdt.locate_2.0_calib ucur\n\n")
  exit(1)
end

vind = 0                                                                      # identify the variable of interest
ARGS[2] == utf8("ucur") && (vind = UCUR)                                      # (note ARGS are UTF8 by default)
ARGS[2] == utf8("vcur") && (vind = VCUR)
if vind == 0
  write(STDERR, "ERROR : didn't recognize $(ARGS[2]) among ucur/vcur\n")
  exit(-1)
end

torus(x) = (x - TIMSTA) / TIMLEN - 0.5                                        # map of daily index to NFFT torus [-0.5,0.5)
dirs   = ["v2.0_global_025_deg_ekman_15m", "v2.0_global_025_deg_ekman_hs", "v2.0_global_025_deg_geostrophic", "v2.0_global_025_deg_total_15m", "v2.0_global_025_deg_total_hs"]
calalp = [                            2.2,                            4.1,                               3.4,                             3.4,                            3.3]
calbet = [                           0.96,                           0.92,                              0.82,                            0.98,                           0.91]
dirn   = length(dirs)

fpa = My.ouvre(ARGS[1],                           "r")
fpb = My.ouvre(ARGS[1]*"."*ARGS[2]*".got$CUTOFF", "w")
fpc = My.ouvre(ARGS[1]*"."*ARGS[2]*".not$CUTOFF", "w")

for line in eachline(fpa)                                                     # for each location, read all analysis data
  data = Array(Float64, TIMLEN, dirn) ; (lat, lon, num) = float(split(line))
  tmp = @sprintf("%9.3f.%9.3f", lat, lon) ; tail = replace(tmp, " ", ".")
  for (a, dir) in enumerate(dirs)
    fpd = My.ouvre("$dir/$dir.$tail", "r", false) ; lines = readlines(fpd) ; close(fpd)
    for b = TIMSTA:TIMSTA + TIMLEN - 1  data[b-TIMSTA+1,a] = float(split(lines[b])[vind])  end
  end

  mask = ones(TIMLEN)                                                         # identify valid days, mapping gaps in
  for a = 1:dirn, b = 1:TIMLEN                                                # any one analysis to all
    if data[b,a] < -333 || data[b,a] > 333
      mask[b] = 0.0
    end
  end

  if sum(mask) < CUTOFF                                                       # augment one of the subset list files
    write(fpc, line)                                                          # and only create a spectra file if all
  else                                                                        # analyses have good temporal coverage
    write(fpb, line)

    half = div(TIMLEN, 2)
    spec = Array(Float64, half + 1, dirn)
    datb = Array(Float64, 0)
#   datc = Array(Complex{Float64}, 0)
    datc = Array(Float64, 0)
    for (a, dir) in enumerate(dirs)
#     if (a == 5 && (vind == AIRT || vind == SSTT)) || (a == 1 &&  vind == LHFX)
#       for b = 1:half + 1  spec[b,a] = MISS  end
#     else
        datb = Array(Float64, 0)                                              # define the NFFT sampled data with
#       datc = Array(Complex{Float64}, 0)                                     # endpoints at 2012-10-15 and 2014-10-15
        datc = Array(Float64, 0)
        for b = TIMSTA:TIMSTA + TIMLEN - 1
#         if -333 < data[b-TIMSTA+1,a] < 333
          if mask[b-TIMSTA+1] == 1.0
            push!(datb, torus(b))
            push!(datc, data[b-TIMSTA+1,a])
          end
        end
#       nums = length(datb)                                                   # get the spectral coefficients and compute the

        nums = length(datb)                                                   # (call out an estimate of each FFT spectrum)
        tmpy = "fft/$(ARGS[2]).$tail.$dir.fftt"
        tmpz = "fft/$(ARGS[2]).$tail.$dir.ffttest"
        fpz = My.ouvre(tmpy, "w", false)
        for z = 1:nums
          linz = @sprintf("%15.8f %15.8f\n", datb[z], datc[z])
          write(fpz, linz)
        end
        close(fpz)

        run(`/home1/homedir1/perso/rdaniels/bin/diag.trajectory.drifters.timeseries.nfft $tmpy`)
        fpz = My.ouvre(tmpz, "r", false) ; linez = readlines(fpz) ; close(fpz)
        for (z, linz) in enumerate(linez)
          spec[z,a] = float(split(linz))[2]
        end
        rm(tmpy) ; rm(tmpz)

#       plan = NFFTPlan(datb, nums)                                           # one-sided spectra for each analysis, where FFT
#       flan = nfft_adjoint(plan, datc)                                       # requires normalization by TIMES^2; NFFT employs
#       absf = abs2(flan)                                                     # TIMES*obsnum to satisfy Parseval's equation
#       for b = 1:half - 1
#         spec[b+1,a] = (absf[half+1+b] + absf[half+1-b]) / float(TIMLEN * nums)
#       end
#       spec[     1,a] = absf[half+1]                     / float(TIMLEN * nums)
#       spec[half+1,a] = absf[     1]                     / float(TIMLEN * nums)
#     end
    end

    if CALIB == 0  tmp = "fft/$(ARGS[2]).$tail.fft"                           # and store the spectra in an fft dir
    else           tmp = "fft/$(ARGS[2]).$tail.fftcal"  end
    fpd = My.ouvre(tmp, "w", false)
    for a = 1:half + 1
      tmp = @sprintf("%15.8f %15.8f %15.8f %15.8f %15.8f %15.8f\n",
        (a - 1.0) / float(TIMLEN), spec[a,1], spec[a,2], spec[a,3], spec[a,4], spec[a,5])
      write(fpd, tmp)
    end
    close(fpd)
  end
end

close(fpa)
close(fpb)
close(fpc)
exit(0)


#=
  for a = 1:TIMLEN
    if mask[a] == 0.0
      for b = 1:dirn
        data[a,b] = MISS
      end
    end
  end

bartl(x) = x < 1736.5 ? (x - TIMSTA) / 1277.5 : (3014 - x) / 1277.5           # map daily index to Bartlett triangular weights
bartl(x) = 1.0

using My, NFFT
const MISS             = -9999.0
const DAYS             = 3745
const TIMES            = 2556
const CALIB            = 0
dirs   = ["cfsr"]
dlen   = length(dirs) + 1
data   = Array(Float64, dlen, DAYS)
torus(x) = (x - TIMSTA) / TIMES - 0.5
bartl(x) = 1.0
tmpa = "cfsr....27.500...-39.000"
fpa = My.ouvre(tmpa, "r", false) ; lines = readlines(fpa) ; close(fpa)
a = 1
for (b, line) in enumerate(lines)
  vals = float(split(line))
  data[a,b] = vals[1]
end
for a = 1:DAYS                                                                # mask days missing any one analysis
  data[dlen,a] = 1
  for b = 1:dlen - 1
    if -333 >= data[b,a] || data[b,a] >= 3333
      data[dlen,a] = 0
    end
  end
end
datb = Array(Float64, 0)
datc = Array(Complex{Float64}, 0)
push!(datb, torus(TIMSTA))
push!(datc, bartl(TIMSTA) * 0.0im)
a = 1
for b = 460:3013
  if data[dlen,b] > 0
    push!(datb, torus(b))
    push!(datc, bartl(b) * data[a,b])
  end
end
push!(datb, torus(3014))
push!(datc, bartl(3014) * 0.0)

nums = length(datb)
p = NFFTPlan(datb, nums)
f = nfft_adjoint(p, datc)
# g = nfft(p, f)
absf = abs2(f)
datd = Array(Float64, 0)
push!(datd, absf[div(TIMES,2)+1])
for a = 1:TIMES/2-1
  push!(datd, absf[TIMES/2+1+a] + absf[TIMES/2+1-a])
end
push!(datd, absf[1])
datd /= TIMES * nums
for a = 1:TIMES/2+1
  @printf("%15.8f %15.8f\n", (float)a / (float)TIMES, datd[a])
end

    for (a = 1; a <= TIMES/2+1; a++) {
      if (a == 1)
        tmp = powf(p.f_hat[TIMES/2+1][0], 2.0) + powf(p.f_hat[TIMES/2+1][1], 2.0);
      else if (a == TIMES/2+1)
        tmp = powf(p.f_hat[1        ][0], 2.0) + powf(p.f_hat[1        ][1], 2.0);
      else
        tmp = powf(p.f_hat[TIMES/2+1+a][0], 2.0) + powf(p.f_hat[TIMES/2+a][1], 2.0) +
              powf(p.f_hat[TIMES/2+1-a][0], 2.0) + powf(p.f_hat[TIMES/2-a][1], 2.0);
      tmp /= (float)TIMES * (float)obsnum;
      fprintf(fpa, "%15.8f %15.8f\n", (float)a / (float)TIMES, tmp);
=#
