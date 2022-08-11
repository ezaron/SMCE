
# Just load the needed components:
include("edznc.jl")
try
    period("M2")
catch
    include("/home/ezaron/NASA-Tides/iho/iho_func.jl")
    using .HAmod
end
using Interpolations
using ColorSchemes
using ElectronDisplay
using CairoMakie
#using GLMakie  ?
#using WGLMakie ?
using Printf

# Input file containing lat, lon, time:
infile = "SWOT_L2_LR_SSH_Expert_018_063_20150404T095346_20150404T104512_DG10_01.nc"
# List the tidal frequencies you want to predict. This is largest set:
#cidvec = ["M2", "S2", "K1", "O1", "MA2", "MB2"]
cidvec = ["M2", "S2", "K1", "O1"] # For checking with precomputed values in the file.

lonVars = ["longitude", "lon"]
latVars = ["latitude", "lat"]
timeVars = ["time"]
lon = ncvarget(infile,lonVars)
lat = ncvarget(infile,latVars)
time = ncvarget(infile,timeVars)
tunits = ncgetatt(infile,timeVars,"units")

if (size(lon) != size(lat))
    println("ERROR: Sizes on longitude and latitude arrays differ.")
end

# Try to figure out the units and reference time from tunits:
tunits = ncgetatt(infile,timeVars,"units")
# Assume the time unit, e.g., seconds, is the first substring in tunits:
m=split(tunits," ")[1]
# Convert time to units of days:
if (m == "seconds")
    time = time/8.64e4
elseif (m == "hours")
    time = time/24.0
elseif (m == "days")
else
    println("ERROR: Could not parse the time units of the input file.")
end
# Get the reference time:
ymd=match(r"(....)-(..)-(..)",tunits)
hms=match(r"(..):(..):(..)",tunits)
if (isnothing(hms))
    hr = 0.0
    mi = 0.0
    se = 0.0
else
    try
        hr = parse(Float64,hms[1])
        mi = parse(Float64,hms[2])
        se = parse(Float64,hms[3])
        global dayfrac = (hr + (mi + se/60.0)/60.0)/24.0
    catch
        println("ERROR: Could not parse the reference time of the input file: bad HH:MM:SS.")
    end
end
try
    global yy = parse(Float64,ymd[1])
    global mm = parse(Float64,ymd[2])
    global dd = parse(Float64,ymd[3])
catch
    println("ERROR: Could not parse the reference time of the input file: bad or missing YYYY-MM-DD.")
end
jd0 = ymd2jd(yy,mm,dd + dayfrac)
# Convert time to decimal Julian Date which begins at midnight rather than noon:
time = jd0 .+ time
# Sanity check:
yy,mm,dd,hr,mi,se = jd2ymdhms(minimum(time))
start_str=@sprintf("First measurement time = %04i-%02i-%02i %02i:%02i:%05.2f",
                   yy,mm,dd,hr,mi,se)
yy,mm,dd,hr,mi,se = jd2ymdhms(maximum(time))
stop_str=@sprintf("Last measurement time = %04i-%02i-%02i %02i:%02i:%05.2f",
                  yy,mm,dd,hr,mi,se)
println(start_str)
println(stop_str)
# This would be a good place to plot the input mesh and indicate
# the start and end times on the plot.

rank=length(size(lon))
f = Figure(resolution=(360,180))
ax = Axis(f[1, 1],
          xlabel="longitude",
          ylabel="latitude",
          title="data locations")
if (rank == 2)
    lines!(ax, lon[1,:], lat[1,:])
    lines!(ax, lon[end,:], lat[end,:])
    lines!(ax, lon[:,1], lat[:,1])
    lines!(ax, lon[:,end], lat[:,end])
    if (prod(size(lon)) < 10000)
        n,m=size(lon)
        for j=1:m
            scatter!(ax,lon[:,j],lat[:,j])
        end
    end
elseif (rank == 1)
    lines!(ax, lon, lat)
else
    println("Not sure how to show (lat,lon). Skipping graphics.")
end
f

# Note how we force time to look like a 1-d vector even if it is
# a matrix:
indu = find(!isnan,time)
iid,F = FMAT(time[indu],cidvec,1,0) # Contains nodal correction
# iid,F = FMAT(time[indu],cidvec,0,0) # Without nodal correction

# Load data for each component frequency one at a time:
fhret = "/home/ezaron/FFTest/HRET8.1/HRET_v8.1.nc"
hlon = ncvarget(fhret,"longitude")
hlat = ncvarget(fhret,"latitude")
hpred = zeros(size(lon))
# Unwrap longitude into [0,360] if needed:
lonx = copy(lon)
indw = find( x -> x < 0.0,lonx)
lonx[indw] = lon[indw] .+ 360.0
for cid = cidvec
    println("Working on cid = ",cid)
    ic = cid * "re"
    is = cid * "im"
    iic = iid[cid * "c"]
    iis = iid[cid * "s"]
    hc = ncvarget(fhret,ic)
    hs = ncvarget(fhret,is)
    ihc = interpolate((hlon,hlat),hc,(Gridded(Linear()),Gridded(Linear())))
    ihs = interpolate((hlon,hlat),hs,(Gridded(Linear()),Gridded(Linear())))
    ehc = extrapolate(ihc,(Periodic(),Flat()))
    ehs = extrapolate(ihs,(Periodic(),Flat()))
    if (size(lon) == size(time))
        xc = ehc.(lonx[indu],lat[indu])
        xs = ehs.(lonx[indu],lat[indu])
        hpred[indu] = hpred[indu] + F[:,iic].*xc + F[:,iis].*xs
        continue
    end
    # Different cases when size(time) != size(lat)
    if (length(size(time)) == 1)
        na,nb=size(lat)
        for k=indu
            if (length(time) == na)
                lon1 = view(lonx,k,:)
                lat1 = view(lat,k,:)
                h1   = view(hpred,k,:)
            else
                lon1 = view(lonx,:,k)
                lat1 = view(lat,:,k)
                h1   = view(hpred,:,k)
            end
            xc = ehc.(lon1,lat1)
            xs = ehs.(lon1,lat1)
            h1[:] = h1[:] + F[k,iic].*xc + F[k,iis].*xs
        end
    end
end
    
checkval = ncvarget(infile,"internal_tide_hret")
using Statistics
std(checkval) # 3mm
# Hmmm. The errors are little larger than I would have expected:
std(checkval - hpred) # 0.1mm
maximum(checkval - hpred) #  0.09mm
minimum(checkval - hpred) # -1.12mm
# Try turning off the nodal correction:
std(checkval - hpred) # 0.03mm
maximum(checkval - hpred) #  0.2mm
minimum(checkval - hpred) #  0.2mm
# Seems like the difference is plausibly related to the quantization
# error in the compressed HRET file Remko used.

tmp = split(infile,".")
tmp[end] = "_hret.nc"
fout = reduce(*,tmp)
run(`rm -f $fout`)

if (length(size(lon)) == 2)
    num_pixels,num_lines = size(lon)
    
    nccreate(fout,"hret",
             "num_pixels",num_pixels,
             "num_lines",num_lines,
             atts=Dict("longname" => "predicted internal tide sea level anomaly",
                       "units" => "meter"),
             gatts=Dict("creator" => "Software written by Edward D. Zaron, 2022-08-10.",
                        "contact" => "edward.d.zaron@oregonstate.edu",
                        "cidvec"  => cidvec,
                        "project" => "NASA SWOT Science Team",
                        "history" => "Based on https://ingria.ceoas.oregonstate.edu/fossil/HA checkout dde20225c3ed7e3364ffc282ff53b4832980c930",
                        "version" => "1.0",
                        "source"  => "https://ingria.ceoas.oregonstate.edu/fossil/SMCE checkout d43216f1d9981b2fd151404ecce38a9220a4f6ce",
                        "infile"  => infile,
                        "outfile" => fout,
                        "fhret"   => fhret
                        )
             )
    nccreate(fout,"longitude",
             "num_pixels",
             "num_lines",
             atts=Dict("longname" => "longitude, copied from infile")
             )
    nccreate(fout,"latitude",
             "num_pixels",
             "num_lines",
             atts=Dict("longname" => "latitude, copied from infile")
             )

    ncwrite(lon,fout,"longitude")
    ncwrite(lat,fout,"latitude")
    ncwrite(hpred,fout,"hret")
    ncsync()
end

if (length(size(1)) == 1)
    nccreate(fout,"hret",
             "time",time,Dict("note" => "copied from input file"),
             atts=Dict("longname" => "predicted internal tide sea level anomaly",
                       "units" => "meter"),
             gatts=Dict("creator" => "Software written by Edward D. Zaron, 2022-08-10.",
                        "contact" => "edward.d.zaron@oregonstate.edu",
                        "cidvec"  => cidvec,
                        "project" => "NASA SWOT Science Team",
                        "history" => "Based on https://ingria.ceoas.oregonstate.edu/fossil/HA checkout dde20225c3ed7e3364ffc282ff53b4832980c930",
                        "version" => "1.0",
                        "source"  => "https://ingria.ceoas.oregonstate.edu/fossil/SMCE checkout d43216f1d9981b2fd151404ecce38a9220a4f6ce",
                        "infile"  => infile,
                        "outfile" => fout,
                        "fhret"   => fhret
                        )
             )
    nccreate(fout,"longitude",
             "time",
             atts=Dict("longname" => "longitude, copied from infile")
             )
    nccreate(fout,"latitude",
             "time",
             atts=Dict("longname" => "latitude, copied from infile")
             )

    ncwrite(lon,fout,"longitude")
    ncwrite(lat,fout,"latitude")
    ncwrite(hpred,fout,"hret")
    ncsync()
end
