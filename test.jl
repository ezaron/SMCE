#
# This program generates baroclinic tide predictions for lat/lon/time
# coordinates provided in netcdf files. The predictions consist of the
# baroclinic sea level anomaly, eastward surface current, and westward
# surface current. The HRET14 empirical tide model is the basis for
# these predictions. This model consists of gridded maps of sea level
# anomaly for the M2, S2, K1, O1, and N2 tides. The surface currents
# are derived from the sea level anomaly using a linear momentum balance.
#
# Sample data files are provided under the ./SWOTdata/, ./RADSdata/, and
# ./GDPdata directories.
#
# Usage:
#    julia test.jl PATHSPEC
#    (0) The parser fails if PATHSPEC begins with ".". Use full path instead.
#    (1) If PATHSPEC ends in .nc (i.e., if it is the full path to a netcdf file)
#        then we compute baroclinic tide SLA predictions and write them in a file 
#        called PATHSPEC_hret.nc.
#    (2) If PATHSPEC ends in .txt then we treat the contents of PATHSPEC
#        as a list of files and generate predictions, as above.
#    (3) If PATHSPEC ends in "/" then we compute tide predictions for all
#        the files with names ending in .nc in this directory, as above.
#
#    This program attempts to read and parse your input files to extract
#    the latitude, longitude, and time at which predictions are requested.
#    It has enough logic to handle cases, such as in the simulated SWOT
#    data files, where the two-dimensional (latitude,longitude) arrays
#    must be matched with a one-dimensional time array.
#
#    The most important thing you must check, as a user, is that the
#    reference time for the time variable has been parsed correctly.
#    The start and end times in YYYY-MM-DD HH:MM:SS format are output
#    to the console and also included on the figures output by this
#    program.
#
# edward.d.zaron@oregonstate.edu
#

# Load the needed components:
include("edznc.jl")
try
    HAmod.isloaded()
catch
    include("./HAmod.jl")
    using .HAmod
end
using Interpolations
#using ColorSchemes
using CairoMakie
using Printf
using Base.GC
using Statistics
if (gethostname() == "ebi")
    using ElectronDisplay
end

# Open window and draw figures or not:
#FIG=true
FIG=false
# Save hard-copies of figures or not:
#SFIG=true
SFIG=false
# Directory for saving figures, if they are saved:
DEST="./Figures/"
if (SFIG) run(`mkdir -p $DEST`) ; end

function main(pathspec)
    
#    println("pathspec = ",pathspec)
    if (string(pathspec[end]) == "/")
        infiles = readdir(pathspec)
        scat(x::String,y::String)=x * y
        infiles = scat.(pathspec,infiles)
    elseif (string(pathspec[end-2:end]) == ".nc")
        infiles = [ pathspec ]
    elseif (string(pathspec[end-4:end]) == ".txt")
        infiles = readlines(pathspec)
    else
        println("Cannot parse commandline argument. See comments at the top of ./test.jl.")
        exit()
    end

    # List the tidal frequencies you want to predict. This is largest set:
    #cidvec = ["M2", "S2", "K1", "O1", "MA2", "MB2"] # HRET8.1
    cidvec = ["M2", "S2", "K1", "O1", "N2"] # HRET14
#    cidvec = ["M2"]
#    cidvec = ["S2"]
#    cidvec = ["K1"]
#    cidvec = ["O1"]
#    cidvec = ["M2", "S2", "K1", "O1"] # HRET8.1

lonVars = ["longitude", "lon", "LON"]
latVars = ["latitude", "lat", "LAT"]
timeVars = ["time", "TIME"]

for infile=infiles
    println("Trying to process infile = ",infile)
    m = match(r".nc",infile)
    if isnothing(m) continue ; end
    m = match(r"_hret.nc",infile)
    if !isnothing(m) continue ; end

lon,lonV = ncvarget(infile,lonVars)
lat,latV = ncvarget(infile,latVars)
time,timeV = ncvarget(infile,timeVars)
tunits = ncgetatt(infile,timeV,"units")

if (size(lon) != size(lat))
    println("ERROR: Sizes on longitude and latitude arrays differ.")
end

# Try to figure out the units and reference time from tunits:
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
    dayfrac = 0.0
else
    try
        hr = parse(Float64,hms[1])
        mi = parse(Float64,hms[2])
        se = parse(Float64,hms[3])
        dayfrac = (hr + (mi + se/60.0)/60.0)/24.0
    catch
        println("ERROR: Could not parse the reference time of the input file: bad HH:MM:SS.")
    end
end
try
    yy = parse(Float64,ymd[1])
    mm = parse(Float64,ymd[2])
    dd = parse(Float64,ymd[3])
catch
    println("ERROR: Could not parse the reference time of the input file: bad or missing YYYY-MM-DD.")
end
jd0 = ymd2jd(yy,mm,dd + dayfrac)
# Convert time to decimal Julian Date which begins at midnight rather than noon:
time = jd0 .+ time
indg = find(!isnan,time)
indgg = find(!isinf,time[indg])
indg = indg[indgg]
# Sanity check:
yy,mm,dd,hr,mi,se = jd2ymdhms(minimum(time[indg]))
start_str=@sprintf("First measurement time = %04i-%02i-%02i %02i:%02i:%05.2f",
                   yy,mm,dd,hr,mi,se)
yy,mm,dd,hr,mi,se = jd2ymdhms(maximum(time[indg]))
stop_str=@sprintf("Last measurement time = %04i-%02i-%02i %02i:%02i:%05.2f",
                  yy,mm,dd,hr,mi,se)
println(start_str)
println(stop_str)
# This would be a good place to plot the input mesh and indicate
# the start and end times on the plot.

rank=length(size(lon))

if (FIG)
f = Figure(resolution=(600,300))
ax = Axis(f[1, 1],
          xlabel="longitude",
          ylabel="latitude",
          title="data locations")
if ( (rank == 2) & (length(time) < 10000) )
    #lines!(ax, lon[1,:], lat[1,:])
    #lines!(ax, lon[end,:], lat[end,:])
    #lines!(ax, lon[:,1], lat[:,1])
    #lines!(ax, lon[:,end], lat[:,end])
    xx = [lon[:,1] ; lon[end,:] ; reverse(lon[:,end]) ; reverse(lon[1,:])]
    yy = [lat[:,1] ; lat[end,:] ; reverse(lat[:,end]) ; reverse(lat[1,:])]
    lines!(ax,xx,yy,color=:black)
    if (prod(size(lon)) < 1000)
        n,m=size(lon)
        for j=1:m
            scatter!(ax,lon[:,j],lat[:,j])
        end
    end
elseif ( (rank == 1) & (length(time) < 10000) )
    lines!(ax, lon, lat)
elseif (rank == 1)
    lines!(ax, lon[1:100:end], lat[1:100:end])    
else
    println("Not sure how to show (lat,lon). Skipping graphics.")
end
text!(ax,start_str,position=(10,60),align=(:left,:center))
text!(ax,stop_str,position=(10,30),align=(:left,:center))
plotcoast(ax)
    display(f)
end
if (SFIG)
    tmp = split(infile,"/")
    tmp = tmp[end]
    tmp = split(tmp,".")
    tmp[end] = ".pdf"
    fout = reduce(*,tmp)
    println("Saving graphics at ",fout)
    save(DEST * fout,f)
end

# Note how we force time to look like a 1-d vector even if it is
# a matrix:
indu = find( x -> (!isnan(x)) & (!isinf(x)),time)
    iid,F = FMAT(time[indu],cidvec,1,0) ; println("nodal modulation is turned ON (most accurate, but worse agreement with checkval)")
    # iid,F = FMAT(time[indu],cidvec,0,0) ;  println("nodal modulation is turned OFF (less accurate, but better agreement with checkval)")

    #        fhret = "/home/ezaron/FFTest/HRET8.1/HRET_v8.1.nc"
    fhret   = "./HRET14/HRET14_SSH_d1.nc"
    fhretuv = nothing
    # Uncomment the next line if you would like to predict tidal currents in addition to SSH:
    fhretuv = "./HRET14/HRET14_UV_d1.nc"

hlon = ncvarget(fhret,"lon")
hlat = ncvarget(fhret,"lat")
hpred = zeros(size(lon))
upred = zeros(size(lon))
vpred = zeros(size(lon))
# Unwrap longitude into [0,360] if needed:
lonx = copy(lon)
indw = find( x -> x < 0.0,lonx)
lonx[indw] = lon[indw] .+ 360.0
# Load data for each component frequency one at a time:    
for cid = cidvec
    println("Working on cid = ",cid)
    ic = cid * "c"
    is = cid * "s"
    iic = iid[cid * "c"]
    iis = iid[cid * "s"]
    # It is possible that loading and building all the interpolators for
    # ssh, u, and v will consume too much memory. Monitor usage during tests!
    hc = ncvarget(fhret,ic)
    hs = ncvarget(fhret,is)
    ihc = interpolate((hlon,hlat),hc,(Gridded(Linear()),Gridded(Linear())))
    ihs = interpolate((hlon,hlat),hs,(Gridded(Linear()),Gridded(Linear())))
    ehc = extrapolate(ihc,(Periodic(),Flat()))
    ehs = extrapolate(ihs,(Periodic(),Flat()))
    if !isnothing(fhretuv)
        uc = ncvarget(fhretuv,"u" * ic)
        us = ncvarget(fhretuv,"u" * is)
        iuc = interpolate((hlon,hlat),uc,(Gridded(Linear()),Gridded(Linear())))
        ius = interpolate((hlon,hlat),us,(Gridded(Linear()),Gridded(Linear())))
        euc = extrapolate(iuc,(Periodic(),Flat()))
        eus = extrapolate(ius,(Periodic(),Flat()))
        vc = ncvarget(fhretuv,"v" * ic)
        vs = ncvarget(fhretuv,"v" * is)
        ivc = interpolate((hlon,hlat),vc,(Gridded(Linear()),Gridded(Linear())))
        ivs = interpolate((hlon,hlat),vs,(Gridded(Linear()),Gridded(Linear())))
        evc = extrapolate(ivc,(Periodic(),Flat()))
        evs = extrapolate(ivs,(Periodic(),Flat()))
    end
    if (size(lon) == size(time))
        xc = ehc.(lonx[indu],lat[indu])
        xs = ehs.(lonx[indu],lat[indu])
        hpred[indu] = hpred[indu] + F[:,iic].*xc + F[:,iis].*xs
        if !isnothing(fhretuv)
            xc = euc.(lonx[indu],lat[indu])
            xs = eus.(lonx[indu],lat[indu])
            upred[indu] = upred[indu] + F[:,iic].*xc + F[:,iis].*xs
            xc = evc.(lonx[indu],lat[indu])
            xs = evs.(lonx[indu],lat[indu])
            vpred[indu] = vpred[indu] + F[:,iic].*xc + F[:,iis].*xs
        end
        continue
    end
    # Different cases when size(time) != size(lat)
    if (length(size(time)) == 1)
        na,nb=size(lat)
        for kk=length(indu)
            k=indu[kk]
            if (length(time) == na)
                lon1 = view(lonx,k,:)
                lat1 = view(lat,k,:)
                h1   = view(hpred,k,:)
                if !isnothing(fhretuv)
                    u1   = view(upred,k,:)
                    v1   = view(vpred,k,:)
                end
            else
                lon1 = view(lonx,:,k)
                lat1 = view(lat,:,k)
                h1   = view(hpred,:,k)
                if !isnothing(fhretuv)
                    u1   = view(upred,:,k)
                    v1   = view(vpred,:,k)
                end
            end
            xc = ehc.(lon1,lat1)
            xs = ehs.(lon1,lat1)
            h1[:] = h1[:] + F[kk,iic].*xc + F[kk,iis].*xs
            if !isnothing(fhretuv)
                xc = euc.(lon1,lat1)
                xs = eus.(lon1,lat1)
                u1[:] = u1[:] + F[kk,iic].*xc + F[kk,iis].*xs
                xc = evc.(lon1,lat1)
                xs = evs.(lon1,lat1)
                v1[:] = v1[:] + F[kk,iic].*xc + F[kk,iis].*xs
            end
        end
    end
end
GC.gc()
println("DONE!")

println("Writing results to file.")

    if (1 == 0)
        # Skip the validation test for HRET14 since there are not corresponding values precomputed in the files.
        checkval = ncvarget(infile,"internal_tide_hret")
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
        println("MAXIMUM error compared to checkval [mm] = ",1e3*maximum(abs.(checkval - hpred)))
        println("This should be less than 0.1 mm when nodal modulation is turned off.")
        println("It could be more than a 1 mm when the nodal modulation is turned on.")
    end
    
    tmp = split(infile,".")
    if (length(cidvec) == 1)
        tmp[end] = @sprintf("_hret_%s.nc",cidvec[1])
    else
        tmp[end] = "_hret.nc"
    end
fout = reduce(*,tmp)
run(`rm -f $fout`)

# Hmmm. For the SWOT examples, the output file is considerably larger than the input file, even though it has
# many fewer variables. Seems like it would be nice to copy the compression and stuff from the input file,
# but maybe it would be easiest to write the new fields into the infile.

if (length(size(lon)) == 2)
    num_pixels,num_lines = size(lon)

    println("Attempting to create this file: ",fout)
    if isnothing(fhretuv)
        fhretuv = "nothing"
    end
    nccreate(fout,"hret",
             "num_pixels",num_pixels,
             "num_lines",num_lines,
             atts=Dict("long_name" => "predicted internal tide sea level anomaly",
                       "units" => "meter"),
             gatts=Dict("creator" => "Software written by Edward D. Zaron, 2023-05-26.",
                        "contact" => "edward.d.zaron@oregonstate.edu",
                        "cidvec"  => cidvec,
                        "project" => "NASA SWOT Science Team",
                        "history" => "Based on https://ingria.ceoas.oregonstate.edu/fossil/HA checkout dde20225c3ed7e3364ffc282ff53b4832980c930",
                        "version" => "1.0",
                        "source"  => "https://ingria.ceoas.oregonstate.edu/fossil/SMCE",
                        "infile"  => infile,
                        "outfile" => fout,
                        "fhret"   => fhret,
                        "fhretuv" => fhretuv
                        )
             )
    if (fhretuv == "nothing")
        fhretuv = nothing
    end

    if !isnothing(fhretuv)
        nccreate(fout,"uhret",
                 "num_pixels",
                 "num_lines",
                 atts=Dict("long_name" => "predicted internal tide eastward surface current",
                           "units" => "meter/second"))
        nccreate(fout,"vhret",
                 "num_pixels",
                 "num_lines",
                 atts=Dict("long_name" => "predicted internal tide northward surface current",
                           "units" => "meter/second"))
    end
    nccreate(fout,lonV,
             "num_pixels",
             "num_lines",
             atts=Dict("long_name" => "longitude, copied from infile")
             )
    nccreate(fout,latV,
             "num_pixels",
             "num_lines",
             atts=Dict("long_name" => "latitude, copied from infile")
             )
    
    ncwrite(lon,fout,lonV)
    ncwrite(lat,fout,latV)
    ncwrite(hpred,fout,"hret")
    if !isnothing(fhretuv)
        ncwrite(upred,fout,"uhret")
        ncwrite(vpred,fout,"vhret")
    end
    ncsync()
end

if (length(size(lon)) == 1)
    nccreate(fout,"hret",
             timeV,time,Dict("note" => "copied from input file"),
             atts=Dict("long_name" => "predicted internal tide sea level anomaly",
                       "units" => "meter"),
             gatts=Dict("creator" => "Software written by Edward D. Zaron, 2023-05-26.",
                        "contact" => "edward.d.zaron@oregonstate.edu",
                        "cidvec"  => cidvec,
                        "project" => "NASA SWOT Science Team",
                        "history" => "Based on https://ingria.ceoas.oregonstate.edu/fossil/HA checkout 1f4512dbcc1140befa49a475d221344a55dc83b0581e5def8006e1ef0c197379",
                        "version" => "1.0",
                        "source"  => "https://ingria.ceoas.oregonstate.edu/fossil/SMCE",
                        "infile"  => infile,
                        "outfile" => fout,
                        "fhret"   => fhret
                        )
             )

    if !isnothing(fhretuv)
        nccreate(fout,"uhret",
                 timeV,
                 atts=Dict("long_name" => "predicted internal tide eastward surface current",
                           "units" => "meter/second"))
        nccreate(fout,"vhret",
                 timeV,
                 atts=Dict("long_name" => "predicted internal tide northward surface current",
                           "units" => "meter/second"))
    end
    nccreate(fout,lonV,
             timeV,
             atts=Dict("long_name" => "longitude, copied from infile")
             )
    nccreate(fout,latV,
             timeV,
             atts=Dict("long_name" => "latitude, copied from infile")
             )
    
    # It would be better if I copy the same name as used in the input file:
    ncwrite(lon,fout,lonV)
    ncwrite(lat,fout,latV)
    ncwrite(hpred,fout,"hret")
    if !isnothing(fhretuv)
        ncwrite(upred,fout,"uhret")
        ncwrite(vpred,fout,"vhret")
    end
    ncsync()
end
GC.gc()
    
end # infile = infiles

end # main()

# Determine if we are running from the Terminal or if we are
# running inside a Notebook
#println("ARGS = ",ARGS)
if (length(ARGS) == 0)
    # We are running inside a Terminal, and no commandline arguments were specified, so we use a sample file:
    pathspec = "/home/jovyan/DEMO_FILES/SWOT_L2_LR_SSH_Expert_001_004_20140412T143420_20140412T152546_DG10_01.nc"   
elseif (length(ARGS) == 1)
    m = match(r"\.json",ARGS[1])
    if (!isnothing(m))
        # We are running inside a Notebook, no commandline arguments are used, we just set the filename as above:
        pathspec = "/home/jovyan/DEMO_FILES/SWOT_L2_LR_SSH_Expert_001_004_20140412T143420_20140412T152546_DG10_01.nc"   
    else
        # We are running inside a Terminal, and a commandline argument was specified:
        pathspec = ARGS[1]
    end
else
    println("Could not determine if we are running in a Terminal or Notebook. See comments at the top of ./test.jl.")
    exit()
end

main(pathspec)
#main("/home/jovyan/DEMO_FILES/")

println("Output files have been written. This script is done!")
if (length(ARGS) > 1)
    exit()
end
