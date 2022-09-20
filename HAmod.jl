
module HAmod
# To use these functions:
#    include("/home/ezaron/NASA-Tides/iho/iho_func.jl")
#    using .HAmod
#
# Remember: You can access the non-exported functions and data
#           structures using their full names, i.e., to get the
#           big list of Doodson numbers, use "HAmod.CID".

# This list of modules could be trimmed for use in the SMCE script if we
# need to speed up the startup time of the script:
using DelimitedFiles
using LinearAlgebra
using Statistics
using Printf
using DSP
using FFTW

try
    function gethostname()
        hostname = Vector{UInt8}(undef, 128)
        ccall((:gethostname, "libc"), Int32,
              (Ptr{UInt8}, Csize_t),
              hostname, sizeof(hostname))
        hostname[end] = 0; # ensure null-termination
        return unsafe_string(pointer(hostname))
    end
    global myhost = gethostname()
catch
    # println("iho_func.jl failed to define gethostname1!")
end

try
    function gethostname()
        hostname = Vector{UInt8}(undef, 128)
        ccall((:gethostname, "/lib/x86_64-linux-gnu/libc-2.27.so"), Int32,
              (Ptr{UInt8}, Csize_t),
              hostname, sizeof(hostname))
        hostname[end] = 0; # ensure null-termination
        return unsafe_string(pointer(hostname))
    end
    global myhost = gethostname()
catch
    # println("iho_func.jl failed to define gethostname2!")
end

function isloaded()
    return true
end

myhost = gethostname()

function loadDoodson(myhost,HAroot)
    try
        X = readdlm(HAroot * "iLongPeriod.txt",',')
        Y = readdlm(HAroot * "iDiurnal.txt",',')
        Z = readdlm(HAroot * "iSemiDiurnal.txt",',')
        P = readdlm(HAroot * "iThirdDiurnal.txt",',')
        Q = readdlm(HAroot * "iFourthDiurnal.txt",',')
        R = readdlm(HAroot * "iSixthDiurnal.txt",',')
        S = readdlm(HAroot * "iEighthDiurnal.txt",',')
        return X,Y,Z,P,Q,R,S
    catch
#        println("iho_func.jl failed to load Doodson!")
        return nothing
    end
end

#HAroots = ["/home/ezaron/NASA-Tides/iho/",
#           "/home/ezaron/HA/",
#           "/home/ezaron/Desktop/NASA-Tides/iho/",
#           "./" ]
# For SMCE:
HAroots = ["./"]

CID = nothing
for HAroot in HAroots
    X = loadDoodson(myhost,HAroot)
    if (X == nothing)
        println("HAmod did not find HA data files here: ",HAroot)
    else
        X,Y,Z,P,Q,R,S = X
        println("HAmod found Doodson number definitions here: ",HAroot)
        global CID = [X; Y; Z; P; Q; R; S]
        break
    end
end
if (isnothing(CID))
    println("HAmod is using a minimal set of Doodson numbers defined inline.")
    CID = [
        "Sa"  , "0.041 069", "0 565 555", "Z ZAZ ZZZ", "z",
        "Xa"  , "0.041 069", "---------", "Z ZAB ZZZ", "z",
        "Ssa" , "0.082 137", "0 575 555", "Z ZBZ ZZZ", "z",
        "Xsa" , "0.082 137", "---------", "Z ZBB ZZZ", "z",
        "Xm"  , "0.471 521", "---------", "Z AVA ZZZ", "x",
        "MSm" , "0.471 521", "0 636 555", "Z AXA ZZZ", "x",
        "Mm"  , "0.544 375", "0 654 555", "Z AZY ZZZ", "y",
        "Xf"  , "1.015 896", "---------", "Z BUZ ZZZ", "b",
        "MSf" , "1.015 896", "0 735 555", "Z BXZ ZZZ", "b",
        "X2"  , "2.031 792", "---------", "Z DTZ ZZZ", "c",
        "2SM" , "2.031 792", "0 915 555", "Z DVZ ZZZ", "c",
        "O1"  ,"13.943 036", "1 455 554", "A YZZ ZZY", "y",
        "K1"  ,"15.041 069", "1 655 556", "A AZZ ZZA", "y",
        "MA2" ,"28.943 036", "2 545 555", "B ZYZ ZZZ", "m",
        "M2"  ,"28.984 104", "2 555 555", "B ZZZ ZZZ", "y",
        "MB2" ,"29.025 173", "2 565 555", "B ZAZ ZZZ", "m",
        "S2"  ,"30.000 000", "2 735 555", "B BXZ ZZZ", "z" ]
    CID = permutedims(reshape(CID,(5,6)),(2,1))
end

NM   = 1
SPD  = 2
XDO  = 3
XDOP = 4
NODE = 5

XDOPm = Dict(
    "N" =>-12. ,
    "O" =>-11. ,
    "P" =>-10. ,
    "Q" => -9. ,
    "R" => -8. ,
    "S" => -7. ,
    "T" => -6. ,
    "U" => -5. ,
    "V" => -4. ,
    "W" => -3. ,
    "X" => -2. ,
    "Y" => -1. ,
    "Z" =>  0. ,
    "A" =>  1. ,
    "B" =>  2. ,
    "C" =>  3. ,
    "D" =>  4. ,
    "E" =>  5. ,
    "F" =>  6. ,
    "G" =>  7. ,
    "H" =>  8. ,
    "I" =>  9. ,
    "J" => 10. ,
    "K" => 11. ,
    "L" => 12. ,
    "M" => 13.
)

nc,nx = size(CID)
nm  = Array{String,1}(undef,nc) ;
spd = Array{Float64,1}(undef,nc) ;
xdo = Array{Float64,2}(undef,nc,7) ;
node= Array{String,1}(undef,nc) ;
for i=1:nc
    nm[i] = CID[i,NM]
    str = CID[i,SPD]
    # I don't trust the speeds transcribed in the tables.
    # For better precision, compute the speeds from the
    # Doodson numbers, below.
    # spd[i] = parse(Float64,join(map(strip,split(str))))
    
    str = CID[i,XDO]
    if (str[1] != '-')
        sxdo = join(strip.(split(str,"")))
        xdo1 = [ parse(Float64,string(sxdo[1])),
                 parse(Float64,string(sxdo[2])) .- 5.0,
                 parse(Float64,string(sxdo[3])) .- 5.0,
                 parse(Float64,string(sxdo[4])) .- 5.0,
                 parse(Float64,string(sxdo[5])) .- 5.0,
                 parse(Float64,string(sxdo[6])) .- 5.0,
                 parse(Float64,string(sxdo[7])) .- 5.0
                 ]
    else
        str = CID[i,XDOP]
        sxdo = join(strip.(split(str,"")))
        xdo1 = Array{Float64,1}(undef,7)
        xdo1 = [XDOPm[string(sxdo[j])] for j=1:7]
    end
    xdo[i,:] = xdo1 ;
    node[i] = CID[i,NODE]
end

# Fundamental angular velocities, degrees per hour, from Simon
# Table 5.4, page 113:
av = [14.492052121,
       0.549016518,
       0.041068640,
       0.004641822,
       0.002206407,
       0.000001962]
# Compute the speed, degrees per hour, instead of reading from
# the tables.
spd = xdo[:,1:6]*av

#
# Use this dict to access cid's in the above array, e.g.,
#  cid = id["S2"]
#
#id = Dict{String,Int}()
id = Dict{String,Int}(" " => 0)
for i=collect(1:nc)
    #    id = merge(id,Dict([(nm[i],i)]))
    merge!(id,Dict([(nm[i],i)]))
end

dtr = pi/180.

#=================================================================================#

# Convenience method: grab period, in days, given cid
export period
function period(cid)
    # cid is the name of the frequency, e.g., "K1"
    spdtide= spd[id[cid]] # degrees per hour
    return 1.0/spdtide*360.0/24. # [days]
end

export freq
function freq(cid)
    # cid is the name of the frequency, e.g., "K1"
    spdtide= spd[id[cid]] # degrees per hour
    return spdtide*pi/180.0/60.0 ; # [rad/sec]
end

export getxdo
function getxdo(cid)
    return convert.(Int,xdo[id[cid],:])
end

export getxdoodson
function getxdoodson(cid)
    dood = getxdo(cid)
    dood[2:end] = dood[2:end] .+ 5
    return dood
end

# Cartwright, Taylor, and Edden format:
export ctestr
function ctestr(cid)
    d = getxdo(cid)
    return @sprintf("%2i%2i%2i%2i%2i%2i%2i",d[1],d[2],d[3],d[4],d[5],d[6],d[7])
end

export doodstr
function doodstr(cid)
    d = getxdoodson(cid)
    ind1 = findall( x -> x == 10, d)
    ind2 = findall( x -> x == 11, d)
    dd = map(x -> @sprintf("%i",x),d)
    if (length(ind1) > 0)
        dd[ind1] .= "X"
    end        
    if (length(ind2) > 0)
        dd[ind2] .= "E"
    end
    return @sprintf("%1s%1s%1s.%1s%1s%1s,%1s",dd[1],dd[2],dd[3],dd[4],dd[5],dd[6],dd[7])
end

# Aliasing formulas:
# From the table, the angular speed is in units of deg per hour
#
export aliasspd1
function aliasspd1(cid,dtau)
    # cid is the name of the frequency, e.g., "K1"
    # dtau is the satellite repeat period in days
    spdsat = 360.0/(dtau*24.)
    spdtide= spd[id[cid]]
    N = convert(Int,floor( spdtide/spdsat ))
    return spdtide .- N*spdsat
end

export aliasspd
function aliasspd(cid,dtau)
    # cid is the name of the frequency, e.g., "K1"
    # dtau is the satellite repeat period in days
    spdtide= spd[id[cid]]
    dpha = dtau*24.0*spdtide
    N = convert(Int,floor( dpha/360. ))
    dpha = dpha .- N*360.
    if (dpha > 180.) dpha = dpha .- 360. ; end
    return dpha/(dtau*24.)
end

export aliasperiod
function aliasperiod(cid,dtau)
    # Given dtau in days, this returns
    # alias period in days.
    return abs.(360.0/aliasspd(cid,dtau)/24.)
end

export days_to_sep
function days_to_sep(cid1,cid2,dtau)
    t1 = aliasperiod(cid1,dtau)
    t2 = aliasperiod(cid2,dtau)
    return 1.0/abs.(1.0/t1 .- 1.0/t2)
end

export aliasf
function aliasf(freq,dt)
    # Freq is in units of cyc/time (not radians/time!)
    # Closest integer multiple:
    if (dt < 0.5./freq) return freq ; end
    nint = round(Int,freq*dt )
    # This is the usual unsigned alias:
    #falias = abs(freq - nint/dt)
    # This is the signed alias, used for interpreting rotary spectra with aliasing!
    falias = freq - nint/dt
    return falias
end
    
# Exact-repeat mission orbit periods [days], relative to a fixed position on earth:
# Uppercase are the main missions. Lower case refer to long-repeat mission phases,
# which are small perturbations of the reference missions. Not sure if these are big
# enough perturbations to be of consequence for tidal aliases.
DTAU = Dict("TXA" => 9.91564280 ,
            "j1c" => 16.9809,
            "j2c" => 8.003378378,
            "E2A" => 35.0,
            "3A"  => 27.0,
            "sab" => 35.00654,
            "GSA" => 23.06888431 ,
            "G1A" => 17.05057808 ,
            "C2A" => 368.239641 ,
            "SWT1" => 0.99349 ,
            "SWT" => 20.86455)

if (1 == 0)
    
    dtau = DTAU["SWT"]
    dtau = DTAU["TXA"]
    dtau = DTAU["G1A"]
    dtau = DTAU["E2A"]
    
    cidvec = ["S1" "Sa" "Ssa" "M2" "S2" "N2" "K2" "K1" "O1" "P1" "Q1"]
    for i=collect(2:length(cidvec))
        for j=collect(1:i-1)
            cid1 = cidvec[i]
            cid2 = cidvec[j]
            dt = days_to_sep(cid1,cid2,dtau)
            println(cid1," and ",cid2," take ",round(dt)," days to separate.")
        end
    end
end

#=================================================================================#
#=
 THE ASTRONOMICAL ARGUMENTS

 The tshpnp function returns the 6 astronomical arguments, as used by Doodson:
     tau : lunar hour angle in degrees, computed from solar time (UT) and s and h.
     s   : mean longitude of the Moon
     h   : mean longitude of the Sun
     p   : mean longitude of the lunar perigee
     N'  : mean longitude of the ascending lunar node
     p1  : mean longitude of the solar perigee
     c90 : one-half pi rad => 90 degrees

 The code is based on RDR's astro5.f code.
 Changes relative to RDR's code:
      Input is Julian Day, not Modified Julian Day.
      [TO-DO:
         My usage of Julian Day is non-standard (incorrect), but internally
         consistent. By definition, the Julian Day (an integer) and Julian
         Date (a real number) are referenced to noon (according to RDR's
         web page on the Modified Julian Date). And, by definition, the
         Modified Julian Date is referenced to midnight.
         Within this code, when Julian Date or Julian Day is an integer,
         it is referenced to noon; when the Julian date is a real number,
         it is referenced to midnight.]
      Return N'=-N rather than N.
      Return lunar hour angle, tau.
      Return 90deg vector for forming full product with xdo.

 I checked the values returned by comparison with Table 1 on p. 163 of
 Schureman. I also implemented a version based on the slightly different
 formulas in Seidelmann, which happen to be quoted earlier in Meeus
 text.
=#

export tshpnp
function tshpnp(TIME)

#= This code is based on R. D. Ray's astro5.f:
*
*---------------------------------------------------------------------
*  Computes the 5 basic astronomical mean longitudes  s, h, p, N, p'.
*
!!! --> EDZ change *  Note that N is N', i.e. N is increasing with time, unlike RDR's version
*
*  TIME is UTC in decimal Julian Day, unlike RDR's version which uses Modified Julian Day (MJD).
*  All longitudes returned in degrees.
*
*  R. D. Ray, NASA/GSFC   August 2003
*
*  Most of the formulae for mean longitudes are extracted from 
*  Jean Meeus, Astronomical Algorithms, 2nd ed., 1998.  
*  Page numbers below refer to this book.
*
*  Note: This routine uses TIME in UT and does not distinguish between
*    the subtle differences of UTC, UT1, etc.  This is more than adequate
*    for the calculation of these arguments, especially in tidal studies.
*
*  Revised 4/15/10 .- minor change to modulo function eliminates some extra code.
*---------------------------------------------------------------------
*
    =#
    CIRCLE=360.0e0
    nt = length(TIME)
    TSHPNP = Array{Float64,2}(undef,nt,7)
    
#    #*     Convert to Julian Day and to Ephemeris Time
#    #*     -------------------------------------------
#    TJD = TIME .+ 2400000.5e0
    TJD = TIME

    #*     UT (SOLAR) TIME IN HOURS
    #*     ------------------------
    jd0 = floor.(TJD .+ 0.5) .- 0.5
    TH = (TJD .- jd0)*24.
    
    #*     Compute time argument in centuries relative to J2000
    #*     ----------------------------------------------------
    T = ( TJD .- 2451545.e0 )/36525.e0
    
    for i=1:nt
        #*     mean longitude of moon (p.338)
        #*     ------------------------------
        TSHPNP[i,2] = (((-1.53388e-8*T[i] .+ 1.855835e-6)*T[i] .- 1.5786e-3)*T[i] .+ 481267.88123421e0)*T[i] .+ 218.3164477e0
        
        #*     mean elongation of moon (p.338)
        #*     -------------------------------
        D = (((-8.8445e-9*T[i] .+ 1.83195e-6)*T[i] .- 1.8819e-3)*T[i] .+ 445267.1114034e0)*T[i] .+ 297.8501921e0
        
        #*     mean longitude of sun
        #*     ---------------------
        TSHPNP[i,3] = TSHPNP[i,2] .- D
        
        #*     mean longitude of lunar perigee (p.343)
        #*     ---------------------------------------
        TSHPNP[i,4] = ((-1.249172e-5*T[i] .- 1.032e-2)*T[i] .+ 4069.0137287e0)*T[i] .+ 83.3532465e0

        #*     mean longitude of ascending lunar node (p.144)
        #*     ----------------------------------------------
        # This is N:
        TSHPNP[i,5] = ((2.22222e-6*T[i] .+ 2.0708e-3)*T[i] .- 1934.136261e0)*T[i] .+ 125.04452e0
        # We use N':
        TSHPNP[i,5] = .- TSHPNP[i,5]

        #*     mean longitude of solar perigee (Simon et al., 1994)  = s .- D .- l'
        #*     ----------------------------------------------------
        TSHPNP[i,6] = 282.93734e0 .+ 1.71953e0 * T[i]

        #*     lunar hour angle: TAU = 15deg/hr*TH .+ h .- s
        #*     TH is time in hours since 2000-01-01-12:00:00, as per Simon p.112.
        #*     ----------------------------------------------------
        TSHPNP[i,1] = 15.e0*TH[i] .+ TSHPNP[i,3] .- TSHPNP[i,2]

        #*     90 degrees, a constant
        #*     ----------------------
        TSHPNP[i,7] = 90.
        
        for j=1:6
            TSHPNP[i,j] = mod( TSHPNP[i,j], CIRCLE )
        end
    end
    return TSHPNP
    # Check values:
    # vs. Meeus p.342:
    # ymd2jd(1992,4,12) = 2448725 at noon ==> jd = 2.4487245e6 at 0hr
    # From this routine, s = -36945.70981828064
    # or s = 134.29018171936332
    # vs. Meeus value  L' = 134.290182
    
    # Noon on 2000-01-01 is jd=2451545, or
    # mjd = 2451545. .-  2400000.5e0 = 51544.5
    # Thus, midnight of 2000-01-01 is mjd=51544.
    # which shall be the check date for the angles.
    # tshpnp(51544.) =
    # tau      s        h        p        N'       p1
    # 68.2452  211.728  279.973  83.2975  234.929  282.937
    # vs. Schureman, p. 163:
    #          s        h        p        N=125.069 p1
    #          211.744  279.973  83.294   234.931  282.940
    #
    # For the date 1900-01-01, 0 hour:
    # tshpnp( ymd2jd(1900.,1.,1.) .- 2400000.5 ) =
    # 3.16772  277.022  280.190  334.385  100.844  281.218
    #          277.026  280.190  334.384  100.844  281.221 
    #
    # 2020-01-10:
    #    Hmmm. The difference in angles, for s, for example, 0.024
    #    is 0.024/360 = 6.7e-5, is almost the same as the relative
    #    fractional difference between TPXO9a and FES14 in deep
    #    water. This is apparently the level of precision of Schureman's
    #    astronomical arguments.
    #    See H.-G. Wenzel (1997) "Tide-Generating Potential for the Earth" in
    #    Tidal Phenomena (H. Wilhelm, W. Zurn, and H.-G. Wenzel, eds.),
    #    Lecture Notes in Earth Sciences, Vol. 66, Springer, Berlin, pages 9 - 26, 398pp.
    #    This article contains a concise tabulation of astronomical argument polynomials,
    #    on p. 19, following Simon et al (1994).
    #    Simon, J. L., P. Bretagnon, J, Chapront, M. Chapront-Touz\'e, G. Francou
    #    and J. Laskar (1994) "Numerical Expressions for precession formulae and mean
    #    elements for the moon and planets" Astron. Astrophys., 282: 663--683.
    #
end

export tshpnpSimon
function tshpnpSimon(TIME)

#= This code is based on Wenzel 1997, Table 4:
*
*---------------------------------------------------------------------
*  Computes the 5 basic astronomical mean longitudes  s, h, p, N, p'.
*
*  Note that N is N' (negative mean longitude of lunar ascending node), 
*  i.e. N is increasing with time, unlike RDR's version
*
*  TIME is UTC in decimal Julian Day, unlike RDR's version which uses Modified Julian Day (MJD).
*  All longitudes returned in degrees.
*
*  Note: This routine uses TIME in UT and does not distinguish between
*    the subtle differences of UTC, UT1, etc.  This is more than adequate
*    for the calculation of these arguments, especially in tidal studies.
*
*---------------------------------------------------------------------
*
* CHECK: Values are very close to Richard's: Largest error is less than 1 part per 10^6 (for s):
*
* julia> map( x -> @sprintf("%e",x) , ( tshpnp(ymd2jd(2020.0,1.0,1.0)) .-  tshpnpSimon(ymd2jd(2020.0,1.0,1.0)) )./360 )
* 1×7 Array{String,2}:
* "4.602395e-08"  "-5.409065e-07"  "-4.948811e-07"  "6.149315e-08"  "9.724262e-08"  "-1.329808e-08"  "0.000000e+00"
    =#
    CIRCLE=360.0e0
    nt = length(TIME)
    TSHPNP = Array{Float64,2}(undef,nt,7)
    
#    #*     Convert to Julian Day and to Ephemeris Time
#    #*     -------------------------------------------
#    TJD = TIME .+ 2400000.5e0
    TJD = TIME

    #*     UT (SOLAR) TIME IN HOURS
    #*     ------------------------
    jd0 = floor.(TJD .+ 0.5) .- 0.5
    TH = (TJD .- jd0)*24.
    
    #*     Compute time argument in 10*centuries (millenia) relative to J2000
    #*     ----------------------------------------------------
    T = ( TJD .- 2451545.0 )/365250.0
    
    for i=1:nt
        #*     mean longitude of moon, s
        #*     -------------------------
        TSHPNP[i,2] = (((-0.00015355*T[i] .+ 0.00185140)*T[i] .- 0.14663889)*T[i] .+ 4812678.81195750)*T[i] .+ 218.31664562999
        
        #*     mean longitude of sun, h
        #*     ------------------------
        TSHPNP[i,3] = (((-0.00006532*T[i] .+ 0.00002000)*T[i] .+ 0.03032222)*T[i] .+ 360007.69748806)*T[i] .+ 280.46645016002
        
        #*     mean longitude of lunar perigee, p
        #*     ----------------------------------
        TSHPNP[i,4] = ((( 0.00052655*T[i] .- 0.01249168)*T[i] .- 1.03217222)*T[i] .+ 40690.13635250)*T[i] .+ 83.35324311998

        #*     NEGATIVE mean longitude of ascending lunar node, N'
        #*     ---------------------------------------------------
        TSHPNP[i,5] = ((( 0.00016501*T[i] .- 0.00213942)*T[i] .- 0.20756111)*T[i] .+ 19341.36261972)*T[i] .+ 234.95544499000

        #*     mean longitude of solar perigee, p_s
        #*     ----------------------------------------------------
        TSHPNP[i,6] = ((( -0.00003323*T[i] .- 0.00001776)*T[i] .+ 0.04568889)*T[i] .+ 17.19457667)*T[i] .+ 282.93734098001

        #*     lunar hour angle, tau
        #*     ---------------------
        # Seems like Richard's expression might be more accurate here since it avoids the big multiplication for the 15deg/hr
        # multiplied by 1000yr!
        TSHPNP[i,1] = (((  0.00008824*T[i] .- 0.00183140)*T[i] .+ 0.17696111)*T[i] .+ 127037328.88553056)*T[i] .+ 242.14980452999

        #*     90 degrees, a constant
        #*     ----------------------
        TSHPNP[i,7] = 90.
        
        for j=1:6
            TSHPNP[i,j] = mod( TSHPNP[i,j], CIRCLE )
        end
    end
    return TSHPNP
end


export shpnOTIS
function shpnOTIS(TIME)

#= This code is based on R. D. Ray's ASTROL routine, part of the OTIS Tidal Prediction Software.
*
*---------------------------------------------------------------------
*  Computes the 4 basic astronomical mean longitudes  s, h, p, N:
*
*  TIME is UTC in decimal Julian Day, unlike RDR's version which uses Modified Julian Day (MJD).
*  All longitudes returned in degrees.
*---------------------------------------------------------------------
*
    =#
    CIRCLE=360.0e0
    nt = length(TIME)
    SHPN = Array{Float64,2}(undef,nt,4)
    
    #*     Convert JD to OTIS time:
    #*     Compute time argument relative to J2000
    #*     ----------------------------------------------------
#    T = TIME .- 51544.4993 <------ Why is it not exactly 51544.5 ? Maybe Cartwright used a different base date and absorbed a coef.
    T = TIME .- ( 51544.4993 + 2400000.5 ) ; # Aha! This was RDR's conversion from UTC to Terrestrial Time. The latter must be used for
    # precise astronomical arguments. See the OTIS software distribution, which has the routine deltat for computing this time-dependent
    # offset. It is not necessary for the nodal corrections which ASTROL is used for within OTIS.
    
    # T = TIME .- ( 51544.5 + 2400000.5 ) ;
    # Noon on 2000-01-01 is jd=2451545
    # Noon on 2000-01-01 is mjd = jd - 2400000.5 = 51544.5
    # It appears that Richards "decimal MJD" is equal to MJD+0.5
    # Or, jd to "decimal MJD" = jd - 2400000
    for i=1:nt
        #*     mean longitude of moon
        #*     ----------------------
        SHPN[i,1] = 218.3164 .+ 13.17639648 * T[i]
        
        #*     mean longitude of sun
        #*     ---------------------
        SHPN[i,2] = 280.4661 .+  0.98564736 * T[i]
        
        #*     mean longitude of lunar perigee
        #*     -------------------------------
        SHPN[i,3] =  83.3535 .+  0.11140353 * T[i]

        #*     mean longitude of ascending lunar node (p.144)
        #*     ----------------------------------------------
        # This is N:
        SHPN[i,4] = 125.0445 .-  0.05295377 * T[i]
        # We use N':
        SHPN[i,4] = .- SHPN[i,4]
        
        for j=1:4
            SHPN[i,j] = mod( SHPN[i,j], CIRCLE )
        end
    end
    return SHPN
#=
CHECK:
julia>  map( x -> @sprintf("%e",x) , ( tshpnp(ymd2jd(2020.0,1.0,1.0))[2:5]' .-  shpnOTIS(ymd2jd(2020.0,1.0,1.0)) )./360 )
1×4 Array{String,2}:
 "-2.577349e-05"  "-1.446836e-06"  "-2.181203e-06"  "-4.932042e-07"
=#
end


#=
 TIME AND DATE CONVERSION FUNCTIONS.
 Be aware of the distinction between integer Julian Day (jd) of a given integer
 date as y-m-d, which corresponds to the jd at noon of the date
  vs.
 jd of a real-valued y-m-d date, which returns the jd at the requested value of
 d, which is referenced to UT=0 hr.
 The latter quantity is 0.5 days or 12 hours prior to noon.
=#

    export jd2ymdhms
    function jd2ymdhms(jd::Float64)
        Zr = trunc(jd .+ 0.5)
        Z = convert(Int,Zr)
        
        F = jd +0.5 .- Zr
        A = Z
        if (Z < 2299161)
        else
            ar = trunc( (Zr .- 1867216.25)/36524.25 )
            a  = convert(Int,ar)
            A = Z .+ 1 .+ a .- convert(Int,trunc(ar/4.0))
        end
        B = A+1524
        C = trunc( (B .- 122.1)/365.25 )
        D = trunc( 365.25*C )
        E = trunc( (B .- D)/30.6001 )
        da = B .- D .- trunc( 30.6001*E ) .+ F
        E = convert(Int,E)
        mo = E .- 1
        if E < 14
        else
            mo = E-13
        end

        C = convert(Int,C)
        yr = C-4716
        if mo > 2
        else
            yr = C-4715
        end
        sec = ( da .- floor(da) )*8.64e4
        hr = floor(Int,sec/3600. )
        sec = sec .- hr*3600.
        mi = floor(Int,sec/60. )
        sec = sec .- mi*60.
        da = floor(Int,da)
        return yr,mo,da,hr,mi,sec
    end

    export ymd2jd
function ymd2jd(y::Int, m::Int, d::Int)
    # y,m,d refers to year, month, day of the Gregorian Calendar.
    # The jd is the Julian day number, which equals the Julian day
    # corresponding to noon on the given date.
    jd::Int = (
    trunc((1461*(y .+ 4800 .+ trunc((m-14)/12)))/4)
    +trunc((367*(m-2-12*trunc((m-14)/12)))/12)
    -trunc((3*(trunc(y+4900+trunc((m-12)/12))/100))/4)
    +d
    -32075
    )
    return jd
    # Check value:
    # ymd2jd(1990,6,25) = 2448068, Seidelmann page 600
end

function ymd2jd(y::Float64,m::Float64,d::Float64)
    # This inputs ymd as a decimal.
    # d=1.0 corresponds to midnight on the given date.
    # Thus, this routine returns whole-number-valued jd
    # for the half-days, i.e.,
    # ymd2jd(1977.0,4.0,26.4) = 2443259.9
    # matches the check value on p. 59 of Meeus.
    if (m < 1. || m > 12.)
        return 0.
    end
    # Algorithm is from Meeus, p.61:
    if (m <= 2)
        m = m .+ 12.
        y = y .- 1.
    end
    A = floor(y/100.)
    B = 2 .- A .+ floor(A/4.0)
    JD = (
        floor(365.25*(y .+ 4716.))
        .+ floor(30.6001*(m+1.))
        .+ d
        .+ B
        .- 1524.5
    )
    return JD
end

# Vector version of above:
function ymd2jd(y::Array{Float64,1},
                m::Array{Float64,1},
                d::Array{Float64,1})
    # This inputs ymd as a decimal.
    # d=1.0 corresponds to midnight on the given date.
    # Thus, this routine returns whole-number-valued jd
    # for the half-days, i.e.,
    # ymd2jd(1977.,4.,26.4) = 2443259.9
    # matches the check value on p. 59 of Meeus.
    out = Array{Float64,1}(undef,length(y))
    for i=1:length(y)
        out[i] = ymd2jd(y[i],m[i],d[i])
    end
    return out
end

# Vector version of jd::Int function.
function ymd2jd(y::Array{Int,1},
                m::Array{Int,1},
                d::Array{Int,1})
    out = Array{Int,1}(undef,length(y))
    for i=1:length(y)
        out[i] = ymd2jd(y[i],m[i],d[i])
    end
    return out
end

export jd2ymd
function jd2ymd(jd::Int)
    # Check value:
    # jd2ymd(2448068) = (1990,6,25), Seidelman page 600
    l::Int = jd .+ 68569
    n::Int = floor((4*l)/146097)
    l = l .- floor((146097*n+3)/4)
    i::Int = floor((4000*(l+1))/1461001)
    l = l .- floor((1461*i)/4) .+ 31
    j::Int = floor((80*l)/2447)
    d::Int = l .- floor((2447*j)/80)
    l = floor(j/11)
    m::Int = j .+ 2 .- 12*l
    y::Int = 100*(n .- 49) .+ i .+ l
    return (y,m,d)
end

    export jd2gmst
function jd2gmst(jd::Float64) # [hr]
    # Convert Julian Date to Greenwich mean sidereal time [hours],
    # according to this US Naval Observatory web page:
    # http://aa.usno.navy.mil/faq/docs/GAST.php
    # Julian dates start at Greenwich mean noon.
    # Convert to time relative to 2000 Jan 1, 12h UT:
    d  = jd  .- 2451545.0
    gmst = 18.697374558 .+ 24.06570982441908*d
    gmst = mod(gmst,24.0)
    return gmst
    # Check, Meeus exaple 12.a, p. 88:
    # GMST on 1987-04-10-00:00:00
    # JD = 2446895.5
    # Meeus value: THETA0=13. .+ 10./60. .+ 46.3668/3600.
    # = 13.179546333333333 hr
    # jd2gmst(2446895.5)
    # = 13.179545921491808
end

    export jd2gmst_meeus
function jd2gmst_meeus(jd::Float64) # [hr] from UT 0 [hr]
    # From Meeus, p. 87, eqn 12.3
    jd0 = floor(jd+0.5) .- 0.5
    T = ( jd0 .- 2451545.e0 )/36525.e0
    THETA0 = 100.46061837 .+ 36000.770053608*T .+ 0.000387933*T*T .- T*T*T/38710000.
    tp = jd .- jd0
    theta0 = tp*360.0*1.00273790935
    gmst = mod(THETA0+theta0,360.0)
    return gmst
    # Check, Meeus exaple 12.a, p. 88:
    # jd2gmst_meeus(2446895.5)/360.*24.
    # = 13.179546340533063
    # vs. Meeus value: THETA0=13. .+ 10./60. .+ 46.3668/3600.
    # = 13.179546333333333
    # Check example 12.b, p. 89:
    # jd2gmst_meeus(2446895.5 .+ (19. .+ 21./60.)/24.)/15.
    # = 8.582524884214278
    # vs. Meeus value: 8. .+ 34./60. .+ 57.0896/3600.
    # = 8.582524888888889
end

function jd2gmst_hp(jd::Float64)
    # High-precision version, better than 0.1 second per century:
    # Convert Julian Date to Greenwich mean sidereal time [hours],
    # according to this US Naval Observatory web page:
    # http://aa.usno.navy.mil/faq/docs/GAST.php
    # Julian dates start at Greenwich mean noon.
    # Get the Julian date of the previous midnight:
    jd0 = floor(jd) .- 0.5
    h = (jd .- jd0)*24.0
    # Convert to time relative to 2000Jan 1, 12h UT:
    d  = jd  .- 2451545.0
    d0 = jd0 .- 2451545.0
    # Number of centuries since 2000:
    t  = d/36525.0
    gmst = 6.697374558 .+ 0.06570982441908*d0 .+ 1.00273790935*h .+ 0.000026*t*t
    gmst = mod(gmst,24.0)
    return gmst
    # Check value:
    # jd2gmst_hp(2.44689630625e6), same as used in Meeus, 12.b, p.89:
    # = 8.58252488633616
    # vs. Meeus value: 8. .+ 34./60. .+ 57.0896/3600.
    # = 8.582524888888889
end


function Nodex(cid1)
#    println("Nodex not implemented -- overtide.")
    return fnodal=Symbol("nodal"*"$cid1")
end
function Nodey(cid1)
    return fnodal=Symbol("nodal"*"$cid1")
end
function Nodez(cid1)
    #    return fnodal=Symbol("nodalz")
    return fnodal=Symbol("nodalz")
end
function Nodea(cid1)
    return fnodal=Symbol("nodalMm")
end
function Nodeb(cid1)
    return fnodal=Symbol("nodal"*"$cid1")
end
function Nodec(cid1)
    return fnodal=Symbol("nodal"*"$cid1")
end
function Noded(cid1)
    return fnodal=Symbol("nodalKQ1")
end
function Nodee(cid1)
    return fnodal=Symbol("nodalK2")
end
function Nodef(cid1)
    return fnodal=Symbol("nodalz")
end
function Nodeg(cid1)
    return fnodal=Symbol("nodal"*"$cid1")
end
function Nodej(cid1)
    return fnodal=Symbol("nodalJ1")
end
function Nodek(cid1)
    return fnodal=Symbol("nodalK1")
end
function Nodem(cid1)
    return fnodal=Symbol("nodalM2")
end
function Nodeo(cid1)
    return fnodal=Symbol("nodalO1")
end

# All the nodal functions return f,u pairs
function nodalz(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    xnt = length(cosp)
    nodeu = zeros(xnt,)
    nodef = ones(xnt,)
    return nodef, nodeu
end

function nodalMm(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

#    xnt = length(cosp)
#    nodeu = Array{Float64,1}(0.,xnt)
    nodeu = zeros(size(cosp))
    nodef = 1. .- 0.1311*cosN .+ 0.0538*cos2p .+ 0.0205*cos2pN
    return nodef, nodeu
end

function nodalMf(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = (-23.7*sinN .+ 2.7*sin2N .- 0.4*sin3N)*dtr
    nodef = 1.084 .+ 0.415*cosN .+ 0.039*cos2N
    return nodef,nodeu
end

function nodalO1(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = (10.80*sinN .- 1.34*sin2N .+ 0.19*sin3N)*dtr
    nodef = 1.0176 .+ 0.1871*cosN .- 0.0147*cos2N
    return nodef,nodeu
end

function nodalK1(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = (-8.86*sinN .+ 0.68*sin2N .- 0.07*sin3N)*dtr
    nodef = 1.0060 .+ 0.1150*cosN .- 0.0088*cos2N .+ 0.0006*cos3N
    return nodef,nodeu
end

function nodalJ1(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = (-12.94*sinN .+ 1.34*sin2N .- 0.19*sin3N)*dtr
    nodef = 1.1029 .+ 0.1676*cosN .- 0.0170*cos2N .+ 0.0016*cos3N
    return nodef,nodeu
end

function nodalM2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = (-2.14*sinN)*dtr
    nodef = 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N
    return nodef,nodeu
end

function nodalMSf(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodef =  nodef1
    nodeu = -nodeu1
    return nodef,nodeu
end

function nodalMSAf(pak)
    return nodalMSf(pak)
end
function nodalMSBf(pak)
    return nodalMSf(pak)
end

function nodalMSm(pak)
    # According to Constituents.pdf, MSm is misnamed, MNum.
    # Both the M2 and the Nu2 pieces have amp modulation like M2,
    # but they subtract to cancel phase modulation.
    nodef1,nodeu1 = nodalM2(pak)
    nodef = nodef1.^2
    nodeu = zeros(size(nodef1))
    return nodef,nodeu
end

function nodalMSAm(pak)
    return nodalMSm(pak)
end
function nodalMSBm(pak)
    return nodalMSm(pak)
end

function nodalMnum(pak)
    # -- copied from MSm
    # According to Constituents.pdf, MSm is misnamed, MNum.
    # Both the M2 and the Nu2 pieces have amp modulation like M2,
    # but they subtract to cancel phase modulation.
    nodef1,nodeu1 = nodalM2(pak)
    nodef = nodef1.^2
    nodeu = zeros(size(nodef1))
    return nodef,nodeu
end

function nodalMnuAm(pak)
    return nodalMnum(pak)
end
function nodalMnuBm(pak)
    return nodalMnum(pak)
end
    
function nodalNO1(pak)
    nodef1,nodeu1 = nodalM2(pak) # <-- nodalN2(pak)
    nodef2,nodeu2 = nodalO1(pak)
    nodef =  nodef1.*nodef2.^2
    nodeu =  nodeu1 .- 2.0*nodeu2
    return nodef,nodeu
end

function nodalMN4(pak)
    nodef1,nodeu1 = nodalM2(pak) # <-- also nodalN2(pak)
    nodef =  2.0*nodef1
    nodeu =  2.0*nodeu1
    return nodef,nodeu
end

function nodalOQ2(pak)
    nodef1,nodeu1 = nodalO1(pak)
    nodef =  nodef1.^2
    nodeu =  2.0*nodeu1
    return nodef,nodeu
end

function nodalKOo(pak)
    nodef1,nodeu1 = nodalK1(pak)
    nodef2,nodeu2 = nodalO1(pak)
    nodef =  nodef1.*nodef2
    nodeu =  nodeu1 .- nodeu2
    return nodef,nodeu
end

function nodalK2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu =         (-17.74*sinN .+   0.68*sin2N .-   0.04*sin3N)*dtr
    nodef = 1.0246 .+ 0.2863*cosN .+ 0.0083*cos2N .- 0.0015*cos3N
    return nodef,nodeu
end


function nodalM1B(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =      2.783*sin2p .+ 0.558*sin2pN .+ 0.184*sinN
    fcosu = 1. .+ 2.783*cos2p .+ 0.558*cos2pN .+ 0.184*cosN
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodalM1(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =      sinp .+ 0.2*sinpN
    fcosu = 2.0*( cosp .+ 0.2*cospN )
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end
    
function nodalM1A(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =     -0.3593*sin2p .- 0.2*sinN .- 0.066*sin2pN
    fcosu = 1. .+ 0.3593*cos2p .+ 0.2*cosN .+ 0.066*cos2pN
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodalgamma2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =      0.147*sin2N2p
    fcosu = 1. .+ 0.147*cos2N2p
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodalalpha2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =    .- 0.0446*sinpp1
    fcosu = 1. .- 0.0446*cospp1
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodaldelta2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =      0.477*sinN
    fcosu = 1. .- 0.477*cosN
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodalxi2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =    .- 0.439*sinN
    fcosu = 1. .+ 0.439*cosN
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodaleta2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    return nodalxi2(pak)
end

function nodalL2(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    fsinu =    .- 0.2505*sin2p .- 0.1102*sin2pN .+ 0.0156*sin2N2p .- 0.037*sinN
    fcosu = 1. .- 0.2505*cos2p .- 0.1102*cos2pN .- 0.0156*cos2N2p .- 0.037*cosN
    nodeu = atan.(fsinu,fcosu)
    nodef = sqrt.( fsinu.^2 .+ fcosu.^2 )
    return nodef,nodeu
end

function nodalMS(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = ( 2.14*sinN)*dtr
    nodef = 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N
    return nodef,nodeu
end

    # Uh-oh. nodalMSo is repeated. Which is correct?
    # What is MS vs MSo? nodalMSo seems to have wrong sign on phase!
function nodalMSo(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodef =  nodef1
    nodeu = -nodeu1
    return nodef,nodeu
end

function nodalMSo(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = ( 2.14*sinN)*dtr
    nodef = 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N
    return nodef,nodeu
end

function nodal2SM(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = 2.0*( 2.14*sinN)*dtr
    nodef = ( 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N ).^2
    return nodef,nodeu
end

function nodal2SMA(pak)
    return nodal2SM(pak)
end
function nodal2SMB(pak)
    return nodal2SM(pak)
end
    
    # Uh-oh. M3 is repeated. Why not just define it in terms of nodalM2?
function nodalM3(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = (-3.21*sinN)*dtr
    nodef = (sqrt.( 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N )).^3
    return nodef,nodeu
end

function nodalM3(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = -3.0*(1.07*sinN)*dtr
    nodef = sqrt.( 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N ).^3
    return nodef,nodeu
end

function nodalM5(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = -5.0*(1.07*sinN)*dtr
    nodef = sqrt.( 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N ).^5
    return nodef,nodeu
end

function nodalM7(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = -7.0*(1.07*sinN)*dtr
    nodef = sqrt.( 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N ).^7
    return nodef,nodeu
end

function nodalM9(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    nodeu = -9.0*(1.07*sinN)*dtr
    nodef = sqrt.( 1.0007 .- 0.0373*cosN .+ 0.0002*cos2N ).^9
    return nodef,nodeu
end

function nodalKQ1(pak)
    (cosp,sinp,
     cosN,sinN,
     cos2p,sin2p,
     cos2N,sin2N,
     cos3N,sin3N,
     cospN,sinpN,
     cospp1,sinpp1,
     cos2pN,sin2pN,
     cos2N2p,sin2N2p) = pak

    #    1 954 556 = 1 4 0 -1 0 0 1 :
    # K2 .- Q1 : 2 755 555 .- 1 356 554 = 2 2 0 0 0 0 0 .- 1 -2 0 1 0 0 -1 = 1 4 0 -1 0 0 1
    nodef1,nodeu1 = nodalK2(pak)
    nodef2,nodeu2 = nodalO1(pak)
    nodeu = nodeu1 .- nodeu2
    nodef = nodef1.*nodef2
    return nodef,nodeu
end

function nodalMK3(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodef2,nodeu2 = nodalK1(pak)
    nodeu = nodeu1 .+ nodeu2
    nodef = nodef1.*nodef2
    return nodef,nodeu
end

function nodal2MK3(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodef2,nodeu2 = nodalK1(pak)
    # Hmmm: not sure about this:
    nodeu = 2*nodeu1 .- nodeu2
    nodef = nodef1.*nodef2
    return nodef,nodeu
end

function nodalMO3(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodef2,nodeu2 = nodalO1(pak)
    nodeu = nodeu1 .+ nodeu2
    nodef = nodef1.*nodef2
    return nodef,nodeu
end

function nodalSK3(pak)
#    nodef1,nodeu1 = nodalS2(pak) This is zero
    nodef2,nodeu2 = nodalK1(pak)
    nodeu = nodeu2
    nodef = nodef2
    return nodef,nodeu
end

function nodalMS4(pak)
    nodef1,nodeu1 = nodalM2(pak)
#    nodef2,nodeu2 = nodalS2(pak) This is zero
    nodeu = nodeu1
    nodef = nodef1
    return nodef,nodeu
end

function nodal2SM2(pak)
    nodef1,nodeu1 = nodalM2(pak)
#    nodef2,nodeu2 = nodalS2(pak) This is zero
    nodeu = 2.0*nodeu1
    nodef = nodef1.*nodef1
    return nodef,nodeu
end

function nodalM4(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodeu = 2.0*nodeu1
    nodef = nodef1.*nodef1
    return nodef,nodeu
end

function nodalM6(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodeu = 3.0*nodeu1
    nodef = nodef1.*nodef1.*nodef1
    return nodef,nodeu
end

function nodalM8(pak)
    nodef1,nodeu1 = nodalM2(pak)
    nodeu = 4.0*nodeu1
    nodef = nodef1.^4
    return nodef,nodeu
end

export HA
function HA(tjd,xh,cidvec,donodal::Int64=1,doz0::Int64=1)

    iid, F = FMAT(tjd,cidvec,donodal,doz0)

    try
        # Solve via normal equations:
        #        x = (F'*F)\ (F'*xh)
        # Solve via QR:
        x = F \ xh
    catch
        println("Harmonic analysis failed:")
        #        U,S,V = svd( F'*F )
        R = svd( F'*F )
        U=R.U
        S=R.S
        V=R.Vt
        #        indx = find([ S[i] <= 1.e-8*S[1] for i=1:length(S)])
        #        indu = find([ S[i] >  1.e-8*S[1] for i=1:length(S)])
        indx = findall( x -> x < 1.e-8*S[1], S)
        indu = findall( x -> x >= 1.e-8*S[1], S)
        Sinv = zeros(size(S))
        Sinv[indu] = 1.0/S[indu]
        x = V*(diagm(0 => Sinv)*(U'*(F'*xh)))
        for ii=1:length(indx)
            println("bad cid = ",indx[ii])
        end
    end

    # Prediction and residual:
    pred = F*x
    res  = xh .- pred

    # Screen for outliers and redo if necessary:
    resstd = median( abs.( res ))*1.4826
    Rcrit  = 4.685
    e = abs.( res/resstd )
    #    indb = find([ e[i] >= Rcrit for i=1:xnt])
    indb = findall( x -> x >= Rcrit, e)

    if (length(indb) > 0)
        #        indg = find([ e[i] < Rcrit for i=1:xnt])
        indg = findall( x -> x < Rcrit, e)
        nt = length(indg)
        F1 = F[indg,:]
        xh1 = xh[indg]
        # via normal equations:
        # x = (F1'*F1)\ (F1'*xh1)
        # via QR:
        x = F1 \ xh1
        # Compute the residual at all the points, even those excluded by the outlier rejection.
        # This will inflate the error bars on the harmonic constants, which is probably okay!
        #pred = F*x
        #res  = xh .- pred
    end
    # Although the harmonic constants are computed after one pass of outlier rejection,
    # the error estimates and residuals are computed based on the original dataset, with
    # outliers included.

    pred = F*x
    res  = xh .- pred

    resvar = var(res)
    
    # My fancy stderr based on pseudo-spectrum for unevenly-spaced data:
    win = hamming(length(res))
    resh2 = abs.(fft(win.*res)).^2 ./mean(win.^2)
    resx2 = copy(resh2)
    resh2[2:end] = 0.5*(resh2[2:end] + resx2[end:-1:2])
    tmp = 0.33*(resh2[1:end-2] + resh2[2:end-1] + resh2[3:end])
    resh2[2:end-1] = tmp ; resh2[1] = resh2[2] ; resh2[end] = resh2[1]
    tmp = 0.33*(resh2[1:end-2] + resh2[2:end-1] + resh2[3:end])
    resh2[2:end-1] = tmp ; resh2[1] = resh2[2] ; resh2[end] = resh2[1]
    # resh2 is now an estimate of the pseudo-spectrum of the residual, i.e., the diagonal error covariance.
    Fh = zeros(Complex{Float64},size(F))
    len,ncid=size(Fh)
    for k=1:ncid
        tmp = view(F,:,k)
        Fh[:,k] = fft(tmp)
    end
    # d = Fx => Fh d = Fh F x => Fh d should have a diagonal covariance.
    # x x' = (Fh F)^-1 Fh d d' Fh' (Fh F)^-'
    # x x' = (Fh F)^-1 Ch ((Fh F)^-1 Ch)'
    # Ch is the square root of the diagonal covariance, C = Fh d d' Fh'
    # With Fh = I, d d'= I*s this would normally collapse to:
    #    x x' = (F)^-1 d d'(F)^-' = (F' F)^-1 s
    # Make an ensemble estimate, so this is just like a colored noise bootstrap except
    # it might be applied to irregularly-spaced data.
    stderrp = zeros(ncid,)
    nens = 10
    resh = sqrt.(resh2)
    for k=1:nens
        e = complex.( randn(len,).*resh, randn(len,).*resh )
        tmp = Fh \ e
        stderrp[:] = stderrp[:] + abs.(tmp).^2
    end
    stderrp = sqrt.( stderrp ./ nens )

#    return iid,x,xstd,res
    return iid,x,stderrp,res

end


export SYN
function SYN(tjd,xvec,cidvec,donodal::Int64=1)
    #
    # Sythesize a time series from the harmonic constants stored in xvec
    # at the requested times tjd.
    # The vector of harmonic constants in xvec must be laid out exactly
    # the same as the x vector output by HA, WITH ONE CHANGE:
    # The Z0 entry of x is omitted.
    # In other words, the iid pointers output by HA must all be decremented
    # by 1 and the xvec must discard the first entry in x.
    #
    # On output we have the iid corresponding to xvec for a check,
    # and the predicted time series.
    #

    iid, F = FMAT(tjd,cidvec,donodal,0)

    # Prediction and residual:
    pred = F*xvec

    # Return both iid and the prediction so that we can double-check the
    # ordering of the real and imaginary parts!
    return iid,pred

end

export FMAT
    function FMAT(tjd,cidvec,donodal::Int64=1,doz0::Int64=1)
        #
        # Just create and return the F-matrix
        # at the requested times tjd.
        #
        
        xnt = length(tjd)
        
        ncvec = length(cidvec)
        idvec = [ id[cidvec[i]] for i=1:ncvec]
        xdovec = zeros(ncvec,7)
        for i=1:ncvec
            for j=1:7
                xdovec[i,j] = xdo[idvec[i],j]
            end
        end
        
        orbel = tshpnp(tjd)
        args = orbel*xdovec'*dtr
    
    if (donodal == 1)
        
        # p and N for the nodal corrections, converted to radians:
        orbp  = orbel[:,4]*dtr
        orbN  = orbel[:,5]*dtr
        orbp1 = orbel[:,6]*dtr
        
        cosp = cos.(orbp)
        sinp = sin.(orbp)
        
        cosN = cos.(orbN)
        sinN = sin.(orbN)
        
        cos2p = cos.(2.0*orbp)
        sin2p = sin.(2.0*orbp)
        
        cos2N = cos.(2.0*orbN)
        sin2N = sin.(2.0*orbN)

        cos3N = cos.(3.0*orbN)
        sin3N = sin.(3.0*orbN)

        cospN = cos.(orbp .- orbN)
        sinpN = sin.(orbp .- orbN)

        cospp1 = cos.(orbp .- orbp1)
        sinpp1 = sin.(orbp .- orbp1)

        cos2pN = cos.(2.0*orbp .- orbN)
        sin2pN = sin.(2.0*orbp .- orbN)

        cos2N2p = cos.(2.0*orbN .- 2.0*orbp)
        sin2N2p = sin.(2.0*orbN .- 2.0*orbp)

        pak = (cosp,sinp,
               cosN,sinN,
               cos2p,sin2p,
               cos2N,sin2N,
               cos3N,sin3N,
               cospN,sinpN,
               cospp1,sinpp1,
               cos2pN,sin2pN,
               cos2N2p,sin2N2p)
        
        nodevec = node[idvec]
        f = zeros(xnt,ncvec)
        u = zeros(xnt,ncvec)
        for i=1:ncvec
            nodallettercode = nodevec[i]
            cid1 = cidvec[i]
            # Find which function will construct the function name:
            fnode = Symbol("Node"*"$nodallettercode")
            # Evaluate to get the function name:
            fnode = eval(fnode)("$cid1")
            # Now evaluate the nodal correction:
            nodef,nodeu = eval(fnode)(pak)
            # Stick these into the array which will be used to create F:
            f[:,i] = nodef[:]
            u[:,i] = nodeu[:]
        end

        if (doz0 == 1)
            F = [ ones(xnt,1) f.*cos.(args .- u) f.*sin.(args .- u) ]
        else
            F = [             f.*cos.(args .- u) f.*sin.(args .- u) ]
        end

    else # donodal == 0

        if (doz0 == 1)
            F = [ ones(xnt,1) cos.(args) sin.(args) ]
        else
            F = [             cos.(args) sin.(args) ]
        end

    end

        # Set up an index to permit grabbing the cos/sin coefficients:
        if (doz0 == 1)
            iid = Dict([("Z0",1)])
        else
            iid = Dict()
        end
        for i=1:ncvec
            s = cidvec[i] * "c"
            iid = merge(iid,Dict([(s,doz0 .+ i)]))
        end
        for i=1:ncvec
            s = cidvec[i] * "s"
            iid = merge(iid,Dict([(s,doz0 .+ ncvec .+ i)]))
        end

    # Return both iid and the prediction so that we can double-check the
    # ordering of the real and imaginary parts!
    return iid,F
end

# More efficient computation of the nodal corrections?
# This is only about 12% faster than FMAT for a single time
# series.
export FMATx
    function FMATx(tjd,cidvec,donodal::Int64=1,doz0::Int64=1)
        #
        # Just create and return the F-matrix
        # at the requested times tjd.
        #
        
        xnt = length(tjd)
        
        ncvec = length(cidvec)
        idvec = [ id[cidvec[i]] for i=1:ncvec]
        xdovec = zeros(ncvec,7)
        for i=1:ncvec
            for j=1:7
                xdovec[i,j] = xdo[idvec[i],j]
            end
        end
        
        orbel = tshpnp(tjd)
        args = orbel*xdovec'*dtr
    
    if (donodal == 1)

        nodevec = node[idvec]
        f = zeros(xnt,ncvec)
        u = zeros(xnt,ncvec)
        nodef = ones(ncvec,)
        nodeu = zeros(ncvec,)

        fnods = Array{Any}(undef,ncvec)
        for i=collect(1:ncvec)
            nodallettercode = nodevec[i]
            cid1 = cidvec[i]
            # Find which function will construct the function name:
            fnode = Symbol("Node"*"$nodallettercode")
            # Evaluate to get the function name:
            fnode = eval(fnode)("$cid1")
            fnods[i] = eval(fnode)
        end
        
        j=0
        dn = true
        tlast = tjd[1]
        while (j < xnt)
            j=j+1
            if (dn)
                # p and N for the nodal corrections, converted to radians:
                orbp  = orbel[j,4]*dtr
                orbN  = orbel[j,5]*dtr
                orbp1 = orbel[j,6]*dtr
                
                cosp = cos.(orbp)
                sinp = sin.(orbp)
                
                cosN = cos.(orbN)
                sinN = sin.(orbN)
                
                cos2p = cos.(2.0*orbp)
                sin2p = sin.(2.0*orbp)
                
                cos2N = cos.(2.0*orbN)
                sin2N = sin.(2.0*orbN)
                
                cos3N = cos.(3.0*orbN)
                sin3N = sin.(3.0*orbN)
                
                cospN = cos.(orbp .- orbN)
                sinpN = sin.(orbp .- orbN)
                
                cospp1 = cos.(orbp .- orbp1)
                sinpp1 = sin.(orbp .- orbp1)
                
                cos2pN = cos.(2.0*orbp .- orbN)
                sin2pN = sin.(2.0*orbp .- orbN)
                
                cos2N2p = cos.(2.0*orbN .- 2.0*orbp)
                sin2N2p = sin.(2.0*orbN .- 2.0*orbp)
                
                pak = (cosp,sinp,
                       cosN,sinN,
                       cos2p,sin2p,
                       cos2N,sin2N,
                       cos3N,sin3N,
                       cospN,sinpN,
                       cospp1,sinpp1,
                       cos2pN,sin2pN,
                       cos2N2p,sin2N2p)

                for i=collect(1:ncvec)
                    fnode = fnods[i]
                    # Now evaluate the nodal correction:
                    ff,uu=fnode(pak)
                    nodef[i]=ff[1]
                    nodeu[i]=uu[1]
                end
                tlast = tjd[j]
            end # if dn
            # Evaluate the nodal correction at most once
            # every 10 days:
            jp1 = min(j+1,xnt)    
            dn = ( (tjd[jp1]-tlast) > 10. )
            # Stick these into the array which will be used to create F:
            f[j,:] = nodef[:]
            u[j,:] = nodeu[:]
        end

        if (doz0 == 1)
            F = [ ones(xnt,1) f.*cos.(args .- u) f.*sin.(args .- u) ]
        else
            F = [             f.*cos.(args .- u) f.*sin.(args .- u) ]
        end

    else # donodal == 0

        if (doz0 == 1)
            F = [ ones(xnt,1) cos.(args) sin.(args) ]
        else
            F = [             cos.(args) sin.(args) ]
        end

    end

        # Set up an index to permit grabbing the cos/sin coefficients:
        if (doz0 == 1)
            iid = Dict([("Z0",1)])
        else
            iid = Dict()
        end
        for i=1:ncvec
            s = cidvec[i] * "c"
            iid = merge(iid,Dict([(s,doz0 .+ i)]))
        end
        for i=1:ncvec
            s = cidvec[i] * "s"
            iid = merge(iid,Dict([(s,doz0 .+ ncvec .+ i)]))
        end

    # Return both iid and the prediction so that we can double-check the
    # ordering of the real and imaginary parts!
    return iid,F
end


end
    
