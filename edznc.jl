using NetCDF
find(a,b)=findall(a,b)

function ncheader(fn)
    # grab the header and global attributes in a structure:
    nci = NetCDF.open(fn)
    NetCDF.close(nci)
    return nci
end
function ncvarget(fn,v;start=[],count=[],stride=[],xend=[])
    if ( (length(count) > 0) & (length(xend) > 0) )
        println("Error: cannot specify both count and xend keyword")
        println("       parameters to ncvarget.")
        return NaN
    end
    if ( (length(start) == 0) & (length(xend) > 0) )
        println("Error: must specify start keyword if xend keyword")
        println("       is given to ncvarget.")
        return NaN
    end
    if ( (length(start) == 0) & (length(count) > 0) )
        println("Error: must specify start keyword if count keyword")
        println("       is given to ncvarget.")
        return NaN
    end
    if ( (length(start) == 0) & (length(stride) > 0) )
        println("Error: must specify start keyword if stride keyword")
        println("       is given to ncvarget.")
        return NaN
    end
    if (length(start) == 0)
        ix=ncread(fn,v)
    else
        if (length(stride) == 0)
            if (length(count) > 0)
                ix=ncread(fn,v,start=start,count=count)
            elseif (length(count) == 0)
                if (length(xend) == 0)
                    ix=ncread(fn,v,start=start)
                else
                    count = xend .- start .+ 1
                    ix=ncread(fn,v,start=start,count=count)
                end
            end
        else
            jj=0
            # Get the dimlen of the second y-dimension:
            # nci=ncinfo(fn)
            # Argh! Latest version of NetCDF does not provide a useful ncinfo command
            # for Julia 1.2!
            nci = NetCDF.open(fn)
            # This is funky. The dot (period) is the fieldname
            # separator. The colon indicates the fieldname.
            nx=Int(nci.:vars[v].:dim[1].:dimlen)
            ny=Int(nci.:vars[v].:dim[2].:dimlen)
            NetCDF.close(nci)
            # Either count or xend or neither are specified, not both:
            if (length(count) > 0)
                xend[1] = start[1] + count[1] - 1
                xend[2] = start[2] + count[2] - 1
                nx = min(xend[1],nx)
                ny = min(xend[2],ny)
            elseif (length(xend) > 0)
                nx = min(xend[1],nx)
                ny = min(xend[2],ny)
            else
                # neither: do nothing.
            end
                
            ix=zeros(length(start[1]:stride[1]:nx),
                     length(start[2]:stride[2]:ny))
            for j=collect(start[2]:stride[2]:ny)
                jj=jj+1
                tmp=ncread(fn,v,start=[start[1],j],count=[-1,1])
                ix[:,jj]=tmp[1:stride[1]:nx+1-start[1]]
            end
        end
    end
    fx=ncgetatt(fn,v,"_FillValue")
    if ( fx == nothing )
        fx=ncgetatt(fn,v,"missing_value")
    end
    if ( fx == nothing )
        
    else
        indx=find(x -> x == fx, ix)
        # We can only handle NaN values if the type is real!
        if (length(indx) > 0)
            ix = convert.(Float64,ix)
        end
    end ;
    sx=ncgetatt(fn,v,"scale_factor")
    if ( sx == nothing )
        x = ix
    else
        x = sx.*ix
    end
    ox=ncgetatt(fn,v,"add_offset")
    if ( ox == nothing )
    else
        x = x .+ ox
    end
    if ( fx == nothing )
    else
        x[indx] .= NaN
    end
    return x
end

# Not sure what is name of variable in the file? Pass a list in and we'll
# search the options.
function ncvarget(fin,vars::Array{String,1})
    for v=vars
        try
            var=ncvarget(fin,v)
            return var
        catch
        end
    end
    println("ERROR: Could not find the requested variable in the input file.")
    println("  fin = ",fin)
    println("  var list = ",vars)
    return nothing
end
import NetCDF.ncgetatt
function ncgetatt(fin,vars::Array{String,1},att)
    for v=vars
        try
            att=ncgetatt(fin,v,att)
            return att
        catch
        end
    end
    println("ERROR: Could not find the requested variable or attribute in the input file.")
    println("  fin = ",fin)
    println("  var list = ",vars)
    println("  att = ",att)
    return nothing
end
