using CSV
using DataFrames

function getBreaks(graphFile::String)
    (sample, amp) = split(replace(graphFile, "_graph.txt" => ""), "_")
    lns = readlines(graphFile)
    ## find the first line that begin with discordant key word
    i = findfirst(x -> startswith(x, "discordant"), lns)
    if isnothing(i)
        return ""
    else
        ## drop out other lines
        newlns = lns[i:end]
        ## parse the breakpoint
        breakpoints = Vector{NamedTuple{(:chrom, :bp, :sample, :amp),Tuple{String,Int64,String,String}}}()
        for line in newlns
            (lft, rgt) = split(split(line, "\t")[2], "->")
            lft = split(rstrip(lft, ['-', '+']), ":")
            rgt = split(rstrip(rgt, ['-', '+']), ":")
            lft_tp = (chrom=lft[1], bp=parse(Int64, lft[2]), sample=sample, amp=amp)
            rgt_tp = (chrom=rgt[1], bp=parse(Int64, rgt[2]), sample=sample, amp=amp)
            if lft_tp.bp != -1
                if !(lft_tp in breakpoints)
                    push!(breakpoints, lft_tp)
                end
                if !(rgt_tp in breakpoints)
                    push!(breakpoints, rgt_tp)
                end
            else
                if !(rgt_tp in breakpoints)
                    push!(breakpoints, rgt_tp)
                end
            end
        end
        return DataFrame(breakpoints)
    end
end

df = getBreaks("13_amplicon10_graph.txt")  ## for example
# println(df)
