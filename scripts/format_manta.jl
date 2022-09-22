using Glob

function BNDdct(fs::String)
    dct = Dict{String,String}()
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            CHROM, POS, ID, _, _, _, _, INFO, _, _, _ = split(line, "\t")
            ID = convert(String, ID)
            svclass = replace(match(r"SVTYPE=\w+", INFO).match, "SVTYPE=" => "")
            if svclass == "BND"
                dct[ID] = line
            end
        end
    end
    dct
end

function getStrand(alt::AbstractString, alt2::AbstractString)
    if occursin("[", alt)
        strand1 = "5"
    else
        strand1 = "3"
    end
    if occursin("[", alt2)
        strand2 = "5"
    else
        strand2 = "3"
    end
    string(strand1, "to", strand2)
end

function parse_line(ln::AbstractString, dct::Dict{String,String})
    shunxu = Dict("chr1" => 1, "chr2" => 2, "chr3" => 3, "chr4" => 4, "chr5" => 5, "chr6" => 6, "chr7" => 7, "chr8" => 8, "chr9" => 9, "chr10" => 10, "chr11" => 11, "chr12" => 12, "chr13" => 13, "chr14" => 14, "chr15" => 15, "chr16" => 16, "chr17" => 17, "chr18" => 18, "chr19" => 19, "chr20" => 20, "chr21" => 21, "chr22" => 22, "chrX" => 23, "chrY" => 24, "chrM" => 25)
    CHROM, POS, ID, _, alt, _, _, INFO, _, _, _ = split(ln, "\t")
    svclass = replace(match(r"SVTYPE=\w+", INFO).match, "SVTYPE=" => "")
    pesupport = 0
    POS = parse(Int64, POS)
    mateid = ""
    if svclass == "BND"
        mateid = string(replace(match(r"MATEID=MantaBND[:\d]+", INFO).match, "MATEID=" => ""))
        CHROM2, POS2, ID2, _, alt2, _, _, INFO2, _, _, _ = split(dct[mateid], "\t")
        POS2 = parse(Int64, POS2)
        if shunxu[CHROM] > shunxu[CHROM2]
            tmp_chrom = CHROM2
            tmp_pos = POS2
            tmp_alt = alt2
            CHROM2 = CHROM
            POS2 = POS
            alt2 = alt
            CHROM = tmp_chrom
            POS = tmp_pos
            alt = tmp_alt
        end
        strand = getStrand(alt2, alt)
    elseif svclass == "INV"
        strand = split(INFO, ";")[end]
        POS2 = replace(match(r"END=\d+", INFO).match, "END=" => "")
        CHROM2 = CHROM
        POS2 = parse(Int64, POS2)
    else
        strand = "NtoN"
        POS2 = replace(match(r"END=\d+", INFO).match, "END=" => "")
        CHROM2 = CHROM
        POS2 = parse(Int64, POS2)
    end
    (CHROM, POS, CHROM2, POS2, strand, pesupport, svclass, ID, mateid)
end

function parse_manta(fs::String, IO)
    mateids = []
    BNDs = BNDdct(fs)
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
    println(IO, HEADER)
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            CHROM, POS, CHROM2, POS2, strand, pesupport, svclass, ID, mateid = parse_line(line, BNDs)
            if svclass == "BND"
                if ID in mateids
                    continue
                else
                    if strand == "5to5"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", "TRA", "\t", "Manta")
                    elseif strand == "3to3"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", "TRA", "\t", "Manta")
                    elseif strand == "5to3"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "TRA", "\t", "Manta")
                    elseif strand == "3to5"
                        println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", "TRA", "\t", "Manta")
                    else
                        @warn "Something was wrong.  You should pay attention"
                    end
                end
                push!(mateids, mateid)
            elseif svclass == "INV"
                if strand == "INV5"
                    println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", svclass, "\t", "Manta")
                elseif strand == "INV3"
                    println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", svclass, "\t", "Manta")
                else
                    @warn "Something was wrong.  You should pay attention"
                end
            elseif svclass == "DEL"
                println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", svclass, "\t", "Manta")
            else
                println(IO, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "INS/DUP", "\t", "Manta")
            end
        end
    end
end


function main()
    for fs in Glob.glob("manta_sv/*_manta.ft.vcf")
        prx = replace(fs, "_manta.ft.vcf" => "", "manta_sv/" => "")
        open("format_manta/format.$(prx).sv.manta.txt", "w+") do io
            parse_manta(fs, io)
        end
    end
end

main()
