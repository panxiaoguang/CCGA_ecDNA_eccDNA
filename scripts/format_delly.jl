using Glob

function parse_line(ln::AbstractString)
    CHROM, POS, ID, _, _, _, _, INFO, _, _, _ = split(ln, "\t")
    svclass = replace(match(r"SVTYPE=\w+", INFO).match, "SVTYPE=" => "")
    pesupport = replace(match(r"PE=\d+", INFO).match, "PE=" => "")
    strand = replace(match(r"CT=\Sto\S", INFO).match, "CT=" => "")
    POS2 = svclass == "BND" ? replace(match(r"POS2=\d+", INFO).match, "POS2=" => "") : replace(match(r"END=\d+", INFO).match, "END=" => "")
    CHROM2 = svclass == "BND" ? replace(match(r"CHR2=chr[\dXYM]{1,2}", INFO).match, "CHR2=" => "") : CHROM
    POS2 = parse(Int64, POS2)
    POS = parse(Int64, POS)
    (CHROM, POS, CHROM2, POS2, strand, pesupport, svclass, ID)
end

function parse_delly(fs::String, OUT)
    shunxu = Dict("chr1" => 1, "chr2" => 2, "chr3" => 3, "chr4" => 4, "chr5" => 5, "chr6" => 6, "chr7" => 7, "chr8" => 8, "chr9" => 9, "chr10" => 10, "chr11" => 11, "chr12" => 12, "chr13" => 13, "chr14" => 14, "chr15" => 15, "chr16" => 16, "chr17" => 17, "chr18" => 18, "chr19" => 19, "chr20" => 20, "chr21" => 21, "chr22" => 22, "chrX" => 23, "chrY" => 24, "chrM" => 25)
    HEADER = "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsv_id\tpe_support\tstrand1\tstrand2\tsvclass\tsvmethod"
    println(OUT, HEADER)
    for line in eachline(fs)
        if startswith(line, "#")
            continue
        else
            CHROM, POS, CHROM2, POS2, strand, pesupport, svclass, ID = parse_line(line)
            if svclass == "BND"
                if strand == "5to5"
                    if shunxu[CHROM] > shunxu[CHROM2]
                        println(OUT, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", "TRA", "\t", "Delly")
                    else
                        println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", "TRA", "\t", "Delly")
                    end
                elseif strand == "3to3"
                    if shunxu[CHROM] > shunxu[CHROM2]
                        println(OUT, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", "TRA", "\t", "Delly")
                    else
                        println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", "TRA", "\t", "Delly")
                    end
                elseif strand == "5to3"
                    if shunxu[CHROM] > shunxu[CHROM2]
                        println(OUT, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", "TRA", "\t", "Delly")
                    else
                        println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "TRA", "\t", "Delly")
                    end
                elseif strand == "3to5"
                    if shunxu[CHROM] > shunxu[CHROM2]
                        println(OUT, CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", CHROM, "\t", POS, "\t", POS + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "TRA", "\t", "Delly")
                    else
                        println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", "TRA", "\t", "Delly")
                    end
                else
                    @warn "Something was wrong.  You should pay attention"
                end
            elseif svclass == "INV"
                if strand == "5to5"
                    println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "-", "\t", svclass, "\t", "Delly")
                elseif strand == "3to3"
                    println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "+", "\t", svclass, "\t", "Delly")
                else
                    @warn "Something was wrong.  You should pay attention"
                end
            elseif svclass == "DEL"
                println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "+", "\t", "-", "\t", svclass, "\t", "Delly")
            else
                println(OUT, CHROM, "\t", POS, "\t", POS + 1, "\t", CHROM2, "\t", POS2, "\t", POS2 + 1, "\t", ID, "\t", pesupport, "\t", "-", "\t", "+", "\t", "INS/DUP", "\t", "Delly")
            end
        end
    end
end

function main()
    for fs in Glob.glob("delly_sv/*_delly.vcf")
        prx = replace(fs, "_delly.vcf" => "", "delly_sv/" => "")
        open("format_delly/format.$(prx).sv.delly.txt", "w+") do io
            parse_delly(fs, io)
        end
    end
end

main()
