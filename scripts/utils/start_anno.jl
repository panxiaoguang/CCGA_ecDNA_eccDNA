using GenomicFeatures

function get_intervals(file::String)
    col = IntervalCollection{Vector{SubString{String}}}()
    for line in eachline(file)
        shuju = split(line, "\t")
        push!(col, Interval(shuju[1], parse(Int64, shuju[2]), parse(Int64, shuju[3]), '.', shuju))
    end
    col
end

function get_overlap(case::IntervalCollection{Vector{SubString{String}}}, db::IntervalCollection{Vector{SubString{String}}}, IO::IOStream)
    for (a, b) in eachoverlap(case, db)
        if leftposition(b) <= leftposition(a) <= rightposition(b)
            println(IO, join(GenomicFeatures.metadata(a), "\t"), "\t", join(GenomicFeatures.metadata(b), "\t"))
        end
    end
end



function Anno(sample::String)
    case = get_intervals("bedFile/randomBEDs/$(sample).random.bed")  ## filtered ecc bed 
    dbs = get_intervals("dbs/hg38.coding.bed")   ## database
    w = open("bedFile/randomBEDs/gene_anno/$(sample).startAnno.bed", "w+")  ### result output
    get_overlap(case, dbs, w)
    close(w)
end

for sample in ["cBca_1N", "cBca_2N", "cBca_3N", "cBca_4N",
    "cBca_5N", "cBca_6N", "cBca_7N", "cBca_8N", "cBca_9N",
    "cBca_10N", "cBca_11N", "cBca_12N", "cBca_13N",
    "cBca_14N", "cBca_15N", "cBca_16N", "cBca_17N",
    "cBca_18N", "cBca_19N", "cBca_20N", "cBca_21N",
    "cBca_23N", "cBca_24N", "cBca_25N", "cBca_26N",
    "cBca_27N", "cBca_28N", "cBca_29N", "cBca_30N",
    "cBca_31N", "cBca_32N", "cBca_33N", "cBca_34N",
    "cBca_35N", "cBca_36N", "cBca_37N", "cBca_38N",
    "cBca_39N", "cBca_40N", "cBca_41N", "cBca_42N",
    "cBca_43N", "cBca_44N", "cBca_45N", "cBca_46N",
    "cBca_47N", "cBca_48N", "cBca_50N", "cBca_51N",
    "cBca_52N", "cBca_54N", "cBca_55N", "cBca_56N",
    "cBca_57N", "cBca_58N", "cBca_59N", "cBca_60N",
    "cBca_62N", "cBca_63N", "cBca_64N", "cBca_65N",
    "cBca_66N", "cBca_67N", "cBca_68N", "cBca_69N",
    "cBca_70N", "cBca_71N", "cBca_72N", "cBca_74N",
    "cBca_75N", "cBca_76N", "cBca_77N", "cBca_78N",
    "cBca_79N", "cBca_80N", "cBca_84N", "cBca_85N",
    "cBca_86N", "cBca_87N", "cBca_88N"]
    Anno(sample)
end

