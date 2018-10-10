using DelimitedFiles
using Plots
using LaTeXStrings

wavesInFile = readdlm("src/wavelist_formated.txt")

path_to_thresholds = "src/thresholds_formated.txt"

thresholds = fill(0.0,size(wavesInFile,1))
let path = path_to_thresholds
    if isfile(path) && size(wavesInFile,1)==88
        v = readdlm(path)
        for i in 1:size(v,1)
            global thresholds[Int64(v[i,1])] = v[i,2]
        end
    else
        warn("Do not consider thresholds!")
    end
end

# gr()
# default(fmt = :png)
let mv = 0.5:0.02:2.5
    plot(mv, [sum((thresholds[i]<m) for i in 1:88) for m in mv],
        xlab = L"M_{3π} bin", ylab="Number of waves")
end
savefig("/tmp/test.pdf")

# plot(rand(10))

# [sum((thresholds[i]<m) for i in 1:88) for m in 0.5:0.02:2.5]
