using DelimitedFiles
using Printf

test_num = "MOSCATOhac"
casename = "salazar"

idarr = []
logtext = readdlm("logMOSCATOhac")
for i = 1:size(logtext,1)
if logtext[i,1]=="Job"
  tmptext = logtext[i,2]
  push!(idarr,parse(Int,tmptext[2:end-1]))
end
end

filelist = []
dirnum = 0.00:0.0045:0.09
for i = 1:length(idarr)
  id = idarr[i]
  idstr = @sprintf("%.4f",dirnum[i])
  push!(filelist,"MOSCATO-hac-$idstr/$casename.$id")
  @show id
end
data = zeros(length(idarr),12)

fc = 1
for file in filelist
run(pipeline(`awk 'x+=/Step 100000/' $file`, stdout = file*"-last"))
text = readdlm(file*"-last")
totlines = size(text,1)

current = 0.0
Vcell = 0.0
ifcurrent = false
ifVcell = false
i = totlines

while ((i >= 1) && !(ifcurrent && ifVcell))
    line = text[i,:]
    if (line[1]=="avg" && line[2]=="anode" && line[3]=="current" && !ifcurrent)
        current = line[5]
        ifcurrent = true
    end
    if (line[1]=="V_cell5:" && !ifVcell)
        Vcell = line[2]
        ifVcell = true
    end
    i = i - 1
end

println("file:$file")
println("current,Vcell:")
println("$current,$Vcell")
data[fc,1] = current
data[fc,2] = Vcell

global fc = fc + 1
end

open("MOSCATO-hac-data/currentanode-$test_num.txt","w") do io 
    writedlm(io,data[:,1])
end
open("MOSCATO-hac-data/Vcell-$test_num.txt","w") do io 
    writedlm(io,data[:,2])
end

