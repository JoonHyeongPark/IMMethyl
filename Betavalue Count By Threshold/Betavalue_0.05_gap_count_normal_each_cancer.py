cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_file2 = []
output_file2 = []

gap = 0.05

probe_count = 485577

threshold_array = []
accumulation_array2 = []

header = "Threshold"

def MakeThresholdArray(gap) :
        
    global threshold_array
    global header

    cutoff = 0.0
    count = 0
        
    while(cutoff <= 1.0) :
        threshold_array.append(cutoff)
        accumulation_array2.append(0)

        count += 1
        cutoff = gap * float(count)

    print(threshold_array)

    for i in range(0, len(threshold_array) - 1) :
        header += "\t%s~%s" % (str(threshold_array[i]), str(threshold_array[i + 1]))

    return

MakeThresholdArray(gap)

header += "\n"

def ReturnRange(value) :

    for i in range(0, len(threshold_array) - 1) :
        if(value >= threshold_array[i] and value <= threshold_array[i + 1]) :
            return i

for i in range(0, len(cancerlist)) :

    input_file2.append(open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r'))
    output_file2.append(open(cancerlist[i] + "." + str(gap) + ".Betavalue_WholeNormalCount.txt", 'w'))
    
    input_file2[i].readline()
    input_file2[i].readline()

    output_file2[i].write(header)

normal_count = []

for i in range(0, probe_count) :

    normal_count.append([])

    for j in range(0, len(cancerlist)) :

        normal_count[i].append([])
        for k in range(0, len(threshold_array)) :
            normal_count[i][j].append(0)

        line2 = input_file2[j].readline().split()
    
        while(len(line2) == 0) :
            line2 = input_file2[j].readline().split()
            print(cancerlist[j] + " normal %d " % i + "has no line value.")

        sample_id = line2.pop(0)
        for normal_value in line2 :
            if(normal_value == "NA") : continue
            normal_count[i][j][ReturnRange(float(normal_value))] += 1

        printline2 = sample_id
        for k in range(0, len(threshold_array) - 1) :
            printline2 += "\t%s" % str(normal_count[i][j][k])
        printline2 += "\n"

        output_file2[j].write(printline2)

    if(i % 10000 == 0) : print("%s completed." % str(i))
