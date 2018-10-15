cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_file1 = []
output_file1 = []

gap = 0.01

probe_count = 485577

threshold_array = []
accumulation_array1 = []

header = "Threshold"

def MakeThresholdArray(gap) :
        
    global threshold_array
    global header

    cutoff = 0.0
    count = 0
        
    while(cutoff <= 1.0) :
        threshold_array.append(cutoff)
        accumulation_array1.append(0)

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

    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
    output_file1.append(open(cancerlist[i] + "." + str(gap) + ".Betavalue_WholeTumorCount.txt", 'w'))
    
    input_file1[i].readline()
    input_file1[i].readline()

    output_file1[i].write(header)

cancer_count = []

for i in range(0, probe_count) :

    cancer_count.append([])

    for j in range(0, len(cancerlist)) :

        cancer_count[i].append([])
        for k in range(0, len(threshold_array)) :
            cancer_count[i][j].append(0)

        line1 = input_file1[j].readline().split()
    
        while(len(line1) == 0) :
            line1 = input_file1[j].readline().split()
            print(cancerlist[j] + " tumor %d " % i + "has no line value.")

        sample_id = line1.pop(0)
        for cancer_value in line1 :
            if(cancer_value == "NA") : continue
            cancer_count[i][j][ReturnRange(float(cancer_value))] += 1

        printline1 = sample_id
        for k in range(0, len(threshold_array) - 1) :
            printline1 += "\t%s" % str(cancer_count[i][j][k])
        printline1 += "\n"

        output_file1[j].write(printline1)

    if(i % 10000 == 0) : print("%s completed." % str(i))
