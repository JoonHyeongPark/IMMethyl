cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_file1 = []

gap = 0.05

output_file1 = open(str(gap) + ".sample." + "Betavalue_WholeTumorCount.csv", 'w')

probe_count = 485577

threshold_array = []

header = "Threshold"

def MakeThresholdArray(gap) :
        
    global threshold_array
    global header
                
    cutoff = 0.0
    threshold_array.append([])

    count = 0
        
    while(cutoff <= 1.0) :    
        threshold_array.append(cutoff)
        count += 1
        cutoff = gap * float(count)
    
    for i in range(0, len(threshold_array) - 1) :
        header += "\t%s~%s" % (str(threshold_array[i]), str(threshold_array[i + 1]))

    return

sample_id = []
sample_count = {}

def GetSample() :
    
    cytoact_file = open("TCGA_methylation_cowork_1.txt", 'r')
    header = cytoact_file.readline().split() # getting header

    id_posit = header.index("id") # sample ID positioning
    cytoact_posit = header.index("CytAct") # CytAct positioning
    cytodata = cytoact_file.readlines() # read data table
    cytoact_file.close()

    count = 0

    global sample_id
    global sample_count
    global threshold_array
    
    for line in cytodata :
        line = line.split()
        sample_id.append(line[id_posit].replace('_', '')) # sample ID extraction

        sample_count[sample_id[count]] = []
        for i in range(len(threshold_array)) : sample_count[sample_id[count]].append(0)

        count += 1
    
    return

MakeThresholdArray(gap)
GetSample()

header += "\n"

output_file1.write(header)

def ReturnRange(value) :

    for i in range(0, len(threshold_array) - 1) :
        if(value >= threshold_array[i] and value <= threshold_array[i + 1]) : return i

for i in range(0, len(cancerlist)) :

    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')) 
    sample_header1 = input_file1[i].readline().split()
    input_file1[i].readline().split()

    del sample_header1[0]; del sample_header1[0]
    for j in range(0, len(sample_header1)) : sample_header1[j] = sample_header1[j][:15].replace('-', '')

    for j in range(0, probe_count) :

        line1 = input_file1[i].readline().split()
        while(len(line1) == 0) : line1 = input_file1[i].readline().split()
        site_id = line1.pop(0)

        for k in range(0, len(line1)) :
            cancer_value = line1[k]
            if(cancer_value == "NA" or sample_header1[k] not in sample_id) : continue
           
            sample_count[sample_header1[k]][ReturnRange(float(cancer_value))] += 1

for name in sample_id :
    printline1 = name
    print_array = sample_count[name].itemes()
    for data in print_array :
        printline1 += "\t%s" % str(data)
    printline1 += "\n"

    output_file1.write(printline1)
