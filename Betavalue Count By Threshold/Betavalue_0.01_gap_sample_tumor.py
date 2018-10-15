cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_file1 = []

gap = 0.01

output_file1 = open(str(gap) + ".sample." + "Betavalue_WholeTumorCount.csv", 'w')

probe_count = 485577

threshold_array = []

header = "Threshold"

def MakeThresholdArray(gap) :
        
    global threshold_array
    global header
                
    cutoff = 0.0
    count = 0
        
    while(cutoff <= 1.0) :    
        threshold_array.append(cutoff)
        count += 1
        cutoff = gap * float(count)
    
    for i in range(0, len(threshold_array) - 1) :
        header += "\t%s~%s" % (str(threshold_array[i]), str(threshold_array[i + 1]))

    return

sample_id = []
sample_index= []
sample_count = []

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
        sample_count.append([])
        for i in range(len(threshold_array)) : sample_count[count].append(0)

        count += 1
    
    return

MakeThresholdArray(gap)

print(threshold_array)

header += "\n"

output_file1.write(header)

def ReturnRange(value) :

    for i in range(0, len(threshold_array) - 1) :
        if(value >= threshold_array[i] and value <= threshold_array[i + 1]) : return i

sample_number = 0          

for i in range(0, len(cancerlist)) :

    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')) 
    sample_header1 = input_file1[i].readline().split()
    input_file1[i].readline().split()

    del sample_header1[0]; del sample_header1[0]

    sample_index = []
    length = len(sample_header1)
    
    for j in range(0, length) :
        sample_header1[j] = sample_header1[j][:15].replace('-', '')

        if(sample_header1[j] not in sample_id) :
            sample_id.append(sample_header1[j])
            sample_count.append([])
            for k in range(len(threshold_array)) : sample_count[sample_number].append(0)

            sample_index.append(sample_number)
            sample_number += 1

        else : sample_index.append(sample_id.index(sample_header1[j]))

    for j in range(0, probe_count) :

        line1 = input_file1[i].readline().split()
        site_id = line1.pop(0)

        for k in range(0, length) :
            cancer_value = line1[k]
            if(cancer_value == "NA") : continue        
            sample_count[sample_index[k]][ReturnRange(float(cancer_value))] += 1

    print(cancerlist[i] + " completed.")
    
for i in range(0, len(sample_id)) :
    printline1 = sample_id[i]
    for j in range(0, len(threshold_array)) :
        printline1 += "\t%s" % str(sample_count[i][j])
    printline1 += "\n"

    output_file1.write(printline1)
