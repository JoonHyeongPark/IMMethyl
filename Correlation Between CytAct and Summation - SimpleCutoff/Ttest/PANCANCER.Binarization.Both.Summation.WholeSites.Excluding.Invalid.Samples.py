from operator import itemgetter

#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]
cancerlist = ["PANCANCER"]

input_file1 = []
input_file2 = []

probe_count = 485577 

sample_id = []
cytoact = []

sample_index = []

def GetSample() :
    
    cytoact_file = open("TCGA_methylation_cowork_1.txt", 'r')
    header = cytoact_file.readline().split() # getting header

    id_posit = header.index("id") # sample ID positioning
    cytoact_posit = header.index("CytAct") # CytAct positioning
    cytodata = cytoact_file.readlines() # read data table
    cytoact_file.close()

    count = 0

    global sample_id
    global cytoact
    
    for line in cytodata :
        line = line.split()
        sample_id.append(line[id_posit].replace('_', '')) # sample ID extraction
        cytoact.append(float(line[cytoact_posit])) # CytAct extraction

        count += 1
    
    return count # Sample number return

sample_number = GetSample()

percentage = [0.005, 0.0125, 0.025, 0.05, 0.0625, 0.075, 0.0875, 0.1, 0.1125, 0.125, 0.1375, 0.15]
print_percentage = [0.01, 0.025, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]

for i in range(0, len(cancerlist)) :

    input_tumor = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')

    sample_header1 = input_tumor.readline().split() # sample line
    input_tumor.readline() # junk line

############################################################################################################################################################################
# make sample index table

    del sample_header1[0]; del sample_header1[0]

    sample_index = []
    sample_binary_table = []

    length = len(sample_header1)

    for j in range(0, length) :
        sample_header1[j] = sample_header1[j][:15].replace('-', '')
        if(sample_header1[j] in sample_id) : sample_index.append(sample_id.index(sample_header1[j]))
        else : sample_index.append(-1)
        
    for k in range(len(percentage)) :
        sample_binary_table.append([])
        for l in range(length) : sample_binary_table[k].append(0)

############################################################################################################################################################################
    
    for j in range(probe_count) :
        line1 = input_tumor.readline().split()
        site_id = line1.pop(0)

############################################################################################################################################################################
# getting betavalue for each cpg site
        
        betavalue_row = []
        new_length = length
        for k in range(0, length) :
            if(line1[k] == "NA" or sample_index[k] == -1) :
                new_length -= 1
                continue
            betavalue_row.append([float(line1[k]), k])

        betavalue_row.sort(key = itemgetter(0))

############################################################################################################################################################################

        if(new_length > 0) :
                    
            for percentage_i in range(len(percentage)) :
                threshold = int(new_length * percentage[percentage_i])
                for l in range(threshold) : sample_binary_table[percentage_i][betavalue_row[new_length - l - 1][1]] += 1
                for l in range(threshold) : sample_binary_table[percentage_i][betavalue_row[l][1]] += 1

        if(j % 10000 == 0) :
            print(cancerlist[i] + " %d completed." % j)

    for k in range(len(percentage)) :
        output_file = open("WholeSites.Percentage." + str(print_percentage[k]) + "." + cancerlist[i] + ".Both.Binarization.Summation.Excluding.Invalid.Samples.txt", 'w')
        for l in range(length) :
            if(sample_index[l] == -1) : continue
            printline = sample_header1[l] + "\t%s\n" % str(sample_binary_table[k][l])
            output_file.write(printline)
