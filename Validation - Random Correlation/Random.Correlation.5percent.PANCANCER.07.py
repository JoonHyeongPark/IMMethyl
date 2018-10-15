import random
import scipy.stats as stat

#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]
cancerlist = ["PANCANCER"]

probe_count = 485577

percentage = 0.05
execution_number = 500

sample_id = []
cytoact = []

ORDER = "RANDOM.07."

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

for i in range(0, len(cancerlist)) :
    
    input_tumor = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')
    
    sample_header1 = input_tumor.readline().split() # sample line
    del sample_header1[0]; del sample_header1[0]
    input_tumor.readline() # junk line
    
    length = len(sample_header1)

    for j in range(length) :
        sample_header1[j] = sample_header1[j][:15].replace('-', '')

    output_file = open(ORDER + cancerlist[i] + ".Correlation.Values." + str(execution_number) + ".Times.txt", 'w')
    
    sample_summation = []
    temp_summation = []
    
    for j in range(length) :
        if(sample_header1[j] not in sample_id) : temp_summation.append(-1)
        else : temp_summation.append(0)
        
    for iteration in range(execution_number) :
        sample_summation.append([])
        for j in range(length) :
            sample_summation[iteration].append(temp_summation[j])
    
    input_index = open(cancerlist[i] + ".Valid.Samples.For.Each.CpGsite.txt", 'r')
        
    for j in range(probe_count) :
          
        line = input_index.readline().split()
        site_id = line.pop(0)

        choosing_number = int(len(line) * percentage)
        
        line = map(int, line)
                
        for iteration in range(execution_number) :
    
            chosen_index = random.sample(line, choosing_number)
            for index in chosen_index : sample_summation[iteration][index] += 1

    cor_cytoact = [] # valid cytact
    for j in range(length) :
        if(temp_summation[j] != -1) : cor_cytoact.append(cytoact[sample_id.index(sample_header1[j])])
        
    for iteration in range(execution_number) :

        cor_summation = []
        for j in range(length) :
            if(sample_summation[iteration][j] != -1) : cor_summation.append(sample_summation[iteration][j])
                
        output_file.write(str(stat.spearmanr(cor_cytoact, cor_summation)[0]) + "\n")
