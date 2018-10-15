cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]

input_file1 = []

probe_count = 485577 

sample_id = []
cytoact = []

sample_index = []

def GetSample() :
    
    cytoact_file = open("PANCANCER.CytAct.txt", 'r')
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
        sample_id.append(line[id_posit].replace('-', '')) # sample ID extraction
        cytoact.append(float(line[cytoact_posit])) # CytAct extraction

        count += 1
    
    return count # Sample number return

sample_number = GetSample()

for i in range(0, len(cancerlist)) :

    input_tumor = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')

    sample_header1 = input_tumor.readline().split() # sample line
    input_tumor.readline() # junk line

    del sample_header1[0]; del sample_header1[0]

    sample_index = []
    sample_binary_table = []

    length = len(sample_header1)
    invalid_sample = []
    
    for j in range(0, length) :
        sample_header1[j] = sample_header1[j][:15].replace('-', '')
        if(sample_header1[j] in sample_id) : sample_index.append(sample_id.index(sample_header1[j]))
        else :
            invalid_sample.append(sample_header1[j])
            sample_index.append(-1)
            
    print(len(invalid_sample))
    
    output_file = open(cancerlist[i] + ".Invalid.Samples.By.CytAct.txt", 'w')
    
    for invalid_sample_id in invalid_sample : output_file.write(invalid_sample_id + "\n") 
