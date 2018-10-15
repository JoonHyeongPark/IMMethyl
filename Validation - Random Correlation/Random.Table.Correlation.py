cancerlist = ["PANCANCER"]

p_threshold = [0.05, 0.005, 0.0005, 0.0000001]
percentage = [0.01, 0.025, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]

input_invalid_samples = open("PANCANCER.Invalid.Samples.By.CytAct.txt", 'r')
invalid_samples = input_invalid_samples.read().splitlines()

sample_id = []
cytoact = []

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

GetSample()

for i in range(len(cancerlist)) :

    for j in range(len(p_threshold)) :

        input_file = open(str(p_threshold[j]) + "." + cancerlist[i] + ".ValidCpGsites.Of.Each.Sample.txt", 'r')

        probability = []
        cytact = []

        lines = input_file.readlines()

        total = 0
        
        for line in lines :
            
            line = line.split()
            sample = line[0][:15].replace('-', '')
            
            if(sample in sample_id) :
                probability.append(int(line[1]))
                cytact.append(cytoact[sample_id.index(sample)])
                total += probability[-1]

        probability = map(lambda x : float(probability[x] / total), range(len(probability)))
               
        print(probability)
        print(cytact)
