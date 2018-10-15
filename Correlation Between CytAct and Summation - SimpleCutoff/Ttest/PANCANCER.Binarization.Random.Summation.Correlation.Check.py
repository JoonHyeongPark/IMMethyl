from scipy import stats
import random

#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]
cancerlist = ["PANCANCER"]

input_file1 = []
input_file2 = []

#probe_count = 485577 
probe_count = 1

p_threshold = [0.05, 0.005, 0.0005, 0.0000001]

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

percentage = [0.01, 0.025, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]
execution_number = 100

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

    for iteration in range(execution_number) :
        sample_binary_table.append([])
        for j in range(len(p_threshold)) :
            sample_binary_table[iteration].append([])
            for k in range(len(percentage)) :
                sample_binary_table[iteration][j].append([])
                for l in range(length) : sample_binary_table[iteration][j][k].append(0)

    print("Default Completed.")

############################################################################################################################################################################

    whole_skew = []
    whole_skew_index = []
        
    for j in range(len(p_threshold)) :
        input_file = open(str(p_threshold[j]) + "." + cancerlist[i] + ".CpGsites.By.TTest.txt", 'r')
        input_file.readline() # junk line

        whole_skew.append([])

        lines = input_file.readlines()
        for line in lines : # Derivation of meaningful CpG sites
                
            line = line.split()
            whole_skew[j].append(line[0])

        whole_skew[j].append("END_POINT")
        whole_kew_index.append(0)

    for j in range(probe_count) :
        line1 = input_tumor.readline().split()
        site_id = line1.pop(0)

############################################################################################################################################################################
# getting betavalue for each cpg site
            
        betavalue_row = []
        new_length = length
        for k in range(0, length) :
            if(line1[k] == "NA") :
                new_length -= 1 
                continue
            betavalue_row.append(k)

############################################################################################################################################################################

        if(new_length > 0) :
                
            for iteration in range(excution_number) :
                    
                for k in range(len(p_threshold)) :
                        
                    if(whole_skew[k][whole_skew_index[k]] == site_id) :

                        for percentage_i in range(len(percentage)) :
                            threshold = int(new_length * percentage[percentage_i])
                            for l in range(threshold) : sample_binary_table[iteration][k][percentage_i][betavalue_row[random.randrange(new_length)]] += 1

                        whole_skew_index[k] += 1
    
        if(j % 10000 == 0) :
            print(cancerlist[i] + " %d completed." % j)

    for j in range(len(p_threshold)) :
        
        for k in range(len(percentage)) :
            
            output_file = open("Pvalue." + str(p_threshold[j]) + ".Percentage." + str(percentage[k]) + ".Correlation.By.Random.Table.txt", 'w')
            
            for iteration in range(execution_number) :
                
                cyt_vector = []
                betavalue_vector = []
                
                for l in range(length) :
                    cyt_vector.append(cytact[sample_index[l]])
                    betavalue_vector.append(sample_binary_table[iteration][j][k][l])
                    
                spearman_pair = stats.spearmanr(cyt_vector, betavalue_vector)
                output_file.write("%dth\t" % (iteration) + str(spearman_cor[0]) + "\n")
