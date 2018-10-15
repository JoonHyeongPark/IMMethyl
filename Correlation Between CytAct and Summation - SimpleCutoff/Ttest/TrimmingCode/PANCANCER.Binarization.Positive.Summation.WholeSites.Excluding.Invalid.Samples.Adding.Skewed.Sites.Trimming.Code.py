from numpy import zeros
from operator import itemgetter

#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]
cancerlist = ["PANCANCER"]

probe_count = 485577 

sample_id = []
cytoact = []

percentage = [0.01, 0.025, 0.05, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3]

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

def GetSampleHeader(sample_header, length) :

    global sample_id

    sample_index = []
    original_sample_header = []
    
    for j in range(0, length) :
    
        original_sample_header.append(sample_header[j])
        sample_header[j] = sample_header[j][:15].replace('-', '')
        
        if(sample_header[j] in sample_id) : sample_index.append(sample_id.index(sample_header[j]))
        else : sample_index.append(-1)
        
    return sample_header, sample_index, original_sample_header

def GetValidBetavalueRowAndSorted(line, length, sample_index) :

    betavalue_row = []
    new_length = length
    
    for k in range(0, length) :
        if(line[k] == "NA" or sample_index[k] == -1) :
                new_length -= 1
                continue
        betavalue_row.append([float(line[k]), k])

    betavalue_row.sort(key = itemgetter(0))

    return betavalue_row, new_length

def MedianInSortedRow(betavalue_row, new_length) :
    
    betavalue_median = betavalue_row[new_length / 2][0]
    if(new_length % 2 == 0) : betavalue_median = (betavalue_median + betavalue_row[new_length / 2 - 1][0]) / 2
    
    return betavalue_median

def GetSampleBinaryTable(cancer) :

    input_tumor = open(cancer + ".humanmethylation450.tumor.txt", 'r')
    
    output_left_skewed = open(cancer + ".Left.Skewed.CpGsites.txt", 'w')
    output_right_skewed = open(cancer + ".Left.Skewed.CpGsites.txt", 'w')

    sample_header1 = input_tumor.readline().split() # sample line
    input_tumor.readline() # junk line

    del sample_header1[0]; del sample_header1[0]

    sample_index = []
    sample_binary_table = []

    length = len(sample_header1)

    sample_header, sample_index, original_sample_header = GetSampleHeader(sample_header1, length)
    sample_binary_table = zeros((len(percentage), length, 2), dtype = int)

    for j in range(probe_count) :
        
        line1 = input_tumor.readline().split()
        site_id = line1.pop(0)

        betavalue_row, new_length = GetValidBetavalueRowAndSorted(line1, length, sample_index) # getting betavalue for each cpg site

        if(new_length > 0) :

            betavalue_median = MedianInSortedRow(betavalue_row, new_length)

            if(betavalue_median > 0.5) :
               for percentage_i in range(len(percentage)) :
                    threshold = int(new_length * percentage[percentage_i])
                    for l in range(threshold) : sample_binary_table[percentage_i][betavalue_row[new_length - l - 1][1]][1] += 1

               output_right_skewed.write(site_id + "\n")
                        
            else :
                for percentage_i in range(len(percentage)) :
                    threshold = int(new_length * percentage[percentage_i])
                    for l in range(threshold) : sample_binary_table[percentage_i][betavalue_row[l][1]][0] += 1

                output_left_skewed.write(site_id + "\n")

        if(j % 10000 == 0) : print(cancer + " %d completed." % j)

    input_tumor.close()
    output_left_skewed.close()
    output_right_skewed.close()

    return cancer, sample_binary_table, original_sample_header, sample_index, length

def PrintSampleBinaryTable(cancer, sample_binary_table, sample_header, sample_index, length) :

    global percentage
    global cytoact

    for k in range(len(percentage)) :
        
        output_file = open("WholeSites.Percentage." + str(percentage[k]) + "." + cancer + ".Positive.Binarization.Summation.Excluding.Invalid.Samples.txt", 'w')
        
        header = "Site\tHypoMethylationBurden\tHyperMethylationBurden\tWholeMethylationBurden\tCytAct\n"
        output_file.write(header)
        
        for l in range(length) :
            if(sample_index[l] == -1) : printline = sample_header[l] + "\tNA\tNA\tNA\tNA\n"
            else : printline = sample_header[l] + "\t%s\t%s\t%s\t%s\n" % str(sample_binary_table[k][l][0]), str(sample_binary_table[k][l][1]), str(sample_binary_table[k][l][0] + sample_binary_table[k][l][1]), str(cytoact[sample_index[l]])
            output_file.write(printline)
            
        output_file.close()
    
    return

def main() :

    sample_number = GetSample()

    for i in range(0, len(cancerlist)) : PrintSampleBinaryTable(GetSampleBinaryTable(cancerlist[i]))

    return

if(__name__ == "__main__") : main()
