import numpy
from scipy.stats import mannwhitneyu

cancerlist = ["OV"]

input_file1 = []
output_file1 = []

probe_count = 1

threshold = 0.2

sample_id = []
cytoact = []

sample_index= []

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
header = "CpGsite\tFC(high_mean/low_mean)\tFC(high_median/low_median)\tP-value\n"

for i in range(0, len(cancerlist)) :

    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))
    output_file1.append(open("Vector.FC.Pvalue." + cancerlist[i] + ".txt", 'w'))

    output_file1[i].write(header)
    
    sample_header1 = input_file1[i].readline().split()
    input_file1[i].readline().split()

    del sample_header1[0]; del sample_header1[0]

    sample_index = []
    length = len(sample_header1)
    
    for j in range(0, length) :
        sample_header1[j] = sample_header1[j][:15].replace('-', '')
        if(sample_header1[j] in sample_id) : sample_index.append(sample_id.index(sample_header1[j]))
        else : sample_index.append(-1)

    for j in range(0, probe_count) :

        line1 = input_file1[i].readline().split()
        site_id = line1.pop(0)

        betavalue_low = []; betavalue_high = []
        cytoact_low = []; cytoact_high = []

        for k in range(0, length) :
            cancer_value = line1[k]
            print(sample_index[k])
            if(cancer_value == "NA" or sample_index[k] == -1) : continue
            
            if(float(cancer_value) <= threshold) :
                betavalue_low.append(float(cancer_value))
                cytoact_low.append(cytoact[sample_index[k]])
            else :
                betavalue_high.append(float(cancer_value))
                cytoact_high.append(cytoact[sample_index[k]])

        high_mean = numpy.mean(cytoact_high)
        low_mean = numpy.mean(cytoact_low)

        high_median = numpy.median(cytoact_high)
        low_median = numpy.median(cytoact_low)
        
        FC1 = high_mean / low_mean
        FC2 = high_median / low_median

        print(cytoact_high)
        print(cytoact_low)

        print(len(cytoact_high), len(cytoact_low))
        manwhitney_pair = mannwhitneyu(cytoact_high, cytoact_low)

        printline = site_id
        printline += "\t%s\t%s\t%s\n" % (str(FC1), str(FC2), str(manwhitney_pair[1]))
        output_file1[i].write(printline)

        if(j % 10000 == 0) : print(cancerlist[i] + " %d completed." % j)

    print(cancerlist[i] + " completed.")
