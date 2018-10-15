cancerlist = ["PANCANCER"]

input_file1 = []
output_file1 = []

threshold = 0.2

probe_count = 485577

for i in range(0, len(cancerlist)) :

    input_file1.append(open(str(threshold) + ".Cutoff.FC.Pvalue." + cancerlist[i] + ".txt", 'r'))
    output_file1.append(open(str(threshold) + ".MeaningfulCpGsitesByPvalue0.05.Without.Extreme." + cancerlist[i] + ".txt", 'w'))
    
    input_file1[i].readline()

    for j in range(0, probe_count) :

        line1 = input_file1[i].readline().split()
        site_id = line1.pop(0)

        if(line1[0] == "NoSamples" or len(line1) == 1) : continue
        if(float(line1[2]) > 0.05) : continue

        printline = site_id
        for k in range(0, len(line1)) :
            printline += "\t" + line1[k]
            
        printline += "\n"
        
        output_file1[i].write(printline)

        if(j % 10000 == 0) : print(cancerlist[i] + " %d completed." % j)
