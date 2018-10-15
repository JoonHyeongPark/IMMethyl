#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]
cancerlist = ["PANCANCER"]

probe_count = 485577 

cpgsites = []

p_threshold = [0.05, 0.005, 0.0005, 0.0000001]

for i in range(0, len(cancerlist)) :

    site_index = []
 
    for j in range(len(p_threshold)) :

        cpgsites.append([])
        site_index.append(0)

        input_cpgsites = open(str(p_threshold[j]) + "." + cancerlist[i] + ".CpGsites.By.TTest.txt", 'r')
        input_cpgsites.readline()
        lines = input_cpgsites.read().splitlines()

        for line in lines : cpgsites[j].append(line.split()[0])
        cpgsites[j].append("END_POINT")

    input_tumor = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')

    sample_header1 = input_tumor.readline().split() # sample line
    del sample_header1[0]; del sample_header1[0]

    print("#Samples : %s" % str(len(sample_header1)))

    valid_cpgsites = []
    for j in range(len(p_threshold)) :
        valid_cpgsites.append([])
        for k in range(len(sample_header1)) : valid_cpgsites[j].append(0)
    
    input_tumor.readline() # junk line
    
    for j in range(probe_count) :
        line1 = input_tumor.readline().split()
        site_id = line1.pop(0)

        for k in range(len(p_threshold)) :

            if(site_id == cpgsites[k][site_index[k]]) :

                for l in range(len(line1)) :
                    if(line1[k] == "NA") : continue
                    valid_cpgsites[k][l] += 1

                site_index[k] += 1

        if(j % 10000 == 0) : print(cancerlist[i] + " %d" % j)

    for j in range(len(p_threshold)) :
        output_file = open(str(p_threshold[j]) + "." + cancerlist[i] + ".Valid.CpGsties.Of.Each.Sample.txt", 'w')
        for k in range(len(sample_header1)) :
            output_file.write(sample_header1[k] + "\t%d\n" % valid_cpgsites[j][k])
