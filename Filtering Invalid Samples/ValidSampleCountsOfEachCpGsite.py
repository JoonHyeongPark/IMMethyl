#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM", "PANCANCER"]
cancerlist = ["PANCANCER"]

probe_count = 485577 

for i in range(0, len(cancerlist)) :

    input_tumor = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')

    sample_header1 = input_tumor.readline().split() # sample line
    input_tumor.readline() # junk line

    output_file = open(cancerlist[i] + ".Valid.Sample.Counts.Of.Each.CpGsite.txt", 'w')
    
    for j in range(probe_count) :
        line1 = input_tumor.readline().split()
        site_id = line1.pop(0)

        new_length = 0
        for data in line1 :
            if(data == "NA") : continue
            new_length += 1

        output_file.write(site_id + "\t%s\n" % str(new_length))
