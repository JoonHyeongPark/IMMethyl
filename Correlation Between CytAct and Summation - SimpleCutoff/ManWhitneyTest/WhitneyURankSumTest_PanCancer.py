from scipy.stats import mannwhitneyu

#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
cancerlist = ["PANCANCER"]

input_file1 = []
input_file2 = []

output_file = []

probe_count = 485577

header = "site\tStatistics\tP-value\n"

cancer_number = len(cancerlist)

for i in range(0, cancer_number) :
    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))

    input_file1[i].readline()
    input_file1[i].readline()

    input_file2.append(open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r'))

    normalsample = input_file2[i].readline().split()
    input_file2[i].readline()

output_file = []
pvalue_threshold = [0.05, 0.005, 0.0005, 0.0000001]

for i in range(len(pvalue_threshold)) :
  output_file.append(open(str(pvalue_threshold[i]) + ".PANCANCER.CpGsites.By.ManWhitneyTest.txt", 'w'))
  output_file[i].write(header)

for i in range(0, probe_count) :

    tumor_value = []
    normal_value = []

    site_id1 = ""
    site_id2 = ""

    for j in range(0, cancer_number) :

        line1 = input_file1[j].readline().split()
        line2 = input_file2[j].readline().split()

        site_id1 = line1.pop(0)
        for value in line1 :
            if(value != "NA") : tumor_value.append(float(value))

        site_id2 = line2.pop(0)
        for value in line2 :
            if(value != "NA") : normal_value.append(float(value))

    if(len(tumor_value) == 0 or len(normal_value) == 0) : continue

    manwhitney_pair = mannwhitneyu(tumor_value, normal_value)
    printline = site_id1 + "\t%f\t%f\n" % (manwhitney_pair[0], manwhitney_pair[1])
    
    for j in range(len(pvalue_threshold)) :
        if(manwhitney_pair[1] < pvalue_threshold[j]) : output_file[j].write(printline)

    if(i % 10000 == 0) : print(str(i) + " completed")
