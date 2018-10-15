from scipy.stats import mannwhitneyu

cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_file1 = []
input_file2 = []

output_file = []

probe_count = 485577

header = "site\tStatistics\tP-value\n"

for i in range(0, len(cancerlist)) :

    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))

    input_file1[i].readline()
    input_file1[i].readline()

    input_file2.append(open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r'))

    normalsample = input_file2[i].readline().split()
    input_file2[i].readline()

    output_file.append(open(cancerlist[i] + "CpGsitesByManWhitneyTest.txt", 'w'))
    output_file[i].write(header)

    if(len(normalsample) < 31) :
        print(cancerlist[i] + " has only %d normal sample." % (len(normalsample) - 2))
        continue

    for j in range(0, probe_count) :

        line1 = input_file1[i].readline().split()
        line2 = input_file2[i].readline().split()

        tumor_value = []
        normal_value = []

        site_id1 = line1.pop(0)
        for value in line1 :
            if(value != "NA") : tumor_value.append(float(value))

        site_id2 = line2.pop(0)
        for value in line2 :
            if(value != "NA") : normal_value.append(float(value))

        manwhitney_pair = mannwhitneyu(tumor_value, normal_value)

        printline = site_id1 + "\t%f\t%f\n" % (manwhitney_pair[0], manwhitney_pair[1])
        output_file[i].write(printline)

    print(cancerlist[i] + " completed.")

    input_file1[i].close()
    input_file2[i].close()

    output_file[i].close()
