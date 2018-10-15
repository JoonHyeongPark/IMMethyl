from scipy import stats

cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
#cancerlist = ["PANCANCER"]

input_file1 = []
input_file2 = []

probe_count = 485577

p_threshold = [0.0000001]
#p_threshold = [0.05, 0.005, 0.0005, 0.0000001]

header = "site\tStatistics\tP-value\n"

for i in range(0, len(cancerlist)) :

    input_file1.append(open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r'))

    input_file1[i].readline()
    input_file1[i].readline()

    input_file2.append(open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r'))

    normalsample = input_file2[i].readline().split()
    input_file2[i].readline()

    output_file = []
    for j in range(len(p_threshold)) :
        output_file.append(open(str(p_threshold[j]) + "." + cancerlist[i] + ".CpGsites.By.TTest.txt", 'w'))
        output_file[j].write(header)

    if(len(normalsample) < 12) :
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

        equal_variance = stats.levene(tumor_value, normal_value)
        if(equal_variance >= 0.05) : ttest_pair = stats.ttest_ind(tumor_value, normal_value, equal_var = True)
        else : ttest_pair = stats.ttest_ind(tumor_value, normal_value, equal_var = False)

        for k in range(len(p_threshold)) :
            if(ttest_pair[1] < p_threshold[k]) : output_file[k].write(site_id1 + "\t%f\t%f\n" % (ttest_pair[0], ttest_pair[1]))

        if(j % 10000 == 0) : print(cancerlist[i] + " %d completed." % j)

    print(cancerlist[i] + " completed.")

    input_file1[i].close()
    input_file2[i].close()
