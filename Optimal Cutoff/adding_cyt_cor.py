from scipy import stats
import matplotlib.pyplot as plt

plt.style.use('ggplot')
probe_count = 485577

cancerlist = ["KIRC", "BRCA", "THCA", "HNSC", "LIHC", "PRAD", "UCEC", "KIRP", "LUSC", "COAD", "LUAD", "PANCANCER"]

sample_table = {}

input_cytact = open("TCGA_methylation_cowork_1.txt",'r')

header = input_cytact.readline().split()
indexing = header.index("CytAct")

whole_table = input_cytact.readlines()

for line in whole_table :
    line = line.split()
    ID = line[0].replace("_", "")

    sample_table[ID] = line[indexing]

input_file = []
output_file = []

for i in range(0, len(cancerlist)) :

    input_file.append(open(cancerlist[i] + ".optimal_cutoff_result.txt", 'r'))
    output_file.append(open(cancerlist[i] + ".optimal_cutoff_including_cytact.txt", 'w'))

    header = "Sample_ID\tNumber_of_sites\tCytAct\n"
    output_file[i].write(header)

    reading_table = input_file[i].readlines()

    num_table = []
    cytact_table = []

    for line in reading_table :

        line = line.split()
        if(line[0] not in sample_table or line[1] == "0") : continue

        line.append(sample_table[line[0]])

        printline = ""
        for j in range(0, len(line)) : printline += "%s\t" % line[j]
        printline += "\n"
        output_file[i].write(printline)

        num_table.append(float(line[1]))
        cytact_table.append(float(line[2]))

    pearson = stats.pearsonr(num_table, cytact_table)
    spearman = stats.spearmanr(num_table, cytact_table)

    last_printline1 = "Pearson : %f\t%f\n" % (pearson[0], pearson[1])
    last_printline2 = "Spearman : %f\t%f\n" % (spearman[0], spearman[1])

    output_file[i].write(last_printline1)
    output_file[i].write(last_printline2)

    plt.scatter(num_table, cytact_table)
    plt.show()
