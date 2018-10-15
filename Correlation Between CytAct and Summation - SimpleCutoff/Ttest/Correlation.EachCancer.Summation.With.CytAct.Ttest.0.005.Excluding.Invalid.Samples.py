from operator import itemgetter
from scipy import stats
import numpy
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot
import matplotlib.pylab

matplotlib.pyplot.style.use('ggplot')

cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]

input_file1 = []
output_file1 = []

probe_count = 485577

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

PN = ["Positive", "Negative", "Both", "Tstat"]
percentage = ["0.01", "0.025", "0.05", "0.1", "0.125", "0.15", "0.175", "0.2", "0.225", "0.25", "0.275", "0.3"]

P_value = 0.005

correlation_table = []

for i in range(0, len(cancerlist)) :

   correlation_table.append([])
    
   for k in range(len(percentage)) : 

        correlation_table[i].append([])

        summation = []
        cytact = []

        for j in range(len(PN)) :
            summation.append([])
            cytact.append([])

        total_count = []
        zero_count = []

        for j in range(len(PN)) :

            total_count.append(0)
            zero_count.append(0)

            correlation_table[i][k].append([])

            input_file1 = open("Pvalue." + str(P_value) + ".Percentage." + percentage[k] + "." + cancerlist[i] + "." + PN[j] + ".Binarization.Summation.Excluding.Invalid.Samples.txt", 'r')
            lines = input_file1.readlines()
    
            for line in lines :
                
                line = line.split()
                if(line[0] in sample_id) :
                    index = sample_id.index(line[0])
                    summation[j].append(float(line[1]))
                    cytact[j].append(cytoact[index])
                    total_count[j] += 1
                    if(int(line[1]) == 0) : zero_count[j] += 1
    
        xlimit = max(map(lambda x : max(summation[x]), range(len(PN))))
        ylimit = max(map(lambda x : max(cytact[x]), range(len(PN))))

        for j in range(len(PN)) :

            if(len(summation[j]) < 1 / float(percentage[k]) or total_count[j] == zero_count[j]) :
                print("Pvalue." + str(P_value) + "." + percentage[j] + ", " + cancerlist[i] + " " + PN[j] + " " + percentage[k] + " has less than %s samples." % str(1 / float(percentage[k])))
                correlation_table[i][k][j].append("NA")
                correlation_table[i][k][j].append("NA")
                continue

            matplotlib.pyplot.hold(False)
            matplotlib.pyplot.title(cancerlist[i] + " " + PN[j] + " " + percentage[k] + " T-test's P-value < " + str(P_value))
            matplotlib.pyplot.xlabel("Summation")
            matplotlib.pyplot.ylabel("CytAct")
            matplotlib.pyplot.legend()
            matplotlib.pyplot.grid(True)
 
    #        pearson_pair = stats.pearsonr(summation[j], cytact[j])
            spearman_pair = stats.spearmanr(summation[j], cytact[j])

            spearman_cor = "{0:.3f}".format(spearman_pair[0])
            spearman_P = "{0:.3f}".format(spearman_pair[1])

            correlation_table[i][k][j].append("{0:.5f}".format(spearman_pair[0]))
            correlation_table[i][k][j].append("{0:.5f}".format(spearman_pair[1]))

            matplotlib.pyplot.scatter(summation[j], cytact[j], marker = '.', c = 'b')
            matplotlib.pyplot.title(cancerlist[i] + " " + PN[j] + "    Spearman Cor = %s, P = %s" % (spearman_cor, spearman_P))
            matplotlib.pyplot.hold(True)
            matplotlib.pylab.plot(summation[j], cytact[j], '.')
            z = numpy.polyfit(summation[j], cytact[j], 1)
            p = numpy.poly1d(z)
            matplotlib.pylab.plot(summation[j], p(summation[j]), "r--")
    
            figure = matplotlib.pyplot.gcf()
            matplotlib.pyplot.xlim(0, xlimit)
            matplotlib.pyplot.ylim(0, ylimit)
            matplotlib.pyplot.show()
            figure.savefig("Pvalue." + str(P_value) + "." + cancerlist[i] + "." + percentage[k] + "." + PN[j] + ".pdf")

output_file = open("Pvalue." + str(P_value) + ".EachCancer.Correlation.CytAct.Between.Meaningful.Betavalue.Count.Of.Each.CpGsite.EveryPercentage.Excluding.Invalid.Samples.txt", 'w')

header2 = "Percentage"
for j in range(len(PN)) : header2 += "\t%s\t%s" % (PN[j] + "(Cor)", PN[j] + "(P)")
header2 += "\n"
output_file.write(header2)
 
for k in range(len(percentage)) :
#    output_file = open(percentage[k] + ".Meaningful.Betavalue.Count.Of.Each.CpGsite.Positive.Negative.Both.Correlation.Table.txt", 'w')
   
    for i in range(0, len(cancerlist)) :
        printline = percentage[k]
        for j in range(len(PN)) : printline += "\t%s\t%s" % (correlation_table[i][k][j][0], correlation_table[i][k][j][1])
        output_file.write(printline + "\n")
