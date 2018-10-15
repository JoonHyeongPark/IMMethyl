import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')

import seaborn as sns

cancerlist = ["PANCANCER"]

execution_number = 10000

target_value = -4.0

for i in range(len(cancerlist)) :
    
    input_file = open("RANDOM.MERGE." + cancerlist[i] + ".Correlation.Values." + str(execution_number) + ".Times.txt", 'r')
    pvalues = map(float, input_file.read().splitlines())

    pvalues.sort(reverse = True)

    rank = 1
    for pvalue in pvalues :
        if(target_value > pvalue) : break
        rank += 1

    sns.set()
    sns.distplot(pvalues, kde = True, hist = True)

    plt.title("#10000 Random Correlation-values' Density Plot")
    plt.xlabel("Random Correlation-value")
#    plt.ylabel("Density")
    plt.grid(True)
    figure = plt.gcf()
    plt.show()
    figure.savefig(cancerlist[i] + ".Random.Correlation.Values.Density.Plot.pdf")
