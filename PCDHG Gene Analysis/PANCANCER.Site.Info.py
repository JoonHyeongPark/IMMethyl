from numpy import mean, median

cancerlist = ["PANCANCER"]

probe_site = 485577

high_threshold = 0.8
low_threshold = 0.2

header = "Site\tTumorMean\tNormalMean\tMeanGap\tTumorMedian\tNormalMedian\tMedianGap\tTumorMax\tNormalMax\tMaxGap\tTumorMin\tNormalMin\tMinGap\t#ValidTumorSamples\t#ValidNormalSamples\n"

for i in range(len(cancerlist)) :

    input_tumor = open(cancerlist[i] + ".humanmethylation450.tumor.txt", 'r')
    input_normal = open(cancerlist[i] + ".humanmethylation450.normal.txt", 'r')

    output = open(cancerlist[i] + ".Extreme.Site.Tumor." + str(high_threshold) + ".Normal." + str(low_threshold) + ".txt", 'w')
    output_site_infor = open(cancerlist[i] + ".Site.INFO.txt", 'w')

    output_site_infor.write(header)

    sample_header1 = input_tumor.readline().split()
    del sample_header1[0]; del sample_header1[0]
    input_tumor.readline()

    tumor_sample_length = len(sample_header1)

    sample_header2 = input_normal.readline().split()
    del sample_header2[0]; del sample_header2[0]
    input_normal.readline()

    normal_sample_length = len(sample_header2)

    for j in range(probe_site) :

        if(j % 10000 == 0) : print(cancerlist[i] + " %d completed." % j)

        tumor_line = input_tumor.readline().split()
        normal_line = input_normal.readline().split()

        site_id1 = tumor_line.pop(0)
        site_id2 = normal_line.pop(0)

        new_length1 = tumor_sample_length
        new_length2 = normal_sample_length

        tumor_beta = []
        normal_beta = []

        for k in range(tumor_sample_length) :

            if(tumor_line[k] == "NA") :
                new_length1 -= 1
                continue

            tumor_beta.append(float(tumor_line[k]))

        for k in range(normal_sample_length) :

            if(normal_line[k] == "NA") :
                new_length2 -= 1
                continue

            normal_beta.append(float(normal_line[k]))

        if(new_length1 == 0 or new_length2 == 0) : continue

        tumor_mean = mean(tumor_beta)
        normal_mean = mean(normal_beta)

        tumor_median = median(tumor_beta)
        normal_median = median(normal_beta)

        tumor_max = max(tumor_beta)
        normal_max = max(normal_beta)

        tumor_min = min(tumor_beta)
        normal_min = min(normal_beta)

        output_site_infor.write("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n" % (site_id1, tumor_mean, normal_mean, abs(tumor_mean - normal_mean), tumor_median, normal_median, abs(tumor_median - normal_median), tumor_max, normal_max, abs(tumor_max - normal_max), tumor_min, normal_min, abs(tumor_min - normal_min), new_length1, new_length2))

        
