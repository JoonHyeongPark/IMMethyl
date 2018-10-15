# -*- coding: utf-8 -*-
#cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "ESCASTAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
cancerlist = ["STAD"]
from operator import itemgetter
#from scipy import stats
#import numpy as np

betavalue_arr = []
cytoact_arr = []

probe_name = []
sample_id = []

def getting_cytoact() :
    cytoact_file = open("TCGA_methylation_cowork_1.txt", 'r')
    header = cytoact_file.readline().split() # header ?쎄린

    id_posit = header.index("id") # sample ID ?꾩튂 ?몃뜳??    cancer_posit= header.index("origin") # cancer type ?꾩튂 ?몃뜳??    cytoact_posit = header.index("CytAct") # cytact ?꾩튂 ?몃뜳??        
    cytoact_posit = header.index("CytAct")
    cytodata = cytoact_file.readlines() # 愿 ??data瑜??듭㎏濡??쎌뼱?ㅼ엫
    cytoact_file.close()
    
    for line in cytodata :
        line = line.split()
        sample_id.append(line[id_posit].replace('_', '')) # sample ID 異붿텧
        
    ############# betavalue 由ъ뒪??珥덇린??###################
    probe_count = 485577; sample_count = len(sample_id)
    #probe_count = 21; sample_count = len(sample_id)
    
    for reset_x in range(0, probe_count) :
        betavalue_arr.append([])
        for reset_y in range(0, len(cancerlist) + 1) :
            betavalue_arr[reset_x].append([])
            for reset_z in range(0, sample_count) : betavalue_arr[reset_x][reset_y].append([None, None])
    ###########################################################
    
    for i in range(0, len(cancerlist)) :
        cytoact_arr.append([]) # 誘몃━ cytoact dictionary?먯꽌 cancer type蹂?cytact data瑜?list ?뺥깭濡?? ?ν븷 寃껋쓣 紐낆떆        
        for j in range(0, sample_count) : cytoact_arr[i].append(None)

    cytoact_arr.append([])

    for line in cytodata :
        line = line.split() # sample 蹂꾨줈 ?뺣━??data瑜??쎌뼱?ㅼ씠怨?遺꾩젅
        
        if(line[cytoact_posit] != "NA") :
            cancer_posit = cancerlist.index(line[cancer_posit].replace('_', ''))
            #cancer_posit = cancerlist.index("STAD")
            sample_posit = sample_id.index(line[id_posit].replace('_', ''))
            
            cytoact_arr[cancer_posit][sample_posit] = float(line[cytoact_posit])

    for i in range(0, len(cancerlist)) :
        if(cytoact_arr[i] != None) : cytoact_arr[len(cancerlist)].extend(cytoact_arr[i]) # pan cancer???곗씠?곕? ?⑹퀜以?        
    return;

getting_cytoact()
first = True

def getting_betavalue(normality_str, name, cancer_number, basic_posit) :
    filename = name + ".humanmethylation450." + normality_str + ".txt" # 媛곴컖??cancer name??? ???뚯씪紐낆쓣 臾몄옄?대줈 ?뺤꽦
    single = open(filename, 'r') # 媛곴컖??cancer??? ???곗씠???뚯씪??single_cancer?쇨퀬 吏 ??    
    sample_name = single.readline().split() # 泥ル쾲吏?以꾩? sample name?대?濡??쎄퀬 遺꾩젅

    del sample_name[0]; del sample_name[0] # sample name?  3踰덉㎏遺 ???쒖옉?섎?濡?1, 2踰덉㎏ 臾몄옄?댁? ?쒓굅
    single.readline() # 2踰덉㎏ line?  ?꾩슂 ?녿뒗 data?대?濡??쎄퀬 踰꾨━湲?    
    probedata = single.readlines() # 媛곴컖??probe??? ??beta-value瑜?2李⑥썝 諛곗뿴濡??듭㎏濡?湲곸뼱??        
    i = 0
    for line in probedata :
        line = line.split() # 諛쏆븘?ㅼ씤 媛믩뱾???쇱감??諛곗뿴 ?뺥깭濡??뺤젣
        if(first) : probe_name.append(line[0])
        del line[0]
        
        for j in range(0, len(line)) : # 2李⑥썝 諛곗뿴 ?뺥깭濡?sample紐낃낵 beta value瑜??숈떆???ｌ뼱以?
            
            sample_name[j] = sample_name[j][:15].replace('-', '')

            if(line[j] != "NA" and sample_name[j] in sample_id) :
                betavalue_arr[i][cancer_number][sample_id.index(sample_name[j])][basic_posit] = float(line[j])
        
        i += 1
        
    single.close()
    return;

import math

# calculate the mean
def mean(x):
    sum = 0.0
    for i in x: sum += i
    return sum / len(x) 

# calculates the sample standard deviation
def sampleStandardDeviation(x):
    sumv = 0.0
    for i in x:
        sumv += (i - mean(x))**2
    return math.sqrt(sumv/(len(x)-1))

# calculates the PCC using both the 2 functions above
def pearson(x, y):
    scorex = []
    scorey = []

    for i in x: 
        scorex.append((i - mean(x))/sampleStandardDeviation(x)) 

    for j in y:
        scorey.append((j - mean(y))/sampleStandardDeviation(y))

    # multiplies both lists together into 1 list (hence zip) and sums the whole list   
    return (sum([i*j for i,j in zip(scorex,scorey)]))/(len(x)-1)

print("CytAct_Completed")

can_num = 0
for cancer_name in cancerlist :
    getting_betavalue("tumor", cancer_name, can_num, 0)
    can_num += 1
    first = False
    
print("Cancer_Completed")
    
can_num = 0
for cancer_name in cancerlist :
    getting_betavalue("normal", cancer_name, can_num, 1)
    can_num += 1
    
print("Normal_Completed")

tumor_pearson_output = open("Tumor_Cor_CpGsite&CytAct_pearson.txt", 'w'); tumor_spearman_output = open("Tumor_Cor_CpGSite&CytAct_spearman.txt", 'w')
minus_pearson_output = open("minus_Cor_CpGsite&CytAct_pearson.txt", 'w'); minus_spearman_output = open("minus_Cor_CpGSite&CytAct_spearman.txt", 'w')
FC_pearson_output = open("FC_Cor_CpGsite&CytAct_pearson.txt", 'w'); FC_spearman_output = open("FC_Cor_CpGSite&CytAct_spearman.txt", 'w')

printline_header = "\t"
for cancer_name in cancerlist : printline_header += "%s\t" % cancer_name
printline_header += "pan_cancer\n"

tumor_pearson_output.write(printline_header); tumor_spearman_output.write(printline_header)
minus_pearson_output.write(printline_header); minus_spearman_output.write(printline_header)
FC_pearson_output.write(printline_header); FC_spearman_output.write(printline_header)

cancerlist.append("pan_cancer") # 遺꾩꽍??pan-cancer???ы븿?섎?濡?cancerlist??pan-cancer瑜??ｌ뼱以? 留??욎뿉 異쒕젰?쒗궗 寃껋씠誘 濡??욎뿉?ㅺ? ?뷀빐以?
cutoff = 0.3

for i in range(0, len(probe_name)) : # 紐⑤뱺 probe瑜?寃 ??   
    printline1 = "%s\t" % probe_name[i]; printline2 = "%s\t" % probe_name[i]; printline3 = "%s\t" % probe_name[i]; printline4 = "%s\t" % probe_name[i]; printline5 = "%s\t" % probe_name[i]; printline6 = "%s\t" % probe_name[i]
    
    #for j in range(0, len(cancerlist) - 1) :
    #    betavalue_arr[i][len(cancerlist) - 1] += betavalue_arr[i][j] # 媛곴컖??CpGSite betavalue??pan-cancer data瑜??삳뒗 怨쇱젙
    
    for j in range(0, len(cancerlist)) :
        
        tumor = [None, None]; correction_FC = [None, None]; correction_minus = [None, None]
        tumor[0] = list(); correction_FC[0] = list(); correction_minus[0] = list(); tumor[1] = list(); correction_FC[1] = list(); correction_minus[1] = list()
        check_tumor = False; check_minus = False; check_FC = False
        
        for k in range(0, len(sample_id)) :
            
            if(betavalue_arr[i][j][k][0] != None and cytoact_arr[j][k] != None) :
                
                betavalue_arr[i][len(cancerlist) - 1] = betavalue_arr[i][j]
                tumor[0].append(betavalue_arr[i][j][k][0]); tumor[1].append(cytoact_arr[j][k])
                check_tumor = True
                
                if(betavalue_arr[i][j][k][1] != None) :
                    if(abs(betavalue_arr[i][j][k][0] - betavalue[i][j][k][1]) > cutoff) : # tumor - normal betavalue > cutoff
                        correction_minus[0].append(betavalue_arr[i][j][k][0]); correction_minus[1].append(cytoact_arr[j][k])
                        check_minus = True
                    
                    if(betavalue_arr[i][j][k][0] > betavalue[i][j][k][1]) :
                        correction_FC[0].append(betavalue_arr[i][j][k][0]); correction_FC[1].append(cytoact_arr[j][k])
                        check_FC = True
        
        if(check_tumor) :
            tumor_pearson_pair = pearson(tumor[0], tumor[1]); printline1 += "%f\t" % tumor_pearson_pair
    #        tumor_pearson_pair = stats.pearsonr(tumor[0], tumor[1]); printline1 += "%f\t" % tumor_pearson_pair[0]
    #        tumor_spearman_pair = stats.spearmanr(tumor[0], tumor[1]); printline2 += "%f\t" % tumor_spearman_pair[0]
        else : printline1 += "NA\t"; printline2 += "NA\t"
            
        if(check_minus) :
            correction_minus_pearson_pair = pearson(correction_minus[0], correction_minus[1]); printline3 += "%f\t" % correction_minus_pearson_pair
    #        correction_minus_pearson_pair = stats.pearsonr(correction_minus[0], correction_minus[1]); printline3 += "%f\t" % correction_minus_pearson_pair[0]
    #        correction_minus_spearman_pair = stats.spearmanr(correction_minus[0], correction_minus[1]); printline4 += "%f\t" % correction_minus_spearman_pair[0]
        else : printline3 += "NA\t"; printline4 += "NA\t"
               
        if(check_FC) :
            correction_FC_pearson_pair = pearson(correction_FC[0], correction_FC[1]); printline5 += "%f\t" % correction_FC_pearson_pair
    #        correction_FC_pearson_pair = stats.pearsonr(correction_FC[0], correction_FC[1]); printline5 += "%f\t" % correction_FC_pearson_pair[0]
    #        correction_FC_spearman_pair = stats.spearmanr(correction_FC[0], correction_FC[1]); printline6 += "%f\t" % correction_FC_spearman_pair[0]
        else : printline5 += "NA\t"; printline6 += "NA\t"
    
    printline1 += "\n"; printline2 += "\n"; printline3 += "\n"; printline4 += "\n"; printline5 += "\n"; printline6 += "\n"
    
    tumor_pearson_output.write(printline1); #tumor_spearman_output.write(printline2)
    minus_pearson_output.write(printline3); #minus_spearman_output.write(printline4)                       
    FC_pearson_output.write(printline5); #FC_spearman_output.write(printline6)

print("END")
        
tumor_pearson_output.close(); tumor_spearman_output.close()     
minus_pearson_output.close(); minus_spearman_output.close()
FC_pearson_output.close(); FC_spearman_output.close()
