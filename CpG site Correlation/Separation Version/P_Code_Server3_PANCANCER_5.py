# -*- coding: utf-8 -*-
cancerlist = ["PANCANCER"] # server1

from operator import itemgetter
from scipy import stats
import numpy as np

betavalue_arr = []
cytoact_arr = []

probe_name = []
sample_id = []

start_number = 150000

SEP_NAME = ".SEP_3.5."

probe_count = 6000
sample_count = 0

probe_separation_number = 1000
probe_iteration = 0

######################################################################################################################################################

def getting_cytoact() :
    cytoact_file = open("TCGA_methylation_cowork_1.txt", 'r')
    header = cytoact_file.readline().split() # header 읽기

    id_posit = header.index("id") # sample ID positioning
    cytoact_posit = header.index("CytAct") # CytAct positioning
    cytodata = cytoact_file.readlines() # 데이터 테이블 통째로 읽어들임
    cytoact_file.close()
    
    for line in cytodata :
        line = line.split()
        sample_id.append(line[id_posit].replace('_', '')) # sample ID 추출 (주형으로 사용할 것)
    
    sample_count = len(sample_id)
    
    for i in range(0, sample_count) : cytoact_arr.append(None) # CytAct value table 초기화

    for line in cytodata :
        line = line.split() # 1 sample data를 분절해서 CytAct value 추출하기 위함
        
        if(line[cytoact_posit] != "NA") : # CytAct value가 결측치가 아니라면
            sample_posit = sample_id.index(line[id_posit].replace('_', ''))
            cytoact_arr[sample_posit] = float(line[cytoact_posit]) # 저장한다
    return;

######################################################################################################################################################

getting_cytoact()
print("CytAct_Completed")

######################################################################################################################################################

def reset_betavalue() :
    
    del betavalue_arr[:]
    for reset_x in range(0, probe_separation_number) : betavalue_arr.append({})
    
    return

######################################################################################################################################################

import math
def mean(x):
    sum = 0.0
    for i in x: sum += i
    return sum / len(x) 

def sampleStandardDeviation(x):
    sumv = 0.0
    for i in x:
        sumv += (i - mean(x)) ** 2
    return math.sqrt(sumv / (len(x) - 1))

def pearson(x, y):
    scorex = []
    scorey = []

    for i in x: scorex.append((i - mean(x)) / sampleStandardDeviation(x)) 
    for j in y: scorey.append((j - mean(y)) / sampleStandardDeviation(y))

    return (sum([i * j for i, j in zip(scorex, scorey)])) / (len(x) - 1)

######################################################################################################################################################

def getting_betavalue(name) :
    
    tumor_pearson_output = open(name + SEP_NAME + "Tumor_Cor_CpGsite&CytAct_pearson.txt", 'w'); tumor_spearman_output = open(name + SEP_NAME + ".Tumor_Cor_CpGSite&CytAct_spearman.txt", 'w')
    tumor_pearson_output.write("CpGsite\t%s\tP_value\n" % name); tumor_spearman_output.write("CpGsite\t%s\tP_value\n" % name)
    
    filename1 = name + ".humanmethylation450.tumor.txt" # cancer name별로 파일명이 다름을 고려해줌
#    filename2 = name + ".humanmethylation450.normal.txt" # cancer name별로 파일명이 다름을 고려해줌
    input_file1 = open(filename1, 'r')
#    input_file2 = open(filename2, 'r')
    
    sample_name1= input_file1.readline().split() # header에 sample ID가 있으므로 별도로 읽어준 후
#    sample_name2 = input_file2.readline().split()

    for separation_count in range(0, start_number) : input_file1.readline()
    
    del sample_name1[0]; del sample_name1[0]
#    del sample_name2[0]; del sample_name2[0] # 쓰레기 데이터를 제외
    
    input_file1.readline() # betavalue임을 명시하는 row를 읽고 폐기
#    input_file2.readline()
   
    probe_separation_number_copy = probe_separation_number

    i = 0
    while i < probe_count :
        
        ##################################################################################################################################################
        del probe_name[:]
        escape = False
        
        for normal_iteration in range(0, probe_separation_number_copy) : 
            line1 = input_file1.readline().split() # 한 probe에 대한 여러 sample ID의 betavalue를 읽어들임
#            line2 = input_file2.readline().split()
            
            if(len(line1) == 0) : # 끝까지 읽은 경우 루프 문을 빠져나감
                escape = True; probe_separation_number_copy = normal_iteration
                break          
            
#            while line1[0] != line2[0] :  # 관심 있는 probe와 같아질 때까지 normal data를 훑어 내려감
#                ############ 관심 없는 probe이므로 출력해준 뒤 Not interested 표시
#                tumor_pearson_output.write("%s\tNot interested\n" % line2[0]); tumor_spearman_output.write("%s\tNot interested\n" % line2[0])
#                line2 = input_file2.readline().split() # 새로운 probe를 읽어들임

            probe_name.append(line1[0]) # CpG site name 추출한 뒤, betavalue가 아니므로 리스트에서 제거
            del line1[0]; #del line2[0]
        
        ##################################################################################################################################################
            
            for j in range(0, len(line1)) : # sample ID의 개수만큼 반복함
                sample_name1[j] = sample_name1[j][:15].replace('-', '') # sample ID의 형식을 통일해줌

                if(line1[j] != "NA" and sample_name1[j] in sample_id) : # 결측치가 아니고, sample id에 포함된 경우
                    betavalue_arr[normal_iteration][sample_name1[j]] = [float(line1[j]), None] # sample ID에 맞는 index에 betavalue 저장
        
        ##################################################################################################################################################
                    
#            for j in range(0, len(line2)) : # sample ID의 개수만큼 반복함
#                sample_name2[j] = sample_name2[j][:15].replace('-', '') # sample ID의 형식을 통일해줌

#                if(line2[j] != "NA" and sample_name2[j] in sample_name1) : # 결측치가 아니고, sample id에 포함된 경우
#                    betavalue_arr[normal_iteration][sample_name2[j]][1] = float(line2[j]) # sample ID에 맞는 index에 betavalue 저장
                
        ##################################################################################################################################################
        
        for j in range(0, probe_separation_number_copy) :
            printline1 = "%s\t" % probe_name[j]; printline2 = "%s\t" % probe_name[j]
        
            tumor = [None, None]; correction_FC = [None, None]; correction_minus = [None, None]
            tumor[0] = list();  tumor[1] = list()
            check_tumor = False

            iteration_number = len(betavalue_arr[j])
            each_site_sample = betavalue_arr[j].items()
            
            for k in range(0, iteration_number) :
                
                if(each_site_sample[k][1][0] != None) :
                
                    tumor[0].append(float(each_site_sample[k][1][0])); tumor[1].append(float(cytoact_arr[sample_id.index(each_site_sample[k][0])]))
                    check_tumor = True
                    
                
            if(check_tumor) :
                tumor_pearson_pair = stats.pearsonr(tumor[0], tumor[1]); tumor_spearman_pair = stats.spearmanr(tumor[0], tumor[1])
                printline1 += "%f\t%.3f" % (tumor_pearson_pair[0], tumor_pearson_pair[1]); printline2 += "%f\t%.3f" % (tumor_spearman_pair[0], tumor_spearman_pair[1])
            
            else : printline1 += "NA\tNA"; printline2 += "NA\tNA"
                                    
            printline1 += "\n"; printline2 += "\n"
        
            tumor_pearson_output.write(printline1)
            tumor_spearman_output.write(printline2)
        
        if(escape) : break
            
        ##################################################################################################################################################
        
        i += probe_separation_number

        print(i, probe_separation_number, "We are processing %s" % SEP_NAME)
        
    input_file1.close(); #input_file2.close()
    
    tumor_pearson_output.close(); tumor_spearman_output.close()
    
    return

def process(cancer_name) :
    
    reset_betavalue()
    getting_betavalue(cancer_name)
    
    return

for cancer_name in cancerlist :
    process(cancer_name); print(cancer_name + " completed")

print("END")
