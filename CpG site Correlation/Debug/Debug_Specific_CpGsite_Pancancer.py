# -*- coding: utf-8 -*-

from operator import itemgetter
from scipy import stats
import numpy as np

betavalue_arr = []
cytoact_arr = []

probe_name = []
sample_id = []


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

output = open("debug2.txt", 'w')
    
filename1 = open("PANCANCER.humanmethylation450.tumor.txt", 'r') # cancer name별로 파일명이 다름을 고려해줌

sample_name = filename1.readline().split(); filename1.readline()
del sample_name[0]; del sample_name[0]

now_target = filename1.readline().split()
probe_name = now_target.pop(0)

output.write("%s\n" % probe_name)

column1 = []
column2 = []

for i in range(0, len(sample_name)) :
    
    sample_name[i] = sample_name[i][:15].replace('-', '')
    
    if(sample_name[i] in sample_id and now_target[i] != "NA") :

        posit = sample_id.index(sample_name[i])
        printline = "%s\t%s\t%s\n" % (sample_name[i], now_target[i], cytoact_arr[posit])
        column1.append(float(now_target[i]))
        column2.append(float(cytoact_arr[posit]))
        output.write(printline)

cor = stats.spearmanr(column1, column2)
lastprint = "%f\t%f\n" % (cor[0], cor[1])

output.write(lastprint)

output.close()

print("END")
