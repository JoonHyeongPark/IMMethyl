#-*- coding: utf-8 -*-
cancerlist = ["ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", "ESCASTAD", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM"]
cancer_number = 0

betavalue_arr = {}
cytoact_arr = {}

for cancer_name in cancerlist :
    filename = cancer_name + ".humanmethylation450.tumor.txt" # 각각의 cancer name에 대한 파일명을 문자열로 형성
    single_cancer = open(filename, 'r') # 각각의 cancer에 대한 데이터 파일을 single_cancer라고 지정
            
    sample_name = single_cancer.readline().split() # 첫번째 줄은 sample name이므로 읽고 분절
    del sample_name[0]; del sample_name[0] # sample name은 3번째부터 시작하므로 1, 2번째 문자열은 제거
    single_cancer.readline() # 2번째 line은 필요 없는 data이므로 읽고 버리기
                      
    probedata = single_cancer.readlines() # 각각의 probe에 대한 beta-value를 2차원 배열로 통째로 긁어냄
                            
    for line in probedata :
        line = line.split() # 받아들인 값들을 일차원 배열 형태로 정제
        probe_name = line.pop(0) # site만 문자열이므로 별도로 저장하고 삭제
                                                    
        betavalue_arr[probe_name] = { cancer_name : list() } # betavalue dictionary에서 각각의 CpGSite에 대한 Cancer type별 data를 list 형태로 넣겠다
                                                            
        for i in range(0, len(line)) : # 2차원 배열 형태로 sample명과 beta value를 동시에 넣어줌
            if(line[i] != "NA") : betavalue_arr[probe_name][cancer_name].append([sample_name[i], float(line[i])])                                                                                       

single_cancer.close()

cytoact_file = open("TCGA_methylation_cowork_1.txt", 'r')

header = cytoact_file.readline().split() # header 읽기

id_posit = header.index("tissue.code") # sample ID 위치 인덱싱
cancer_posit= header.index("origin") # cancer type 위치 인덱싱
cytoact_posit = header.index("CytAct") # cytact 위치 인덱싱

cytodata = cytoact_file.readlines() # 관련 data를 통째로 읽어들임
cytoact_file.close()

for cancer_name in cancerlist : cytoact_arr[cancer_name] = list() # 미리 cytoact dictionary에서 cancer type별 cytact data를 list 형태로 저장할 것을 명시
cytoact_arr["pan_cancer"] = list() # pan-cancer 형태도 추가해줌

for line in cytodata :
    line = line.split() # sample 별로 정리된 data를 읽어들이고 분절
            
    if(line[cytoact_posit] != "NA") :
        sample_name = line[id_posit] # id 위치로부터 sample 명을 읽어들인 후 통일된 양식으로 가공
        sample_name = sample_name.replace('.', '-')
                                
        cancer_name = line[cancer_posit] # cancer type 위치로부터 sample 명을 읽어들인 후 통일된 양식으로 가공
        cancer_name = cancer_name.replace('_', '')
                                                    
        cytoact_arr[cancer_name].append([sample_name, float(line[cytoact_posit])])

from operator import itemgetter
from scipy import stats
import pandas as pd

for cancer_name in cancerlist :
    if(cytoact_arr[cancer_name] != None) : cytoact_arr["pan_cancer"] += cytoact_arr[cancer_name] # pan cancer에 데이터를 합쳐줌

pearson_output = open("Cor_CpGsite&CytAct_pearson.txt", 'w')
spearman_output = open("Cor_CpGSite&CytAct_spearman.txt", 'w')

allprobe = betavalue_arr.keys() # CpgSite의 이름만 별도로 추출해서 배열 형태로 저장
allprobe.sort() # CpgSite를 이름 순서대로 정렬

printline_header = "\tpan_cancer\t"
for cancer_name in cancerlist : printline_header += "%s\t" % cancer_name

printline_header += "\n"
                                                   
pearson_output.write(printline_header)
spearman_output.write(printline_header)

cancerlist.insert(0, "pan_cancer") # 분석에 pan-cancer도 포함되므로 cancerlist에 pan-cancer를 넣어줌, 맨 앞에 출력시킬 것이므로 앞에다가 더해줌

for probe_name in allprobe : # 모든 probe를 검색    
    printline1 = printline2 = "%s\t" % probe_name

    betavalue_arr[probe_name]["pan_cancer"] = list() # betavalue의 pan-cancer data를 합치기 위한 예비작업
                
    for cancer_name in cancerlist : 
        if(cancer_name != "pan_cancer") :
            betavalue_arr[probe_name]["pan_cancer"] += betavalue_arr[probe_name][cancer_name] # 각각의 CpGSite betavalue의 pan-cancer data를 얻는 과정
                                        
    for cancer_name in cancerlist :
        if(len(betavalue_arr[probe_name][cancer_name]) > 0 and len(cytoact_arr[cancer_name]) > 0) :
            sample1, beta = zip(*betavalue_arr[probe_name][cancer_name])
            sample2, cyt = zip(*cytoact_arr[cancer_name])
                                                                        
            beta_table = pd.DataFrame({"ID" : sample1, "value" : beta})
            cyt_table = pd.DataFrame({"ID" : sample2, "value" : cyt})
                                                                                                
            total = pd.merge(beta_table, cyt_table, on = "ID", how = "inner")
                                                                                                            
            pearson_pair = stats.pearsonr(total["value_x"], total["value_y"]); printline1 += "%f\t" % pearson_pair[0]
            spearman_pair = stats.spearmanr(total["value_x"], total["value_y"]); printline2 += "%f\t" % spearman_pair[0]
        else : printline1 += "NA\t"; printline2 += "NA\t"
                                                                                                                                        
    printline1 += "\n"; printline2 += "\n"
                                                                                                                        
    pearson_output.write(printline1)
    spearman_output.write(printline2)

pearson_output.close()
spearman_output.close()
