GeneMethyl 1.0.0
Written by Joon-Hyeong Park, Seoul National University, Applied Biology and Chemistry, clearclouds@snu.ac.kr.



A. Goal

	DNA methylation, adding methyl group to 5' carbon of a cytosine pyrimidine ring, is considered as a important biomarker for some diseases.
	And CpG sites, regions of a single DNA strand where a cytosine nucleotide is followed by a guanine nucleotide, are known as easily methylated in human body.
	Methylation of CpG sites can affect specific genes' expression levels and in some diseases like cancer, CpG sites' hypermethylation silence tumor suppressor genes' activities and hypomethylation promote retrotransposons' activities like LINE-1 to make chromosome instability.
	So, numerous clinical researches about methylation of CpG sites were already completed and still be underway.
	To make easily analyze relationship between DNA methylation status and target genes' activities, I made the packaged named GeneMethyl.
	I hope this package will contribute to many researches. And I plan to add more functions to solve more complicated problems.



	
B. Input data format

	Required Files : DNA methylation data, RNAseq data(raw counts or normalized by RSEM) of samples about a target disease.

	1. File name : [TargetDiseaseName].DNA_methylation_450K.tsv (Separated by tab)

		Column : Each sample
		Row : Each CpG site

		Missing values are denoted as nothing. (Just separted by tab)
		We recognize sample length up to 12. (excluding '-', or '_')

	Illumina's Infinium HumanMethylation450K Beadchip made us easily get DNA methylation data for selected CpG sites about 480,000.
	So, this package use beta-value data table generated from HumanMethylation450K.
	And the format of a DNA methylation dataset is the same with a TCGA's pan-cancer atlas methylation dataset.

		

	2. File name : [TargetDiseaseName].PANCANCER.RNAseq.tsv (Separated by tab)

		Column : Each sample
		Row : Each gene

		Missing values are denoted as nothing. (Just separted by tab)
		We recognize sample length up to 12. (excluding '-', or '_')

	Your target gene expression levels are automatically calculated from RNAseq raw counts data or normalized data by RSEM package.
	And the format of a RNAseq dataset is the same with a TCGA's pan-cancer atlas RNAseq dataset.





C. Package functions

	1. BetavalueDistribution

		(1) BetavalueDistribution.Draw("TargetDiseaseName", Cutoff, WhetherHistogram) : Drawing beta-value distribution by density plot based on histogram.

			Input parameters : 

				"TargetDiseaseName" : Target disease name
				Cutoff : Cutoff must be in [0, 1] and sections are divided by this cutoff.
					ex) Cutoff : 0.1 	-> 	Your Sections : 0~0.1, 0.1~0.2, 0.2~0.3, 0.3~0.4, 0.4~0.5, 0.5~0.6, 0.6~0.7, 0.7~0.8, 0.8~0.9, 0.9~1.0
				WhetherHistogram : True or False. If you choose True, histogram of beta-values divided by cutoff is also shown.


			Description : 

				DNA methylation data is too large to easily handle, so receiving whole DNA methylation data and drawing density plot at once is not effective from the perspective of memory.
				So, I received DNA methylation data line by line and approximated DNA methylation beta-values into several sections divided by specific cutoff. (You can choose cutoff)
				Then by using histogram, I made density plot of beta-values' distribution.
	

			Output file :

				/Result/DistributionPlot/[TargetDiseaseName].Betavalue.Distribution.Plot.pdf





	2. SimpleCutoff

		(1) SimpleCutoff.TargetGeneActivity("TargetDiseaseName", [TargetGenesList]) : Calculating target genes' activities.

			Input parameters : 

				"TargetDiseaseName" : Target disease name
				[TargetGenesList] : Target genes' list. You can choose multiple genes.

		
			Description : 

				You can simply calculate target genes' activities.
				I decided representative target genes' activity by using logarithm.
					base : a number of gene
					anti-logarithm : geometric mean of target genes' RNAseq data added to pseudocount(1).
					cf) I added whole RNAseq data to pseudocount(1) to prevent from minus value of representative target genes' activity.


			Output file :
				
				/Result/SimpleCutoff/[TargetDiseaseName].TargetGeneActivity.txt




		(2) SimpleCutoff.View_Correlation_AND_ScatterPlot("TargetDiseaseName", [TargetGenesList], [Cutoff], Type, WhetherFoldChange) : Calculating sperman's correlation between representative target genes' activities and summations of whole samples. Drawing scatter plots of this correlation.

			Input parameters : 

				"TargetDiseaseName" : Target disease name
				[TargetGenesList] : Target genes' list. You can choose multiple genes.
				Cutoff, Type : Calculating summations by using Cutoff depending on Type.
					Type : "Lower", "Higher", "Both", "All"
						"Lower" : If beta-value is lower than cutoff, sample's beta-value is converted into 1. Or, into 0. Then summate this values to each sample. -> You can get the number of beta-values lower than cutoff to each sample. (in [0, cutoff])
						"Higher" : If beta-value is higher than cutoff, sample's beta-value is converted into 1. Or, into 0. Then summate this values to each sample. -> You can get the number of beta-values higher than cutoff to each sample. (in [cutoff, 1])
						"Both" : If beta-value is lower than cutoff/2 or higher than 1 - cutoff/2, sample's beta-value is converted into 1. Or, into 0. Then summate this values to each sample. -> You can get the number of beta-values in [0, cutoff/2] or [1 - cutoff/2, 1] to each sample.
						"All" : Doing all of theses types respectively.
					Cutoff : Cutoff must be in [0, 1] and it determine the section. You can choose multiple cutoffes.
				WhetherFoldChange : Calculating the fold change
					Fold change = Mean or Median of the representative target genes' activity of samples included in the section by cutoff / Mean or Median of the representative target genes' activity of samples NOT included in the section by cutoff

			Output file :

				/Result/SimpleCutoff/FC_CpGsites/WholeSites.Cutoff.[Cutoff].[TargetDiseaseName].[Type].FC.CpGsites.txt
				/Result/SimpleCutoff/Summation/WholeSites.Cutoff.[Cutoff].[TargetDiseaseName].[Type].Binarization.Summation.txt
				/Result/SimpleCutoff/Correlation/WholeSites.[TargetDiseaseName].[Type].Correlation.Summation.And.TargetGeneActivity.txt -> including whole cutoffes to compare easily
				/Result/SimpleCutoff/Correlation/WholeSites.[TargetDiseaseName].CompareAll.Correlation.Summation.And.TargetGeneActivity.txt -> only emerging if Type is All to compare easily
				/Result/SimpleCutoff/ScatterPlot/WholeSites.Cutoff.[Cutoff].[TargetDiseaseName].[Type].ScatterPlot.pdf
	




	3. TopPercentageCutoff

		(1) TopPercentageCutoff.TargetGeneActivity("TargetDiseaseName", [TargetGenesList]) : Calculating target genes' activities.

			Input parameters : 

				"TargetDiseaseName" : Target disease name
				[TargetGenesList] : Target genes' list. You can choose multiple genes.


			Description : 

				You can simply calculate target genes' activities.
				I decided representative target genes' activity by using logarithm.
					base : a number of gene
					anti-logarithm : geometric mean of target genes' RNAseq data added to pseudocount(1).
					cf) I added whole RNAseq data to pseudocount(1) to prevent from minus value of representative target genes' activity.
	

			Output file :
				
				/Result/TopNpercentageCutoff/[TargetDiseaseName].TargetGeneActivity.txt




		(2) TopPercentageCutoff.View_Correlation_AND_ScatterPlot("TargetDiseaseName", [TargetGenesList], [Percentage], Type, WhetherMeanMethod, WhetherFoldChange) : Calculating sperman's correlation between representative target genes' activities and summations of whole samples. Drawing scatter plots of this correlation.

			Background : 

				Before explaining TopNpercentageCutoff, we have to determine what percentage means in this method.

				In some diseases like cancer, beta-value distribution of each CpG site is look like roughly 2-peaked graph.
				So, I roughly classified beta-value distribution of each CpG site as 2 categories, left-skewed and right-skewed.
				To determine skewedness of CpG sites, we compare median of each CpG site with mean or 0.5( = an exact half of [0, 1]).
					If you want to use 'mean method', you need to make WhetherMeanMethod True.
					If you want to use '0.5 method', you need to make WhetherMeanMethod False.

				After determining skewedness, we can classify types as "Positive", "Negative", "Both".
					If type is "Positive" and CpG site's beta-value distribution is right-skewed, top N% of sample beta-value is converted into 1. Or, into 0.
					If type is "Positive" and CpG site's beta-value distribution is left-skewed, bottom N% of sample beta-value is converted into 1. Or into 0.
					Thus, "Positive" type means we count the number of beta-values following the tendency of each CpG site's distribution for each sample.

					If type is "Negative" and CpG site's beta-value distribution is right-skewed, bottom N% of sample beta-value is converted into 1. Or, into 0.
					If type is "Negative" and CpG site's beta-value distribution is left-skewed, top N% of sample beta-value is converted into 1. Or into 0.
					Thus, "Negative" type means we count the number of beta-values not following the tendency of each CpG site's distribution for each sample.

					If type is "Both", regardless of CpG site's skewedness top N/2% and bottom N/2% of sample beta-value is converted into 1. Or, into 0.
					Thus, "Both" type means we count the number of strongly methylated or demethylated CpG sites for each sample.

				By using this method, we can simply count the number of hypermethylated or hypomethylated CpG sites for the specific situation.
				To explain this, let's take an example of cancer.
					Globally, CpG sites are hypomethylated to increase chromosomal instability, by expressing retrotranspons(just one example).
					But, CpG sites are hypermethylated near the promoter regions to silence life-critical genes.
					So, If we use TopNpercentageCutoff Positive type method, we can count the number of these hypermethylated or hypomethylated CpG sites.
					Then, we can correlate this values with the representative target genes' activity.

				In short, we can analyze the epigenetic impact of diseases.


			Input parameters :

				"TargetDiseaseName" : Target disease name
				[TargetGenesList] : Target genes' list. You can choose multiple genes.
				Percentage, Type : Calculating summations by using Percentage depending on Type.
					Type : "Positive", "Negative", "Both", "All"
					Percentage : What percentage do you want?
				WhetherFoldChange : Calculating the fold change
					Fold change = Mean or Median of the representative target genes' activity of samples included in the section by percentage / Mean or Median of the representative target genes' activity of samples NOT included in the section by percentage.
				WhetherMeanMethod : Choosing the method of determining skewedness


			Output file :

				/Result/TopPercentageCutoff/FC_CpGsites/WholeSites.Percentage.[Percentage].[TargetDiseaseName].[Type].FC.CpGsites.txt
				/Result/TopPercentageCutoff/Summation/WholeSites.Percentage.[Percentage].[TargetDiseaseName].[Type].Binarization.Summation.txt
				/Result/TopPercentageCutoff/Skewed/WholeSites.[TargetDiseaseName].Left.Skewed.CpGsites.txt
				/Result/TopPercentageCutoff/Skewed/WholeSites.[TargetDiseaseName].Right.Skewed.CpGsites.txt
				/Result/TopPercentageCutoff/Correlation/WholeSites.[TargetDiseaseName].[Type].Correlation.Summation.And.TargetGeneActivity.txt -> including whole percentages to compare easily
				/Result/TopPercentageCutoff/Correlation/WholeSites.[TargetDiseaseName].CompareAll.Correlation.Summation.And.TargetGeneActivity.txt -> only emerging if Type is All to compare easily
				/Result/TopPercentageCutoff/ScatterPlot/WholeSites.Percentage.[Percentage].[TargetDiseaseName].[Type].ScatterPlot.pdf





D. Usage Example

	from GeneMethyl import *

	BetavalueDistribution.DrawDensityPlot("PANCANCER", 0.001, True)
	
	SimpleCutoff.TargetGeneActivity("PANCANCER", ["GZMA", "PRF1"])
	SimpleCutoff.View_Correlation_AND_ScatterPlot("PANCANCER", ["GZMA", "PRF1"], [0.1, 0.2], "All", True)
		cf) Even before not calling TargetGeneActivity, it will automatically executed in this function.
	
	TopPercentageCutoff.TargetGeneActivity("PANCANCER", ["GZMA", "PRF1"])
	TopPercentageCutoff.View_Correlation_AND_ScatterPlot("PANCANCER", ["GZMA", "PRF1"], [0.05, 0.1, 0.15, 0.2], "All", False, True)
		cf) Even before not calling TargetGeneActivity, it will automatically executed in this function.
