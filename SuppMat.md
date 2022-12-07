# Supplementary material

1. [Assessment of depth and missingness in the final dataset](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/blob/main/SuppMat.md#1-assessment-of-depth-and-missingness-in-the-final-dataset)

2. [Is there a difference in SNPs sequenced per population or per library?](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/blob/main/SuppMat.md#2-is-there-a-difference-in-snps-sequenced-per-population-or-per-library)

3. [Does a higher genotyping rate change the results and conclusions in our MS?](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/blob/main/SuppMat.md#3-analysis-of-smaller-dataset-with-higher-genotyping-rate-95)


## 1. Assessment of depth and missingness in the final dataset

The final filtered dataset comprised 251 individuals genotyped at 61210 loci. 

Estimate depth and missingness per individual and plot per population: 
```
export PATH=/share/apps/genomics/vcftools-0.1.16/bin:$PATH

vcftools --vcf AA251.FINAL.MAF0.01.missing0.5perpop.vcf --depth
vcftools --vcf AA251.FINAL.MAF0.01.missing0.5perpop.vcf --missing-indv
```

Read into R and determine median and range of depth. Plot per population ranges. 
```
#R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
#Copyright (C) 2022 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin17.0 (64-bit)

library(ggplot2)

depth <- read.table("out.idepth", header=T)
miss <- read.table("out.imiss", header=T)

head(depth)
         INDV N_SITES MEAN_DEPTH
1 BAR_10_2013   61210    410.177
2 BAR_11_2014   61210    246.913
3 BAR_12_2013   61210    238.969
4 BAR_13_2014   61210    255.450
5 BAR_14_2013   61210    147.284
6 BAR_14_2014   61210    382.788

head(miss)
         INDV N_DATA N_GENOTYPES_FILTERED N_MISS    F_MISS
1 BAR_10_2013  61210                    0   3000 0.0490116
2 BAR_11_2014  61210                    0   4033 0.0658879
3 BAR_12_2013  61210                    0   2520 0.0411697
4 BAR_13_2014  61210                    0   2992 0.0488809
5 BAR_14_2013  61210                    0   3552 0.0580297
6 BAR_14_2014  61210                    0   1270 0.0207482

summary(depth)
     INDV              N_SITES        MEAN_DEPTH     
 Length:251         Min.   :61210   Min.   :  9.604  
 Class :character   1st Qu.:61210   1st Qu.:161.687  
 Mode  :character   Median :61210   Median :215.247  
                    Mean   :61210   Mean   :235.247  
                    3rd Qu.:61210   3rd Qu.:263.539  
                    Max.   :61210   Max.   :812.052
		    
#There are some low depth individuals. Let's see who they are: 		    
depth[depth$MEAN_DEPTH<20,]

INDV N_SITES MEAN_DEPTH
146 LYD_34_2014   61210   15.49270
211  SWD_7_2014   61210    9.60392
221 WIS_20_2014   61210   13.29600

#We decided to keep SWD_7_2014 because it's mean depth is close to 10x, our lower threshold. The rest of the individuals all have a mean depth >20x
#We can visualise this: 
#First get pop info from the out.idepth file in linux
awk -F "_" '{print $1}' > pop

#Read into R
pop <- read.table("pop", header=F)
depth$pop <- pop$V1

#plot
pdf("AA251.depth_perpop.pdf")
ggplot(depth, aes(y=MEAN_DEPTH, x=pop))+geom_boxplot()+geom_point(color="black", size=1, alpha=0.9)
dev.off()

```

### Mean depth per individual, with boxplots grouped per population. 

![alt_txt][out.idepth]

[out.idepth]:https://user-images.githubusercontent.com/12142475/206146870-26fb5731-e093-4697-820d-044c896284f6.png



```
##Missingness per population

summary(miss)
     INDV               N_DATA      N_GENOTYPES_FILTERED     N_MISS     
 Length:251         Min.   :61210   Min.   :0            Min.   :  562  
 Class :character   1st Qu.:61210   1st Qu.:0            1st Qu.: 2410  
 Mode  :character   Median :61210   Median :0            Median : 3187  
                    Mean   :61210   Mean   :0            Mean   : 4781  
                    3rd Qu.:61210   3rd Qu.:0            3rd Qu.: 4055  
                    Max.   :61210   Max.   :0            Max.   :32982  
     F_MISS        
 Min.   :0.009182  
 1st Qu.:0.039381  
 Median :0.052067  
 Mean   :0.078108  
 3rd Qu.:0.066247  
 Max.   :0.538834  


## We can see that the majority of the samples have a very low missingness rate, although a few samples have >10% missingness.
## Let's identify all individuals with >20% missingness (the usual allowed missingness in pop gen papers): 

dim(miss[miss$F_MISS>0.2,])
[1] 19  5

miss[miss$F_MISS>0.2,]

INDV N_DATA N_GENOTYPES_FILTERED N_MISS   F_MISS
15  BAR_22_2014  61210                    0  32982 0.538834
16  BAR_24_2014  61210                    0  13409 0.219066
25   BAR_6_2014  61210                    0  15772 0.257670
71  BRO_18_2013  61210                    0  12464 0.203627
72  BRO_19_2013  61210                    0  21061 0.344078
91  FOR_22_2014  61210                    0  30230 0.493874
93  FOR_24_2014  61210                    0  16141 0.263699
94  FOR_25_2014  61210                    0  28421 0.464320
95  FOR_26_2014  61210                    0  30387 0.496438
96  FOR_29_2014  61210                    0  14856 0.242705
97  FOR_30_2014  61210                    0  28172 0.460252
109 HOD_20_2014  61210                    0  21605 0.352965
116 HOD_34_2014  61210                    0  12270 0.200457
119 HOD_38_2014  61210                    0  14439 0.235893
120 HOD_39_2014  61210                    0  15041 0.245728
139 LYD_25_2014  61210                    0  18262 0.298350
143 LYD_31_2014  61210                    0  27881 0.455497
146 LYD_34_2014  61210                    0  29949 0.489283
211  SWD_7_2014  61210                    0  31076 0.507695

## LYD_34_2014 and SWD_7_2014 are again the worst two samples. 

## Let's plot this per population to see if any particular population is worse off: 

miss$pop <- pop$V1  #add population information

#plot
pdf("AA251.miss_perpop.pdf")
ggplot(miss, aes(y=F_MISS, x=pop))+geom_boxplot()+geom_point(color="black", size=1, alpha=0.9)
dev.off()

```
### The proportion of missing data per individual, grouped by population. 

![alt_txt][miss]

[miss]:https://user-images.githubusercontent.com/12142475/206148905-8bfe17bf-fda3-43a2-91b4-e803240f5d0c.png





# 2. Is there a difference in SNPs sequenced per population or per library?

First, we have a look at the number of private alleles sequenced in each library that are present in our final dataset. ie, does "Library" explain a difference in the variants present in the final dataset: 

Split the vcf files to include only individuals sequenced in each library. The pop files we used below can be found [here](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/tree/main/Files/PopFiles_Splitvcfs)
```
export PATH=/share/apps/genomics/vcftools-0.1.16/bin:$PATH

for i in $(ls Lib*names); do vcftools --vcf AA251.FINAL.MAF0.01.missing0.5perpop.vcf --keep $i --recode --recode-INFO-all --out $i.subset; done
```

And create vcf files that exclude each of these populations in turn (i.e., leave-one-out datasets)
```
for i in $(ls Lib*names); do vcftools --vcf AA251.FINAL.MAF0.01.missing0.5perpop.vcf --remove $i --recode --recode-INFO-all --out No_$i.subset; done
```


Next, we'd like to determine the number of private alleles in each of these libraries. We'll compare the loci called in a library that isn't found in the rest of the dataset: 
```
export PATH=/share/apps/genomics/bcftools-1.15/bin:$PATH
#bcf requires bgzipped and indexed vcf files: 
for i in $(ls Lib*vcf); do bcftools view $i -Oz -o $i.gz; done
for i in $(ls Lib*vcf.gz); do bcftools index $i; done

for i in $(ls No*vcf); do bcftools view $i -Oz -o $i.gz; done
for i in $(ls No*vcf.gz); do bcftools index $i; done

#Create a file containing Lib1-Lib6, one per line
seq -f "Lib%g" 1 6 > Lib.names

#And run each in turn: 
for i in $(cat Lib.names); do bcftools isec -p ${i}_pvtalleles No_${i}.names.subset.recode.vcf.gz ${i}.names.subset.recode.vcf.gz; done 

#Check which file we're interested in: 
cat Lib1_pvtalleles/README.txt

This file was produced by vcfisec.
The command line was:	bcftools isec  -p Lib1_pvtalleles No_Lib1.names.subset.recode.vcf.gz Lib1.names.subset.recode.vcf.gz

Using the following file names:
Lib1_pvtalleles/0000.vcf	for records private to	No_Lib1.names.subset.recode.vcf.gz
Lib1_pvtalleles/0001.vcf	for records private to	Lib1.names.subset.recode.vcf.gz
Lib1_pvtalleles/0002.vcf	for records from No_Lib1.names.subset.recode.vcf.gz shared by both	No_Lib1.names.subset.recode.vcf.gz Lib1.names.subset.recode.vcf.gz
Lib1_pvtalleles/0003.vcf	for records from Lib1.names.subset.recode.vcf.gz shared by both	No_Lib1.names.subset.recode.vcf.gz Lib1.names.subset.recode.vcf.gz
```


Count the number of private alleles in each of the libraries: 
```
#We're interested in private alleles within each library. These are written to 0001.vcf in each folder: 

for i in $(cat Lib.names); do vcftools --vcf ${i}_pvtalleles/0001.vcf; done

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib1_pvtalleles/0001.vcf

After filtering, kept 45 out of 45 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib2_pvtalleles/0001.vcf

After filtering, kept 44 out of 44 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib3_pvtalleles/0001.vcf

After filtering, kept 43 out of 43 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib4_pvtalleles/0001.vcf

After filtering, kept 42 out of 42 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib5_pvtalleles/0001.vcf

After filtering, kept 41 out of 41 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib6_pvtalleles/0001.vcf

After filtering, kept 36 out of 36 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds


### We're also interested in whether there are any loci that don't appear in a particular library. These are written to 0000.vcf
for i in $(cat Lib.names); do vcftools --vcf ${i}_pvtalleles/0000.vcf; done

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib1_pvtalleles/0000.vcf

After filtering, kept 206 out of 206 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib2_pvtalleles/0000.vcf

After filtering, kept 207 out of 207 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib3_pvtalleles/0000.vcf

After filtering, kept 208 out of 208 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib4_pvtalleles/0000.vcf

After filtering, kept 209 out of 209 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib5_pvtalleles/0000.vcf

After filtering, kept 210 out of 210 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib6_pvtalleles/0000.vcf

After filtering, kept 215 out of 215 Individuals
After filtering, kept 0 out of a possible 0 Sites
File does not contain any sites
Run Time = 0.00 seconds


### For completeness we'll check how many loci are shared between the datasets in each case: 0002.vcf
for i in $(cat Lib.names); do vcftools --vcf ${i}_pvtalleles/0002.vcf; done

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib1_pvtalleles/0002.vcf

After filtering, kept 206 out of 206 Individuals
After filtering, kept 61210 out of a possible 61210 Sites
Run Time = 1.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib2_pvtalleles/0002.vcf

After filtering, kept 207 out of 207 Individuals
After filtering, kept 61210 out of a possible 61210 Sites
Run Time = 1.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib3_pvtalleles/0002.vcf

After filtering, kept 208 out of 208 Individuals
After filtering, kept 61210 out of a possible 61210 Sites
Run Time = 0.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib4_pvtalleles/0002.vcf

After filtering, kept 209 out of 209 Individuals
After filtering, kept 61210 out of a possible 61210 Sites
Run Time = 1.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib5_pvtalleles/0002.vcf

After filtering, kept 210 out of 210 Individuals
After filtering, kept 61210 out of a possible 61210 Sites
Run Time = 1.00 seconds

VCFtools - 0.1.16
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf Lib6_pvtalleles/0002.vcf

After filtering, kept 215 out of 215 Individuals
After filtering, kept 61210 out of a possible 61210 Sites
Run Time = 1.00 seconds
```


2. Next, we'll use a linear model to see if Library or Population explain individual missingness in the dataset

Use vcftools to calculate the missingness per individual
```

```

Use R to run the model. The file pop_lib contains population and library information for each individual and can be found [pop_lib]([pop_lib](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/blob/main/Files/pop_lib))
```
#R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
#Copyright (C) 2022 The R Foundation for Statistical Computing
#Platform: x86_64-apple-darwin17.0 (64-bit)

library("lme4")

miss <- read.table("out.imiss", header = T)
pop_lib <- read.table("pop_lib", header=T)
summary(pop_lib)
    INDV               POP            Library_nr
 Length:251         Length:251         1:45      
 Class :character   Class :character   2:44      
 Mode  :character   Mode  :character   3:43      
                                       4:42      
                                       5:41      
                                       6:36    



miss$pop <- as.factor(pop_lib$POP)
miss$lib <- as.factor(pop_lib$Library_nr)
summary(miss)
    INDV               N_DATA      N_GENOTYPES_FILTERED     N_MISS          F_MISS         lib         pop    
 Length:251         Min.   :61210   Min.   :0            Min.   :  562   Min.   :0.009182   1:45   BCH    :38  
 Class :character   1st Qu.:61210   1st Qu.:0            1st Qu.: 2410   1st Qu.:0.039381   2:44   SWD    :38  
 Mode  :character   Median :61210   Median :0            Median : 3187   Median :0.052067   3:43   WIS    :38  
                    Mean   :61210   Mean   :0            Mean   : 4781   Mean   :0.078108   4:42   HOD    :29  
                    3rd Qu.:61210   3rd Qu.:0            3rd Qu.: 4055   3rd Qu.:0.066247   5:41   BAR    :28  
                    Max.   :61210   Max.   :0            Max.   :32982   Max.   :0.538834   6:36   MOF    :25  
                                                                                                   (Other):55  




#linear model of missingness with library and pop as explanatory variables
fulllm <- lm(miss$F_MISS~miss$lib+miss$pop)
summary(fulllm)



Call:
lm(formula = miss$F_MISS ~ miss$lib + miss$pop)

Residuals:
     Min       1Q   Median       3Q      Max 
-0.14052 -0.03867 -0.01058  0.01299  0.46607 

Coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept)  0.093277   0.020470   4.557 8.32e-06 ***
miss$lib2   -0.006784   0.018057  -0.376 0.707470    
miss$lib3   -0.020511   0.018267  -1.123 0.262631    
miss$lib4    0.034245   0.019677   1.740 0.083097 .  
miss$lib5   -0.021216   0.019602  -1.082 0.280212    
miss$lib6   -0.029843   0.021000  -1.421 0.156600    
miss$popBCH -0.041405   0.021693  -1.909 0.057508 .  
miss$popBRO -0.018603   0.029651  -0.627 0.531005    
miss$popFOR  0.094334   0.025483   3.702 0.000266 ***
miss$popHOD  0.006415   0.023617   0.272 0.786141    
miss$popLYD  0.019090   0.025293   0.755 0.451151    
miss$popMOF -0.044540   0.023292  -1.912 0.057051 .  
miss$popSWD -0.015049   0.022843  -0.659 0.510670    
miss$popWIS -0.027311   0.020793  -1.313 0.190288    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.08334 on 237 degrees of freedom
Multiple R-squared:  0.2085,	Adjusted R-squared:  0.1651 
F-statistic: 4.802 on 13 and 237 DF,  p-value: 2.053e-07
```

### Conclusion: 

1. There are no loci that were sequenced only in a particular library. 
2. No library significantly affected the missingness rate (inverse of the genotyping rate) in the dataset


# 3. Analysis of smaller dataset with higher genotyping rate (95%)




# Analysis of ddRAD loci using a de novo assembled genome

A prior analysis of these data was conducted before the *A.agestis* was available. This resulted in a much smaller dataset, but the results matched those reported in this manuscript. Here we present the electronic lab book for the previous 1) data assessment, and 2) population structure analysis. 


## 1. Raw data to variants

# Raw data to variants

Call variants on ddRAD data. 

1. demultiplex and remove adapter using ipyRAD

2. map to reference genome with bwa mem

3. call variants with samtools mpileup and bcftools call

4. SNP filters


## Demultiplex and remove adapter

Data:

Paired-end ddRAD; 2 files per library

6 libraries in total, sequenced in 2014 (4 libraries on HD - Elements. 2 on fgcz server)

Barcodes and indiv names in BrownArgus_Barcodes.xlsx (on Mac)

Restriction enzymes

```
PstI TGCAG - 3'

EcoRI AATT - 3'
```

For previous datasets I've used process_radtags (from Stacks) and Trimmomatic to preprocess the reads before SNP calling. The main motivation being that I can specify a longer recognition sequence for each individual by concatenating the barcode and restriction site.

Example barcode file (here lib5). All files are here: alexjvr@fgcz-c-047:/srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/barcodes

```
AACCAACGAATT    AACCAACGTGCAG   WIS_3
AACCAGAAAATT    AACCAGAATGCAG   WIS_5
AACCATGCAATT    AACCATGCTGCAG   WIS_7
AACCGAACAATT    AACCGAACTGCAG   WIS_8
AACCGGTTAATT    AACCGGTTTGCAG   WIS_9
AACCTAGAAATT    AACCTAGATGCAG   WIS_10
ACCATAGGAATT    ACCATAGGTGCAG   WIS_13
ACTGACCAAATT    ACTGACCATGCAG   WIS_15
ACTTCTAAAATT    ACTTCTAATGCAG   WIS_16
AGGTTATGAATT    AGGTTATGTGCAG   WIS_17
AGTTATGCAATT    AGTTATGCTGCAG   WIS_19
ATATTACGAATT    ATATTACGTGCAG   WIS_21
CAACCTCTAATT    CAACCTCTTGCAG   WIS_22
CAGCGGACAATT    CAGCGGACTGCAG   WIS_23
```


Demultiplex using the output from the HiSeq. --disable_rad_check will allow the inclusion of the restriction site in the barcode

-r : rescue barcodes and radtags

-c : clean data, remove and read with an uncalled base

-q : discard reads with low quality scores (Q20)

-D : capture discarded reads to a file

(I didn't use -c -q filters, since data will be cleaned downstream):


```
/usr/local/ngseq/stow/stacks-1.28/bin/process_radtags -i gzfastq -1 rawData/library5/160325_D00248_0159_AC8TGPANXX_2_1.sanfastq.gz  -2 rawData/library5/160325_D00248_0159_AC8TGPANXX_2_2.sanfastq.gz -o ./demultiplexed/ -y fastq -b barcodes/barcodes_lib5 --inline_inline --disable_rad_check -r -D
```

When I tried this for the BA data, these strings can't be found. I can find the individual barcodes, but not with the restriction site afterwards. Why would this be?


Start ipyrad in Jupyter notebook
Information on how to start a Jupyter notebook that runs on a remote server can be found at these links:

http://kawahara.ca/how-to-run-an-ipythonjupyter-notebook-on-a-remote-machine/

http://ipyrad.readthedocs.io/HPC_Tunnel.html

Basically,

1. Set a password for the jupyter notebook
at alexjvr@fgcz-c-047:/srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads:

Password can be set up beforehand. The following code will prompt you to enter a password. If an instance is started without this command, a password will be automatically generated and can be found in the jupyter-notebook.config file, or on the server as the "token=xxx" script.

jupyter-notebook password
2. Start a jupyter notebook instance on the remote server. Do this in screen in case the connection times out.
at alexjvr@fgcz-c-047:/srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads

jupyter notebook --no-browser --port=8898
The port can be whatever you want it to be. And if a port is still in use, it can be closed. See troubleshooting in the first link.

3. Start notebook on local computer
Now open the notebook in a browser from the local computer. On terminal:

ssh -N -f -L 127.0.0.1:8898:127.0.0.1:8898 fgcz47
And then in an web browser address bar:

http://127.0.0.1:8898/
4. Section off a subset of nodes on the server for pyrad to use
On fgcz (remember to start these in screen in case the connection is lost).

ipcluster start --n=20
ipyRAD step 1: Demultiplex comparisons
Compare number of reads obtained from Stacks, and ipyrad 0 mismatches or 1 mismatch.

This is all run in the jupyter notebook: ipyRADopt.ipynb

This test showed that allowing 1 mismatch in the barcode allows more sequences to be assigned. Barcodes are 2+ bases different from each other, so a single mismatch won't mis-assign sequences to individuals.




## 2. Map to Brown Argus genome

Romain sent me the first draft assembly of the Aricia agestis genome (romain.villoutreix@gmail.com)

Data is from a paired-end HiSeq run. Initial assembly with Discovar (denovo). These data are currently on Blue Crystal:

bluecp3

Aricia_agestis/02b_DDN_assembly_ada09/a.final/a.lines.fasta

Strategy:
Check whether mapping to the draft genome is working

Optimise mapping and SNP calling using GATK best practices

2.1. Demultiplex all data (test 0-1 mismatches)

2.2. Remove all adapter dimer

2.3. Map using BWA-mem

      2.3.1 Optimise insert size
       
      2.3.2 Optimise mismatch and gap penalties
      
      2.3.3 Include RAD specific filter
      
2.4. Remove unmapped reads using Samtools

2.5. Sort sam files using Picard tools

2.6. Add read groups and remove PCR duplicates

2.7 HaplotypeCaller to call variants

2.8 Include only biallelic SNPs (filter using BCFtools)

Mapping and variant calling on the full dataset.

## Discovar de novo stats

### assembly statistics
### please see also frags.dist.png in parent directory

contig line N50: 8,315
scaffold line N50: 8,808
total bases in 1 kb+ scaffolds: 619,253,409
total bases in 10 kb+ scaffolds: 274,647,024
There are 167,594,314 reads of mean length 249.6 and mean base quality 35.9.
MPL1 = mean length of first read in pair up to first error = 222
(normal range is 175-225 for 250 base reads)
Estimated chimera rate in read pairs (including mismapping) = 0.15%.
genomic read coverage, using 1 kb+ scaffolds for genome size estimate: 67.6
1. Test map to genome
I'm copying the draft genome to the fgcz server, since all the raw data is there.

Checked with md5sum that everything copied over correctly

1.1 Index the genome
/usr/local/ngseq/packages/Aligner/BWA/0.7.15/bin/bwa index a.lines.fasta

## Real time: 1459.824 sec; CPU: 1460.068 sec  ## ~24min
1.2 test alignment
Map WIS_35_2014 demultiplexed using ipyrad (i.e. no .rem files). Using multithreading (-t 20) on the fgcz server.

/usr/local/ngseq/packages/Aligner/BWA/0.7.15/bin/bwa mem -t 20 a.lines.fasta /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6/WIS_35_2014_R1_.fastq.gz /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6/WIS_35_2014_R2_.fastq.gz > aln-pe.WIS35.sam

##Real time: 235.844 sec; CPU: 3329.964 sec ## ~2min
Processing of sequences

##index the draft genome
/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools faidx a.lines.fasta


##convert sam to bam
/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools import a.lines.fasta.fai aln-pe.WIS35.sam aln-pe.WIS35.bam

##sort bam file
/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools/1.5/bin/samtools sort aln-pe.WIS35.bam -o WIS35.bam.sorted
Get stats

/usr/local/ngseq/packages/Tools/samtools/1.5/bin/samtools flagstat WIS35.bam.sorted

3412445 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
266029 + 0 supplementary
0 + 0 duplicates
3396463 + 0 mapped (99.53% : N/A)
3146416 + 0 paired in sequencing
1573208 + 0 read1
1573208 + 0 read2
2615844 + 0 properly paired (83.14% : N/A)
3122622 + 0 with itself and mate mapped
7812 + 0 singletons (0.25% : N/A)
502576 + 0 with mate mapped to a different chr
219198 + 0 with mate mapped to a different chr (mapQ>=5)
Or with samstats for an html of the mapping statistics per sample. The command is simply to call samstat and then name all the files to be analysed.

/usr/local/ngseq/packages/QC/SAMStat/1.5.1/bin/samstat *sam    
2. Optimise mapping and variant calling
55 indivs in test run

2.1. Optimise demutliplexing
I used pyrad to demultiplex all the samples and tested whether 0, 1 or 2 mismatches were best. https://github.com/alexjvr1/Butterflies/blob/master/2.ipyRADopt.ipynb

Results can all be found here: /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexingOpt.

I trimmed Illumina adapters with pyRAD, but when I look at the log files, it looks like cutadapt failed because the files were not fastq. I'm not sure why I got this error, or why pyRAD didn't flag this while I was running the analysis.

I will rerun adapter trimming using Trimmomatic.

/srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6

for i in *_R1_.fastq.gz; do java -jar /usr/local/ngseq/packages/QC/Trimmomatic/0.36/trimmomatic-0.36.jar PE -threads 10 -trimlog test.log -basein $i -baseout $i ILLUMINACLIP:/usr/local/ngseq/packages/QC/Trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36; done

2.1 Optimise mapping
There are several mappers available. I'm using Teaser to compare them for this dataset. Paper.

Teaser
This can be downloaded in the working directory.

I'm using BWA-mem for mapping. I will run a test on a subset of individuals that should represent the total genetic diversity:

First the simplest possible method

for i in $_R1_.fastq.gz
/usr/local/ngseq/packages/Aligner/BWA/0.7.15/bin/bwa mem -t 20 a.lines.fasta /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6/WIS_35_2014_R1_.fastq.gz /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6/WIS_35_2014_R2_.fastq.gz > aln-pe.WIS35.sam

Jonathan Purtiz makes some suggestsions:

https://github.com/jpuritz/Winter.School2017/blob/master/Exercises/Day%201/Mapping%20Exercise.md

I (insert size): This is estimated by bwa from the mapping distributions. But it can be forced with e.g. -I [mean, SD]. Vilhelmiina tried this manually with some WPA data, but there was not much difference in the number of mapped reads.

Mismatch and opening gap penalties

RAD specific filter. (This is only useful for traditional RAD. The sequence is expected to have better quality on the 5' end where the restriction enzyme cut, so the sequences can be cut based on the decrease in quality).

for i in *fastq.gz; do 
/usr/local/ngseq/packages/Aligner/BWA/0.7.15/bin/bwa mem -t 20 a.lines.fasta /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6/test/$iR1* /srv/kenlab/alexjvr_p1795/Butterflies/BrownArgus_ddRAD_raw_reads/demultiplexed.ipyrad_lib1.6/test/*R2* > aln-pe.WIS35.sam


## 3. Mapping and variant calling: full dataset

Final mapping and variant calling was performed on bluecrystal. 

Mapping with BWA mem using the script in 2_MapwithBWAmem.ARRAY.sh

Variants were called using mpileup and bcftools call. Script: call_SNVs_bluecrystal.pl

The only filters were mpileup -p 0.05, and bcftools call -e FMT/DP=0


Final raw variant file can be found on bluecrystal here: /panfs/panasas01/bisc/aj18951/1a_Aricia_agestis_PopGenomics/03_variants/variants.raw.vcf



## 4. SNP filters

Assess final dataset on bluecrystal

```
module avail

module load apps/vcftools-0.1.12b

vcftools --vcf variants.raw.vcf

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf variants.raw.vcf

After filtering, kept 286 out of 286 Individuals
After filtering, kept 282092 out of a possible 282092 Sites
Run Time = 105.00 seconds

```

Remove the two populations that we won't be working with (TP and GTT)

```
#load bcftools module
module load apps/bcftools-1.8
#use bcftools to extract all sample names from vcf file
bcftools query -l variants.raw.vcf > raw.sample.names  

#extract TP and GTT names from this list
grep TP raw.sample.names >> removeTP.GTT.names
grep GTT raw.sample.names >> removeTP.GTT.names 
#Check the list
cat removeTP.GTT.names 

#remove these samples from the raw vcf file
vcftools --vcf variants.raw.vcf --remove removeTP.GTT.names --recode --recode-INFO-all --out AAgestis.276.raw

After filtering, kept 276 out of 286 Individuals
Outputting VCF file...
After filtering, kept 282092 out of a possible 282092 Sites
Run Time = 270.00 seconds

```

Look at the distribution of missingness across the dataset: 

```
vcftools --vcf variants.raw.vcf --missing-indv


```

## 2. Variant filtering

# Filtering raw SNP file 

Filtering variants from the raw vcf output to produce two datasets: 

1) Population structure analyses (little missing data)

2) Signals of selection (include as many loci as possible). 



#Filtering has been redone using reads mapped to the Sanger genome. See new filter sets [here](https://github.com/alexjvr1/AriciaAgestis_GWASMS/blob/main/2.VariantFiltering.md)



Final dataset for PopGen statistics (MAF 1%):

251 individuals

84117 variants

10097 loci



Final dataset for Outlier analyses (MAF 5%):

251 individuals

35463 variants

6785 loci






# See below for filtering on old dataset using deNovo Reference genome  ######################

### Filters to apply

1) Remove loci with QUAL < 100 (i.e. Phred confidence in variant site)

2) Minimum mean depth of 10 (i.e. remove loci with lower than mean 6x depth)

3) Max mean depth of mean + 2xSD of meanDP (here 255)

4) remove all multi-allelic SNPs

5) Remove all loci genotyped in <30% of individuals

6) Remove individuals with >60% missingness

Finally 

7) Remove loci with <1% MAF (to exclude any errors).

8) Filter for <5% MAF when doing outlier detection. 



##### Dataset1
```
#find meanDP with --depth flag in vcftools

vcftools --vcf AAgestis.276.raw.recode.vcf --minQ 100 --min-meanDP 10 --max-meanDP 255 --max-alleles 2 --max-missing 0.3 --recode --recode-INFO-all --out AAgestis_FINAL

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AAgestis.276.raw.recode.vcf
	--recode-INFO-all
	--max-alleles 2
	--max-meanDP 255
	--min-meanDP 10
	--minQ 100
	--max-missing 0.3
	--out AAgestis_FINAL
	--recode

After filtering, kept 276 out of 276 Individuals
Outputting VCF file...
After filtering, kept 63054 out of a possible 282092 Sites
Run Time = 151.00 seconds

```

Remove individuals with a high proportion of missingness. 

```
vcftools --vcf AAgestis_FINAL.recode.vcf --missing-indv

awk '!/IN/' out.imiss | cut -f5 > totalmissing

gnuplot << \EOF 
set terminal dumb size 120, 30
set autoscale 
unset label
set title "Histogram of % missing data per individual"
set ylabel "Number of Occurrences"
set xlabel "% of missing data"
#set yr [0:100000]
binwidth=0.01
bin(x,width)=width*floor(x/width) + binwidth/2.0
plot 'totalmissing' using (bin( $1,binwidth)):(1.0) smooth freq with boxes
pause -1
EOF

```

The vast majority of the samples have <60% missing data. 


![alt_txt][missing1]

[missing1]:https://user-images.githubusercontent.com/12142475/52953422-57e50200-337f-11e9-82ba-f5acce21e50f.png



Remove individuals with >60% missingness
```
awk '$5>0.6 {print $1}' out.imiss > indivstoremove

cat indivstoremove

INDV
BAR_1_2013_R1_.fastq.gz
BAR_23_2014_R1_.fastq.gz
BCH_38_2013_R1_.fastq.gz
BRO_4_2014_R1_.fastq.gz
HOD_13_2014_R1_.fastq.gz
HOD_18_2014_R1_.fastq.gz
HOD_6_2014_R1_.fastq.gz
HOD_8_2014_R1_.fastq.gz
LYD_20_2014_R1_.fastq.gz
LYD_28_2014_R1_.fastq.gz
MOF_42_2014_R1_.fastq.gz
SWD_17_2013_R1_.fastq.gz
```


Remove samples
```
vcftools --vcf AAgestis_FINAL.recode.vcf --remove indivstoremove --recode --recode-INFO-all --out AAgestis.264_FINAL

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AAgestis_FINAL.recode.vcf
	--exclude indivstoremove
	--recode-INFO-all
	--out AAgestis.264_FINAL
	--recode

Excluding individuals in 'exclude' list
After filtering, kept 264 out of 276 Individuals
Outputting VCF file...
After filtering, kept 63054 out of a possible 63054 Sites
Run Time = 47.00 seconds

```


And check if the filters remove any more loci: 
```
vcftools --vcf AAgestis.264_FINAL.recode.vcf --minQ 100 --min-meanDP 6 --max-missing 0.3 --recode --recode-INFO-all --out AAgestis.264.63054_FINAL

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AAgestis.264_FINAL.recode.vcf
	--recode-INFO-all
	--min-meanDP 6
	--minQ 100
	--max-missing 0.3
	--out AAgestis.264.63054_FINAL
	--recode

After filtering, kept 264 out of 264 Individuals
Outputting VCF file...
After filtering, kept 63054 out of a possible 63054 Sites
Run Time = 51.00 seconds

```

Rename the individuals in this file
```
module load apps/bcftools-1.8

bcftools query -l AAgestis.264_FINAL.recode.vcf > A264.oldnames

sed 's/_R1_fastq.gz//g' A264.oldnames > A264.newnames

bcftools reheader AAgestis.264_FINAL.recode.vcf -s A264.newnames -o AAgestis.264_FINAL.newnames.vcf

```




#### Assess the data

These stats are based on the dataset before filtering for MAF 1% (66054 loci)

*Mean depth per site per population* 

Individual mean depth vs number of genotyped loci. Coloured per population. 

Coverage per individual and per pop


vcftools finds these statistics. 

[OUTPUT DEPTH STATISTICS](http://vcftools.sourceforge.net/man_latest.html)

--depth

Generates a file containing the mean depth per individual. This file has the suffix ".idepth".

--site-depth

Generates a file containing the depth per site summed across all individuals. This output file has the suffix ".ldepth".

--site-mean-depth

Generates a file containing the mean depth per site averaged across all individuals. This output file has the suffix ".ldepth.mean".


1. Plot of mean depth per individual grouped by population: raw data first, and then filtered


```
/Users/alexjvr/2018.postdoc/BrownArgus_2018/201902_DataAnalysis/StatsReseqData

R version 3.5.0 (2018-04-23) -- "Joy in Playing"
Copyright (C) 2018 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin15.6.0 (64-bit)

DP.indivs <- read.table("out.idepth_vcf.raw", header=T)
head(DP.indivs)
DP.indivs$pop <- gsub("_.*gz", "", DP.indivs$INDV)
DP.indivs$dataset <- "raw"
head(DP.indivs)

DP.indivs.filtered <- read.table("out.idepth", header=T)
head(DP.indivs.filtered)
DP.indivs.filtered$pop <- gsub("_.*_.*", "", DP.indivs.filtered$INDV)
DP.indivs.filtered$dataset <- "filtered"
head(DP.indivs.filtered)

Tot.DP.indivs <- rbind(DP.indivs, DP.indivs.filtered)

pdf("1a_AriciaAgestis_MeanDP.pdf")
ggplot(Tot.DP.indivs, aes(x=pop, y=MEAN_DEPTH, colour=dataset)) + geom_boxplot()
dev.off()
```

![alt_txt][meanDP]

[meanDP]:https://user-images.githubusercontent.com/12142475/53504646-3b954380-3aaa-11e9-9255-b6e5503f1744.png




2. Mean depth per site in raw & filtered data


```
DP.sites <- read.table("out.ldepth_vcf.raw", header=T)
DP.sites.filtered <- read.table("out.ldepth", header=T)
DP.sites.filtered$dataset <- "filtered"
DP.sites$dataset <- "raw"
TotDP.sites <- rbind(DP.sites, DP.sites.filtered)

pdf("a1_AriciaAgestis_DepthPerSite.pdf")
ggplot(TotDP.sites, aes(SUM_DEPTH, colour=dataset)) + geom_histogram()
dev.off()


summary(DP.sites.filtered)
         CHROM            POS           SUM_DEPTH      SUMSQ_DEPTH      
 contig_365  :   82   Min.   :     4   Min.   : 2518   Min.   :   59081  
 m_scaff_3898:   80   1st Qu.:  2275   1st Qu.:13229   1st Qu.: 1643958  
 m_scaff_714 :   80   Median :  5135   Median :29046   Median : 5223654  
 contig_2942 :   79   Mean   :  7790   Mean   :30056   Mean   : 6068789  
 m_scaff_5830:   78   3rd Qu.: 10551   3rd Qu.:45418   3rd Qu.: 9723112  
 contig_16698:   77   Max.   :147561   Max.   :68598   Max.   :18003826  
 (Other)     :62578                                                      
   dataset         
 Length:63054      
 Class :character  
 Mode  :character  
                   
                   
summary(DP.sites)
          CHROM             POS           SUM_DEPTH      SUMSQ_DEPTH      
 m_scaff_6129:   176   Min.   :     4   Min.   :    1   Min.   :       1  
 m_scaff_5252:   168   1st Qu.:  2197   1st Qu.:    4   1st Qu.:       9  
 m_scaff_6232:   153   Median :  5058   Median :   39   Median :     226  
 contig_2942 :   152   Mean   :  7780   Mean   : 7724   Mean   : 1528150  
 m_scaff_910 :   150   3rd Qu.: 10497   3rd Qu.: 4038   3rd Qu.:  390659  
 m_scaff_714 :   144   Max.   :147699   Max.   :71626   Max.   :18744992  
 (Other)     :281149                                                      

```

![alt_txt][meanDP_sites]

[meanDP_sites]:https://user-images.githubusercontent.com/12142475/53504859-a9416f80-3aaa-11e9-92d6-73c1f1969d9d.png




*Number of variants per population (with individual variation)*


```
out.imiss <- read.table("out.imiss", header=T)
out.imiss$N_Variants <- (out.imiss$N_DATA-out.imiss$N_MISS)
out.imiss$pop <- gsub("_.*_.*", "", out.imiss$INDV)

pdf("1a_AriciaAgestis_nrVariantsPerPop.pdf")
ggplot(out.imiss, aes(x=pop, y=N_Variants)) + geom_boxplot()
dev.off()
```

![alt_txt][Nvariants]

[Nvariants]:https://user-images.githubusercontent.com/12142475/53504951-d8f07780-3aaa-11e9-9efc-ad74416d245a.png


There's some variance in the number of variants between populations. Most concerning is FOR which seems to have several indivdiduals with a low number of variants. I need to check whether this is a problem with the data before I continue. 

An obvious problem could be coverage per sample: 


*Nr variants vs Depth*

```
Tot.filtered <- merge(out.imiss, DP.indivs.filtered, by="INDV")
head(Tot.filtered)


pdf("1a_AriciaAgestis_variantsVsDepth.pdf")
ggplot(Tot.filtered, aes(x=N_Variants, y=MEAN_DEPTH, colour=pop.x)) + geom_point()
dev.off()
```

![alt_txt][variants_DP]

[variants_DP]:https://user-images.githubusercontent.com/12142475/53505386-9ed3a580-3aab-11e9-8c81-ba6671317704.png



And depth vs het
```
het.264 <- read.table("out.het", header=T)
het.264$het <- (het.264$N_SITES-het.264$O.HOM.)/het.264$N_SITES
het.264$pop <- gsub("_.*_.*", "", het.264$INDV)
head(het.264)
Tot.filtered$het <- het.264$het


```

![atl_txt][Tot.filtered_DPvshet]

[Tot.filtered_DPvshet]:https://user-images.githubusercontent.com/12142475/53509175-2f61b400-3ab3-11e9-9465-5577737ba756.png



There's a worrying correlation between MEAN_DEPTH and number of variants called in an individual. However, the mean depth is really high (>50x), and I've filtered for minQ (variant quality) PHRED 100. 
However, I can't filter individual genotypes without GQ (Genotype Quality) flag in the vcf file. 


When I filter the dataset to include only unlinked loci (assuming a max ddRAD locus of 600bp based on the insert size in the library prep): 

```
vcftools --vcf AAgestis.264.thin600.recode.vcf --het

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AAgestis.264.thin600.recode.vcf
	--het

After filtering, kept 264 out of 264 Individuals
Outputting Individual Heterozygosity
After filtering, kept 5625 out of a possible 5625 Sites
Run Time = 1.00 seconds

```

And then look at depth vs N_variants: 
```
out.het.unlinked <- read.table("AAgestis.264.thin600/out.het", header=T)
out.het.unlinked$het <- (out.het.unlinked$N_SITES-out.het.unlinked$O.HOM.)/out.het.unlinked$N_SITES
out.het.unlinked$pop <- gsub("_.*_.*", "", out.het.unlinked$INDV)

ggplot(out.het.unlinked, aes(N_SITES, MEAN_DEPTH)) + geom_point()

```

![alt_txt][unlinked_DPTHvsN_VAR]

[unlinked_DPTHvsN_VAR]:https://user-images.githubusercontent.com/12142475/53508333-7189f600-3ab1-11e9-94ce-197347ce24bc.png


and N-variants (unlinked) vs het:
```
ggplot(out.het.unlinked, aes(N_SITES, het)) + geom_point()

```

![alt_txt][unlinked_het]

[unlinked_het]:https://user-images.githubusercontent.com/12142475/53509025-db56cf80-3ab2-11e9-8dba-659e659ffde4.png


But the het doesn't increase with MEAN_DEPTH
```
ggplot(out.het.unlinked, aes(MEAN_DEPTH, het)) + geom_point()

```

![alt_txt][unlinked_meandpvshet]

[unlinked_meandpvshet]:https://user-images.githubusercontent.com/12142475/53508240-430c1b00-3ab1-11e9-94a8-555ce8c7b0c7.png




F (inbreeding) vs Number of variants. 

![alt_txt][meanHet]

[meanHet]:https://user-images.githubusercontent.com/12142475/53207386-7c650680-362b-11e9-94e7-56c321d5d8b5.png








#### Population genetics statistics

Before continuing with these analyses I will filter the dataset to include only loci genotyped in at least 50% of individuals in each population. 
I'm also filtering for MAF 1% as these loci are more likely to include errors introduced during PCR or sequencing. 

I am also filtering out individuals that are related more than would be expected under random mating. 

##### 1. MAF 1%

Filter for MAF 1%
```
vcftools --vcf AAgestis.264_FINAL.newnames.vcf --maf 0.01 --recode --recode-INFO-all --out AA.264.41508_FINAL

VCFtools - v0.1.12b
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AAgestis.264_FINAL.newnames.vcf
	--recode-INFO-all
	--maf 0.01
	--out AA.264.41508_FINAL
	--recode

After filtering, kept 264 out of 264 Individuals
Outputting VCF file...
After filtering, kept 41508 out of a possible 63054 Sites
Run Time = 35.00 seconds

```



##### 2. Keep only loci genotyped across all populations

To compare between populations I will keep only loci that are genotyped in at least 50% of individuals across all populations.
	
First split the vcf file by population
```
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep BCH.names --recode --recode-INFO-all --out BCH
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep LYD.names --recode --recode-INFO-all --out LYD
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep SWD.names --recode --recode-INFO-all --out SWD
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep BAR.names --recode --recode-INFO-all --out BAR
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep HOD.names --recode --recode-INFO-all --out HOD
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep MOF.names --recode --recode-INFO-all --out MOF
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep WIS.names --recode --recode-INFO-all --out WIS
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep BRO.names --recode --recode-INFO-all --out BRO
vcftools --vcf AA.264.41508_FINAL.recode.vcf --keep FOR.names --recode --recode-INFO-all --out FOR

```


Then filter each of these for 50% missingness
```
/newhome/aj18951/1a_Aricia_agestis_PopGenomics/03_variants/PerPopVCFfile
for i in $(ls *vcf); do vcftools --vcf $i --max-missing 0.5 --recode --recode-INFO-all --out $i.maxmiss

```

bgzip all of these (this was done on the mac) and find the intersection using bcftools. This keeps only loci that are genotyped in all of the files. Various options available. See [here](https://samtools.github.io/bcftools/bcftools.html#isec)
```
/Users/alexjvr/2018.postdoc/BrownArgus_2018/201902_DataAnalysis/DiversityStats_AA264.41508_FINAL

for i in $(ls *vcf); do bgzip $i; done
for i in $(ls test/*gz); do tabix $i; done

bcftools isec -n 9 BAR.recode.vcf.maxmiss.recode.vcf.gz BCH.recode.vcf.maxmiss.recode.vcf.gz BRO.recode.vcf.maxmiss.recode.vcf.gz FOR.recode.vcf.maxmiss.recode.vcf.gz HOD.recode.vcf.maxmiss.recode.vcf.gz LYD.recode.vcf.maxmiss.recode.vcf.gz MOF.recode.vcf.maxmiss.recode.vcf.gz SWD.recode.vcf.maxmiss.recode.vcf.gz WIS.recode.vcf.maxmiss.recode.vcf.gz -p test
```

And merge all these files using bcftools.
```

for i in $(ls test/*gz); do tabix $i; done

bcftools merge test/0000.vcf.gz test/0001.vcf.gz test/0002.vcf.gz test/0003.vcf.gz test/0004.vcf.gz test/0005.vcf.gz test/0006.vcf.gz test/0007.vcf.gz test/0008.vcf.gz -O v > test/AA264.merge0.5missing.vcf


vcftools --vcf test/AA264.merge0.5missing.vcf 

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf test/AA264.merge0.5missing.vcf

After filtering, kept 264 out of 264 Individuals
After filtering, kept 31522 out of a possible 31522 Sites
Run Time = 3.00 seconds

```


This still leaves a really big file with 31k loci! 



##### 3. Relatedness filter

Estimate relatedness using the [KING method](https://academic.oup.com/bioinformatics/article/26/22/2867/228512) implemented in vcftools.  

```
/Users/alexjvr/2018.postdoc/BrownArgus_2018/201902_DataAnalysis/DiversityStats_AA264.41508_FINAL/

#first thin to include only one variant per locus

vcftools --vcf AA264.merge0.5missing.vcf --thin 600 --recode --recode-INFO-all --out
AA264.merge0.5missing.thin600

vcftools --vcf AA264.merge0.5missing.thin600.recode.vcf --relatedness2
mv out.relatedness AA264.merged.thinned.relatedness2
```

Read into R and plot
```
AA264.thin.relatedness2 <- read.table("AA264.merged.thinned.relatedness2", header=T)
head(AA264.thin.relatedness2)
AA264.thin.relatedness2$pop2 <- gsub("_.*_.*", "", AA264.thin.relatedness2$INDV2)
AA264.thin.relatedness2$pop1 <- gsub("_.*_.*", "", AA264.thin.relatedness2$INDV1)

pdf("AA264.thin.relatedness2.pdf")
ggplot(AA264.thin.relatedness2[which(AA264.thin.relatedness2$pop1==AA264.thin.relatedness2$pop2),], aes(x=pop1, y=RELATEDNESS_PHI)) + geom_boxplot()
dev.off()
```

![alt_txt][AA264.relatedness2]

[AA264.relatedness2]:https://user-images.githubusercontent.com/12142475/53795582-30b73480-3f2a-11e9-92ae-5acb79cf5390.png




Figure1 from [Manichaikul et al. 2010](https://academic.oup.com/bioinformatics/article/26/22/2867/228512) shows that there us a lot of variance in the expected measure of relatedness when using only a few markers. As this method assumes unlinked markers, the thinned dataset includes only 3831 variants. 

![alt_txt][KING.fig1]

[KING.fig1]:https://user-images.githubusercontent.com/12142475/53795669-5a705b80-3f2a-11e9-8c6a-5ada97bc8836.png


Given the small number of markers I can only accurately identify first and second order relatives. Third order relatives can't be distinguished from unrelated individuals. 
Using a (second order relative) cut-off of 0.05 I identified a few related individuals, mostly within the WIS pop: 

```
AA264.thin.relatedness2[which(AA264.thin.relatedness2$pop1==AA264.thin.relatedness2$pop2 & AA264.thin.relatedness2$RELATEDNESS_PHI>0.05 & AA264.thin.relatedness2$RELATEDNESS_PHI<0.5),]
            INDV1       INDV2 N_AaAa N_AAaa N1_Aa N2_Aa RELATEDNESS_PHI pop1
20144 BRO_24_2013  BRO_6_2013    264    101   472   435       0.0683572  BRO
58832 WIS_13_2013 WIS_15_2013    241     75   460   432       0.1020180  WIS
58840 WIS_13_2013 WIS_22_2013    260     76   460   444       0.1194690  WIS
58849 WIS_13_2013 WIS_28_2013    267    104   460   479       0.0628328  WIS
59372 WIS_16_2013 WIS_24_2013    248     86   415   471       0.0857788  WIS
```

I'm using vcftools to check the missingness within each of these individuals to choose which to remove: 
```
vcftools --vcf AA264.merge0.5missing.vcf --missing-indv

for i in $(cat related.indivs.toremove); do grep $i out.imiss ; done
BRO_24_2013	31522	0	1615	0.0512341
BRO_6_2013	31522	0	2239	0.0710298
WIS_13_2013	31522	0	1682	0.0533596
WIS_15_2013	31522	0	1444	0.0458093
WIS_22_2013	31522	0	1867	0.0592285
WIS_28_2013	31522	0	1037	0.0328977
WIS_16_2013	31522	0	3662	0.116173
WIS_24_2013	31522	0	997	0.0316287
```

I'll remove: 

BRO_6_2013, WIS_13_2013, and WIS_16_2013


And apply all the previous filters: 
```

vcftools --vcf AA.264.41508_FINAL.recode.vcf --remove related.indivs.toremove --recode --recode-INFO-all --out AA261

vcftools --vcf AA261.recode.vcf --minQ 100 --min-meanDP 6 --max-missing 0.3 --maf 0.01 --recode --recode-INFO-all --out AA261.FINAL

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AA261.recode.vcf
	--recode-INFO-all
	--maf 0.01
	--min-meanDP 6
	--minQ 100
	--max-missing 0.3
	--out AA261.FINAL
	--recode

After filtering, kept 261 out of 261 Individuals
Outputting VCF file...
After filtering, kept 41462 out of a possible 41508 Sites
Run Time = 30.00 seconds

```

And filter to keep only loci genotyped in 50% of individuals across all species: 
First split the vcf file by population
```
vcftools --vcf AA261.FINAL.recode.vcf --keep BCH.names --recode --recode-INFO-all --out BCH
vcftools --vcf AA261.FINAL.recode.vcf --keep LYD.names --recode --recode-INFO-all --out LYD
vcftools --vcf AA261.FINAL.recode.vcf --keep SWD.names --recode --recode-INFO-all --out SWD
vcftools --vcf AA261.FINAL.recode.vcf --keep BAR.names --recode --recode-INFO-all --out BAR
vcftools --vcf AA261.FINAL.recode.vcf --keep HOD.names --recode --recode-INFO-all --out HOD
vcftools --vcf AA261.FINAL.recode.vcf --keep MOF.names --recode --recode-INFO-all --out MOF
vcftools --vcf AA261.FINAL.recode.vcf --keep WIS.names --recode --recode-INFO-all --out WIS
vcftools --vcf AA261.FINAL.recode.vcf --keep BRO.names --recode --recode-INFO-all --out BRO
vcftools --vcf AA261.FINAL.recode.vcf --keep FOR.names --recode --recode-INFO-all --out FOR
```


Then filter each of these for 50% missingness
```
/newhome/aj18951/1a_Aricia_agestis_PopGenomics/03_variants/PerPopVCFfile_AA261/
for i in $(ls *vcf); do vcftools --vcf $i --max-missing 0.5 --recode --recode-INFO-all --out $i.maxmiss; done

```

bgzip all of these (this was done on the mac) and find the intersection using bcftools. This keeps only loci that are genotyped in all of the files. Various options available. See [here](https://samtools.github.io/bcftools/bcftools.html#isec)
```
/Users/alexjvr/2018.postdoc/BrownArgus_2018/201902_DataAnalysis/DiversityStats_AA264.41508_FINAL

for i in $(ls *maxmiss.recode.vcf); do bgzip $i; done
for i in $(ls *gz); do tabix $i; done


bcftools isec -n 9 BAR.recode.vcf.maxmiss.recode.vcf.gz BCH.recode.vcf.maxmiss.recode.vcf.gz BRO.recode.vcf.maxmiss.recode.vcf.gz FOR.recode.vcf.maxmiss.recode.vcf.gz HOD.recode.vcf.maxmiss.recode.vcf.gz LYD.recode.vcf.maxmiss.recode.vcf.gz MOF.recode.vcf.maxmiss.recode.vcf.gz SWD.recode.vcf.maxmiss.recode.vcf.gz WIS.recode.vcf.maxmiss.recode.vcf.gz -p AA261.0.5miss.9popsMerged
```

And merge all these files using bcftools.
```
for i in $(ls AA261.0.5miss.9popsMerged/*.vcf); do bgzip $i; done
for i in $(ls AA261.0.5miss.9popsMerged/*gz); do tabix $i; done

bcftools merge AA261.0.5miss.9popsMerged/0000.vcf.gz AA261.0.5miss.9popsMerged/0001.vcf.gz AA261.0.5miss.9popsMerged/0002.vcf.gz AA261.0.5miss.9popsMerged/0003.vcf.gz AA261.0.5miss.9popsMerged/0004.vcf.gz AA261.0.5miss.9popsMerged/0005.vcf.gz AA261.0.5miss.9popsMerged/0006.vcf.gz AA261.0.5miss.9popsMerged/0007.vcf.gz AA261.0.5miss.9popsMerged/0008.vcf.gz -O v > AA261.0.5miss.9popsMerged/AA261.0.5miss.9popsMerged.vcf


vcftools --vcf AA261.0.5miss.9popsMerged/AA261.0.5miss.9popsMerged

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf AA261.0.5miss.9popsMerged/AA261.0.5miss.9popsMerged.vcf

After filtering, kept 261 out of 261 Individuals
After filtering, kept 31381 out of a possible 31381 Sites
Run Time = 3.00 seconds
```

Final dataset for PopGen statistics (MAF 1%): 

261 individuals

31381 variants

3819 loci



Final dataset for Outlier analyses (MAF 5%):

261 individuals

18338 variants

3559 loci 



## 3. Population structure

# Population structure

Initial tests are run on the SE1.s3.recode.vcf dataset. 

179 individuals

6352 loci 

(only ~3000 unique loci, but I can't get rid of the duplicates). 


## Planned analyses

1. Fst table and heatmap

2. Isolation by distance

3. PCAdapt to look at initial population structure

4. fastStructure

5. chromopainter/fineStructure


### 1. Fst heatmap

Convert plink to structure (not fastStructure) input using pgdspider

```
/Users/alexjvr/2018.postdoc/BrownArgus_2018/BA_DataAnalysis_20180628/SumStats

##in R

library(adegenet)
library(hierfstat)
library(reshape)

BA179 <- read.structure("BA179.str")

How many genotypes are there? 179

 How many markers are there? 6352

 Which column contains labels for genotypes ('0' if absent)? 1

 Which column contains the population factor ('0' if absent)? 0

 Which other optional columns should be read (press 'return' when done)? 1: 

 Which row contains the marker names ('0' if absent)? 1

 Are genotypes coded by a single row (y/n)? n

 Converting data from a STRUCTURE .stru file to a genind object... 

Warning message:
In df2genind(X = X, pop = pop, ploidy = 2, sep = sep, ncode = ncode) :
  duplicate labels detected for some loci; using generic labels


pop.BA179 <- read.table("BA179.popnames", header=T) #read in pop names
pop.factor <- as.factor(pop.BA179$pop)
BA179@pop <- pop.factor #assign to genind

hier.BA179 <- genind2hierfstat(BA179) # convert to hierfstat format

BA179.fst <- pairwise.fst(BA179, pop=NULL, res.type=c("dist", "matrix"))  #calculate fst

m <- as.matrix(BA179.fst)
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2)<- c("c1","c2", "distance")

library(gplots)

shadesOfGrey <- colorRampPalette(c("grey100", "grey0"))  ##define the colourpalette. 

Dend <- read.table("heatmap.popcolours", header=T)  ##list of colour names for each population based on R colour palatte. In alphabetical order (as in genind file): Samples automatically get reordered by hierfstat
Dend.Colours <- as.character(Dend$colours.pop)

##Plot 1: Basic
par(oma=c(1,1,2,1))
heatmap.2(as.matrix(BA179.fst), trace="none", RowSideColors=Dend.Colours, ColSideColors=Dend.Colours, col=shadesOfGrey, labRow=F, labCol=F, key.ylab=NA, key.xlab=NA, key.title="Fst Colour Key", keysize=0.9, main="Pairwise Fst of BA: 179 indivs, 9 pops, 6352loci")  ##RowSideColors is for the dendrogram on the row, ColSideColors for the upper dendrogram. Colour order should be the same as the input. The pop order is alphabetical in the output. 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

popnames.all <- as.character(Dend$pop)
legend("bottom", popnames.all, xpd = TRUE, horiz = TRUE, inset = c(0, 0), bty="o", pch=15, col=Dend.Colours, title="Side Dendrogram:Region")


##Plot 2: half a matrix
BA.popnames = c("1.BCH", "2.LYD", "3.SWD", "4.BHH", "5.HOD", "6.MOF", "7.WIS", "8.BRO", "9.FOR")
BA179.fst.modified <- m
BA179.fst.modified[BA179.fst.modified == 0] <- NA  ##To remove the 0 inflation in the heatmap    
#BA.popnames <- rownames(BA179.fst)

pdf("Fig2.BA179.Fst.July2018.pdf")
par(oma=c(1,1,2,1))
heatmap.2(as.matrix(BA179.fst.modified), na.rm=T, trace="none", RowSideColors=Dend.Colours, ColSideColors=Dend.Colours, col=shadesOfGrey, labRow=BA.popnames, labCol=F, key.ylab=NA, key.xlab=NA, key.title="Fst Colour Key", keysize=0.9, main="Pairwise Fst, 6352loci")  ##RowSideColors is for the dendrogram on the row, ColSideColors for the upper dendrogram. Colour order should be the same as the input. The pop order is alphabetical in the output. 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
dev.off()
```

![alt_txt][Fst.heatmap]

[Fst.heatmap]:https://user-images.githubusercontent.com/12142475/42625408-9eabfe8c-85bf-11e8-8581-c7e0d141c5f1.png




### 2. Isolation by distance


In the same directory as before: /Users/alexjvr/2018.postdoc/BrownArgus_2018/BA_DataAnalysis_20180628/SumStats

```
##Fst 9 pops -> Fst/(1-Fst)

library(reshape)
library(fields)


m <- as.matrix(BA179.fst)
m
m2 <- melt(m)[melt(upper.tri(m))$value,]
names(m2) <- c("c1", "c2", "distance")
m2
m2$IBD <- m2$distance/(1-m2$distance)


BA.pop.coords <- read.table("BA.pop.coords", header=T)
BApop_lon.lat <- cbind(BA.pop.coords$Long, BA.pop.coords$Lat)
distance.matrix.BApop <- rdist.earth(BApop_lon.lat, miles=F)  ##great circle dist based on the coordinates
m.dist <- as.matrix(distance.matrix.BApop)
summary(m.dist)

m2.dist <- melt(m.dist)[melt(upper.tri(m.dist))$value,]
names(m2.dist) <- c("c1", "c2", "distance")
summary(m2.dist)
m2.dist$log.km <- log(m2.dist$distance)


library(MASS)
#dens <- kde2d(m2$IBD,m2.dist$log.km, n=10)
#myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
plot(m2$IBD~m2.dist$log.km, pch=20,cex=.5, xlab="log Geographic distance (km)", ylab="Fst/(1-Fst)")
#image(dens, col=transp(myPal(10),.7), add=TRUE)
abline(fit <- lm(m2$IBD~m2.dist$log.km))
legend("bottomright", bty="n", legend=paste("R2 =", format(summary(fit)$adj.r.squared, digits=4)))  ##and paste R2
title("Isolation by distance plot - BA179")


```

![alt_txt][IBD]

[IBD]:https://user-images.githubusercontent.com/12142475/42630324-3dbc6ca6-85ce-11e8-9272-4a27f0d1177a.png


Is this significant? Test with AMOVA 

```
library(adegenet)
library(poppr)

popInfo <- read.table("BA179.sampleInfo", header=T)

BA179.strata <- popInfo[,1:4]  ##from text file. each column has one hierarchy level specified for all individuals. (indiv, pop, region)
BA179@other <- BA179.strata 

strata(BA179) <- other(BA179)

BA179.genclone <- as.genclone(BA179)

BA179.amova <- poppr.amova(BA179.genclone, ~host/pop)

BA179.amova 
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                            Df    Sum Sq   Mean Sq
Between host                 1  1082.165 1082.1654
Between pop Within host      7  3885.544  555.0777
Between samples Within pop 170 51996.027  305.8590
Within samples             179 40071.685  223.8642
Total                      357 97035.421  271.8079

$componentsofcovariance
                                            Sigma          %
Variations  Between host                 2.939043   1.072296
Variations  Between pop Within host      6.288167   2.294208
Variations  Between samples Within pop  40.997411  14.957712
Variations  Within samples             223.864161  81.675784
Total variations                       274.088783 100.000000

$statphi
                         Phi
Phi-samples-total 0.18324216
Phi-samples-pop   0.15478807
Phi-pop-host      0.02319076
Phi-host-total    0.01072296


BA179.amovatest <- randtest(BA179.amova, nrepet=999)

BA179.amovatest

class: krandtest lightkrandtest 
Monte-Carlo tests
Call: randtest.amova(xtest = BA179.amova, nrepet = 999)

Number of tests:   4 

Adjustment method for multiple comparisons:   none 
Permutation number:   999 
                        Test        Obs    Std.Obs   Alter Pvalue
1  Variations within samples 223.864161 -17.840353    less  0.001
2 Variations between samples  40.997411  14.224792 greater  0.001
3     Variations between pop   6.288167  27.803661 greater  0.001
4    Variations between host   2.939043   2.268624 greater  0.047

plot(BA179.amovatest)
```

![alt_txt][AMOVA]

[AMOVA]:https://user-images.githubusercontent.com/12142475/42632923-c475e01c-85d6-11e8-8b38-feb1be4e0178.png


Randomised test
```
##Randomised test
BA179.new <- BA179.genclone
set.seed(9001)
strata(BA179.new) <- strata(BA179.genclone)[sample(nInd(BA179.genclone)), -1]
head(strata(BA179.new))
head(strata(BA179.genclone))
BA179.new.amova <- poppr.amova(BA179.new, ~host/pop,within=F, quiet=T)

BA179.new.amova   ##now all the variation is within samples and within populations. So no population structure evident. 

$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                             Df    Sum Sq  Mean Sq
Between host                  1   164.086 164.0860
Between samples Within host   7  1136.912 162.4160
Within samples              170 27186.423 159.9201
Total                       178 28487.421 160.0417

$componentsofcovariance
                                               Sigma            %
Variations  Between host                  0.01826172   0.01140898
Variations  Between samples Within host   0.12594976   0.07868695
Variations  Within samples              159.92013587  99.90990406
Total variations                        160.06434734 100.00000000

$statphi
                           Phi
Phi-samples-total 0.0009009594
Phi-samples-host  0.0007869593
Phi-host-total    0.0001140898


BA179.new.amova.test<- randtest(BA179.new.amova, nrepet=999) 

BA179.new.amova.test

class: krandtest lightkrandtest 
Monte-Carlo tests
Call: randtest.amova(xtest = BA179.new.amova, nrepet = 999)

Number of tests:   3 

Adjustment method for multiple comparisons:   none 
Permutation number:   999 
                        Test          Obs    Std.Obs   Alter Pvalue
1  Variations within samples 159.92013587 -0.5969837    less  0.266
2 Variations between samples   0.12594976  0.5372204 greater  0.277
3    Variations between host   0.01826172  0.1657525 greater  0.460

pdf("BA179.AMOVA.randomised.pdf")
plot(BA179.new.amova.test)
dev.off()
```


![alt_txt][AMOVA.randtest]

[AMOVA.randtest]:https://user-images.githubusercontent.com/12142475/42633022-12dd6004-85d7-11e8-905d-366c8c8bdfb3.png



Second AMOVA test between for significant difference between pop and history

```
library(adegenet)
library(poppr)

popInfo <- read.table("BA179.sampleInfo", header=T)

BA179.strata <- popInfo[,1:4]  ##from text file. each column has one hierarchy level specified for all individuals. (indiv, pop, region)
BA179@other <- BA179.strata 

strata(BA179) <- other(BA179)

BA179.genclone <- as.genclone(BA179)

BA179.amova <- poppr.amova(BA179.genclone, ~history/pop)

BA179.amova 
$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                            Df    Sum Sq   Mean Sq
Between history              1  1275.014 1275.0140
Between pop Within history   7  3692.695  527.5279
Between samples Within pop 170 51996.027  305.8590
Within samples             179 40071.685  223.8642
Total                      357 97035.421  271.8079

$componentsofcovariance
                                            Sigma          %
Variations  Between history              4.270607   1.554568
Variations  Between pop Within history   5.581293   2.031678
Variations  Between samples Within pop  40.997411  14.923699
Variations  Within samples             223.864161  81.490055
Total variations                       274.713473 100.000000

$statphi
                         Phi
Phi-samples-total 0.18509945
Phi-samples-pop   0.15478807
Phi-pop-history   0.02063761
Phi-history-total 0.01554568


BA179.amovatest <- randtest(BA179.amova, nrepet=999)

BA179.amovatest

class: krandtest lightkrandtest 
Monte-Carlo tests
Call: randtest.amova(xtest = BA179.amova, nrepet = 999)

Number of tests:   4 

Adjustment method for multiple comparisons:   none 
Permutation number:   999 
                        Test        Obs    Std.Obs   Alter Pvalue
1  Variations within samples 223.864161 -17.358869    less  0.001
2 Variations between samples  40.997411  13.600294 greater  0.001
3     Variations between pop   5.581293  26.055948 greater  0.001
4 Variations between history   4.270607   3.218949 greater  0.034

plot(BA179.amovatest)
```

![alt_txt][AMOVA.history]

[AMOVA.history]:https://user-images.githubusercontent.com/12142475/42633489-92d8c716-85d8-11e8-9813-2f94ab970b14.png


Randomised test
```
##Randomised test
BA179.new <- BA179.genclone
set.seed(9001)
strata(BA179.new) <- strata(BA179.genclone)[sample(nInd(BA179.genclone)), -1]
head(strata(BA179.new))
head(strata(BA179.genclone))
BA179.new.amova <- poppr.amova(BA179.new, ~history/pop,within=F, quiet=T)

BA179.new.amova   ##now all the variation is within samples and within populations. So no population structure evident. 

$call
ade4::amova(samples = xtab, distances = xdist, structures = xstruct)

$results
                                Df     Sum Sq  Mean Sq
Between history                  1   158.2274 158.2274
Between samples Within history   7  1142.7707 163.2530
Within samples                 170 27186.4231 159.9201
Total                          178 28487.4212 160.0417

$componentsofcovariance
                                                  Sigma            %
Variations  Between history                 -0.05753439  -0.03595215
Variations  Between samples Within history   0.16783085   0.10487433
Variations  Within samples                 159.92013587  99.93107782
Total variations                           160.03043233 100.00000000

$statphi
                              Phi
Phi-samples-total    0.0006892218
Phi-samples-history  0.0010483664
Phi-history-total   -0.0003595215



BA179.new.amova.test<- randtest(BA179.new.amova, nrepet=999) 

BA179.new.amova.test

class: krandtest lightkrandtest 
Monte-Carlo tests
Call: randtest.amova(xtest = BA179.new.amova, nrepet = 999)

Number of tests:   3 

Adjustment method for multiple comparisons:   none 
Permutation number:   999 
                        Test          Obs    Std.Obs   Alter Pvalue
1  Variations within samples 159.92013587 -0.6348411    less  0.250
2 Variations between samples   0.16783085  0.7095399 greater  0.230
3 Variations between history  -0.05753439 -0.3921211 greater  0.636


pdf("BA179.AMOVA.randomised.pdf")
plot(BA179.new.amova.test)
dev.off()
```


![alt_txt][AMOVA.randtest]

[AMOVA.randtest]:https://user-images.githubusercontent.com/12142475/42633792-6fd70844-85d9-11e8-8f39-dff7fbcf0fa4.png




### 3. PCAdapt


input is .ped from plink
```
vcftools --vcf SE1.s3.recode.vcf --plink --out SE1.s3.plink

VCFtools - v0.1.14
(C) Adam Auton and Anthony Marcketta 2009

Parameters as interpreted:
	--vcf SE1.s3.recode.vcf
	--out SE1.s3.plink
	--plink

After filtering, kept 179 out of 179 Individuals
Writing PLINK PED and MAP files ... 

Unrecognized values used for CHROM: un - Replacing with 0.
Done.
After filtering, kept 6352 out of a possible 6352 Sites
Run Time = 1.00 seconds
```


import to R
```
R version 3.5.0

library(pcadapt)

BA179 <- read.pcadapt("SE1.s3.plink.ped", type="ped")
x.BA179 <- pcadapt(BA179, K=20)
plot(x.BA179, option="screeplot")
```

![alt_txt][screeplot]

[screeplot]:https://user-images.githubusercontent.com/12142475/42581924-c5acdf32-8525-11e8-93c3-8fe137853b9d.png

The PC components that explain population structure should be above the straight line. Here it looks like K=2. Possibly K=3 or 4. 

This is the plot

```
pop.BA179 <- read.table("BA179.popnames", header=T)
head(pop.BA179)

       indiv pop history     host
1 BCH_12_SE1 BCH     old Rockrose
2 BCH_14_SE1 BCH     old Rockrose
3 BCH_15_SE1 BCH     old Rockrose
4 BCH_19_SE1 BCH     old Rockrose
5  BCH_1_SE1 BCH     old Rockrose
6 BCH_20_SE1 BCH     old Rockrose

poplist.history <- as.character(pop.BA179[,3])

pdf("BA179.pca.pdf")
plot(x.BA179, option="scores", pop=poplist.host)
dev.off()
```

![alt_txt][pca]

[pca]:https://user-images.githubusercontent.com/12142475/42584953-485fe360-852c-11e8-99ea-45e37ab9563e.png



### 4. fastStructure

** Mac Housekeeping

After upgrading to MacOS HighSierra my dynamic libraries wouldn't load (dylib). The solution that finally worked was just to copy the new 
dylib to the name called by the package. I'm sure there's a better way of doing this, but I couldn't get any of the other solutions to work. It seems this is a problem with how homebrew installs some packages, similar to gcc being installed in an unexpected folder or with unexpected name. (brew doctor is useful for checking for problems with any of the homebrew packages). 
```
cp /usr/local/opt/gsl/lib/libgsl.23.dylib /usr/local/opt/gsl/lib/libgsl.0.dylib
```


Structure with no population or geographic prior

This can be run locally

```
pwd

python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 2 --format=str --input=BA179.fast --output=BA179_K2.1
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 2 --format=str --input=BA179.fast --output=BA179_K2.2
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 2 --format=str --input=BA179.fast --output=BA179_K2.3
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 2 --format=str --input=BA179.fast --output=BA179_K2.4
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 2 --format=str --input=BA179.fast --output=BA179_K2.5

python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 3 --format=str --input=BA179.fast --output=BA179_K3.1
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 3 --format=str --input=BA179.fast --output=BA179_K3.2
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 3 --format=str --input=BA179.fast --output=BA179_K3.3
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 3 --format=str --input=BA179.fast --output=BA179_K3.4
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 3 --format=str --input=BA179.fast --output=BA179_K3.5

python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 4 --format=str --input=BA179.fast --output=BA179_K4.1
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 4 --format=str --input=BA179.fast --output=BA179_K4.2
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 4 --format=str --input=BA179.fast --output=BA179_K4.3
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 4 --format=str --input=BA179.fast --output=BA179_K4.4
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 4 --format=str --input=BA179.fast --output=BA179_K4.5

python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 5 --format=str --input=BA179.fast --output=BA179_K5.1
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 5 --format=str --input=BA179.fast --output=BA179_K5.2
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 5 --format=str --input=BA179.fast --output=BA179_K5.3
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 5 --format=str --input=BA179.fast --output=BA179_K5.4
python /Users/alexjvr/Applications/fastStructure/fastStructure/structure.py -K 5 --format=str --input=BA179.fast --output=BA179_K5.5

... etc
```

I ran each K 5 times, K1:5. 

Choose the best K

```
python /Users/alexjvr/Applications/fastStructure/fastStructure/chooseK.py --input=BA179_K*
Model complexity that maximizes marginal likelihood = 2
Model components used to explain structure in data = 2
```

Check that the meanQ files for K=2 are similar to each other

```

```


And plot K=2. Samples are initially in alphabetical order

```
data1 <- read.table("BA179_K2.1.2.meanQ", header=F)
data1$indiv <- sample.info$indiv
data1$pop <- sample.info$pop
data1.sorted <- data1[order(pop),]
attach(data1)
data1.sorted <- data1[order(pop),]
data1.sorted


pdf("BA179.fastStructure.K2.pdf")

popsizes = tapply(data1.sorted$indiv, data1.sorted$pop, length)  ##find pop length
numpops = length(popsizes)
poploc=0
k=2
klast=k+5
popbegins=rep(NA,k)
popends=rep(NA,k)
poplabels = c("BCH", "LYD", "SWD", "BHH", "HOD", "MOF", "WIS", "BRO", "FOR")

bp <- barplot(t(data1.sorted[,1:2]), space=0, axes=F)
for(x in 1:numpops){
popbegins[x] = poploc
popends[x] = poploc + popsizes[x]
poploc = poploc + popsizes[x]
}
popmidpoints = (popbegins+popends)/2
#Puts a dark line between each orignial population
abline(v=c(0,popends),lwd=3,xpd=T)
#Label the populations underneath, at the midpoint of each population.
for(x in 1:numpops){
mtext(poplabels[x],side=1,at=popmidpoints[x],padj=2,cex=1.5,font=4)
}
dev.off()
```

![alt_txt][BA179.structure]

[BA179.structure]:https://user-images.githubusercontent.com/12142475/42640345-02de9728-85ea-11e8-8d6d-9163cd6b3c5c.png



Location of python and fastStructure on FGCZ47
```
/usr/bin/python2.7 
/usr/local/ngseq/src/fastStructure-1.0_20160929/structure.py
```

