Data analysis and scripts for the manuscript entitled: 

### Rapid evolution of novel biotic interactions in the UK Brown Argus butterfly uses genomic variation from across its geographical range


Running title:

Rapid adaptation in the UK Brown Argus butterfly uses widespread genome variation  
 

Maaike de Jong1,5*, Alexandra Jansen van Rensburg1,4*, Samuel Whiteford2, Carl J. Yung2, Mark Beaumont1, Chris Jiggins3, Jon Bridle1,4+

*these authors contributed equally to this publication 
+ To whom correspondence should be addressed: j.bridle@ucl.ac.uk



## 0.Map with sample sites

Use [0_SampleMap.r](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/blob/main/Scripts/0_SampleMap.r) with [BA_SiteInfo](https://github.com/alexjvr1/BrownArgus_PopGenMS_MolEcol/blob/main/Files/BA_SiteInfo) as input file to generate the sample map in Figure 1. 

## 1.RawData_to_Variants.md

Processing of raw Illumina data to a variant file. 

## 2.SNPfiltering.md

Filtering of SNPs for 1) population genetic analyses, and 2) outlier/Environmental association analyses. 

## 3.GeneticDiversity&Structure.md

1) Estimates of genetic diversity

2) Genetic distance between populations (Fst)

3) Estimate of Isolation by distance

4) PCA

5) fineRADstructure

## 4.ColonisationHistory

1. How did A.agestis colonise the UK? 

 - fastSimCoal analyses

2. Where were the new sites colonised from? 

 - fastSimCoal analyses 

## 5.SpecialisationAtRangeMargin.md

Have population at the range margin become specialised? What are the main drivers of genetic divergence here? 

   - Redundancy and partial Redundancy analyses


## 6.AdaptationtoHostPlant.md

Has A.agestis specialised on new host plant at range margins? 

   - Fst outlier analysis (Bayescan)
      
   - Identify outliers 

   - Find annotations for outliers


## 7.HaplotypeNetwork.md


      
      

