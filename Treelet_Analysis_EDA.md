---
title: "Treelet Transform: Identifing clusters of ICD-9 Diagnoses in a Boston Trauma Center"
subtitle: "Exploratory Data Analysis & Visualization"
author: "Dominic DiSanto"
date: "Updated 9/20/2020"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_depth: '5'
    code_folding: show
---

## Preparation


## Libraries

```r
library(magrittr) # Ceci n'est pas une %>% 
library(dplyr) # General data management, cleaning (admittedly I switch between Base R and tidyverse as I code, somewhat stream-of-consciousness ly)
library(ggplot2) # Visualization
library(comorbidity) # Used to easily generate Elixhauser comorbdity grouping/categorization [8/23/2020 Note: may be excluded if Elixhauser or Charlson not used]
library(tidyr) # pivot functions for transposing data to/from long and wide
library(icd) # used in validity check of diagnoses codes
# library(lubridate) # used in working with datetime variables, primarily calculating length of stay 
```

## File Paths

```r
lib <- "C:/Users/Dominic DiSanto/Documents/MIMIC-III/mimic-iii-clinical-database-1.4/" # General location of all MIMIC files
icd_fp <- paste0(lib, "Diagnoses_ICD.csv")
```

### Diagnoses Codes

#### Import & Validity Check


```r
icd_raw <- read.csv(icd_fp, stringsAsFactors = F) %>% select(-ROW_ID)

icd_raw %>% nrow()
```

```
## [1] 651047
```

```r
icd_raw %>% distinct(SUBJECT_ID) %>% nrow()
```

```
## [1] 46520
```

```r
icd_raw %>% distinct(ICD9_CODE) %>% nrow()
```

```
## [1] 6985
```


In our very introductory exploration of the data, we see that the MIMIC-III diagnoses code data contains 6,985 unique ICD-9 codes among 46,520 patients. The diagnoses table in total includes 651,047 rows, with each row being a unique diagnoses for a given patient and visit. Before working with the diagnosis data, we want to do some preliminary cleaning. I will first remove any duplicated diagnoses codes within a patient *and* visit. I will also remove the V and E codes which correspond to Health Services/Factors and Causes of Injury/Illness respecively, separate from diagnoses. 


```r
icd_precln <- icd_raw %>% distinct(SUBJECT_ID, HADM_ID, ICD9_CODE, .keep_all = T) %>% 
  filter(substring(ICD9_CODE, 1, 1)!="E" & substring(ICD9_CODE, 1, 1)!="V") #removing V and E codes

icd_precln %>% head()
```

```
##   SUBJECT_ID HADM_ID SEQ_NUM ICD9_CODE
## 1        109  172335       1     40301
## 2        109  172335       2       486
## 3        109  172335       3     58281
## 4        109  172335       4      5855
## 5        109  172335       5      4254
## 6        109  172335       6      2762
```

Using the `icd` package's built-in `is_defined` function, which tests whether a given input value 1) follows valid formatting for an ICD-9 code (5 or less characters, numeric or alphanumeric for V, E codes (which we've excluded)) and 2) is defined using a call to CMS, which keeps a list of what the package refers to as "canonical" ICD-9 codes:



```r
icd_precln %>% mutate(valid = is_defined(ICD9_CODE)) %>% filter(valid==F & ICD9_CODE!="")
```

```
## [1] SUBJECT_ID HADM_ID    SEQ_NUM    ICD9_CODE  valid     
## <0 rows> (or 0-length row.names)
```

Thankfully all our codes are valid/defined! We also know that within this data, the same patient (i.e. same `SUBJECT_ID`) may have multiple admissions (i.e. multiple `HADM_ID`'s). I will look now at the patients for which this happens. My current thought (and "solution") is to simply retain the most recent visit, believing that this most accurately reflects a patient's disease progression and current diagnoses. 



```r
icd_precln %>% distinct(SUBJECT_ID, HADM_ID, .keep_all = T) %>% arrange(SUBJECT_ID) %>% group_by(SUBJECT_ID) %>% count() %>% filter(n>1)
```

```
## # A tibble: 7,479 x 2
## # Groups:   SUBJECT_ID [7,479]
##    SUBJECT_ID     n
##         <int> <int>
##  1         17     2
##  2         21     2
##  3         23     2
##  4         34     2
##  5         36     3
##  6         61     2
##  7         67     2
##  8         68     2
##  9         84     2
## 10         85     2
## # ... with 7,469 more rows
```

We see that a sizeable proportion of this patient population (7,479 patients of the ~46,000 we first identified) have multiple `HADM_ID` values. I will import the `ADMISSIONS` table to pull in the `ADMITTIME` variable.

###### Identifying Patients with Multiple Stays


```r
admit <- read.csv(paste0(lib, "ADMISSIONS.csv"))

icd_raw_nodups <- admit %>% select(HADM_ID, ADMITTIME) %>%  merge(., icd_precln, by="HADM_ID", all.y=T) %>%  group_by(SUBJECT_ID) %>% 
  filter(ADMITTIME==max(ADMITTIME)) %>% ungroup()

icd_raw_nodups %>% distinct(SUBJECT_ID, HADM_ID, .keep_all = T) %>% arrange(SUBJECT_ID) %>% group_by(SUBJECT_ID) %>% count() %>% filter(n>1)
```

```
## # A tibble: 0 x 2
## # Groups:   SUBJECT_ID [0]
## # ... with 2 variables: SUBJECT_ID <int>, n <int>
```

```r
data.frame(
Data_Frame=c("icd_raw_nodups", "icd_precln"),

Unique_Patients = c(icd_raw_nodups %>% distinct(SUBJECT_ID) %>% nrow(),
                    icd_precln %>% distinct(SUBJECT_ID) %>% nrow()),
Unique_Visits = c(icd_raw_nodups %>% distinct(HADM_ID) %>% nrow(),
                  icd_precln %>% distinct(HADM_ID) %>% nrow()),
Rows=c(icd_raw_nodups  %>% nrow(),
       icd_precln %>% nrow())
)
```

```
##       Data_Frame Unique_Patients Unique_Visits   Rows
## 1 icd_raw_nodups           44319         44319 417679
## 2     icd_precln           44319         56716 553742
```

We can now see we've eliminated duplicate stays by filtering such that the visit's admittance time `ADMITTIME` must be the maximum within a given patient (`SUBJECT_ID`). I also check/demonstrate in the final two lines of the above chunk that the same number of patients are retained in this "no dups" (for no duplicates) data frame. But a much smaller number of rows and unique visits are retained.


#### EDA of Diagnoses Codes 

##### Diagnosis Frequency

Examining simply the most common diagnoses codes, arbitrarily picking the top 15 for legibility of plots. 


```r
icd_counts <- icd_raw_nodups %>% group_by(ICD9_CODE) %>% count() %>% arrange(desc(n))
```


##### Graphing Top 15 Most Common Diagnoses


```r
icd_descr <- read.csv(paste0(lib, "D_ICD_DIAGNOSES.csv"))
top_codes <- merge(icd_counts, icd_descr, by="ICD9_CODE", all.x=T) %>% arrange(desc(n)) %>% ungroup() %>%  filter(row_number()<=15)


ggplot(top_codes, aes(x=reorder(SHORT_TITLE, -n), y=n)) +
  geom_bar(stat="identity", fill="navyblue", alpha=0.65) +
  ggtitle("Frequency of the 15 Most Common Diagnoses", 
          subtitle = paste("n=", nrow(icd_raw_nodups), 
                           " unique diagnoses, among", 
                           length(unique(icd_raw_nodups$SUBJECT_ID)), "patients")) + 
  ylab("Frequency") + xlab("ICD-9 Code") + theme_minimal() +
  theme(axis.text.x=element_text(angle=80, vjust=0.9, hjust=0.8))
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
    # this is a really unfortunately x-axis, couldn't find a better angle or adjustment for the x-axis unfortunately 
```



##### Correlation Matrix Among Top Diagnoses

Looking at the correlation matrix of these most common diagnoses codes (as the correlation of diagnoses is what will determine the hierarchy of clustering in the treelet method):


```r
wide <- icd_raw_nodups %>% filter(ICD9_CODE %in% top_codes[1][[1]]) %>% mutate(values=1) %>% pivot_wider(id_cols="SUBJECT_ID", names_from="ICD9_CODE", values_from="values")
wide[is.na(wide)] <- 0

wide %>% select(-SUBJECT_ID) %>% cor() %>%  corrplot::corrplot(type="upper", diag=F, order="hclust")
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-10-1.png)<!-- -->


I had initially hoped to explore the most highly correlated diagnoses codes among all diagnoses. However the size of this data frame became fairly intractable for my computer, and I was not able to estimate a correlation matrix for all diagnoses codes. I've then been exploring additional ways to aggregate/subset diagnoses codes prior to dimension reduction based on either  
  1. Prevalence in the data set (with a low threshold to inclusion)  
  2. Previously published grouping algorithms (i.e. Elixhauser)  
  3. ICD-9's somewhat implicit structure of coding  

#### Prevalence of >1%

Here I simply slice the diagnoses codes to only include diagnoses present in at least 1% of our patients:


```r
pct_cutoff <- length(unique(icd_raw_nodups$SUBJECT_ID))*0.01
pct_cutoff_codes <- merge(icd_counts, icd_descr, by="ICD9_CODE", all.x=T) %>% arrange(desc(n)) %>% ungroup() %>%  filter(n>=pct_cutoff)
```


```r
prev_wide <- icd_raw_nodups %>% filter(ICD9_CODE %in% pct_cutoff_codes[1][[1]]) %>% mutate(values=1) %>% pivot_wider(id_cols="SUBJECT_ID", names_from="ICD9_CODE", values_from="values")
prev_wide[is.na(prev_wide)] <- 0

corr_mat <- cor(prev_wide)
corr_mat[abs(corr_mat)>0.3 & abs(corr_mat)!=1] # printing only a subset to see if there are any highly correlated variables
```

```
##  [1] 0.3230168 0.3230168 0.4911920 0.6132118 0.3418081 0.6065727 0.5288702
##  [8] 0.4046031 0.4353799 0.4222282 0.5670920 0.3374261 0.3418081 0.3144627
## [15] 0.4474047 0.7011283 0.3457130 0.7011283 0.3361047 0.3144627 0.5725083
## [22] 0.5725083 0.6065727 0.4252590 0.3748477 0.4932475 0.5288702 0.4252590
## [29] 0.4218611 0.4888356 0.3617305 0.6512039 0.3617305 0.4092391 0.7690218
## [36] 0.7690218 0.3886728 0.3361047 0.6512039 0.3886728 0.4911920 0.7447449
## [43] 0.6132118 0.7447449 0.3402869 0.4046031 0.3529470 0.3046853 0.3595060
## [50] 0.3457130 0.4353799 0.4474047 0.4222282 0.3748477 0.4218611 0.6050777
## [57] 0.5670920 0.4932475 0.4888356 0.3529470 0.6050777 0.3196800 0.3077870
## [64] 0.3196800 0.3005606 0.3077870 0.3005606 0.3595060 0.4092391 0.3402869
## [71] 0.3046853 0.7900441 0.3374261 0.7900441
```
  

#### Generating Elixhauser Group

Here I simply use the `comorbidity` function from the eponymous package to generate Elixhauser categories. The function generates 30 Elixhauser membership variables and additional variables related to a scoring or indexing system based upon these Elixhauser category memberships. Values of `1` represent membership in an Elixhauser category (i.e. a present ICD-9 diagnosis code for that categorie's relevant codes) and 0 representing no membership (i.e. no relevant codes recorded for a given patient). 



```r
icd_elix <- comorbidity(icd_raw_nodups, id = "SUBJECT_ID", code = "ICD9_CODE", score = "elixhauser", icd = "icd9", assign0 = F) %>% 
  select(-c(index, score, starts_with("wscore"), starts_with("windex"))) 
# removing aggregate variables of Charlson scoring that this function generates based on several different methods/publications 
icd_elix %>% head()
```

```
##   SUBJECT_ID chf carit valv pcd pvd hypunc hypc para ond cpd diabunc diabc
## 1          3   1     0    0   0   0      0    0    0   0   0       0     0
## 2          4   0     0    0   0   0      0    0    0   0   0       0     0
## 3          6   0     0    0   0   0      0    1    0   0   0       0     0
## 4          8   0     0    0   0   0      0    0    0   0   0       0     0
## 5          9   1     0    0   0   0      1    0    0   0   0       0     0
## 6         10   0     0    0   0   0      0    0    0   0   0       0     0
##   hypothy rf ld pud aids lymph metacanc solidtum rheumd coag obes wloss fed
## 1       0  0  0   0    0     0        0        0      0    0    0     1   0
## 2       0  0  1   0    1     0        0        0      0    0    0     1   1
## 3       0  1  0   0    0     0        0        0      0    0    0     0   1
## 4       0  0  0   0    0     0        0        0      0    0    0     0   0
## 5       0  0  0   0    0     0        0        0      0    0    0     0   1
## 6       0  0  0   0    0     0        0        0      0    0    0     0   0
##   blane dane alcohol drug psycho depre
## 1     0    0       0    0      0     0
## 2     0    0       0    0      0     0
## 3     0    0       0    0      0     0
## 4     0    0       0    0      0     0
## 5     0    0       0    0      0     0
## 6     0    0       0    0      0     0
```

```r
cat("Rows in Elixhauser output: ", nrow(icd_elix), " // Unique subject ID's in input data:", length(unique(icd_raw_nodups$SUBJECT_ID)),
    "\nNote: Values should be equal")
```

```
## Rows in Elixhauser output:  44319  // Unique subject ID's in input data: 44319 
## Note: Values should be equal
```


###### Graphing Group Frequency


```r
summary <- data.frame()

for (i in 2:ncol(icd_elix)){
  summary[i-1, "Freq"] <- sum(icd_elix[,i])  
  summary[i-1, "Group"] <- colnames(icd_elix)[i]
}


ggplot(summary, aes(x=reorder(Group, -Freq), y=Freq)) +
  geom_bar(stat="identity", fill="navyblue", alpha=0.65) +
  ggtitle("Frequency of the 30 Elixhauser Comorbidity Groups") + 
  ylab("Frequency") + xlab("Elixhauser Comorbidity Group") + theme_minimal() + theme(axis.text.x=element_text(angle=80, vjust=0.955))
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
    # this is a really unfortunately x-axis, couldn't find a better angle or adjustment for the x-axis unfortunately 
```


###### Elixhauser Correlation Matrix


```r
icd_elix %>% select(-SUBJECT_ID) %>% cor() %>% corrplot::corrplot(type="upper", diag=F, order="hclust")
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


### Using ICD-9's Implicit Structure

Coding my own "grouping" of codes based on the ICD-9's existing coding of data. 


```r
icd_dom_recode <-  icd_raw_nodups %>% mutate(Group=
                            case_when(
                              ICD9_CODE!="" & ICD9_CODE<"14000" ~ "Infections",
                              ICD9_CODE>="14000" & ICD9_CODE<"24000" ~ "Neoplasms",
                              ICD9_CODE>="24000" & ICD9_CODE<"28000" ~ "Endocrine",
                              ICD9_CODE>="28000" & ICD9_CODE<"29000" ~ "Blood Disorders",
                              ICD9_CODE>="29000" & ICD9_CODE<"32000" ~ "Mental Disorders",
                              ICD9_CODE>="32000" & ICD9_CODE<"39000" ~ "Nervous System",
                              ICD9_CODE>="39000" & ICD9_CODE<"46000" ~ "Circulatory",
                              ICD9_CODE>="46000" & ICD9_CODE<"52000" ~ "Respiratory",
                              ICD9_CODE>="52000" & ICD9_CODE<"58000" ~ "Digestive",
                              ICD9_CODE>="58000" & ICD9_CODE<"63000" ~ "Geniotourinary",
                              ICD9_CODE>="63000" & ICD9_CODE<"68000" ~ "Pregnancy and Childbirth",
                              ICD9_CODE>="68000" & ICD9_CODE<"71000" ~ "Skin and Tissue",
                              ICD9_CODE>="71000" & ICD9_CODE<"74000" ~ "Musculoskeletal",
                              ICD9_CODE>="74000" & ICD9_CODE<"76000" ~ "Congenital",
                              ICD9_CODE>="76000" & ICD9_CODE<"78000" ~ "Perinatal Conditions",
                              ICD9_CODE>="78000" & ICD9_CODE<"80000" ~ "Ill-Defined",
                              ICD9_CODE>="80000" & ICD9_CODE<"99999" ~ "Injury and Poisoning",
                            ))

icd_dom_recode %>% filter(!is.na(Group)) %>% count(Group) %>% arrange(desc(n))
```

```
## # A tibble: 17 x 2
##    Group                         n
##    <chr>                     <int>
##  1 Circulatory              104783
##  2 Endocrine                 51539
##  3 Injury and Poisoning      34459
##  4 Respiratory               33041
##  5 Digestive                 28107
##  6 Geniotourinary            23479
##  7 Ill-Defined               22590
##  8 Perinatal Conditions      19482
##  9 Mental Disorders          18969
## 10 Blood Disorders           17778
## 11 Nervous System            17283
## 12 Infections                14532
## 13 Neoplasms                 11163
## 14 Musculoskeletal            9784
## 15 Skin and Tissue            6434
## 16 Congenital                 3604
## 17 Pregnancy and Childbirth    612
```

```r
# icd_dom_recode_wide <-

icd_dom_recode_wide <- icd_dom_recode %>% filter(!is.na(Group)) %>% group_by(SUBJECT_ID, Group) %>% distinct(SUBJECT_ID, Group, .keep_all = T) %>% ungroup() %>% mutate(values=1) %>% select(SUBJECT_ID, Group, values) %>% pivot_wider(id_cols="SUBJECT_ID", names_from="Group", values_from="values")

icd_dom_recode_wide[is.na(icd_dom_recode_wide)] <- 0
```


```r
summary_dom <- data.frame()

for (i in 2:ncol(icd_dom_recode_wide)){
  summary_dom[i-1, "Freq"] <- sum(icd_dom_recode_wide[,i])  
  summary_dom[i-1, "Group"] <- colnames(icd_dom_recode_wide)[i]
}


ggplot(summary_dom, aes(x=reorder(Group, -Freq), y=Freq)) +
  geom_bar(stat="identity", fill="navyblue", alpha=0.65) +
  ggtitle("Frequency of the \"Implicit\" ICD-9 Groups ") + 
  ylab("Frequency") + xlab("\"Custom\" Comorbidity Group") + theme_minimal() + theme(axis.text.x=element_text(angle=80, vjust=0.62))
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


```r
icd_dom_recode_wide %>% select(-SUBJECT_ID) %>% cor() %>% corrplot::corrplot(type="upper", diag=F, order="hclust")
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-18-1.png)<!-- -->




### Patient Level Data

In addition to the exploration of the diagnoses codes, I began to aggregate and explore some of what I refer t o as the "patient-level data". That is, data elements corresponding to characteristics of each patient. 


```r
pts_raw <- read.csv(paste0(lib, "PATIENTS.csv"))
```


```r
pts_raw %>% count(GENDER) %>% ggplot(aes(x=reorder(GENDER, -n), y=n)) + 
  geom_bar(stat="identity") + ylab("Frequency") + xlab("Gender") +
  ggtitle("Frequency of Gender in MIMIC Data")
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

```r
pts_raw %>% mutate(MortalityType=
                     case_when(EXPIRE_FLAG==0 ~ "Survival",
                               DOD_HOSP!="" ~ "In-Hospital",
                               TRUE ~ "Other (SSN)")
                   )  %>% ggplot(aes(x=as.factor(EXPIRE_FLAG), fill=MortalityType)) + 
  geom_bar() + ylab("Frequency") + xlab("Mortality Status") + scale_fill_brewer(palette=5, type = "qual") +
  ggtitle("Frequency of Mortality Status (including Data Source)")
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-20-2.png)<!-- -->


```r
los_graph <- admit %>% group_by(SUBJECT_ID) %>% filter(ADMITTIME==max(ADMITTIME)) %>% ungroup() %>% 
  mutate(GenLOS=difftime(DISCHTIME, ADMITTIME, units = "days") %>% as.numeric()) %>% select(SUBJECT_ID, GenLOS, DISCHTIME, ADMITTIME) # %>% 
  
los_graph %>%   ggplot(aes(x=GenLOS)) + geom_density() + theme_minimal() +
  xlab("General Hospital Length of Stay") + ylab("Density") + 
  ggtitle("Distribution of General Hospital Length of Stay (Days)") +
  annotate(geom="text", x=100, y=0.05, label=paste0("Length of stay values ranged from ", min(los_graph$GenLOS), " to ", max(los_graph$GenLOS), "days.")) +
  annotate(geom="text", x=114, y=0.04, label=paste0("Our length of stay values have a mean of ", round(mean(los_graph$GenLOS), 2), " and variance \n of ",
                                                    round(var(los_graph$GenLOS),2), ", suggesting overdispersion of this variable."))
```

![](Treelet_Analysis_EDA_files/figure-html/unnamed-chunk-21-1.png)<!-- -->

