# Synopsis

Abstract—Breast cancer treatments using methylation inhibitors are an effective new therapeutic option for breast cancer
patients. This repository has the implementation of different regression models, i.e. proportional hazards regression, elastic net regression, ridge regression, and lasso regression, to identify genes of which methylation rates are strongly correlated to breast cancer survival. With each of the regression models we identified genes of which high methylation rates are strongly favorably correlated with breast cancer survival, and genes of which high methylation rates are strongly adversely correlated with breast cancer survival. A better understanding of the relationship between DNA methylation rates and breast cancer survival can assist in the development of patient-tailored therapy strategies, and the discovery of therapeutic targets.

# Task Definition and Input/Ouput

We used two datasets from TCGA (The Cancer Genomic Atlas); one contains survival data of cancer patients, the other contains genomic data, copy number variation (CNV) data, and methylation data of cancer patients. The survival dataset contains the type of cancer (11 different cancer types in total), the ”time to last contact or event,” and whether the event occurred (1: death) or not (0: no death at time of last contact) for 8089 different patients. Many of the survival times are censored, i.e. the time of observation was cut off before death occurred; this indicates that the patient either was still alive at the end of the study or that the patient withdrew from the study before the end of the study. The dataset with genomic data, copy number variation (CNV) data, and methylation data contains methylation data for more than 16,500 different genes of over 1,000 different cancer patients.

Input: survival data of breast cancer patients including the ”time to last contact or event,” and whether the event occurred (1: death) or not (0: no death at time of last contact), and the methylation data of these breast cancer
patients. 
Output: a set of genes of which the methylation rates are strong predictors of breast cancer survival.
