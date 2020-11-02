library(tidyverse)
library(here)
library(gghighlight)
library(treelet)

cohort_full <-  read.csv(here("/Grad School/Masters/Thesis/Analysis/Data", "cohort_full.csv"))
colnames(cohort_full) <- cohort_full %>% colnames() %>% gsub(pattern = "X", "", x = .)

set.seed(2824)

hold_out_pts <- sample(1:nrow(cohort_full), size=nrow(cohort_full)/5, replace = F)

holdout_test <- cohort_full[hold_out_pts,]
# nrow(holdout_test)

cv_data <- cohort_full[setdiff(1:nrow(cohort_full), hold_out_pts),]
# nrow(cv_data)

(nrow(holdout_test) + nrow(cv_data)) == nrow(cohort_full)

cv_data$fold <- sample(c(rep(1, ceiling(nrow(cv_data)/5)),
                         rep(2, ceiling(nrow(cv_data)/5)),
                         rep(3, ceiling(nrow(cv_data)/5)),
                         rep(4, ceiling(nrow(cv_data)/5)),
                         rep(5, ceiling(nrow(cv_data)/5))
                         ),
                       size=nrow(cv_data), replace=F
                       )

table(cv_data$fold)

treelet_process <- function(x_mat, cov_mat){

  tt_results <- tt_results <- treelet::Run_JTree(cov_mat, nrow(cov_mat)-1, 1:nrow(cov_mat)-1)
  energy <- list()

      for(L in 1:length(tt_results$basis)) {
        
        basisk <- tt_results$basis[[L]]
        w_x <- t(basisk) %*% t(x_mat)
    
          num_vec <- rowSums(abs(w_x)^2)
          den_vec <- x_mat^2 %>% colSums()
              names(num_vec) <- NULL
              names(den_vec) <- NULL
    
        energy[[L]] <- matrix(c(1:ncol(x_mat), num_vec / den_vec), ncol=2, dimnames = list(NULL, c("W_i", "Energy")))
    }

  # Creating blank objects  
    optimal_L <- matrix(c(1:length(energy), rep(NA, length(energy))), nrow=length(energy), dimnames = list(NULL, c("K", "Optimal L")))
    retained_fts <- rep(list(rep(list(rep(NA, length(energy))), length(energy))), length(energy))

  # Reordering the energy matrices in descending order of normed energy score
    energy_ordered <- lapply(1:length(energy), function(L) energy[[L]][energy[[L]][,2] %>% order(decreasing = T),])

  # Identifying optimal L
    optimal_L <- matrix(c(1:length(energy_ordered),
                      sapply(1:length(energy_ordered), 
                             function(K) which.max(sapply(1:length(energy), 
                                                          function(x) sum(energy_ordered[[x]][1:K,2])
                                                          )
                                                   ))),
                      ncol=2, dimnames=list(NULL, c("GivenK", "OptimalBasis_L")))

  # And retained fts
    retained_fts <- lapply(1:length(energy_ordered),
                            function(x) energy_ordered[[optimal_L[x,2]]][optimal_L[1:x,1], 1])

  return(list(basis_mats=tt_results$basis,
              optimal_params=optimal_L,
              retained_fts=retained_fts))
}


#### Mortality
mortality_performance <- read.csv(here("Grad School/Masters/Thesis/Treelet/Results/Treelet_KLOpt_WithinCVLoop/MortalityModel_CVPerformance_NoLOS_NewKLCode.csv"))

k_1sd <- mortality_performance[mortality_performance$BS_TestAvg<=(min(mortality_performance$BS_TestAvg) + sd(mortality_performance$BS_TestAvg)), ] %>% .[1,1]

mortality_performance <- mortality_performance %>% mutate(ParamFlag = 
                                   case_when(
                                     BS_TestAvg==min(BS_TestAvg) ~ "Minimizes Briers Score",
                                     K==k_1sd ~ "More Sparse Parameter",
                                     TRUE ~ NA_character_
                                   )) %>% ungroup()

mortality_performance[!is.na(mortality_performance$ParamFlag),]

cv_xmat <- cv_data  %>% select(matches("0|1|2|4|5|6|9")) %>% select(-Yr1Readmit) %>% as.matrix()
cv_cov <- cov(cv_xmat)

tt_fnc_mortality <- treelet_process(cv_xmat, cv_cov)

tt_fnc_mortality$optimal_params[c(123, 174),]

final_basis_mortality <- tt_fnc_mortality$basis_mats[[57]][,1:123] %>% 
	as.data.frame() %>% mutate(LabelIndex = row_number(),
                                 RowMissCount = rowSums(.==0)) %>%
			filter(RowMissCount<123)

dim(final_basis_mortality)


labels_df <- cv_readmit_cov %>% colnames() %>% 
		data.frame(code = ., label=1:ncol(cv_readmit_cov))

loading_mat_mortality <- merge(final_basis_mortality, labels_df,
						 all.x=T, by.y="label", by.x="LabelIndex")

write.csv(loading_mat_mortality, 
here("Grad School/Masters/Thesis/Treelet/Results/Treelet_KLOpt_WithinCVLoop/LoadingMatrix_Mortality.csv"))



ggplot(mortality_performance, aes(x=K, y=BS_TestAvg, color = as.factor(ParamFlag))) +
  geom_line(lwd=1.1, alpha=0.6) + geom_point(size=2.5) +
  theme_minimal() + ggtitle("In-Hospital Mortality Model") + 
  xlab("Value of Parameter K") + ylab("Average Briers Score (Across 5 Test Folds)") + 
  gghighlight(ParamFlag!=0) + labs(color="Optimal Parameters") +
  scale_color_brewer(type = "qual", palette = 6) + 
  theme(legend.position=c(0.75, 0.75), text = element_text(size=13.5))


### Readmission
readmit_performance <- read.csv(here("Grad School/Masters/Thesis/Treelet/Results/Treelet_KLOpt_WithinCVLoop/ReadmissionModel_CVPerformance_NewKLCode.csv"))


readmit_performance <- readmit_performance %>% mutate(BS_TestAvg = 
                            rowMeans(select(readmit_performance, starts_with("BS"))),
                          AUC_TestAvg =
                            rowMeans(select(readmit_performance, starts_with("AUC"))))

k_1sd_readmit <- readmit_performance[readmit_performance$BS_TestAvg<=(min(readmit_performance$BS_TestAvg) + sd(readmit_performance$BS_TestAvg)), ] %>% .[1,1]


readmit_performance <- readmit_performance %>% mutate(ParamFlag = 
                                   case_when(
                                     BS_TestAvg==min(BS_TestAvg) ~ "Minimizes Briers Score",
                                     K==k_1sd_readmit ~ "More Sparse Parameter",
                                     TRUE ~ NA_character_
                                   )) %>% ungroup()


readmit_performance [!is.na(readmit_performance $ParamFlag),]

cv_readmit_xmat <- cv_data  %>% select(matches("0|1|2|4|5|6|9")) %>% select(-Yr1Readmit) %>% as.matrix()
cv_readmit_cov <- cov(cv_readmit_xmat)

tt_fnc_readmit <- treelet_process(cv_readmit_xmat, cv_readmit_cov)

tt_fnc_readmit$optimal_params[c(5, 30),]

# Matrix of loadings 
final_basis <- tt_fnc_readmit$basis_mats[[52]][,1:5] %>% 
	as.data.frame() %>% mutate(LabelIndex = row_number()) %>%
			filter(V1!=0 | V2!=0 | V3!=0 | V4!=0 | V5!=0)

labels_df <- cv_readmit_cov %>% colnames() %>% 
		data.frame(code = ., label=1:ncol(cv_readmit_cov))

loading_mat_readmit <- merge(final_basis, labels_df, all.x=T, by.y="label", by.x="LabelIndex")
write.csv(loading_mat_readmit, 
here("Grad School/Masters/Thesis/Treelet/Results/Treelet_KLOpt_WithinCVLoop/LoadingMatrix_Readmit.csv"))


ggplot(readmit_performance, aes(x=K, y=BS_TestAvg, color=as.factor(ParamFlag))) +
  geom_line(lwd=1.1, alpha=0.6) + geom_point(size=2.5) +
  theme_minimal() + ggtitle("Hospital Readmission Model") + 
  xlab("Value of Parameter K") + ylab("Average Briers Score (Across 5 Test Folds)") + 
  gghighlight(ParamFlag!=0) + labs(color="Optimal Parameters") +
  scale_color_brewer(type = "qual", palette = 6) + 
  theme(legend.position=c(0.5, 0.8), text = element_text(size=13.5))
