# install.packages('modelsummary')
library(modelsummary)

load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/dat_matched_PM25_bis.RData')
# load('/Users/alicesommer/Desktop/Bureau/DOCTORATE/data_pipeline_microbiome/dat_matched_smoke_bis.RData')

matched_df <- matched_df[,c("W","u3talteru", "u3csex", "u3tbmi",  
                            "u3talkkon", "u3tedyrs",
                            "u3tcigsmk", 
                            "u3tdiabet","u3tphys")]

matched_df$u3csex[matched_df$u3csex == 0] <- "Male"
matched_df$u3csex[matched_df$u3csex == 1] <- "Female"

matched_df$u3tphys[matched_df$u3tphys == 0] <- "No"
matched_df$u3tphys[matched_df$u3tphys == 1] <- "Yes"

matched_df$u3tcigsmk[matched_df$u3tcigsmk == 1] <- "Smoker"
matched_df$u3tcigsmk[matched_df$u3tcigsmk == 2] <-  "Ex-Smoker"
matched_df$u3tcigsmk[matched_df$u3tcigsmk == 3] <-  "Never-Smoker"

matched_df$u3tdiabet[matched_df$u3tdiabet == 0] <- "No"
matched_df$u3tdiabet[matched_df$u3tdiabet == 1] <- "Yes"

datasummary_balance(~W,
                    data = matched_df,
                    output = 'latex')
