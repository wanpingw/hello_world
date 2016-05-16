
# Sophia He, Kay Wang
library(dplyr)
library(RPostgreSQL)
library(lubridate)
library(tidyr)
library(PhViD)
library(magrittr)
library(plyr)

hcopen <- src_postgres(host = "shiny.hc.local", user = "hcreader", dbname = "hcopen", password = "canada1")

# read tables, select columns
cv_reports <- as.data.frame(tbl(hcopen, sql("SELECT * FROM cv_reports")),n=-1) %>% dplyr:: select(REPORT_ID, DATRECEIVED_CLEAN) %>% mutate(month = floor_date(DATRECEIVED_CLEAN, "month")) 
cv_reports <- filter(cv_reports, ymd(month) >=ymd("2015-01-01"))
cv_report_drug <- as.data.frame(tbl(hcopen, sql("SELECT * FROM cv_report_drug")),n=-1)
cv_drug_product_ingredients <- as.data.frame(tbl(hcopen, sql("SELECT * FROM cv_drug_product_ingredients")),n=-1)
cv_reactions <- as.data.frame(tbl(hcopen, sql("SELECT * FROM cv_reactions")),n=-1) %>% select(REPORT_ID, PT_NAME_ENG)


# sort tables
cv_drug_product_ingredients_sorted <- arrange(cv_drug_product_ingredients,DRUG_PRODUCT_ID)
cv_report_drug_sorted <- arrange(cv_report_drug,REPORT_ID)


# 1.select & merge only common varaibles (DRUG_PRODUCT_ID,REPORT_ID & ACTIVE_INGREDIENT_NAME) 
#   from cv_report, cv_report_drug & cv_drug_product_ingredients
# join ingred names and report_id using drug_product_ID
top25_S1 <- dplyr:: select(cv_drug_product_ingredients_sorted, DRUG_PRODUCT_ID, ACTIVE_INGREDIENT_NAME) 
top25_S2 <- dplyr:: select(cv_report_drug_sorted,DRUG_PRODUCT_ID, REPORT_ID)
top25_S3 <- top25_S1 %>% left_join(top25_S2)  %>% left_join(cv_reports)
top25_S3_sorted <- arrange(top25_S3,REPORT_ID) %>% subset(!is.na(month))


# there are many ingredients that have the same frequency
# the top 25 reported ingredients  
top25_ingd <-as.data.frame(table(top25_S3_sorted$ACTIVE_INGREDIENT_NAME)) # obtain the frequency of ingredients
top25_ingd_sorted <- top25_ingd[order(-top25_ingd$Freq),]
top25_ingd_sorted <- rename(top25_ingd_sorted, ACTIVE_INGREDIENT_NAME= Var1) # ingred/frequency table
#join top 25 ingredients with reactions
top_ingrd_final <- top25_ingd_sorted[1:25,] %>% left_join(top25_S3_sorted) %>% left_join(cv_reactions) %>% dplyr:: select(ACTIVE_INGREDIENT_NAME, PT_NAME_ENG)
# frequency table for all pairs
combo_freq <- data.frame(table(top_ingrd_final$ACTIVE_INGREDIENT_NAME, top_ingrd_final$PT_NAME_ENG)) %>% subset( Freq > 0) 
combo_freq <- rename(combo_freq, ACTIVE_INGREDIENT_NAME= Var1)
combo_freq <- rename(combo_freq, PT_NAME_ENG= Var2)

bayes_test <- subset.data.frame(combo_freq, ACTIVE_INGREDIENT_NAME=="rituximab" )

#bayesian analysis using package PhVid
bayes_table <- as.PhViD(bayes_test, MARGIN.THRES = 1)
PRR(bayes_table, RR0 = 1, MIN.n11 = 3, DECISION = 1,DECISION.THRES = 0.05, RANKSTAT = 1)
bayes_result <- BCPNN(bayes_table, RR0 = 1, MIN.n11 = 3, DECISION = 1,DECISION.THRES = 0.05, RANKSTAT = 2, MC=FALSE)
gps_result <- GPS(bayes_table, RR0 = 1, MIN.n11 = 3, DECISION = 1, DECISION.THRES = 0.05,RANKSTAT = 1, TRONC = FALSE, TRONC.THRES = 1,PRIOR.INIT = c(alpha1 = 0.2, beta1 = 0.06, alpha2 = 1.4,beta2 = 1.8, w = 0.1), PRIOR.PARAM = NULL)

PhViD.search(bayes_result, DRUG= "prednisone", EVENT= "Arterial haemorrhage")
PhViD.search(bayes_table, DRUG= "acetaminophen", EVENT= "Pediatric ingestion errors")
