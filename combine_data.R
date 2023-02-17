library(tidyverse)
library(covidData)
library(MMWRweek)
library(lubridate)
library(padr)

hosp_dat <- readr::read_csv("./data/covid_data_cached_2023-02-13.csv")
clade_dat <- readr::read_csv("./data/clade_data_cached_2023-02-15.csv")

# padding clade data
date_padded_dat <- clade_dat %>%
    padr::pad() 
unique_clades <- unique(clade_dat$Nextstrain_clade)
unique_dates <- unique(date_padded_dat$date)
for(i in 1:length(unique_dates)){
    date_dat <- date_padded_dat %>%
        filter(date == unique_dates[i]) 
    existing_clades <- unique(date_dat$Nextstrain_clade)
    new_dat <- date_dat %>%
        add_row(date=unique_dates[i],
                Nextstrain_clade=setdiff(unique_clades,existing_clades),
                total_samples=date_dat$total_samples[1]) 
    date_padded_dat <- rbind(date_padded_dat,new_dat) 
} 
# clean up and fill in columns
date_padded_dat_final <- date_padded_dat %>%
    distinct() %>%
    filter(!is.na(Nextstrain_clade)) %>%
    arrange(date,Nextstrain_clade) %>%
    select(-c("diff_log","nonvar_diff_log")) %>%
    mutate(clade_samples = ifelse(is.na(clade_samples),0,clade_samples),
           clade_pct=ifelse(is.na(clade_pct),1/(10^5),clade_pct))
readr::write_csv(date_padded_dat_final, "data/padded_clade_data_cached_2023-02-15.csv")
