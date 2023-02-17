library(tidyverse)
library(covidData)
library(MMWRweek)
library(lubridate)

for (pathogen in c("covid")) {
    data <- covidData::load_data(
        spatial_resolution = c("state"),
        temporal_resolution = "daily",
        location_code = "25",
        drop_last_date = TRUE,
        measure = ifelse(pathogen == "covid",
                         "hospitalizations",
                         "flu hospitalizations")
    ) %>%
        dplyr::group_by(location) %>%
        dplyr::filter(date >= "2020-03-01", date < max(date)) %>%
        dplyr::transmute(location, date, hosps = inc) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(location, date)
    
    data <- data %>%
        dplyr::left_join(
            covidData::fips_codes %>%
                dplyr::select(location, location_name, population),
            by = "location") %>%
        dplyr::mutate(
            pop100k = population / 100000,
            hosp_rate = hosps / pop100k
        )

    readr::write_csv(data, paste0("data/", pathogen, "_data_cached_", Sys.Date(), ".csv"))
}

## clean genomic data
genbank_global <- read_tsv("data/metadata.tsv")
MA_dat <- genbank_global %>%
    filter(country=="USA",
           division == "Massachusetts",
           host == "Homo sapiens") %>%
    mutate(date = ymd(date),
           date_submitted = ymd(date_submitted),
           reporting_lag = as.numeric(date_submitted - date)) %>%
    select(strain, virus, pango_lineage, Nextstrain_clade, GISAID_clade,## info on the virus
           region, country, division, location, ## info on location
           host, sampling_strategy, ## info about the sample
           date, date_submitted, reporting_lag) %>%
    # formatting
    filter(!is.na(date),
           !is.na(pango_lineage),
           pango_lineage!="unclassifiable",
           !is.na(Nextstrain_clade),
           Nextstrain_clade!="?") %>%
    mutate(year = MMWRweek(date)$MMWRyear,
           epiweek = MMWRweek(date)$MMWRweek,
           location_name=division) %>% ## using this year function as it handles the dates at shoulders better
    group_by(date) %>%
    mutate(total_samples = n()) %>%
    ungroup() %>%
    group_by(date, pango_lineage) %>%
    mutate(lineage_samples = n(),
           lineage_pct = lineage_samples/total_samples,
           sampling_rate = total_samples/data$pop100k[1]) %>%
    ungroup() %>%
    group_by(date, Nextstrain_clade) %>%
    mutate(clade_samples = n(),
           clade_pct = clade_samples/total_samples) %>%
    ungroup() %>%
    select("pango_lineage", "Nextstrain_clade","GISAID_clade",
           "location_name", "date", "reporting_lag", "year", 
           "epiweek", "total_samples","sampling_rate",
           "lineage_samples", "lineage_pct", "clade_samples","clade_pct") 
readr::write_csv(MA_dat, paste0("data/genomic_data_cached_", Sys.Date(), ".csv"))

## clean genomic data
MA_dat <- readr::read_csv("data/genomic_data_cached_2023-02-15.csv")
clade_dat <- MA_dat %>%
    select(date,clade_samples,total_samples,Nextstrain_clade,clade_pct) %>%
    distinct() %>%
    dplyr::group_by(Nextstrain_clade) %>%
    dplyr::arrange(date, .by_group = TRUE) %>%
    dplyr::mutate(diff_log=log(clade_samples)-log(lag(clade_samples))) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(date) %>%
    # add average diff log
    dplyr::group_by(date) %>%
    dplyr::mutate(nonvar_diff_log = (sum(diff_log, na.rm = TRUE)-diff_log)/(n()-1)) %>%
    dplyr::ungroup() 
readr::write_csv(clade_dat, "data/clade_data_cached_2023-02-15.csv")
#rm("genbank_global")
# #file of national level human genomic data samples. Filtered in the file "filtering genomic file.RMD"
# MA_data_raw <- MA_dat %>% 
#     filter(!is.na(date)) %>% #filter out rows with no date 
#     mutate(year = MMWRweek(date)$MMWRyear) %>% ## using this year function as it handles the dates at shoulders better
#     mutate(epiweek = MMWRweek(date)$MMWRweek) 
# write_csv(MA_data_raw, file="./data/MA_raw_data.csv")

# d2 <- MA_data_raw %>% 
#     dplyr::filter(!is.na(pango_lineage)) %>%
#     group_by(year, epiweek, pango_lineage) %>% 
#     summarize(lineage_samples = n()) %>% 
#     ungroup() %>% 
#     group_by(year, epiweek) %>% 
#     mutate(total_samples = sum(lineage_samples)) %>% 
#     ungroup() %>% 
#     mutate(epiweek_year = paste(year, formatC(epiweek, width=2, flag="0"), sep="_"),
#            epidate = MMWRweek2Date(year, epiweek),
#            lineage_pct = lineage_samples/total_samples) 

# write_csv(d2, file="./data/MA_modeling_data.csv")
# if (FALSE) {
#     deaths <- covidData::load_data(
#         spatial_resolution = c("state", "national"),
#         temporal_resolution = "daily",
#         measure = "deaths"
#     ) %>%
#         dplyr::filter(date >= "2020-10-01") %>%
#         dplyr::rename(deaths = inc)
#     data <- hosps %>%
#         dplyr::left_join(deaths %>% dplyr::select(location, date, deaths),
#                          by = c("location", "date"))
# } else {
#     data <- hosps
# }
