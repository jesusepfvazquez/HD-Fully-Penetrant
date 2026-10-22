################################################################################
## Extract 2019 single-age death rates for every country from the United Nations
## World Population Prospects 2024 (WPP 2024), for the global supplementary
## analysis in global_penetrance.R.
##
## Inputs (R package wpp2024, version 1.1-3, United Nations Population Division),
## downloaded from https://github.com/PPgp/wpp2024/tree/main/data :
##   mx1dt.rda       single-year age-specific death rates, ages 0-99 and 100+
##   popAge1dt.rda   single-year population by age (thousands)
##   e01dt.rda       life expectancy at birth
##   UNlocations.txt location names, regions, and location type
## Accessed 2026-10-21. The WHO GHO API used for the main analysis was not
## reachable from this session, and WPP 2024 gives one consistent source for
## every country.
##
## Run from the code/ folder with the four files in raw/:
##   Rscript extract_wpp2024.R
## Output: data/wpp2024_mx_2019.csv (one row per country and single year of age)
################################################################################

raw <- "raw"
load(file.path(raw, "mx1dt.rda")); load(file.path(raw, "popAge1dt.rda")); load(file.path(raw, "e01dt.rda"))
mx1dt <- as.data.frame(mx1dt); popAge1dt <- as.data.frame(popAge1dt); e01dt <- as.data.frame(e01dt)
loc <- read.delim(file.path(raw, "UNlocations.txt"), stringsAsFactors = FALSE)

## Countries and areas only (location_type 4), with at least 90,000 people in 2019,
## the population threshold the WHO uses for its life tables.
loc <- loc[loc$location_type == 4, c("country_code", "name", "reg_name", "area_name")]
pop <- aggregate(pop ~ country_code, data = popAge1dt[popAge1dt$year == 2019, ], FUN = sum)
loc <- merge(loc, pop, by = "country_code")
loc <- loc[loc$pop >= 90, ]

mx <- mx1dt[mx1dt$year == 2019 & mx1dt$country_code %in% loc$country_code,
            c("country_code", "age", "mxB")]
e0 <- e01dt[e01dt$year == 2019, c("country_code", "e0B")]
out <- merge(merge(loc, e0, by = "country_code"), mx, by = "country_code")
out <- out[order(out$country_code, out$age), ]
names(out) <- c("country_code", "name", "region", "area", "pop2019_thousands",
                "e0_wpp", "age", "mx")
dir.create("data", showWarnings = FALSE)
write.csv(out, "data/wpp2024_mx_2019.csv", row.names = FALSE)
cat(length(unique(out$country_code)), "countries written\n")
