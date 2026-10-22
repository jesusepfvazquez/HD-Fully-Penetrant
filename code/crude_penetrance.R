################################################################################
## Crude versus net penetrance of Huntington's disease (HD) alleles under
## country-specific background mortality
##
## Produces every computed number in main.tex:
##   output/table1.csv          Table 1 (net and crude penetrance)
##   output/validation_*.csv    checks of the transcribed WHO data
##   output/sensitivity.csv     sensitivity analyses
##   output/checks.txt          text claims checked against the data
##   output/key-numbers.tex     LaTeX macros for the numbers quoted in the text
##   ../figures/figure1-crude-penetrance.{pdf,png}
##
## Run from the code/ folder:  Rscript crude_penetrance.R
## Base R only; no packages required.
################################################################################

dir.create("output", showWarnings = FALSE)

#-------------------------------------------------------------------------------
# 1. Inputs
#-------------------------------------------------------------------------------

## 1a. WHO Global Health Observatory life tables, both sexes.
##     nmx = age-specific death rate (indicator LIFE_0000000029)
##     ex  = life expectancy at the start of the age group (indicator LIFE_0000000035)
##     Source: https://ghoapi.azureedge.net/api/LIFE_0000000029 and .../LIFE_0000000035
##     filtered by SpatialDim (ISO3 country), TimeDim (year), Dim1 = SEX_BTSX.
##     Accessed 2026-10-09. SOURCES.md lists the exact query URLs.
lt <- read.csv("data/who_life_tables.csv", stringsAsFactors = FALSE)
lt <- lt[order(lt$country, lt$year, lt$age_start), ]

## Table order: high-income comparators, South America, Mexico and Central America
groups <- list("High-income comparators" = c("ESP", "SWE", "GBR", "CAN", "USA"),
               "South America" = c("CHL", "PER", "COL", "ARG", "BRA"),
               "Mexico and Central America" = c("MEX", "SLV", "HND", "GTM"))
countries <- unlist(groups, use.names = FALSE)
central_am <- c("SLV", "HND", "GTM")
cnames <- c(MEX = "Mexico", ESP = "Spain", SWE = "Sweden", GBR = "United Kingdom", CAN = "Canada", USA = "United States",
            ARG = "Argentina", CHL = "Chile", COL = "Colombia", PER = "Peru",
            BRA = "Brazil", SLV = "El Salvador", HND = "Honduras", GTM = "Guatemala")

## 1b. Age-at-onset model of Langbehn et al. (2004), Clin Genet 65:267-277.
##     Onset age ~ logistic(mean = mu(CAG), variance = sigma^2(CAG)):
##       mu(CAG)      = 21.54 + exp(9.556 - 0.1460 * CAG)
##       sigma^2(CAG) = 35.55 + exp(17.72 - 0.3269 * CAG)
##     The mean constants are confirmed in published secondary sources
##     (SOURCES.md). The variance constants MUST be checked against the
##     original paper and its 2004 erratum; Section 6 shows how much the
##     results move if sigma is 20% smaller or larger.
langbehn <- function(cag, sd_scale = 1) {
  mu <- 21.54 + exp(9.556 - 0.1460 * cag)
  sigma <- sqrt(35.55 + exp(17.72 - 0.3269 * cag)) * sd_scale
  s <- sigma * sqrt(3) / pi          # logistic scale parameter
  list(mu = mu, sigma = sigma,
       cdf = function(t) plogis(t, location = mu, scale = s),
       pdf = function(t) dlogis(t, location = mu, scale = s))
}

#-------------------------------------------------------------------------------
# 2. Background survival from a piecewise-constant hazard
#-------------------------------------------------------------------------------

## Within each WHO age group the death rate nmx is treated as a constant hazard.
## The open interval 85+ uses the WHO 85+ rate (main analysis) or 1/e85
## (sensitivity analysis, Section 6).
breaks <- c(0, 1, seq(5, 85, 5))
make_hazard <- function(cc, yr = 2019, open = "nmx") {
  d <- lt[lt$country == cc & lt$year == yr, ]
  h <- d$nmx
  if (open == "inv_e85") h[length(h)] <- 1 / d$ex[nrow(d)]
  function(t) h[findInterval(t, breaks)]
}

## Survival from age a0 to age t (t >= a0): exp(-integral of the hazard from a0 to t),
## computed exactly for a piecewise-constant hazard.
make_surv <- function(haz) {
  function(t, a0) {
    sapply(t, function(tt) {
      if (tt <= a0) return(1)
      knots <- sort(unique(c(a0, breaks[breaks > a0 & breaks < tt], tt)))
      exp(-sum(haz(knots[-length(knots)]) * diff(knots)))
    })
  }
}

## Life expectancy at age a: integral of S(t | a) from a to infinity
## (upper limit 130, where survival is effectively zero).
life_exp <- function(cc, a, yr = 2019, open = "nmx") {
  S <- make_surv(make_hazard(cc, yr, open))
  integrate(function(t) S(t, a), a, 130, subdivisions = 2000)$value
}

#-------------------------------------------------------------------------------
# 3. Validation of the transcribed WHO data
#-------------------------------------------------------------------------------

## 3a. Age-group check. With deaths spread evenly within a group (a_x = n/2),
##     the WHO ex values imply a survival probability p for each group:
##       e_x = n (1 + p) / 2 + p e_{x+n}   =>   p = (e_x - n/2) / (n/2 + e_{x+n})
##     and an implied death rate m = 2 (1 - p) / (n (1 + p)).
##     Comparing implied m with the transcribed nmx group by group detects a
##     transcription error in any single age group. Ages 0 and 1-4 are skipped
##     because WHO uses a_x well below n/2 at those ages.
val_group <- do.call(rbind, lapply(countries, function(cc) {
  d <- lt[lt$country == cc & lt$year == 2019, ]
  do.call(rbind, lapply(3:(nrow(d) - 1), function(i) {
    n <- d$age_start[i + 1] - d$age_start[i]
    p <- (d$ex[i] - n / 2) / (n / 2 + d$ex[i + 1])
    m_impl <- 2 * (1 - p) / (n * (1 + p))
    data.frame(country = cc, age_group = d$age_group[i], nmx_who = d$nmx[i],
               nmx_implied = m_impl, rel_diff = m_impl / d$nmx[i] - 1)
  }))
}))
write.csv(val_group, "output/validation_by_age_group.csv", row.names = FALSE)

## 3b. Life-expectancy check: e0, e40, e60 recomputed from nmx vs WHO ex.
who_ex <- function(cc, a, yr = 2019) lt$ex[lt$country == cc & lt$year == yr & lt$age_start == a]
val_ex <- do.call(rbind, lapply(countries, function(cc) {
  data.frame(country = cc,
             e0_who = who_ex(cc, 0),   e0_ours = life_exp(cc, 0),
             e40_who = who_ex(cc, 40), e40_ours = life_exp(cc, 40),
             e60_who = who_ex(cc, 60), e60_ours = life_exp(cc, 60))
}))
write.csv(val_ex, "output/validation_life_expectancy.csv", row.names = FALSE)

#-------------------------------------------------------------------------------
# 4. Net and crude penetrance
#-------------------------------------------------------------------------------

## Net penetrance by age t, among carriers without onset at age a0:
##   [F(t) - F(a0)] / [1 - F(a0)]
## Crude (lifetime) penetrance, among carriers alive and without onset at a0:
##   integral from a0 to infinity of f(u) S(u | a0) du / [1 - F(a0)]
## (Equation 1 in main.tex, conditioned on age a0.)
net_pen <- function(cag, t, a0 = 18, sd_scale = 1) {
  m <- langbehn(cag, sd_scale)
  (m$cdf(t) - m$cdf(a0)) / (1 - m$cdf(a0))
}
crude_pen <- function(cag, cc, a0 = 18, yr = 2019, open = "nmx", sd_scale = 1) {
  m <- langbehn(cag, sd_scale)
  S <- make_surv(make_hazard(cc, yr, open))
  num <- integrate(function(u) m$pdf(u) * S(u, a0), a0, 130, subdivisions = 2000)$value
  num / (1 - m$cdf(a0))
}

cags <- 40:45
tab <- do.call(rbind, lapply(cags, function(k) {
  m <- langbehn(k)
  crude <- sapply(countries, function(cc) crude_pen(k, cc))
  data.frame(cag = k, median_onset = m$mu, sd_onset = m$sigma,
             net65 = net_pen(k, 65), net75 = net_pen(k, 75), net85 = net_pen(k, 85),
             t(crude), gap_ESP_MEX = unname(crude["ESP"] - crude["MEX"]))
}))
write.csv(tab, "output/table1.csv", row.names = FALSE)
print(round(tab, 3))

#-------------------------------------------------------------------------------
# 5. Text claims checked against the data
#-------------------------------------------------------------------------------

sink("output/checks.txt")
cat("Claim: death rates in Mexico exceed those in Spain in every age group (2019)\n")
dm <- lt[lt$country == "MEX" & lt$year == 2019, ]
ds <- lt[lt$country == "ESP" & lt$year == 2019, ]
print(data.frame(age = dm$age_group, MEX = dm$nmx, ESP = ds$nmx, ratio = round(dm$nmx / ds$nmx, 2)))
cat("All ages higher in Mexico:", all(dm$nmx > ds$nmx), "\n\n")

cat("Claim: crude penetrance is lowest in Guatemala, Honduras, or El Salvador at every CAG,\n")
cat("and all three Central American countries are below Mexico at every CAG\n")
print(setNames(sapply(seq_along(cags), function(i) names(which.min(unlist(tab[i, countries])))), cags))
cat("All Central American < Mexico at every CAG:",
    all(sapply(seq_along(cags), function(i) all(unlist(tab[i, central_am]) < tab$MEX[i]))), "\n")
cat("Mexico rank (1 = lowest) by CAG:",
    sapply(seq_along(cags), function(i) rank(unlist(tab[i, countries]))[["MEX"]]), "\n\n")

cat("Claim: crude penetrance at 40 CAG tracks remaining life expectancy at 40\n")
e40 <- sapply(countries, function(cc) who_ex(cc, 40))
c40 <- unlist(tab[tab$cag == 40, countries])
print(data.frame(e40 = round(e40, 1), crude40 = round(c40, 3))[order(e40), ])
cat("Spearman correlation:", round(cor(e40, c40, method = "spearman"), 3), "\n\n")

cat("Claim: net penetrance by 65 for 40 and 50 repeats\n")
cat("40 CAG:", round(net_pen(40, 65), 3), "  50 CAG:", round(net_pen(50, 65), 4), "\n\n")

cat("Claim: share of carriers of 40 repeats dying before onset (1 - crude)\n")
for (cc in c("MEX", "ESP", "GTM")) {
  q <- 1 - tab[tab$cag == 40, cc]
  cat(cc, ": ", round(q, 3), " = about 1 in ", round(1 / q, 1), "\n", sep = "")
}
sink()

#-------------------------------------------------------------------------------
# 6. Sensitivity analyses
#-------------------------------------------------------------------------------

sens <- function(label, f) {
  do.call(rbind, lapply(c(40, 42, 45), function(k) {
    v <- sapply(c("MEX", "ESP", "GBR"), function(cc) f(k, cc))
    data.frame(analysis = label, cag = k, MEX = v[1], ESP = v[2], GBR = v[3],
               gap_ESP_MEX = v[2] - v[1], row.names = NULL)
  }))
}
## brute-force grid sum, to check integrate()
grid_crude <- function(k, cc, a0 = 18, dt = 0.001) {
  m <- langbehn(k); h <- make_hazard(cc)
  u <- seq(a0, 130, by = dt)
  S <- exp(-cumsum(c(0, h(u[-length(u)]) * dt)))
  sum(m$pdf(u) * S) * dt / (1 - m$cdf(a0))
}
sens_tab <- rbind(
  sens("main (2019, start 18, 85+ = WHO nmx)", function(k, cc) crude_pen(k, cc)),
  sens("85+ hazard = 1/e85",                   function(k, cc) crude_pen(k, cc, open = "inv_e85")),
  sens("2023 life tables",                     function(k, cc) crude_pen(k, cc, yr = 2023)),
  sens("start age 30",                         function(k, cc) crude_pen(k, cc, a0 = 30)),
  sens("onset SD x 0.8",                       function(k, cc) crude_pen(k, cc, sd_scale = 0.8)),
  sens("onset SD x 1.2",                       function(k, cc) crude_pen(k, cc, sd_scale = 1.2)),
  sens("grid sum, dt = 0.001",                 grid_crude)
)
write.csv(sens_tab, "output/sensitivity.csv", row.names = FALSE)
print(sens_tab, digits = 3)

#-------------------------------------------------------------------------------
# 7. LaTeX macros for every number quoted in main.tex
#-------------------------------------------------------------------------------

pct <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
yrs <- function(x) formatC(x, format = "f", digits = 1)
r <- function(k) tab[tab$cag == k, ]
macros <- c(
  "% Generated by code/crude_penetrance.R; do not edit by hand",
  sprintf("\\newcommand{\\eZeroMEX}{%s}", yrs(who_ex("MEX", 0))),
  sprintf("\\newcommand{\\eZeroESP}{%s}", yrs(who_ex("ESP", 0))),
  sprintf("\\newcommand{\\eZeroGBR}{%s}", yrs(who_ex("GBR", 0))),
  sprintf("\\newcommand{\\eFortyMEX}{%s}", yrs(who_ex("MEX", 40))),
  sprintf("\\newcommand{\\eFortyESP}{%s}", yrs(who_ex("ESP", 40))),
  sprintf("\\newcommand{\\medianFortyCAG}{%s}", formatC(r(40)$median_onset, format = "f", digits = 0)),
  sprintf("\\newcommand{\\netSixtyFiveFortyCAG}{%s}", pct(r(40)$net65)),
  sprintf("\\newcommand{\\netSeventyFiveFortyCAG}{%s}", pct(r(40)$net75)),
  sprintf("\\newcommand{\\crudeFortyMEX}{%s}", pct(r(40)$MEX)),
  sprintf("\\newcommand{\\crudeFortyESP}{%s}", pct(r(40)$ESP)),
  sprintf("\\newcommand{\\crudeFortySWE}{%s}", pct(r(40)$SWE)),
  sprintf("\\newcommand{\\crudeFortyCAN}{%s}", pct(r(40)$CAN)),
  sprintf("\\newcommand{\\crudeFortyGBR}{%s}", pct(r(40)$GBR)),
  sprintf("\\newcommand{\\gapForty}{%s}", pct(r(40)$gap_ESP_MEX)),
  sprintf("\\newcommand{\\gapFortyTwo}{%s}", pct(r(42)$gap_ESP_MEX)),
  sprintf("\\newcommand{\\gapFortyFive}{%s}", pct(r(45)$gap_ESP_MEX)),
  sprintf("\\newcommand{\\crudeFortyBRA}{%s}", pct(r(40)$BRA)),
  sprintf("\\newcommand{\\crudeFortySLV}{%s}", pct(r(40)$SLV)),
  sprintf("\\newcommand{\\crudeFortyHND}{%s}", pct(r(40)$HND)),
  sprintf("\\newcommand{\\crudeFortyGTM}{%s}", pct(r(40)$GTM)),
  sprintf("\\newcommand{\\eFortyGTM}{%s}", yrs(who_ex("GTM", 40))),
  sprintf("\\newcommand{\\gapFortyGTM}{%s}", pct(r(40)$ESP - r(40)$GTM)),
  sprintf("\\newcommand{\\spearmanForty}{%s}",
          formatC(cor(sapply(countries, function(cc) who_ex(cc, 40)),
                      unlist(r(40)[countries]), method = "spearman"), format = "f", digits = 2))
)
writeLines(macros, "output/key-numbers.tex")

#-------------------------------------------------------------------------------
# 8. Figure 1
#-------------------------------------------------------------------------------

## Median age at death among people alive at 18: the age t with S(t | 18) = 0.5.
## Shown as dotted vertical lines (labelled) in Figure 1, 40 CAG only.
median_death18 <- function(cc) {
  S <- make_surv(make_hazard(cc))
  uniroot(function(t) S(t, 18) - 0.5, c(18, 120))$root
}
fig_cc <- c("ESP", "GBR", "BRA", "MEX", "GTM")
med_death <- sapply(fig_cc, median_death18)
write.csv(data.frame(country = fig_cc, median_age_at_death_from18 = med_death),
          "output/median_death_from18.csv", row.names = FALSE)

## Macros for the worked example in Supplementary Methods S1 (40 CAG; Mexico and Spain),
## appended to output/key-numbers.tex
m40 <- langbehn(40)
S_mex <- make_surv(make_hazard("MEX")); S_esp <- make_surv(make_hazard("ESP"))
nmx_at <- function(cc, a) lt$nmx[lt$country == cc & lt$year == 2019 & lt$age_start == a]
dot_at <- function(cc) {
  S <- make_surv(make_hazard(cc))
  integrate(function(u) m40$pdf(u) * S(u, 18), 18, med_death[[cc]])$value / (1 - m40$cdf(18))
}
grid40 <- sens_tab[sens_tab$analysis == "grid sum, dt = 0.001" & sens_tab$cag == 40, ]
main40 <- sens_tab[grepl("^main", sens_tab$analysis) & sens_tab$cag == 40, ]
f4 <- function(x) formatC(x, format = "f", digits = 4)
write(c(
  "% Worked example, Supplementary Methods S1",
  sprintf("\\newcommand{\\muFortyCAG}{%s}", yrs(m40$mu)),
  sprintf("\\newcommand{\\sdFortyCAG}{%s}", yrs(m40$sigma)),
  sprintf("\\newcommand{\\scaleFortyCAG}{%s}", formatC(m40$sigma * sqrt(3) / pi, format = "f", digits = 2)),
  sprintf("\\newcommand{\\FeighteenFortyCAG}{%s}", formatC(m40$cdf(18), format = "f", digits = 4)),
  sprintf("\\newcommand{\\hTwentyMEX}{%s}", f4(nmx_at("MEX", 20))),
  sprintf("\\newcommand{\\hSixtyMEX}{%s}", f4(nmx_at("MEX", 60))),
  sprintf("\\newcommand{\\hEightyMEX}{%s}", f4(nmx_at("MEX", 80))),
  sprintf("\\newcommand{\\hFifteenMEX}{%s}", formatC(nmx_at("MEX", 15), format = "f", digits = 6)),
  sprintf("\\newcommand{\\survSixtyMEX}{%s}", pct(S_mex(60, 18))),
  sprintf("\\newcommand{\\survEightyMEX}{%s}", pct(S_mex(80, 18))),
  sprintf("\\newcommand{\\survSixtyESP}{%s}", pct(S_esp(60, 18))),
  sprintf("\\newcommand{\\survEightyESP}{%s}", pct(S_esp(80, 18))),
  sprintf("\\newcommand{\\medDeathMEX}{%s}", yrs(med_death[["MEX"]])),
  sprintf("\\newcommand{\\medDeathESP}{%s}", yrs(med_death[["ESP"]])),
  sprintf("\\newcommand{\\medDeathGTM}{%s}", yrs(med_death[["GTM"]])),
  sprintf("\\newcommand{\\dotMEX}{%s}", pct(dot_at("MEX"))),
  sprintf("\\newcommand{\\gridDiffMax}{%s}",
          formatC(100 * max(abs(unlist(grid40[c("MEX", "ESP", "GBR")]) -
                                unlist(main40[c("MEX", "ESP", "GBR")]))), format = "f", digits = 3))
), "output/key-numbers.tex", append = TRUE)

draw_fig <- function() {
  par(mar = c(4.2, 4.2, 1, 0.5))
  cols <- c(ESP = "#2c3e50", GBR = "#7f8c8d", BRA = "#27ae60", MEX = "#c0392b", GTM = "#e67e22")
  ages <- seq(30, 95, by = 0.5)
  k <- 40
  m <- langbehn(k)
  plot(ages, net_pen(k, ages), type = "n", ylim = c(0, 1), las = 1, yaxt = "n",
       xlab = "Age (years)", ylab = "Cumulative probability of onset")
  axis(2, at = seq(0, 1, 0.2), labels = paste0(seq(0, 100, 20), "%"), las = 1)
  abline(h = seq(0, 1, 0.2), v = seq(30, 90, 10), col = "grey90")
  ## Dotted lines at the median age at death from 18, labelled vertically in small italics
  ## Each dotted line stops where it meets its own country's crude curve
  y_stop <- sapply(names(cols), function(cc) {
    S <- make_surv(make_hazard(cc))
    integrate(function(u) m$pdf(u) * S(u, 18), 18, med_death[cc])$value / (1 - m$cdf(18))
  })
  segments(med_death[names(cols)], -0.05, med_death[names(cols)], y_stop,
           col = cols, lty = 3, lwd = 1.5)
  points(med_death[names(cols)], y_stop, pch = 19, cex = 0.6, col = cols)
  ## Crude penetrance reached by the median age at death, printed at each dot
  ## Spain above the dot; Brazil above-right; Mexico above-left; Guatemala below-left;
  ## UK below and to the right of the Spain line, clear of lines and each other
  dx <- c(ESP = 0, GBR = 1.4, BRA = 0.3, MEX = -0.3, GTM = -0.4)
  dy <- c(ESP = 0.02, GBR = -0.012, BRA = 0.015, MEX = 0.015, GTM = -0.035)
  adj_lab <- list(ESP = c(0.5, 0), GBR = c(0, 1), BRA = c(0, 0), MEX = c(1, 0), GTM = c(1, 1))
  for (cc in names(cols)) {
    lab <- sprintf("%.1f%%", 100 * y_stop[cc])
    aj <- adj_lab[[cc]]
    x0 <- med_death[cc] + dx[cc]; y0 <- y_stop[cc] + dy[cc]
    w <- strwidth(lab, cex = 0.7, font = 3); h <- strheight(lab, cex = 0.7, font = 3)
    rect(x0 - aj[1] * w - 0.15, y0 - aj[2] * h - 0.008, x0 + (1 - aj[1]) * w + 0.15,
         y0 + (1 - aj[2]) * h + 0.008, col = "white", border = NA)
    text(x0, y0, lab, adj = aj, font = 3, cex = 0.7, col = cols[cc])
  }
  write.csv(data.frame(country = names(cols), median_age_at_death_from18 = med_death[names(cols)],
                       crude_penetrance_by_that_age = y_stop),
            "output/median_death_from18.csv", row.names = FALSE)
  lab_nm <- c(ESP = "Spain", GBR = "UK", BRA = "Brazil", MEX = "Mexico", GTM = "Guatemala")
  side <- c(ESP = 1, GBR = -1, BRA = 1, MEX = -1, GTM = -1)   # Spain right of its line to avoid the UK label
  for (cc in names(cols))
    text(med_death[cc] + 0.3 * side[cc], 0.04, sprintf("%s, %.1f years", lab_nm[cc], med_death[cc]),
         srt = 90, font = 3, cex = 0.7, col = cols[cc], adj = c(0, if (side[cc] > 0) 1 else 0))
  lines(ages, net_pen(k, ages), lty = 2, lwd = 1.3)
  for (cc in names(cols)) {
    S <- make_surv(make_hazard(cc))
    ci <- sapply(ages, function(a)
      integrate(function(u) m$pdf(u) * S(u, 18), 18, a)$value) / (1 - m$cdf(18))
    lines(ages, ci, col = cols[cc], lwd = 2)
  }
  legend("topleft", bty = "n", cex = 0.85,
         legend = c("Net (no other deaths)", cnames[names(cols)]),
         lty = c(2, rep(1, length(cols))), lwd = c(1.3, rep(2, length(cols))), col = c("black", cols))
}
draw_fig()
pdf("../figures/figure1-crude-penetrance.pdf", width = 9, height = 4.2); draw_fig(); invisible(dev.off())
png("../figures/figure1-crude-penetrance.png", width = 9, height = 4.2, units = "in", res = 200); draw_fig(); invisible(dev.off())

#-------------------------------------------------------------------------------
# 9. Table 1 as LaTeX (../tables/table1-crude-penetrance.tex)
#-------------------------------------------------------------------------------

f3 <- function(x) formatC(x, format = "f", digits = 3)
## Rows: countries (grouped); columns: life expectancy at 0 and 40, crude penetrance for 40 to 45 CAG
row_for <- function(cc) paste(c(cnames[[cc]], yrs(who_ex(cc, 0)), yrs(who_ex(cc, 40)),
                                f3(tab[[cc]])), collapse = " & ")
body <- unlist(lapply(names(groups), function(g)
  c(sprintf("\\multicolumn{9}{l}{\\textit{%s}} \\\\", g),
    paste0("\\quad ", sapply(groups[[g]], row_for), " \\\\"))))
tex <- c(
  "% Generated by code/crude_penetrance.R; do not edit by hand",
  "\\begin{table}[h!]",
  "\\centering",
  "\\caption{\\textbf{Net and crude probability of Huntington's disease onset by CAG repeat length and country.} Net probabilities assume no death from other causes; crude probabilities give the probability of onset before death from other causes under 2019 national mortality. All probabilities are conditional on survival without onset to age 18.}",
  "\\label{tab:crude}",
  "\\small",
  "\\begin{tabular}{lcccccccc}",
  "\\toprule",
  " & \\multicolumn{2}{c}{Life expectancy} & \\multicolumn{6}{c}{CAG repeat length} \\\\",
  "\\cmidrule(lr){2-3} \\cmidrule(lr){4-9}",
  paste0(" & At birth & At 40 & ", paste(cags, collapse = " & "), " \\\\"),
  "\\midrule",
  paste0("Median age at onset & & & ", paste(formatC(tab$median_onset, format = "f", digits = 1), collapse = " & "), " \\\\"),
  paste0("Net penetrance by age 75 & & & ", paste(f3(tab$net75), collapse = " & "), " \\\\"),
  "\\midrule",
  "\\multicolumn{9}{l}{Crude lifetime penetrance} \\\\",
  body,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{0.5em}",
  "\\footnotesize",
  "Note: The onset distribution follows the logistic model of \\citet{langbehn2004new}, with median onset (years) equal to $21.54 + \\exp(9.556 - 0.146 \\times \\text{CAG})$. Background mortality uses age-specific death rates for both sexes from the World Health Organization Global Health Observatory life tables for 2019 \\citep{who2019lifetables}, the last year before the COVID-19 pandemic, treated as constant within each age group. Life expectancies (years) are the values published by the World Health Organization for 2019.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tex, "../tables/table1-crude-penetrance.tex")

#-------------------------------------------------------------------------------
# 10. Supplementary Table S1 as LaTeX (../tables/tableS1-sensitivity.tex)
#-------------------------------------------------------------------------------

labs <- c("main (2019, start 18, 85+ = WHO nmx)" = "Main analysis",
          "85+ hazard = 1/e85" = "Death rate at 85+ set to $1/e_{85}$",
          "2023 life tables" = "2023 life tables",
          "start age 30" = "Conditional on age 30 instead of 18",
          "onset SD x 0.8" = "Onset standard deviation $\\times$ 0.8",
          "onset SD x 1.2" = "Onset standard deviation $\\times$ 1.2",
          "grid sum, dt = 0.001" = "Grid sum, step 0.001 years")
rows <- sapply(seq_len(nrow(sens_tab)), function(i) with(sens_tab[i, ],
  paste(labs[analysis], cag, f3(MEX), f3(ESP), f3(GBR),
        formatC(100 * gap_ESP_MEX, format = "f", digits = 1), sep = " & ")))
texS <- c(
  "% Generated by code/crude_penetrance.R; do not edit by hand",
  "\\begin{table}[h!]",
  "\\centering",
  "\\caption{\\textbf{Sensitivity of crude lifetime penetrance to analytic choices.}}",
  "\\label{tab:sensitivity}",
  "\\small",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Analysis & CAG & Mexico & Spain & UK & Gap (points) \\\\",
  "\\midrule",
  paste0(rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{0.5em}",
  "\\footnotesize",
  "Note: Gap is crude penetrance in Spain minus crude penetrance in Mexico, in percentage points. The standard deviation rows rescale the onset standard deviation of \\citet{langbehn2004new} at every repeat length. UK, United Kingdom.",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(texS, "../tables/tableS1-sensitivity.tex")

sessionInfo()

#-------------------------------------------------------------------------------
# 11. Supplementary Table S2: survival from age 18 in Mexico, term by term
#     (worked example for Step 4 of Supplementary Methods S1)
#-------------------------------------------------------------------------------

mx_mex <- lt[lt$country == "MEX" & lt$year == 2019 & lt$age_start >= 15 & lt$age_start < 80, ]
yrs_in <- pmin(mx_mex$age_start + 5, 80) - pmax(mx_mex$age_start, 18)   # years of each group between 18 and 80
contrib <- mx_mex$nmx * yrs_in
cumH <- cumsum(contrib)
end_age <- mx_mex$age_start + 5
S2 <- data.frame(group = sprintf("%d--%d", mx_mex$age_start, mx_mex$age_start + 4),
                 years = yrs_in, nmx = mx_mex$nmx, contrib = contrib, cumH = cumH,
                 S = exp(-cumH), end_age = end_age)
stopifnot(abs(exp(-cumH[end_age == 60]) - S_mex(60, 18)) < 1e-12,
          abs(exp(-cumH[end_age == 80]) - S_mex(80, 18)) < 1e-12)
f6 <- function(x) formatC(x, format = "f", digits = 6)
rowsS2 <- with(S2, sprintf("%s & %d & %s & %s & %s & %s & %s \\\\", group, years, f6(nmx),
                           f6(contrib), f6(cumH), end_age, formatC(100 * S, format = "f", digits = 1)))
texS2 <- c(
  "% Generated by code/crude_penetrance.R; do not edit by hand",
  "\\begin{table}[h!]",
  "\\centering",
  "\\caption{\\textbf{Survival from age 18 in Mexico, 2019, computed term by term.} Each row adds the death rate of one World Health Organization age group times the number of years of that group between age 18 and the end of the group. The cumulative hazard $H$ is the running sum of the contributions, and survival from age 18 to the end of the group is $\\exp(-H)$.}",
  "\\label{tab:S2}",
  "\\footnotesize",
  "\\begin{tabular}{lrcccrc}",
  "\\toprule",
  "Age & Years & Death rate & & Cumulative & & Survival \\\\",
  "group & after 18 & ${}_n m_x$ & Contribution & hazard $H$ & To age & (\\%) \\\\",
  "\\midrule",
  rowsS2,
  "\\bottomrule",
  "\\end{tabular}",
  "\\end{table}"
)
writeLines(texS2, "../tables/tableS2-survival-mex.tex")
write(c(
  sprintf("\\newcommand{\\cumHazSixtyMEX}{%s}", formatC(cumH[end_age == 60], format = "f", digits = 4)),
  sprintf("\\newcommand{\\cumHazEightyMEX}{%s}", formatC(cumH[end_age == 80], format = "f", digits = 4)),
  sprintf("\\newcommand{\\cumHazSixtyEightyMEX}{%s}", formatC(cumH[end_age == 80] - cumH[end_age == 60], format = "f", digits = 4))
), "output/key-numbers.tex", append = TRUE)
## Rates written out for the explicit sums in Step 4
g2055 <- mx_mex$age_start >= 20 & mx_mex$age_start <= 55
g6075 <- mx_mex$age_start >= 60 & mx_mex$age_start <= 75
write(c(
  sprintf("\\newcommand{\\ratesTwentyThirtyFiveMEX}{%s}", paste(f6(mx_mex$nmx[g2055][1:4]), collapse = " + ")),
  sprintf("\\newcommand{\\ratesFortyFiftyFiveMEX}{%s}", paste(f6(mx_mex$nmx[g2055][5:8]), collapse = " + ")),
  sprintf("\\newcommand{\\sumRatesTwentyFiftyFiveMEX}{%s}", f6(sum(mx_mex$nmx[g2055]))),
  sprintf("\\newcommand{\\ratesSixtySeventyFiveMEX}{%s}", paste(f6(mx_mex$nmx[g6075]), collapse = " + ")),
  sprintf("\\newcommand{\\sumRatesSixtySeventyFiveMEX}{%s}", f6(sum(mx_mex$nmx[g6075]))),
  sprintf("\\newcommand{\\hFifteenMEXsix}{%s}", f6(mx_mex$nmx[mx_mex$age_start == 15]))
), "output/key-numbers.tex", append = TRUE)
