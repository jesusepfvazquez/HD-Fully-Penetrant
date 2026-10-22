################################################################################
## Global supplementary analysis: crude lifetime penetrance of HD alleles under
## the 2019 background mortality of every country with at least 90,000 people.
##
## Uses the same onset model and the same formula as crude_penetrance.R
## (Equation 3 of main.tex, conditioned on age 18), with single-year death rates
## from WPP 2024 (data/wpp2024_mx_2019.csv, made by extract_wpp2024.R) in place
## of the WHO abridged life tables.
##
## Run from the code/ folder:  Rscript global_penetrance.R
## Outputs:
##   output/global_penetrance.csv           Supplementary Table S2
##   output/global_vs_who.csv               the 12 main-text countries, WPP vs WHO
##   output/key-numbers-global.tex          LaTeX macros quoted in the text
##   ../figures/figureS1-global.{pdf,png}   Supplementary Figure S1
################################################################################

dir.create("output", showWarnings = FALSE)
d <- read.csv("data/wpp2024_mx_2019.csv", stringsAsFactors = FALSE)

## Onset model of Langbehn et al. (2004); identical to crude_penetrance.R
langbehn <- function(cag) {
  mu <- 21.54 + exp(9.556 - 0.1460 * cag)
  s <- sqrt(35.55 + exp(17.72 - 0.3269 * cag)) * sqrt(3) / pi
  list(cdf = function(t) plogis(t, mu, s), pdf = function(t) dlogis(t, mu, s))
}

## Survival from age a0 with a constant hazard within each single year of age;
## the WPP rate for ages 100+ applies to every age from 100 onward.
make_surv <- function(mx, a0 = 18) {
  breaks <- 0:100
  function(t) sapply(t, function(tt) {
    if (tt <= a0) return(1)
    knots <- sort(unique(c(a0, breaks[breaks > a0 & breaks < tt], tt)))
    exp(-sum(mx[findInterval(knots[-length(knots)], breaks)] * diff(knots)))
  })
}

crude_pen <- function(cag, S, a0 = 18) {
  m <- langbehn(cag)
  integrate(function(u) m$pdf(u) * S(u), a0, 130, subdivisions = 2000)$value / (1 - m$cdf(a0))
}

ccs <- unique(d$country_code)
res <- do.call(rbind, lapply(ccs, function(k) {
  dk <- d[d$country_code == k, ]
  S18 <- make_surv(dk$mx, 18); S40 <- make_surv(dk$mx, 40)
  data.frame(country_code = k, name = dk$name[1], region = dk$region[1], area = dk$area[1],
             e0 = dk$e0_wpp[1],
             e40 = integrate(S40, 40, 130, subdivisions = 2000)$value,
             median_death_from18 = uniroot(function(t) S18(t) - 0.5, c(19, 129))$root,
             crude40 = crude_pen(40, S18), crude42 = crude_pen(42, S18),
             crude45 = crude_pen(45, S18))
}))
res <- res[order(res$crude40), ]
write.csv(res, "output/global_penetrance.csv", row.names = FALSE)

## Check against the WHO-based main analysis for the 12 main-text countries
main12 <- c(ESP = 724, SWE = 752, GBR = 826, CAN = 124, USA = 840, CHL = 152, PER = 604, COL = 170,
            ARG = 32, BRA = 76, MEX = 484, SLV = 222, HND = 340, GTM = 320)
who <- read.csv("output/table1.csv")
cmp <- data.frame(iso3 = names(main12),
                  wpp = res$crude40[match(main12, res$country_code)],
                  who = unlist(who[who$cag == 40, names(main12)]))
cmp$diff_pts <- 100 * (cmp$wpp - cmp$who)
write.csv(cmp, "output/global_vs_who.csv", row.names = FALSE)
print(cmp)

## Summaries
lac <- res$area == "Latin America and the Caribbean"
rk <- function(code) which(res$country_code[order(res$crude40)] == code)
pct <- function(x) formatC(100 * x, format = "f", digits = 1)
q <- quantile(res$crude40, c(0, .5, 1))
asia <- c(CHN = 156, JPN = 392, IND = 356, KOR = 410)
macros <- c(
  sprintf("\\newcommand{\\nGlobal}{%d}", nrow(res)),
  sprintf("\\newcommand{\\globalMin}{%s}", pct(q[1])),
  sprintf("\\newcommand{\\globalMinName}{%s}", res$name[1]),
  sprintf("\\newcommand{\\globalMax}{%s}", pct(q[3])),
  sprintf("\\newcommand{\\globalMaxName}{%s}", res$name[nrow(res)]),
  sprintf("\\newcommand{\\globalMedian}{%s}", pct(q[2])),
  sprintf("\\newcommand{\\globalSpearman}{%s}",
          formatC(cor(res$e40, res$crude40, method = "spearman"), format = "f", digits = 2)),
  sprintf("\\newcommand{\\globalRankMEX}{%d}", rk(484)),
  sprintf("\\newcommand{\\globalWhoMaxDiff}{%s}", formatC(max(abs(cmp$diff_pts)), format = "f", digits = 1)),
  sprintf("\\newcommand{\\globalLACmin}{%s}", pct(min(res$crude40[lac]))),
  sprintf("\\newcommand{\\globalLACmax}{%s}", pct(max(res$crude40[lac]))),
  sprintf("\\newcommand{\\globalWPP%s}{%s}", names(asia), pct(res$crude40[match(asia, res$country_code)])),
  sprintf("\\newcommand{\\globalWPP%s}{%s}", c("MEX", "ESP", "GTM"),
          pct(res$crude40[match(c(484, 724, 320), res$country_code)]))
)
writeLines(macros, "output/key-numbers-global.tex")
cat(macros, sep = "\n")

## Supplementary Figure S1: crude lifetime penetrance at 40 CAG against
## remaining life expectancy at 40, every country, colored by continent
areas <- c("Africa", "Asia", "Europe", "Latin America and the Caribbean",
           "Northern America", "Oceania")
cols <- c("#e67e22", "#8e44ad", "#2c3e50", "#c0392b", "#16a085", "#7f8c8d")
names(cols) <- areas
lab <- c(MEX = 484, GTM = 320, ESP = 724, CHN = 156, JPN = 392, IND = 356, CAF = 140,
         USA = 840, BRA = 76, NGA = 566)
labnames <- c(MEX = "Mexico", GTM = "Guatemala", ESP = "Spain", CHN = "China", JPN = "Japan",
              IND = "India", CAF = "Central African Republic", USA = "United States", BRA = "Brazil", NGA = "Nigeria")
draw <- function() {
  par(mar = c(4.2, 4.6, 1, 1), las = 1)
  plot(res$e40, 100 * res$crude40, pch = 21, bg = adjustcolor(cols[res$area], 0.75),
       col = "white", cex = 1.3, xlab = "Remaining life expectancy at age 40 (years)",
       ylab = "", yaxt = "n", ylim = c(20, 95))
  axis(2, at = seq(20, 90, 10), labels = paste0(seq(20, 90, 10), "%"))
  mtext("Crude lifetime penetrance, 40 CAG repeats", side = 2, line = 3.4, las = 0)
  i <- match(lab, res$country_code)
  points(res$e40[i], 100 * res$crude40[i], pch = 21, bg = cols[res$area[i]], col = "black", cex = 1.3)
  text(res$e40[i], 100 * res$crude40[i], labnames[names(lab)], pos = c(4, 4, 1, 3, 3, 1, 4, 4, 2, 4),
       cex = 0.75, font = 3)
  legend("bottomright", legend = areas, pt.bg = cols, pch = 21, col = "white", pt.cex = 1.3,
         bty = "n", cex = 0.8)
}
dir.create("../figures", showWarnings = FALSE)
pdf("../figures/figureS1-global.pdf", width = 7.5, height = 5); draw(); dev.off()
png("../figures/figureS1-global.png", width = 7.5, height = 5, units = "in", res = 200); draw(); dev.off()
