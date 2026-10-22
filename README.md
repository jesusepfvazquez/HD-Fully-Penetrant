# HD-Fully-Penetrant

Code and tables for the paper *Fully Penetrant Within Whose Lifetime? Background Mortality and the Penetrance of Huntington's Disease Alleles in Mexico and Central America*.

The paper combines the age-at-onset model of Langbehn et al. (2004) with national life tables to compare **net penetrance** (the probability of onset if death from other causes were impossible) with **crude penetrance** (the probability of onset before death from other causes) for carriers of 40 to 45 CAG repeats.

## Contents

| Path | What it is |
|---|---|
| `analysis.Rmd` / `analysis.html` | Step-by-step walkthrough of every calculation, with the worked example for 40 repeats in Mexico, checked against the script outputs |
| `code/crude_penetrance.R` | Main analysis: Table 1, Figure 1, Supplementary Tables S1 and S2, validation, and the LaTeX macros for every number quoted in the paper |
| `code/global_penetrance.R` | Global analysis for 201 countries (Supplementary Figure S1) |
| `code/extract_wpp2024.R` | Extracts 2019 single-year death rates from the UN World Population Prospects 2024 (`wpp2024` R package data) |
| `code/data/who_life_tables.csv` | WHO Global Health Observatory life tables, both sexes, 2019 (14 countries) and 2023 (Mexico, Spain, United Kingdom) |
| `code/data/wpp2024_mx_2019.csv` | UN WPP 2024 single-year death rates, 2019, 201 countries |
| `code/output/` | CSV outputs, validation checks, and LaTeX macros (`key-numbers.tex`, `key-numbers-global.tex`) |
| `tables/` | LaTeX tables as they appear in the paper |
| `figures/` | Figure 1 and Supplementary Figure S1 (PDF and PNG) |

## Reproducing the results

All scripts use base R (version 4.3.3); the R Markdown file also needs `rmarkdown` and `knitr`.

```sh
cd code
Rscript crude_penetrance.R     # Table 1, Figure 1, Tables S1-S2, macros
Rscript global_penetrance.R    # Supplementary Figure S1 (needs output/table1.csv from the line above)
cd ..
Rscript -e 'rmarkdown::render("analysis.Rmd")'
```

`extract_wpp2024.R` needs the raw `wpp2024` data files (`mx1dt.rda`, `popAge1dt.rda`, `e01dt.rda`, `UNlocations.txt`) from <https://github.com/PPgp/wpp2024/tree/main/data> in `code/raw/`; its output is already included in `code/data/`.

## Data sources

- World Health Organization. Global Health Observatory life tables, indicators `LIFE_0000000029` (age-specific death rates) and `LIFE_0000000035` (life expectancy), <https://www.who.int/data/gho>.
- United Nations, Department of Economic and Social Affairs, Population Division. World Population Prospects 2024, <https://population.un.org/wpp/>.
- Langbehn DR, Brinkman RR, Falush D, Paulsen JS, Hayden MR. A new model for prediction of the age of onset and penetrance for Huntington's disease based on CAG length. *Clin Genet* 2004;65:267-277.
