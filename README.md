# 📦 TOTKA Shiny Application

**TOTKA** is a comprehensive **R-Shiny application** for **integrated AOP‑anchored framework that combines transcriptomic dose-response modeling with physiologically based kinetic (PBK) simulation**.  

It integrates following core R packages into a single graphical interface, enabling streamlined **physiologically based kinetic (PBK) simulation** and **adverse outcome pathway (AOP) mapping**:

- [PfasDosim package](https://github.com/ntua-unit-of-control-and-informatics/PfasDosim) 
- [AOPfingerprintR package](https://github.com/fhaive/AOPfingerprintR)  

AOPfingerprint APIs are also hosted [here](https://test.enaloscloud.novamechanics.com/insight/aop_fingerprint/__docs__/ )  

---

## Docker

- TOTKA is also available as a docker at [docker hub](https://hub.docker.com/r/fhaive/totka)

How to use: 

```
docker pull fhaive/totka:latest
docker run -p 80:80 totka:latest /bin/bash -c "cd /TOTKA/ && R -e \" options('shiny.port'=80,shiny.host='0.0.0.0');library(TOTKA);setwd('/TOTKA/');TOTKA::run_app() \""
```
Open your browser and visit: 127.0.0.1:80

- AOPFingerprint APIs are also available as a docker at [docker hub](https://hub.docker.com/r/fhaive/aopfingerprint_api)

---

## 🚀 Quick Start (Demo)

You can try out **TOTKA** immediately with the demo dataset included in the repository under the sample_data/liver_spheroids/ folder.

This is omics data generated from PFOS and PFOA exposure in in-vitro liver spheroids across multiple dose levels and four time points (1, 4, 10 and 14 day) (Rowan-Carroll et al., 2021).
[Click here for the original article describing the dataset](https://doi.org/10.1093/toxsci/kfab039)  

## Launch the app locally

```r
remotes::install_github("https://github.com/fhaive/totka")
library(TOTKA)
TOTKA::run_app()
```

---

## 📁 Folder Structure

```plaintext
.
├── data/                         # Data used by the app
├── sample_data/                  # Sample data to run the app
├── PIPELINE_SCRIPT/              # Scripts to execute the analysis pipeline
├── manual/                       # Documentation or user guides
├── server.R                      # Main Shiny server function
└── ui.R                          # Main Shiny UI definition
```

---


## 📑 Description

- This app provides a user-friendly interface for navigating integrated OMICS and PBK analysis outputs.
- Assessment in the context of Adverse Outcome Pathway (AOP).

---

## 📦 Requirements

- R ≥ 4.0
- R packages are managed via `renv` (see `renv.lock`)

---

## 📬 Contact

For questions, suggestions, or contributions, please contact [Angela Serra] at [angela.serra@tuni.fi].
