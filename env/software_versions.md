# Software versions

Written automatically by `scripts/dev/run_all.sh`. These are the versions that
produced the published results; do not edit by hand.

| Item | Value |
|---|---|
| Run started | 2026-11-22T18:09:10Z |
| Commit | `6be6fb29c4896ca6aee7f1212133c74bf12176e1` |
| Data release | `v0.3-data` |
| Operating system | Ubuntu 20.04.4 LTS |
| PLINK 2 | PLINK v2.0.0-a.7LM AVX2 Intel (7 Jul 2025) |
| PLINK 1.9 | PLINK v1.90b6.25 64-bit (5 Mar 2022) |
| REGENIE | v4.1.2.gz |
| R | R version 4.5.1 (2025-06-13) |

## Executables

| Tool | Executable | MD5 | Note |
|---|---|---|---|
| plink2 | `/usr/bin/plink2` | `47a67647af6c9160c9443f8cbd2df017` |  |
| plink | `/usr/local/bin/plink19` | `d5f949f706d8ab017cf9a05426e18b3b` | PLINK 1.9 |
| regenie | `tools/micromamba-root/envs/regenie_env/bin/regenie` | `cca7fc19f2b6aa83106dd96c983df34b` | binary run by the tools/bin/regenie wrapper |
| metal | not found | | not used by the canonical run |
| Rscript | `/usr/bin/Rscript` | `2e423d4fe873986b71010b707d38b2d0` |  |

## R packages

| Package | Version |
|---|---|
| data.table | 1.17.8 |
| optparse | 1.8.2 |
| ggplot2 | 4.0.2 |
| survival | 3.5.7 |
