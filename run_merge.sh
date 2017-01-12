#!/bin/bash
#
#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem 18000
#SBATCH -t 0-2:00
# --mail-type=END,FAIL
# --mail-user=EMAILHERE@uab.edu

rm -R output; ./merge_kraken_reports.py input output kraken_ffn_fulldb.report
