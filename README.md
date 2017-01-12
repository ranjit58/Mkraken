## Merge Kraken
##### This script will merge reports from the kraken pipeline into a single report and can filter and sort results as .OTU.
* Usage: ./merge_kraken_reports.py input output file_format [-s -f -a -v -h]
* Enter ./merge_kraken_reports.py -h for usage help

Use ./run_merge.sh to run the merge_kraken_reports.py script on the SLURM scheduler

Note: The file_format file is missing as it is larger than 100MB. Need to implement git-lfs on the cluster.
