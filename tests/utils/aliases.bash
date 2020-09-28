alias nextflowtest="nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report /artifacts/report.html"

alias assert_first_csv_row="python tests/utils/assert_first_csv_row.py"
