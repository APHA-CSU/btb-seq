nextflow run bTB-WGS_process.nf \
--outdir "/results/" \
--reads "/reads/*_{S*_R1,S*_R2}*.fastq.gz" \
--lowmem '"--memory-map"' \
-with-report /artifacts/report.html