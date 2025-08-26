process READSTATS {
    errorStrategy 'finish'
    tag "$pair_id"
    maxForks 2

    input:
        tuple val(pair_id), path(raw1), path(raw2), path(uniq1), path(uniq2), path(trim1), path(trim2), path(bam)
        val (rm_inter)

    output:
        tuple val(pair_id), path("${pair_id}_stats.csv"), emit: stats
        tuple path("${pair_id}_stats.csv"), emit: stats_table
        tuple val(pair_id), path("${pair_id}_outcome.txt"), emit: outcome

    script:

        """
        # Error handling
        set -eo pipefail

        # Assign the pair_id variable
        pair_id="${pair_id}"

        # Count reads in each category; in fastq files each read consists of four lines
        raw_R1=\$(( \$(zcat ${raw1} | wc -l) / 4 ))
        uniq_R1=\$(( \$(cat ${uniq1} | wc -l) / 4 ))
        trim_R1=\$(( \$(cat ${trim1} | wc -l) / 4 ))
        num_map=\$(samtools view -c ${bam})
        
        #Add ifelse to prevent empty bam errors stopping run
        samtools depth -a ${bam} > depth.txt
        sites=\$(wc -l < depth.txt || echo 0)
        if [ "\$sites" -gt 0 ]; then
            avg_depth=\$(awk '{sum+=\$3} END { printf "%.3f", (NR>0)?sum/NR:0 }' depth.txt)
            zero_cov=\$(awk '\$3<1 {c++} END {print (c+0)}' depth.txt)
        else
            avg_depth=0
            zero_cov=0
        fi
        
        # Cleanup deduplicated reads 
        rm depth.txt
        if [ ${rm_inter} = 'true' ]; then
            rm -f -- "\$(readlink -f -- \"${uniq1}\")" "\$(readlink -f -- \"${uniq2}\")"
        fi

        # Calculate values and percentages
        num_raw=\$((\$raw_R1*2))
        num_uniq=\$((\$uniq_R1*2))
        num_trim=\$((\$trim_R1*2))
        pc_aft_dedup=\$(echo "scale=2; (\$num_uniq*100/\$num_raw)" | bc)
        pc_aft_trim=\$(echo "scale=2; (\$num_trim*100/\$num_uniq)" | bc)
        pc_mapped=\$(echo "scale=2; (\$num_map*100/\$num_trim)" | bc)
        genome_cov=\$(echo "scale=2; (100-(\$zero_cov*100/\$sites))" | bc)

        # Define thresholds for flag assignment
        mindepth=10
        minpc=60
        minreads=600000
        minafttrim=60

        # Assign flags based on thresholds
        if [[ \${pc_aft_trim%%.*} -lt \$minafttrim ]]; then
            flag="LowQualData"
        elif [[ \${avg_depth%%.*} -ge \$mindepth ]] && [[ \${pc_mapped%%.*} -gt \$minpc ]]; then
            flag="Pass"
        elif [[ \${avg_depth%%.*} -lt \$mindepth ]] && [[ \${pc_mapped%%.*} -lt \$minpc ]] && [[ \$num_trim -gt \$minreads ]]; then
            flag="Contaminated"
        elif [[ \${avg_depth%%.*} -lt \$mindepth ]] && [[ \$num_trim -lt \$minreads ]]; then
            flag="InsufficientData"
        else
            flag="CheckRequired"
        fi

        # Write values to CSV file
        echo "Sample,NumRawReads,NumDedupReads,%afterDedup,NumTrimReads,%afterTrim,NumMappedReads,%Mapped,MeanDepth,GenomeCov,Outcome" > "${pair_id}_stats.csv"
        echo "${pair_id},\${num_raw},\${num_uniq},\${pc_aft_dedup},\${num_trim},\${pc_aft_trim},\${num_map},\${pc_mapped},\${avg_depth},\${genome_cov},\${flag}" >> "${pair_id}_stats.csv"
        echo "\${flag}" > "${pair_id}_outcome.txt"
        """
}