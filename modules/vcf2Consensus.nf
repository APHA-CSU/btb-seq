process VCF2CONSENSUS {
	errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/consensus/", mode: 'copy', pattern: '*.fas'
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/snpTables/", mode: 'copy', pattern: '*.tab'
	maxForks 2
	
	input:
		path ref
		tuple val(pair_id), path (mask_bed), path (nonmasked_bed), path (vcf), path (csi)
		val min_freq_alt
        val outdir
        val today
        path DataDir
	
	output:
		tuple val (pair_id), path("${pair_id}_consensus.fas"), path("${pair_id}_snps.tab"), emit: consensus
		path("${pair_id}_ncount.csv"), emit: nCount
	
	script:
    """
    if [ ! -s ${nonmasked_bed} ]; then
        echo "LT708304-Mycobacteriumbovis-AF2122-97	-1	-1" > ${nonmasked_bed}
    fi

    bcftools filter -i "ALT!='.' && INFO/AD[1]/(INFO/AD[0]+INFO/AD[1]) >= ${min_freq_alt}" ${vcf} -Ob -o ${pair_id}_filtered.bcf
    bcftools index ${pair_id}_filtered.bcf

    bcftools consensus -f ${ref} -e 'TYPE="indel"' -m ${mask_bed} ${pair_id}_filtered.bcf |
    sed "/^>/ s/.*/>${pair_id}/" > ${pair_id}_consensus.fas 

    echo -e "Sample,Ncount,ResultLoc" > ${pair_id}_ncount.csv
    N_count=\$(grep -o "N" ${pair_id}_consensus.fas  | wc -l)
    echo -e "${pair_id},\${N_count},$outdir/Results_${DataDir}_${today}/" >> ${pair_id}_ncount.csv

    echo -e 'CHROM\tPOS\tTYPE\tREF\tALT\tEVIDENCE' > ${pair_id}_snps.tab
    bcftools query -R ${nonmasked_bed} -e 'TYPE="REF"' -f '%CHROM,%POS,%TYPE,%REF,%ALT,%DP4\n' ${pair_id}_filtered.bcf |  
    awk -F, '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5"\t"\$5":"\$8+\$9" "\$4":"\$6+\$7}' >> ${pair_id}_snps.tab
    """
}