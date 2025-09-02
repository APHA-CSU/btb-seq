process IDNONBOVIS {
	errorStrategy 'finish'
    tag "$pair_id"
	publishDir "$params.outdir/Results_${params.DataDir}_${params.today}/NonBovID", mode: 'copy', pattern: '*tab'
	maxForks 1
	
	input:
		tuple val(pair_id), path(outcome), path(trim1), path(trim2)
		val params.kraken2db
		val lowmem
        val rm_inter
	
	output:
		path("${pair_id}_bovis.csv"), emit: queryBovis
		tuple val(pair_id), path("${pair_id}_*_brackensort.tab"), path("${pair_id}_*_kraken2.tab"), optional: true, emit: krakenOut
	
	script:
	"""
    # Error handling
    set -eo pipefail

    # Read the outcome
    outcome=\$(cat ${outcome})
    if [[ "\$outcome" != "Pass" ]]; then
        # Run Kraken2
        kraken2 --threads 2 --quick ${lowmem} --db ${params.kraken2db} --output - --report interim_kraken2.tab --paired ${trim1} ${trim2}

        # Process Kraken2 output
        sed 's/\\t\\t1\\troot/\\tR\\t1\\troot/g' interim_kraken2.tab |
        sed 's/\\t1\\t131567\\t/\\tR1\\t131567\\t/g' |
        sed 's/\\t\\t2\\t/\\tD\\t2\\t/g' > ${pair_id}_"\$outcome"_kraken2.tab

        # HACK: Ignore Bracken errors when read counts are low
        set +e

        # Run Bracken for level S
        bracken -d ${params.kraken2db} -r 150 -l S -t 10 -i ${pair_id}_"\$outcome"_kraken2.tab -o ${pair_id}_"\$outcome"_bracken.out
        sed 1d ${pair_id}_"\$outcome"_bracken.out | sort -t \$'\\t' -k7,7 -nr - | head -20 > ${pair_id}_"\$outcome"_brackensort.tab

        # Run Bracken for level S1
        bracken -d ${params.kraken2db} -r 150 -l S1 -i ${pair_id}_"\$outcome"_kraken2.tab -o ${pair_id}_"\$outcome"-S1_bracken.out
        ( sed -u 1q; sort -t \$'\\t' -k7,7 -nr ) < ${pair_id}_"\$outcome"-S1_bracken.out > ${pair_id}_"\$outcome"-S1_brackensort.tab

        set -e

        # Extract Mycobacterium bovis information
        BovPos=\$(grep -m 1 'Mycobacterium' ${pair_id}_"\$outcome"-S1_brackensort.tab | grep 'variant bovis' |
        awk '{printf \$1" "\$2" "\$3" "\$4","\$9"," "%.3f", (\$10*100)}' || true)

        # Write bovis information to CSV
        echo "Sample,ID,TotalReads,Abundance" > "${pair_id}_bovis.csv"
        echo "${pair_id},\$BovPos" >> "${pair_id}_bovis.csv"
    else
        # If outcome is "Pass", no further processing is required
        echo "ID not required"
        echo "Sample,ID,TotalReads,Abundance" > "${pair_id}_bovis.csv"
    fi

    if [ ${rm_inter} = 'true' ]; then
        rm -f -- "\$(readlink -f -- \"${trim1}\")" "\$(readlink -f -- \"${trim2}\")"
    fi

	"""
}