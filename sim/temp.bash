FASTA=./Mycobacterium_bovis_AF212297_LT78304.fa
OUT=results/simulated
dwgsim \
    -e 0.001-0.01 \
    -i \
    -d 330 \
    -N $N_read_pairs \
    -r 0 \
    -R 0 \
    -X 0 \
    -y 0.01 \
    -H \
    -z $seed_number \
    $FASTA $OUT 