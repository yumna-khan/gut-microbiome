DB=~/scratch/kraken2_db
INPUT_DIR=~/scratch/asg3/kraken2_output
OUTPUT_DIR=~/scratch/asg3/bracken_output

mkdir -p ${OUTPUT_DIR}

for REPORT in ${INPUT_DIR}/*.report; do

    SAMPLE=$(basename ${REPORT} .report)

    bracken \
        -d ${DB} \
        -i ${REPORT} \
        -o ${OUTPUT_DIR}/${SAMPLE}.bracken \
        -w ${OUTPUT_DIR}/${SAMPLE}_bracken.report \
        -r 150 \
        -l S \
        -t 10

done

echo "All samples complete!"
