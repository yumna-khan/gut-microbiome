# Create a loop
for R1 in *_1.fastq; do

    # Get sample name
    sample=$(basename ${R1} _1.fastq)
    echo "Processing: ${sample}"

    fastp \
        --in1 ${sample}_1.fastq \
        --in2 ${sample}_2.fastq \
        --out1 trimmed_reads/${sample}_1_trimmed.fastq \
        --out2 trimmed_reads/${sample}_2_trimmed.fastq \
        --html fastp_logs/${sample}_report.html \
        --json fastp_logs/${sample}_report.json \
        --detect_adapter_for_pe \
	--qualified_quality_phred 20 \
	--length_required 50 \
        --thread 4

    echo "Done: ${sample}"

done

echo "All samples complete!"
