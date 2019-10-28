
params.samples = 'input/*_R{1,2}_001.fastq.gz'
params.bwa_threads = 1
params.reference_genome =  '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.fasta'
params.temp_dir = '/tmp/.temp'


reference_genome = file(params.reference_genome)

Channel
  .fromFilePairs(params.samples, flat: true) 
  .set { reads }



process trim_reads_with_fastp{
    input:
    set val(id), file(read1), file(read2) from reads 

    output:
    set val(id), file("${id}_R1_trimmed.fastq.gz") into trimmed_reads_channel_r1
    set val(id), file("${id}_R2_trimmed.fastq.gz") into trimmed_reads_channel_r2
    file "${id}_trimmed.html" into fastp_html
    file "${id}_trimmed.json" into fastp_json

	"""
	fastp -i $read1 -I $read2 -o ${id}_R1_trimmed.fastq.gz -O ${id}_R2_trimmed.fastq.gz -h ${id}_trimmed.html -j ${id}_trimmed.json

	"""
}


process align_reads_with_bwa{
    input:
    set val(id1), file(read1) from trimmed_reads_channel_r1
    set val(id2), file(read2) from trimmed_reads_channel_r2

    output:
    file "${id1}.bam" into inital_bam_channel
    file "${id1}.bam.bai" into inital_bam_indexes_channel

    """
	bwa mem -t $params.bwa_threads -M -R '@RG\\tID:test.test\\tCN:test\\tSM:${id1}\\tLB:test\\tPL:ILLUMINA' $reference_genome $read1 $read2 | samtools view -Sb - | samtools sort -T $params.temp_dir -O bam > ${id1}.bam
    samtools index ${id1}.bam
	"""
}

process merge_bams_and_remove_duplicates{
    input:
        set val(key), file(bams) from inital_bam_channel.map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple()
        set val(key2), file(bam_indexes) from inital_bam_indexes_channel.map { file ->
        def key = file.name.toString().tokenize('_').get(0)
        return tuple(key, file)
     }
    .groupTuple()

    output:
        file("test.bam") into final_bam_channel
        file("test_MarkDuplicatesMetrics.txt")

    script:
    """
    picard MarkDuplicates {bams.collect { "I=$it " }.join()} O=test.bam METRICS_FILE=test_MarkDuplicatesMetrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT MAX_RECORDS_IN_RAM=200000 TMP_DIR=/tmp/
    """



}