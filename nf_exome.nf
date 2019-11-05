#!/usr/bin/env nextflow

/*
========================================================================================
Nextflow pipeline for germline variant calling
========================================================================================

Author: Joseph Halstead

*/


/*
========================================================================================
Set up Parameters
========================================================================================
*/

params.samples = 'input/*/*_R{1,2}_001.fastq.gz'
params.fastp_threads = 1
params.bwa_threads = 1
params.reference_genome =  '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.fasta'
params.sequence_dict = '/media/joseph/Storage/genomic_resources/bwa/human_g1k_v37.dict'
params.bwa_temp_dir = '/tmp/'
params.sequencing_run = 'test'
params.picard_temp_dir = '/tmp/'
params.picard_max_ram_records = 3000000
params.capture_bed = 'config/AgilentOGTFH_ROI_b37.bed'
params.platypus_threads = 1
params.chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y", "MT"]
params.gnomad_exomes = '/media/joseph/Storage/genomic_resources/gnomad/gnomad.exomes.r2.0.1.sites.vcf.gz'
params.gnomad_genomes = '/media/joseph/Storage/genomic_resources/gnomad/gnomad.genomes.r2.0.1.sites.noVEP.vcf.gz'
params.vep_threads = 1
params.vep_cache = '/media/joseph/Storage/genomic_resources/vep_caches/vep'
params.sequencing_centre = 'AWMGL'
params.platypus_min_reads = 5
params.platypus_hap_score_threshold = 20
params.platypus_min_var_freq = 0.05
params.gnomad_max_pop = 0.01


/*
========================================================================================
Define initial files
========================================================================================
*/

reference_genome = file(params.reference_genome)
capture_bed = file(params.capture_bed)
sequence_dict = file(params.sequence_dict)
gnomad_exomes = file(params.gnomad_exomes)
gnomad_genomes = file(params.gnomad_genomes)
vep_cache = file(params.vep_cache)


/*
========================================================================================
Define initial channels
========================================================================================
*/

Channel
  .fromFilePairs(params.samples, flat: true) 
  .set { reads }

chromosomes_ch = Channel
    .from( params.chromosomes)


/*
========================================================================================
Main pipeline
========================================================================================
*/


// Trim reads with fastp
process trim_reads_with_fastp{

    input:
    set val(id), file(read1), file(read2) from reads 

    output:
    set val(id), file("${id}_R1_trimmed.fastq.gz") into trimmed_reads_channel_r1
    set val(id), file("${id}_R2_trimmed.fastq.gz") into trimmed_reads_channel_r2
    file "${id}_trimmed.html" into fastp_html
    file "${id}_trimmed.json" into fastp_json

	"""
	fastp \
    -i $read1 \
    -I $read2 \
    -o ${id}_R1_trimmed.fastq.gz \
    -O ${id}_R2_trimmed.fastq.gz \
    -h ${id}_trimmed.html \
    -j ${id}_trimmed.json \
    -w $params.fastp_threads
	"""
}

// Align reads with bwa
process align_reads_with_bwa{

    input:
    set val(id1), file(read1) from trimmed_reads_channel_r1
    set val(id2), file(read2) from trimmed_reads_channel_r2

    output:
    file "${id1}.bam" into inital_bam_channel
    file "${id1}.bam.bai" into inital_bam_indexes_channel

    script:
    lane = read1.name.toString().tokenize('_').get(2)
    sample_id = read1.name.toString().tokenize('_').get(0)

    """
	bwa mem \
    -t $params.bwa_threads \
    -M \
    -R '@RG\\tID:${params.sequencing_run}.${lane}\\tCN:${params.sequencing_centre}\\tSM:${sample_id}\\tLB:${params.sequencing_run}\\tPL:ILLUMINA' \
    $reference_genome \
    $read1 \
    $read2 | samtools view -Sb - | \
    samtools sort -T $params.bwa_temp_dir -O bam > ${id1}.bam

    samtools index ${id1}.bam
	"""
}

// Merge bams from same sample (different lanes) and remove duplicates
process merge_bams_and_remove_duplicates{

        publishDir "results/final_bams/"

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
    file("${key}_markduplicates.bam") into mark_duplicates_bam_channel
    file("${key}_markduplicates.bai") into mark_duplicates_bam_index_channel
    file("${key}_markduplicate_metrics.txt")

    """
    picard MarkDuplicates \
    ${bams.collect { "I=$it " }.join()} \
    O=${key}_markduplicates.bam \
    METRICS_FILE=${key}_markduplicate_metrics.txt \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT \
    MAX_RECORDS_IN_RAM=${params.picard_max_ram_records} \
    TMP_DIR=${params.picard_temp_dir}
    """
}

// Ensure capture bed is correctly sorted
process sort_capture_bed{

    input:
    file(capture_bed)

    output:
    file("${capture_bed.baseName}_sorted.bed") into sorted_capture_bed

    """
    sort-bed $capture_bed > ${capture_bed.baseName}_sorted.bed
    """
}

// Split capture bed so that we can parallelise variant calling
process split_bed_by_chromosome{

    input:
    val(chr) from chromosomes_ch
    file(sorted_capture_bed) from sorted_capture_bed

    output:
    file("${chr}.bed") into per_chromosome_bed_channel

    """
    bedextract ${chr} $sorted_capture_bed > ${chr}.bed
    """
}

// Call variants with platypus
process call_variants_with_platypus{
    conda 'env/platypus.yaml'

    input:
    file (bams) from mark_duplicates_bam_channel.collect()
    file (bam_indexes) from mark_duplicates_bam_index_channel.collect()
    each file(chr) from per_chromosome_bed_channel

    output:
    file ("${params.sequencing_run}_${chr.baseName}.vcf") into per_chromsome_vcf_channel
    file ("${params.sequencing_run}_${chr.baseName}.log") into per_chromsome_vcf_log_channel

    """
    platypus \
    callVariants \
    --bamFiles=${bams.collect { "$it" }.join(",")} \
    --refFile=$reference_genome \
    --output=${params.sequencing_run}_${chr.baseName}.vcf \
    --regions=${chr} \
    --nCPU=${params.platypus_threads} \
    --minReads=${params.platypus_min_reads} \
    --hapScoreThreshold=${params.platypus_hap_score_threshold} \
    --minVarFreq=${params.platypus_min_var_freq} \
    --logFileName=${params.sequencing_run}_${chr.baseName}.log
    """
}

// PLatypus VCF has strange header so add a correct sequence dict
process fix_platypus_headers{

    input:
    file(per_chr_vcfs) from per_chromsome_vcf_channel

    output:
    file ("${params.sequencing_run}_${chr}_sd.vcf") into per_chromsome_vcf_sd_channel

    script:
    chr = per_chr_vcfs.baseName.toString().tokenize('_').get(1)

    """
    picard UpdateVcfSequenceDictionary \
    I=${per_chr_vcfs} \
    O=${params.sequencing_run}_${chr}_sd.vcf \
    SD=${sequence_dict}
    """
}

// Collect all vcfs from the different chromosomes and edit headers then normalise using vt
process collect_per_chromosome_vcfs_split_and_normalise{

    input:
    file (per_chr_vcfs) from per_chromsome_vcf_sd_channel.collect()

    output:
    file("${params.sequencing_run}_merged_normalised_roi.vcf") into merged_and_normalised_vcf_channel

    """
    bcftools concat ${per_chr_vcfs.collect { "$it " }.join()} | bcftools sort | \
    sed 's/##FORMAT=<ID=GQ,Number=1/##FORMAT=<ID=GQ,Number=./' | \
    sed 's/##FORMAT=<ID=NR,Number=./##FORMAT=<ID=NR,Number=A/' | \
    sed 's/##FORMAT=<ID=NV,Number=./##FORMAT=<ID=NV,Number=A/' | \
    vt decompose -s - | vt normalize -r $reference_genome - > ${params.sequencing_run}_merged_normalised.vcf

    bgzip ${params.sequencing_run}_merged_normalised.vcf
    tabix ${params.sequencing_run}_merged_normalised.vcf.gz

    bcftools view -R $capture_bed ${params.sequencing_run}_merged_normalised.vcf.gz > "${params.sequencing_run}_merged_normalised_roi.vcf"

    """
}

// Annotate using VEP
process annotate_with_vep{

    publishDir "results/annotated_vcf/"

    input:
    file(normalised_vcf) from merged_and_normalised_vcf_channel

    output:
    file("${params.sequencing_run}_merged_normalised_roi_annotated.vcf") into annotated_vcf

    """
    vep \
    --verbose \
    --format vcf \
    --everything \
    --fork $params.vep_threads \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file $normalised_vcf \
    --output_file ${params.sequencing_run}_merged_normalised_roi_annotated.vcf \
    --force_overwrite \
    --cache \
    --dir  $vep_cache \
    --fasta $reference_genome \
    --offline \
    --cache_version 94 \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq \
    --flag_pick \
    --pick_order biotype,canonical,appris,tsl,ccds,rank,length \
    --exclude_predicted \
    --custom ${gnomad_genomes},gnomADg,vcf,exact,0,AF_POPMAX \
    --custom ${gnomad_exomes},gnomADe,vcf,exact,0,AF_POPMAX 
    """
}

// Duplicate annotated VCF channel as multiple downstream processes need it
annotated_vcf.into {
  annotated_vcf_report
  annotated_vcf_sample_names
}

// Extract sample names from vcf
process get_sample_names_from_vcf{

    input:
    file(annotated_vcf) from annotated_vcf_sample_names

    output:
    file("${params.sequencing_run}_sample_names.txt") into sample_names_from_vcf

    """
    bcftools query -l $annotated_vcf > ${params.sequencing_run}_sample_names.txt
    """

}

// Create a channel with sample names as the values
// This is used to create a sample report
sample_names_from_vcf.splitCsv(header:['col1']).map{ row-> tuple(row.col1)}.set { samples_ch }

// Use a custom python script to create csv with variants in
process create_snp_indel_variant_report{

    publishDir "results/variant_reports/"

    input:
    file(annotated_vcf) from annotated_vcf_report
    val sample_names from samples_ch

    output:
    file("${params.sequencing_run}_${sample_names[0]}_filtered.csv") into variant_report_filtered
    file("${params.sequencing_run}_${sample_names[0]}_unfiltered.csv") into variant_report_unfiltered

    """
    create_variant_report.py \
    --vcf $annotated_vcf \
    --sample_id ${sample_names[0]} \
    --worklist_id $params.sequencing_run \
    --output ${params.sequencing_run}_${sample_names[0]} \
    --gnomad_max_af $params.gnomad_max_pop
    """

}