#!/usr/bin/env nextflow 
/*
    Authors:
    - Katie Evans <kathrynevans2015@u.northwestern.edu>
*/

nextflow.preview.dsl=2
date = new Date().format( 'yyyyMMdd' )

contigs = Channel.from("I","II","III","IV","V","X")
params.cores = 4

// imputation parameters
// try same window for all chromosomes, if errors might need to adjust. ce had different param for X and for autosomes
// "window,overlap"
params.chrI = "5,2"
params.chrII = "5,2"
params.chrIII = "5,2"
params.chrIV = "5,2"
params.chrV = "5,2"
params.chrX = "5,2"


// set debug mode
if(params.debug) {
    println """

    *** Using debug mode ***

    """
    params.vcf = "${workflow.projectDir}/test_data/WI.20201230.hard-filter.vcf.gz"
    params.out = "impute-${date}-debug"
} else {
    params.out = "impute-${date}"

    // Check inputs
    if ( params.vcf==null ) error "Parameter --vcf is required. Specify path to the full vcf."
    if ( params.species==null) error "Parameter --species is required. Please provide species c_elegans, c_tropicalis or c_briggsae"
}


def log_summary() {
/*
    Generates a log
*/

out = '''

-----------------
--- IMPUTE-NF ---
-----------------

Subset isotype reference strains from hard-filter vcf, create SNV-only VCF and impute VCF.

''' + """

nextflow main.nf --debug=true

nextflow main.nf --vcf=hard-filtered.vcf --sample_sheet=sample_sheet.tsv --species=c_elegans 

    parameters                 description                                              Set/Default
    ==========                 ===========                                              ========================
    --debug                    Set to 'true' to test                                    ${params.debug}
    --species                  Species: 'c_elegans', 'c_tropicalis' or 'c_briggsae'     ${params.species}
    --vcf                      hard filtered vcf to calculate variant density           ${params.vcf}
    --out                      (Optional) output folder name                            ${params.out}
    --chrI | chrII | chrIII,   Window and overlap for each chromosome                   I:${params.chrI} | II:${params.chrII} | III:${params.chrIII} | IV:${params.chrIV} | V:${params.chrV} | X:${params.chrX}
     chrIV | chrV | chrX
    }
 
    username                                                                      ${"whoami".execute().in.text}

    HELP: http://andersenlab.org/dry-guide/pipeline-postGATK   
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId] 
"""
out
}


log.info(log_summary())


if (params.help) {
    exit 1
}


// prepare inputs
input_vcf = Channel.fromPath("${params.vcf}")
input_vcf_index = Channel.fromPath("${params.vcf}.tbi")


workflow { 

    
    // impute vcf
    input_vcf.combine(input_vcf_index) | subset_snv
    contigs.combine(subset_snv.out) | imputation
    imputation.out 
        .flatten()
        .toSortedList() | concat_imputed

}




/*==============================================
~ ~ ~ > *   Impute hard filter VCF     * < ~ ~ ~
==============================================*/

process subset_snv {

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    publishDir "${params.out}/variation", mode: 'copy'

    input:
        tuple file(hardvcf), file(hardvcf_index)

    output:
        tuple file("WI.${date}.hard-filter.isotype.SNV.vcf.gz"), file("WI.${date}.hard-filter.isotype.SNV.vcf.gz.tbi")

    """
    bcftools view -O u ${hardvcf} | \
    bcftools view -O v --types snps --min-af 0.000001 --max-af 0.999999 | \
    vcffixup - | \
    bcftools view --threads=3 -O z > WI.${date}.hard-filter.isotype.SNV.vcf.gz
    
    bcftools index --tbi WI.${date}.hard-filter.isotype.SNV.vcf.gz

    """

}

// going to need different chr maps for different species.
process imputation { 

    //errorStrategy 'ignore'

    tag {CHROM} 
    cpus params.cores 
    // memory { 25.GB + 10.GB * task.attempt }
    memory '35 GB'


    conda "/projects/b1059/software/conda_envs/beagle"

    input:
        tuple val(CHROM), file(hardvcf), file(hardvcf_index)

    output:
        tuple file("${CHROM}.b5.vcf.gz"), file("${CHROM}.b5.vcf.gz.csi")

    """

    # get window and overlap parameters
    if [ ${CHROM} == "I" ]
    then
        win=`echo "${params.chrI}" | cut -d ',' -f 1`
        ov=`echo "${params.chrI}" | cut -d ',' -f 2`
    elif [ ${CHROM} == "II" ]
    then
        win=`echo "${params.chrII}" | cut -d ',' -f 1`
        ov=`echo "${params.chrII}" | cut -d ',' -f 2`
    elif [ ${CHROM} == "III" ]
    then
        win=`echo "${params.chrIII}" | cut -d ',' -f 1`
        ov=`echo "${params.chrIII}" | cut -d ',' -f 2`
    elif [ ${CHROM} == "IV" ]
    then
        win=`echo "${params.chrIV}" | cut -d ',' -f 1`
        ov=`echo "${params.chrIV}" | cut -d ',' -f 2`
    elif [ ${CHROM} == "V" ]
    then
        win=`echo "${params.chrV}" | cut -d ',' -f 1`
        ov=`echo "${params.chrV}" | cut -d ',' -f 2`
    else 
        win=`echo "${params.chrX}" | cut -d ',' -f 1`
        ov=`echo "${params.chrX}" | cut -d ',' -f 2`
    fi

    beagle -Xmx98g chrom=${CHROM} window=\$win overlap=\$ov impute=true ne=100000 nthreads=1 imp-segment=0.5 imp-step=0.01 cluster=0.0005 gt=${hardvcf} map=${workflow.projectDir}/bin/${params.species}/chr${CHROM}.map out=${CHROM}.b5

    bcftools index ${CHROM}.b5.vcf.gz

    """

}


process concat_imputed { 

    //errorStrategy 'ignore'

    publishDir "${params.out}/variation/", mode: 'copy'

    conda "/projects/b1059/software/conda_envs/popgen-nf_env"

    memory '64 GB'
    cpus 20

    input:
        file("*")

    output:
        tuple file("WI.${date}.impute.isotype.vcf.gz"), file("WI.${date}.impute.isotype.vcf.gz.tbi")

    """
    bcftools concat *.b5.vcf.gz > WI.${date}.impute.isotype.vcf

    bgzip WI.${date}.impute.isotype.vcf

    bcftools index -t WI.${date}.impute.isotype.vcf.gz

    bcftools stats --verbose WI.${date}.impute.isotype.vcf.gz > WI.${date}.impute.isotype.stats.txt

    """
}