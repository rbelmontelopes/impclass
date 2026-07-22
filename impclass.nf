#!/usr/bin/env nextflow

// do not run if the parameter busco_db is not indicated
if (!params.busco_db) {
    error "Missing required parameter: --busco_db"
}

// do not run if there is no taxon indicated to create the database and there is no prebuild orthodb passed 
if (!params.orthodb && !params.taxon) {
    error "When --orthodb is not provided, --taxon is required to build the reference database"
}

//check if a valid BUSCO predictor is indicated
if (!params.busco_predictor in ['metaeuk', 'miniprot']) {
    error "busco_predictor must be either 'metaeuk' or 'miniprot'"
}


// if no busco_downloads path is passed, create a folder in the baseDir
if (!params.busco_downloads) {
    params.busco_downloads = "${baseDir}/busco_downloads"
}


// check if directory to identify exists and contain files with the extension .fna
if (params.identify_dir) {

    identify_path = file(params.identify_dir)

    if (!identify_path.exists()) {
        error "identify_dir does not exist: ${params.identify_dir}"
    }

    identify_files = file("${params.identify_dir}/*.fna")

    if (identify_files.size() == 0) {
        error "No .fna files found in identify_dir: ${params.identify_dir}"
    }
}

// check if the orthodb passed exists
if (params.orthodb) {

    orthodb_path = file(params.orthodb)

    if (!orthodb_path.exists()) {
        error "Provided OrthoFinder database does not exist: ${params.orthodb}"
    }

}


// download the genomes using NCBI datasets. Will download only the reference genomes for a given taxon
process DOWNLOAD_DEHYDRATED {

    tag "${params.taxon}"
    publishDir "${params.outdir}", mode: 'copy'
    
    conda "${baseDir}/envs/ncbi.yml"

    output:
    path "ncbi_dataset.zip"

    script:
    """
    reference_flag=""

    if [ "${params.reference_only}" = "true" ]; then
        reference_flag="--reference"
    fi

    datasets download genome taxon "${params.taxon}" \
        --assembly-level ${params.assembly_level} \
        \$reference_flag \
        --dehydrated \
        --filename ncbi_dataset.zip

    """
}

// rehydrate the downloaded genomes
process REHYDRATE {

    //publishDir "${params.outdir}/rehydrated", mode: 'copy'
    
    conda "${baseDir}/envs/ncbi.yml"

    input:
    path zipfile

    output:	 	
    path "ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl"
    path "ncbi_dataset/ncbi_dataset/data/*/*.fna"
    

    script:
    """
    unzip -q ${zipfile} -d ncbi_dataset
    datasets rehydrate --directory ncbi_dataset
    """
}

// rename the files downloaded from the NCBI to include the species name and accession. If needed filter duplicates if the download was done without the '--referece' flag
process RENAME_FASTA {

    //publishDir "${params.outdir}", mode: 'copy'

    input:
    path report
    path fastas

    output:
    path "renamed/*.fna", emit: fasta

    script:
    """
    mkdir -p renamed

    jq -r '
      .accession as \$acc |
      .organism.organismName as \$org |
      (\$org | split(" ") | .[0:2] | join("_")) + "_" + \$acc
    ' "$report" > mapping.tsv


    for fasta in ${fastas}; do

        acc=\$(basename "\$fasta" | cut -d'_' -f1-2)

        # remove duplicated GCA assemblies if a GCF exists
        if [ "${params.reference_only}" = "false" ]; then

            if [[ "\$acc" == GCA_* ]]; then

                gcf_acc=\$(echo "\$acc" | sed 's/GCA_/GCF_/')

                if grep -q "\$gcf_acc" mapping.tsv; then
                    echo "Skipping duplicate \$acc (reference \$gcf_acc exists)"
                    continue
                fi

            fi
        fi


        name=\$(grep "\$acc" mapping.tsv | cut -f1)

        [ -z "\$name" ] && name="\$acc"

        cp "\$fasta" "renamed/\${name}.fna"

    done
    """
}
// pre-download busco database to avoid crashing during parallel execution in the next stages

process PREPARE_BUSCO_DB {

    tag "busco_db_${params.busco_db}"

    conda "${baseDir}/envs/busco.yml"

    output:
    val params.busco_downloads

    script:
    """
    mkdir -p ${params.busco_downloads}

    if [ ! -d "${params.busco_downloads}/lineages/${params.busco_db}" ]; then

        busco \
            --download_path ${params.busco_downloads} \
            --download ${params.busco_db}

    fi
    """
}
// run busco to get the single copy busco proteins to construct a database
process RUN_BUSCO {
    tag "${busco_label}"
    
    maxForks params.busco_parallel_reference

    //publishDir "${params.outdir}/${busco_label}", mode: 'copy'
    
    conda "${baseDir}/envs/busco.yml"
    
    input:
    path fasta
    val busco_db
    val busco_label


    output:
    path "single_copy_buscos/*.faa"
    path "single_copy_buscos/*.txt"

    script:
    """
    mkdir -p single_copy_buscos
    
 
        busco -i ${fasta} \
            -m genome \
            -l ${params.busco_db} \
            -c ${task.cpus} \
            -o ${fasta.simpleName}_busco_genome_${params.busco_db} \
            --${params.busco_predictor} \
            -f \
            --download_path ${params.busco_downloads} \
            --offline \

        cat ${fasta.simpleName}_busco_genome_${params.busco_db}/run_${params.busco_db}/busco_sequences/single_copy_busco_sequences/*.faa \
            > single_copy_buscos/${fasta.simpleName}.faa
        cp ${fasta.simpleName}_busco_genome_${params.busco_db}/run_${params.busco_db}/short_summary.txt single_copy_buscos/${fasta.simpleName}_${params.busco_db}.summary.txt    
    
    """
}

// run busco to get the single copy busco proteins from the genomes to identify
process RUN_BUSCO_IDENTIFY {
    tag "${busco_label}"
    
    maxForks params.busco_parallel_identify

    //publishDir "${params.outdir}/${busco_label}", mode: 'copy'
    
    conda "${baseDir}/envs/busco.yml"
    
    input:
    path fasta
    val busco_db
    val busco_label


    output:
    path "single_copy_buscos/*.faa"
    path "single_copy_buscos/*.txt"

    script:
    """
    mkdir -p single_copy_buscos
    
    busco -i ${fasta} \
            -m genome \
            -l ${params.busco_db} \
            -c ${task.cpus} \
            -o ${fasta.simpleName}_busco_genome_${params.busco_db} \
            --${params.busco_predictor} \
            -f \
            --download_path ${params.busco_downloads} \
            --offline \

        cat ${fasta.simpleName}_busco_genome_${params.busco_db}/run_${params.busco_db}/busco_sequences/single_copy_busco_sequences/*.faa \
            > single_copy_buscos/${fasta.simpleName}.faa
        cp ${fasta.simpleName}_busco_genome_${params.busco_db}/run_${params.busco_db}/short_summary.txt single_copy_buscos/${fasta.simpleName}_${params.busco_db}.summary.txt    
    """
}

// create a orthofinder DB from the reference single copy buscos
process CREATE_ORTHOFINDER_DB {
    
    publishDir "${params.outdir}", mode: 'copy'
    
    conda "${baseDir}/envs/orthofinder.yml"

    input:
    path proteins
    
    output:
    path "orthofinder_db/Results_*"
    
    script:
    """
    
    orthofinder \
    -t ${task.cpus} \
    -a ${params.orthofinder_alg} \
    -M msa \
    -o orthofinder_db \
    -f ${proteins}

    """
}    

// add the genomes to identify to the orthofinder analysis
process RUN_ASSIGN_CORE {

    publishDir "${params.outdir}/orthofinder_assign", mode: 'copy'

    conda "${baseDir}/envs/orthofinder.yml"

    input:
    path orthofinder_results
    path to_identify

    output:
    path "Results_*", emit: assigned_results

    script:
    """

    ln -s ${orthofinder_results} reference_results

    orthofinder \
        -t ${task.cpus} \
        -a ${params.orthofinder_alg} \
        --assign \
        ${to_identify} \
        --core \
        reference_results

    """
}

// use if there is an existing orthofinder db
process USE_EXISTING_ORTHODB {

    input:
    path db

    output:
    path "orthofinder_db/Results_*"

    script:
    """
    mkdir -p orthofinder_db
    cp -r ${db} orthofinder_db/
    """
}

// put all proteins in the same folder to run orthofinder
process COLLECT_BUSCO_IDENTIFY_PROTEINS {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path faa_files
    path txt_files

    output:
    path "identify_buscos"

    script:
    """
    mkdir -p identify_buscos

    cp ${faa_files} identify_buscos/
    cp ${txt_files} identify_buscos/
    """
}

// put all proteins in the same folder to run orthofinder
process COLLECT_BUSCO_PROTEINS {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    path faa_files
    path txt_files

    output:
    path "reference_buscos"

    script:
    """
    mkdir -p reference_buscos

    cp ${faa_files} reference_buscos/
    cp ${txt_files} reference_buscos/
    """
}




workflow {

  
    busco_db = PREPARE_BUSCO_DB()
   
    
    //if no existing database folder is passed create a database
    
    if (params.orthodb) {

    orthodb = USE_EXISTING_ORTHODB(
        file(params.orthodb)
    )

    } else {
    
    zip = DOWNLOAD_DEHYDRATED()
    data = REHYDRATE(zip)

    renamed=RENAME_FASTA(
        data[0],   // assembly_data_report.jsonl
        data[1]    // FASTA files
    ).fasta
    
    renamed = renamed.flatten()
    
    reference_proteins = RUN_BUSCO(
    renamed,
    busco_db,
    "reference_buscos"
    )
 

    reference_proteins = COLLECT_BUSCO_PROTEINS(
    reference_proteins[0].collect(),
    reference_proteins[1].collect()
    )

    orthodb = CREATE_ORTHOFINDER_DB(reference_proteins[0])

    }
    
    // only run if there the directory with the genomes to identify is passed in the script call otherwise just create the database
    
    if (params.identify_dir) {

        identify_genomes = Channel
            .fromPath("${params.identify_dir}/*.fna")


        identify_proteins = RUN_BUSCO_IDENTIFY(
            identify_genomes,
            busco_db,
            "identify_buscos"
        )


        identify_proteins = COLLECT_BUSCO_IDENTIFY_PROTEINS(
        identify_proteins[0].collect(),
        identify_proteins[1].collect()
        )


        RUN_ASSIGN_CORE(
            orthodb,
            identify_proteins[0]
        )

    }
    
   
}
