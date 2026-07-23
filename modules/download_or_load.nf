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

//load dataset manually downloaded from the NCBI from the zip file or from the ncbi_dataset folder already extracted
process LOAD_MANUAL_DATASET {

    input:
    path dataset

    output:
    path "ncbi_dataset/data/assembly_data_report.jsonl"
    path "ncbi_dataset/data/*/*.fna"

    script:
    """
    if [[ "${dataset.name}" == *.zip ]]; then

        unzip -q ${dataset} -d .

    else

        cp -r ${dataset} ncbi_dataset

    fi
    """
}
