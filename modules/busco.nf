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
    path "single_copy_buscos/*.FAILED", optional: true

    script:
	"""
	mkdir -p single_copy_buscos

	busco \
    	-i ${fasta} \
    	-m genome \
    	-l ${params.busco_db} \
    	-c ${task.cpus} \
    	-o ${fasta.simpleName}_busco_genome_${params.busco_db} \
    	--${params.busco_predictor} \
    	-f \
    	--download_path ${params.busco_downloads} \
    	--offline

	BUSCO_DIR="${fasta.simpleName}_busco_genome_${params.busco_db}/run_${params.busco_db}"
	FAA_DIR="\$BUSCO_DIR/busco_sequences/single_copy_busco_sequences"

	if ls "\$FAA_DIR"/*.faa >/dev/null 2>&1; then

    		cat "\$FAA_DIR"/*.faa \
        	> single_copy_buscos/${fasta.simpleName}.faa

    		cp "\$BUSCO_DIR/short_summary.txt" \
        	single_copy_buscos/${fasta.simpleName}_${params.busco_db}.summary.txt

	else

    	{
        	echo "Genome: ${fasta.simpleName}"
        	echo "Reason: BUSCO completed but no single-copy BUSCO proteins were generated."
        	echo
        	grep -A8 "^C:" "\$BUSCO_DIR/short_summary.txt"
    	} > single_copy_buscos/${fasta.simpleName}.FAILED

    	cp "\$BUSCO_DIR/short_summary.txt" \
        	single_copy_buscos/${fasta.simpleName}_${params.busco_db}.summary.txt \
        	2>/dev/null || true

	fi
	"""
	
}	

// report any BUSCO run that failed	
process REPORT_FAILED_BUSCOS {

    	publishDir "${params.outdir}", mode: 'copy'

    	input:
    	path failed

    	output:
    	path "BUSCO_failed_genomes.txt"

    	script:
    	"""
    	cat ${failed} > BUSCO_failed_genomes.txt
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
    path "single_copy_buscos/*.FAILED", optional: true

    script:
    """
    mkdir -p single_copy_buscos

	busco \
    	-i ${fasta} \
    	-m genome \
    	-l ${params.busco_db} \
    	-c ${task.cpus} \
    	-o ${fasta.simpleName}_busco_genome_${params.busco_db} \
    	--${params.busco_predictor} \
    	-f \
    	--download_path ${params.busco_downloads} \
    	--offline

	BUSCO_DIR="${fasta.simpleName}_busco_genome_${params.busco_db}/run_${params.busco_db}"
	FAA_DIR="\$BUSCO_DIR/busco_sequences/single_copy_busco_sequences"

	if ls "\$FAA_DIR"/*.faa >/dev/null 2>&1; then

    		cat "\$FAA_DIR"/*.faa \
        	> single_copy_buscos/${fasta.simpleName}.faa

    		cp "\$BUSCO_DIR/short_summary.txt" \
        	single_copy_buscos/${fasta.simpleName}_${params.busco_db}.summary.txt

	else

    	{
        	echo "Genome: ${fasta.simpleName}"
        	echo "Reason: BUSCO completed but no single-copy BUSCO proteins were generated."
        	echo
        	grep -A8 "^C:" "\$BUSCO_DIR/short_summary.txt"
    	} > single_copy_buscos/${fasta.simpleName}.FAILED

    	cp "\$BUSCO_DIR/short_summary.txt" \
        	single_copy_buscos/${fasta.simpleName}_${params.busco_db}.summary.txt \
        	2>/dev/null || true

	fi
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
