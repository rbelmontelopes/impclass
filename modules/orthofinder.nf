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
