#!/usr/bin/env nextflow


include { DOWNLOAD_DEHYDRATED; REHYDRATE; LOAD_MANUAL_DATASET } from './modules/download_or_load'
include { RENAME_FASTA } from './modules/rename_fasta'

include {
    PREPARE_BUSCO_DB
    RUN_BUSCO
    RUN_BUSCO_IDENTIFY
    COLLECT_BUSCO_PROTEINS
    COLLECT_BUSCO_IDENTIFY_PROTEINS
    REPORT_FAILED_BUSCOS
} from './modules/busco'

include {
    CREATE_ORTHOFINDER_DB
    RUN_ASSIGN_CORE
    USE_EXISTING_ORTHODB
} from './modules/orthofinder'


workflow {

  // do not run if the parameter busco_db is not indicated
	if (!params.busco_db) {
   	 	error "Missing required parameter: --busco_db"
		}

// One of these must be provided:
//   --orthodb
//   --manual_dataset
//   --taxon
if (!params.orthodb && !params.manual_dataset && !params.taxon) {
    error """
You must provide one of:

  --orthodb <Results_* directory>
  --manual_dataset <NCBI zip or extracted ncbi_dataset folder>
  --taxon <taxon name>
"""
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
	
// check manual NCBI dataset
if (params.manual_dataset) {

    manual_path = file(params.manual_dataset)

    if (!manual_path.exists()) {
        error "manual_dataset does not exist: ${params.manual_dataset}"
    }

}
  
  
  
    busco_db = PREPARE_BUSCO_DB()
   
    
    //if no existing database folder is passed create a database
    
    if (params.orthodb) {

    orthodb = USE_EXISTING_ORTHODB(
        file(params.orthodb)
    )

} else {

    if (params.manual_dataset) {

        data = LOAD_MANUAL_DATASET(
            file(params.manual_dataset)
        )

    } else {

        zip = DOWNLOAD_DEHYDRATED()

        data = REHYDRATE(zip)

    }

    renamed = RENAME_FASTA(
        data[0],      // assembly_data_report.jsonl
        data[1]       // FASTA files
    ).fasta
    
    renamed = renamed.flatten()
    
    reference_proteins = RUN_BUSCO(
    renamed,
    busco_db,
    "reference_buscos"
    )
 
    failed = reference_proteins[2].collect()
    REPORT_FAILED_BUSCOS(failed)

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
