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
