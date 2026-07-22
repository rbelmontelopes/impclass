# Objectives of this pipeline

First of all, ***this pipeline is not a primary classification tool*** (although you could use it as one if you have enough computer resources), and is instead intended to improve pre-existing classification. The motivation behind this pipeline was the classification of metagenome assembled genomes from fungi and algae from Antarctic cryptoendolithic communities, for which BUSCO (Tegenfeldt et al. 2025) autolineage was returning classification at Phyla, Division or Order levels, and for which Eukcc2 (Saary et al 2020) was sometimes able to arrive to Family or Genera level, with the caveat that for lichenized fungi it was classifing most of them as Leotiomyceta or Parmeliaceae, but no Parmeliaceae was ever reported from this communities using culturomics or targeted metagenomics. Aiming to confirm or improve the given classifications, the pipeline downloads genomes from the NCBI, uses BUSCO to generate sets of single copy markers from the selected BUSCO database, constructs a core database with Orthofinder3 (Emms et al. 2026), and then add the the genomes to this core database. The use of the BUSCO databases to generate the sets of markers was selected as it can generate hundreds to thousands of markers specific for a selected database, while at the same time avoiding to use all the predicted proteins of the given genomes, allowing for a speed-up in the Orthofinder steps. Additionally, although the pipeline uses only BUSCO's predicted as single copy, these are not always single copy in comparison between a large number of species, but the use of the Orthofinder3 assign to a core database option allows the use of duplicated genes as it infers a species tree using Astral-Pro (Zhang et al 2025), which was designed to handle multi-copy genes. 

In a certain sense, the present pipeline follows a similar logic used in Buscogeny (Webster & Chapman 2026) with the use of BUSCO markers to infer trees, which also also allows the use of nucleotide sequences, while the present pipeline is restricted to proteins. The two methods differ in their functionalities for downloading genomes from the NCBI (integrated in the present pipeline and with an auxiliary script in Buscogeny), and specially in the way the trees are inferred, with Buscogeny using a concatenated supermatrix and Maximum Likelihood (ML), while the present pipeline use the inference of ML gene trees for the production of a species tree based on the multispecies coalescent process. The present pipeline also allows the addition of new genomes to previously computed Orthofinder3 runs using its assign function, what avoids to running the analysis from the start every time there is need to add a new genome.

# Installation

The easiest way is simple clone the git repository, but is also possible to download everything manually.

```bash
git clone https://github.com/rbelmontelopes/impclass.git
```

# Usage

You can edit the parameters in the params.yml file, or overide in the command line. A example run that will create a reference orthofinder database and run the genomes to identify would be like. The *--busco_db* parameter is mandatory in any case, while the *--taxon* parameter is needed in any run that will construct a database (as the taxon name indicates which genomes to download). Below are some examples (assuming running it from the folder of the cloned repository).

```bash
nextflow run impclass.nf -params-file params.yml --taxon Ascomycota
```
OR to construct a database and identify the genomes
```bash
nextflow run impclass.nf -params-file params.yml --taxon Ascomycota --identify_dir PATH/to_identify  --busco_db ascomycota_odb10 
```
A example run with a prebuild database
```bash
nextflow run impclass.nf -params-file params.yml --orthodb PATH/orthofinder_db/Results_30jun --identify_dir PATH/to_identify  --busco_db ascomycota_odb10 
```
A example run for just constructing a database
```bash
nextflow run impclass.nf -params-file params.yml --taxon Ascomycota  --busco_db ascomycota_odb10 --busco_downloads my_folder/busco_downloads
```
The first time the pipeline is run, it will create the conda enviroments that are needed to run all the steps, what could take some time. The pipeline is enabled to use mamba instead of conda to be faster, but you need to have mamba installed, otherwise will fallback to using conda to create the environments. If you already have all the dependencies installed, an option is to change in the *impclass.nf* the lines *conda "${baseDir}/envs/ncbi.yml"*, *conda "${baseDir}/envs/busco.yml"*, and *conda "${baseDir}/envs/orthofinder.yml"* to the respective paths were the conda enviroments are (i.e. ~/anaconda/envs/orthofinder) or simply remove these lines if everything is already reachable from the base terminal.

## Additional parameters and others

The pipeline have a lot of parameters (see the *params.yml* file for more details), of which the only obligatory in all cases is the busco database to be used (i.e. *--busco_db ascomycota_odb10*). You can change these directly in the *params.yml* file or pass the parameters using the call in the command line.

Note that the database name needs to be recognized by BUSCO as a valid database, with normally have the taxa names all in lowercase followed by the database version (i.e. *_odb10* or *_odb12*). Check the BUSCO site (https://busco-data.ezlab.org/v6/data/lineages/) to know which databases are avaliable (or use the command *busco --list-datasets8*). If you already have the BUSCO databases downloaded is highly indicate to set the parameter *--busco_downloads* to the directory where the databases are, otherwise it will perform the download inside the folder from where the pipeline is running. 

BUSCO is set to use miniprot as default gene predictor, but sometimes it can crash if running on very fragmented genomes, and in this cases metaeuk seems to perform better (but runtimes are longer). You can change the gene predictor used with the parameter *--busco_predictor* for any of the options BUSCO accepts (*augustus*, *metaeuk*, or *miniprot* in the version used in the pipeline)

If you do not pass a previously build orthofinder database folder (i.e.  *--orthodb ./orthofinder_db/Results_30jun*, always use the Results_XXXX folder), is obligatory to indicate the taxon that will be used for downloading the genomes and constructing the database (i.e. *--taxon Ascomycota*). By default the pipeline will download only reference genomes with assembly level indicated as *complete,chromosome,scaffold* (controled by the *--assembly_level* parameter). 

Is highly indicated that you set the number of cpus for BUSCO (i.e. *--cpus 8*) as well as number of parallel BUSCO runs you can do in your system for the reference genomes (i.e *--busco_parallel_reference 1*) and for the genomes that you want to identify (i.e *--busco_parallel_identify 1*). Multiple BUSCO runs at the same time can slow this step, so set the number of parallel tasks accordingly to the capacity of your system. On a laptop with 20 cores one parallel task for each seems to be the best option.

You can also adjust the number of threads for Orthofinder steps with (i.e *--orthofinder_threads 12 --orthofinder_alg 4*).

You can change the default output dir (results) with *--outdir NEW_DIR_NAME*

# Expected results

The basic structure of the results directory is show below.

```bash
.
├── busco_proteins
├── identify_buscos
├── ncbi_dataset.zip
├── orthofinder_assign
└── orthofinder_db
```

-**busco_proteins** contains the single copy BUSCOs predicted for the reference genomes and the summary of the BUSCO runs for each genome

-**identify_buscos** contains the single copy BUSCOs predicted for the genomes to identify and the summary of the BUSCO runs for each genome

-**ncbi_dataset.zip** contains the dehidrated genomes downloaded from the NCBI

-**orthofinder_assign** contains the results of the assign step of Orthofinder3, which is where you will find the species tree with the genomes to identify inside the **Species_Tree folder**. ***The SpeciesTree_rooted.txt*** file inside this folder is the tree you need to look at to see in which clade your genomes ended, as **currently there is no automatic output of the classification** of the genomes to identify. The folder includes also many other results files that are explained in details in the Orthofinder sites (https://github.com/OrthoFinder/OrthoFinder) and (https://orthofinder.org/). 

-**orthofinder_db** contains the core database used to assign the new genomes. This can be re-used for further analysis (i.e. if you find that just the NCBI reference genomes do not give you enough resolution, add additional genomes to the directory with genomes to identify an run the pipeline indicating the folder with the results of the core database)

### Technical stuff

The pipeline was build on Ubuntu 22.04 and nextflow 24.04.4, and also tested in Ubuntu 26 with nextflow 26.04.06, and is expected to work in other Unix based systems, but is unsure if it will work on Windows based machines due to the use of Unix command line tools.

It relies on Conda/Mamba to install the needed enviroments (*ncbi_datasets_cli*,*BUSCO*, and *Orthofinder3*), and it also uses Unix command line tools as *cat*, *grep*, *jq*, and *mkdir*.


# References

Emms, D.M., Liu, Y., Belcher, L. et al. OrthoFinder: improved phylogenetic orthology inference with enhanced accuracy and scalability. Nat Methods (2026). https://doi.org/10.1038/s41592-026-03126-6

Saary P., Mitchell A.L., Finn R.D. Estimating the quality of eukaryotic genomes recovered from metagenomic analysis with EukCC. Genome Biol 21, 244 (2020). https://doi.org/10.1186/s13059-020-02155-4

Tegenfeldt F., Kuznetsov D., Manni M., Berkeley M., Zdobnov E.M., Kriventseva E.V. OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. Nucleic Acids Research, Volume 53, Issue D1, 6 January 2025, Pages D516–D522, https://doi.org/10.1093/nar/gkae987

Webster J., Chapman T.A. Buscogeny: A BUSCO leveraged phylogenomic tree builder. Int Microbiol 29, 255–267 (2026). https://doi.org/10.1007/s10123-025-00752-6

Zhang C., Nielsen R., Mirarab S. ASTER: A Package for Large-Scale Phylogenomic Reconstructions, Molecular Biology and Evolution, Volume 42, Issue 8, August 2025, msaf172, https://doi.org/10.1093/molbev/msaf172


