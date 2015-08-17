#!/usr/bin/python

##########
# miRNA discovery (in carrot)
# Arlo Hoogeveen & Sam Nooij
##########

'''
One file to rule them all.
This snakefile wraps up the whole miRNA analysis that we have
done in carrot. It starts after the QC and ends with a
complete list of the miRNAs found with miRBase and miRDeep-P.
It goes through these steps:
== preparation ==
1. Find spike-ins
2. Analyse spike-ins
3. Remove spike-ins and size restriction ## -here already?
                                ## what about tRNA and rRNA?
4. Align reads to miRBase (selected part)
5. Align reads to tRNA database
6. Align reads to rRNA database
7. Label fastq file with alignment scores
8. Discard tRNA and rRNA reads
== miRBase miRNA identification ==
9. Generate list of miRNAs?
== miRDeep ==
10. Align reads to genome
11. Excise candidate hairpins from genome
12. Fold hairpins with RNAfold
13. Align reads to hairpins
14. miRDeep predictions
15. Filter predictions based on redundancy and plant miRNA
    criteria by Meyers et al. (2008)
16. Check miRDeep output with miRBase
17. Check miRDeep output with pakibase
18. Make a table of the identified miRNAs

In the end, it should return comprehensive RESULTS/lists of
the miRNAs we have identified and are confident are trully
miRNAs. Besides, it gives the analysis results from the
spike-ins, which gives an indication of how well the
sequencing has gone.

_________________________________________________________
! Be careful with running this snakefile, the output is
    probably *dozens of gigabytes* large. Be prepared and 
    don't say I didn't warn you!
_________________________________________________________

-- Sam Nooij, 24 April 2015

CHANGELOG:
24-04-2015: Story/comments written
28-04-2015: added variables for databases and samples
            filled in rules for miRDeep
            testing/debugging of bottom 4 rules
29-04-2015  removed QC
            filled in rules for comparing miRDeep-output and
             making RESULTS and for taking the miRNAs found
             with miRBase
            (only the rules for labelling and filtering tRNA/
             rRNA are left now)
            changed tRNA database to "TAIR_ath" (was gtRNAdb
             "dicots")
            updated (bottom) rules to work with multiple samples
12-05-2015  added and edited pakibase comparisons
            edited comparisons of miRDeep output
'''
import os


DATABASES = "./Databases/"
FILTERED_FASTQ = "./Filtered_fastq/"
MIRDEEP_RESULTS = "./miRDeep/"
BOWTIE_OUTPUT = "./Mapping_results/"
REPORTS = "./Report/"
TABLES = "./Report/Tables/"
SCRIPTS = "./Scripts/"
MIRDEEP_SCRIPTS = SCRIPTS + "miRDP/"
RESULTS = "./Results/"

# databases
FASTQ_FOLDER = "./raw/"
SPIKEDB = DATABASES + "SpikeIns/SpikeIns"
RRNADB = DATABASES + "rRNA/ath_silva"
TRNADB = DATABASES + "tRNA/ath"
MIRBASEDB = DATABASES + "miRNA/mature_cca-ath"
MIRBASE_HAIRPIN = DATABASES + "miRNA/hairpins_cca-ath"
PAKIBASE_M = DATABASES + "Barozai/mature"
PAKIBASE_P = DATABASES + "Barozai/EST"
CARROT_GENOME = DATABASES + "Genome/71.chromosomes.fna" 
CARROT_GENOME_BOWTIE = DATABASES + "Genome/Carrot71ch"

# samples
SAMPLES =  [s[:-6] for s in os.listdir(FASTQ_FOLDER) if s.endswith(".fastq")]

#bowtie 2 options
SPIKE_BOWTIE_PARAMS_LIST = [
    "-L 6 ",                     # seed length
    "-i S,0,0.5",                # interval between extracted seeds
    "--norc",                    # do not align to reverse strand
    "--score-min L,-1,-0.6",     # -1-0.6*read_length -- 10% mismatches allowed
    "--no-unal",                 # no unaligned
    "-p 16"                      # number of threads
]
SPIKE_BOWTIE_PARAMS = " ".join(SPIKE_BOWTIE_PARAMS_LIST)

OTHER_RNA_BOWTIE_PARAMS_LIST = [
    "--norc",                    # do not align to reverse strand
    "--no-unal",                 # no unaligned
    "-p 16"                      # number of threads
]
OTHER_BOWTIE_PARAMS = " ".join(OTHER_RNA_BOWTIE_PARAMS_LIST)

MIRNA_BOWTIE_PARAMS_LIST = [
    "-L 6 "                      # seed length
    "-i S,0,0.5",                # interval between extracted seeds
    "--norc",                    # do not align to reverse strand
    "--no-unal",                 # no unaligned
    "--rdg 1,2",                 # read open and extend penalties (default 5,3)
    "--score-min L,-1,-0.6",     # -1-0.6*read_length -- 10% mismatches allowed
    "-p 16"                      # number of threads
]
MIRNA_BOWTIE_PARAMS = " ".join(MIRNA_BOWTIE_PARAMS_LIST)

###
#graphs?
###

rule requirement:
    input:  
            #expand(BOWTIE_OUTPUT + "Sample_{sample}_miRDP_hairpin-miRBase.sam", sample = SAMPLES),
            #RESULTS + "venn_diagram.pdf",
            expand(MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-hairpin.fasta", sample = SAMPLES),
            expand(MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-mature.fasta", sample = SAMPLES),
            expand(BOWTIE_OUTPUT + "Sample_{sample}_miRDP_mature-miRBase.sam", sample = SAMPLES),
            expand(BOWTIE_OUTPUT + "Sample_{sample}_miRDPtoPakibase-hairpin.sam", sample = SAMPLES),
            expand(BOWTIE_OUTPUT + "Sample_{sample}_miRDPtoPakibase-mature.sam", sample = SAMPLES),
            #expand(RESULTS + "Sample_{sample}_mature_notinmiRBase.csv", sample = SAMPLES),
            #expand(RESULTS + "Sample_{sample}_mature-miRBase.csv", sample = SAMPLES),
            expand(RESULTS + "miRBase_found_{sample}.csv", sample = SAMPLES),
            expand(FILTERED_FASTQ + "sample_{sample}_filtered.fastq", sample = SAMPLES),
            expand(REPORTS + "spike_in_analysis{sample}.csv", sample = SAMPLES)

#####
# values for venn diagram
rule venn_counts:
    input:      lit = "Barozai_miRNAs.txt",
                mirbase = RESULTS + "miRBase_found-all.fasta",
                mirdeep = TABLES + "miRDP_found-all.fasta"
    output:     RESULTS + "venn_values.csv"
    message:    "counts of miRNA discovery"
    shell:
        """
        python Scripts/venn-diagram.py {input.lit} {input.mirbase} {input.mirdeep} {output}
        """

# making venn diagram
rule venn_make:
    input:      RESULTS + "venn_values.csv"
    output:     RESULTS + "venn_diagram.pdf"
    message:    "creating venn diagram with R"
    shell:
        """
        Rscript {SCRIPTS}make_venn.R {input} {output}
        """

#####
# make RESULTS out of miRDeep output
rule miRDeep_table_notinmiRBase:
    input:      infile = BOWTIE_OUTPUT + "Sample_{sample}_miRDP_mature-miRBase.sam",
                hairpins = MIRDEEP_RESULTS + "sample_{sample}_filter_P_prediction.txt",
                predictions = MIRDEEP_RESULTS + "sample_{sample}_predictions.txt"
    output:     RESULTS + "Sample_{sample}_mature_notinmiRBase.csv"
    params:     "n"
    message:    "Making a table of the mature miRNAs that have been predicted by miRDeep and do not exist in miRBase."
    shell:
        """
        python2 {SCRIPTS}"miRDP_miRNA_table+MFEI&AMFE.py" {params} {input.infile} {output} {input.hairpins} {input.predictions}
        """
    
## this one may exist in several versions - which miRNAs do we use? cca and ath?
rule miRDeep_table_miRBase:
    input:      infile = BOWTIE_OUTPUT + "Sample_{sample}_miRDP_mature-miRBase.sam",
                hairpins = MIRDEEP_RESULTS + "sample_{sample}_filter_P_prediction.txt",
                predictions = MIRDEEP_RESULTS + "sample_{sample}_predictions.txt"
    output:     RESULTS + "Sample_{sample}_mature-miRBase.csv"
    params:     "y"
    message:    "Making a table of the mature miRNAs that have been predicted by miRDeep and exist in miRBase."
    shell:
        """
        python {SCRIPTS}"miRDP_miRNA_table+MFEI&AMFE.py" {params} {input.infile} {output} {input.hairpins} {input.predictions}
        """

#####    
#compare miRDeep output to known miRNAs
# 3 pakibase
rule miRDeep_to_pakibase_hairpin:
    input:  MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-hairpin.fasta"
    output: BOWTIE_OUTPUT + "Sample_{sample}_miRDPtoPakibase-hairpin.sam"
    params: MIRNA_BOWTIE_PARAMS + "--no-sq -f"
    message:    "Aligning miRDeep predicted pre-miRNA (hairpin) sequences to the ESTs used by Barozai et al. (2013, Pakistan)."
    shell:
        """
        bowtie2 {params} -x {PAKIBASE_P} -U {input} -S {output}
        """

## does still need to have the output checked! the first time I did it by hand!
rule miRDeep_to_pakibase:
    input:  MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-mature.fasta"
    output: BOWTIE_OUTPUT + "Sample_{sample}_miRDPtoPakibase-mature.sam"
    params: MIRNA_BOWTIE_PARAMS + "--no-sq -f"
    message:    "Aligning predicted mature miRNA to those found by Barozai et al. (2013, Pakistan)."
    shell:
        """
        bowtie2 {params} -x {PAKIBASE_M} -U {input} -S {output}
        """
'''
# 2 miRBase hairpin/precursor
rule miRDeep_to_miRBase_hairpin:
    input:  MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-hairpin.fasta"
    output: BOWTIE_OUTPUT + "Sample_{sample}_miRDP_hairpin-miRBase.sam"
    params: MIRNA_BOWTIE_PARAMS + "--no-sq -f"
    message:    "Aligning hairpin (pre-miRNA) sequences from miRDeep to miRBase."
    shell:
        """
        bowtie2 {params} -x {MIRBASE_HAIRPIN} -U {input} -S {output}
        """
'''  
# 1 miRBase mature
rule miRDeep_to_miRBase_mature:
    input:  MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-mature.fasta"
    output: BOWTIE_OUTPUT + "Sample_{sample}_miRDP_mature-miRBase.sam"
    params:     "-L 6 -i S,0,0.5 --norc --rdg 1,2 --score-min L,-1,-0.6 -p 16 --no-sq -f"
    message:    "Aligning mature miRNA sequences from miRDeep to miRBase."
    shell:
        """
        bowtie2 {params} -x {MIRBASEDB} -U {input} -S {output}
        """

#####
rule miRDeep_to_fasta_hairpin:
    input:      MIRDEEP_RESULTS + "sample_{sample}_filter_P_prediction.txt"
    output:     MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-hairpin.fasta"
    message:    "Making fasta file of pre-miRNAs (hairpins) predicted by miRDeep."
    shell:
        """
        python {SCRIPTS}"miRDP_to_fasta_hairpin.py" {input} {output}
        """

rule miRDeep_to_fasta_mature:
    input:      MIRDEEP_RESULTS + "sample_{sample}_filter_P_prediction.txt"
    output:     MIRDEEP_RESULTS + "Sample_{sample}_filter_P_prediction-mature.fasta"
    message:    "Making fasta file of mature miRNAs predicted by miRDeep."
    shell:
        """
        python {SCRIPTS}"miRDP_to_fasta_mature.py" {input} {output}
        """

#####
#miRDeep
rule miRDeep_redundancy_filter_plant_miRNA_criteria:
    input:      script = MIRDEEP_SCRIPTS + "rm_redundant_meet_plant.pl",
                chromosome = MIRDEEP_RESULTS + "chromosome_length.txt",
                precursors = MIRDEEP_RESULTS + "sample_{sample}_precursors.fa",
                predictions = MIRDEEP_RESULTS + "sample_{sample}_predictions.txt"
    output:     no_redundant = MIRDEEP_RESULTS + "sample_{sample}_nr_prediction.txt",
                plant_criteria = MIRDEEP_RESULTS + "sample_{sample}_filter_P_prediction.txt"
    message:    "Removing redundant miRNAs and filtering miRNAs by plant miRNA criteria."
    shell:
        """
        {input.script} {input.chromosome} {input.precursors} {input.predictions} {output.no_redundant} {output.plant_criteria}
        """

rule make_chromosome_length_file:
    input:      genome = CARROT_GENOME
    output:     MIRDEEP_RESULTS + "chromosome_length.txt"
    message:    "creating chromosome_length file."
    shell:
        """
        Scripts/chromosome_length.py {input} {output}
        """

rule miRDeep_predictions:
    input:      script = MIRDEEP_SCRIPTS + "miRDP.pl",
                signatures = MIRDEEP_RESULTS + "sample_{sample}_signatures.txt",
                structures = MIRDEEP_RESULTS + "sample_{sample}_structures.txt"
    output:     MIRDEEP_RESULTS + "sample_{sample}_predictions.txt"
    message:    "miRDeep predictions are made!"
    shell:
        """
        {input.script} {input.signatures} {input.structures} >{output}
        """

rule miRDeep_sort_predictions:
    input:      MIRDEEP_RESULTS + "sample_{sample}_precursors.bst"
    output:     MIRDEEP_RESULTS + "sample_{sample}_signatures.txt"
    message:    "Sorting files."
    shell:
        """
        sort +3 -25 {input} > {output}
        """

rule miRDeep_convert_to_blast2:
    input:      script = MIRDEEP_SCRIPTS + "convert_bowtie_to_blast.pl",
                bowtie = MIRDEEP_RESULTS + "sample_{sample}_precursors.aln",
                fasta_filtered = MIRDEEP_RESULTS + "sample_{sample}_filtered.fa",
                fasta_precursor = MIRDEEP_RESULTS + "sample_{sample}_precursors.fa"
    output:     MIRDEEP_RESULTS + "sample_{sample}_precursors.bst"
    message:    "Converting bowtie-output to BLAST."
    shell:
        """
        {input.script} {input.bowtie} {input.fasta_filtered} {input.fasta_precursor} > {output}
        """

rule miRDeep_align_reads_to_hairpins:
    input:      fasta = MIRDEEP_RESULTS + "sample_{sample}_filtered.fa",
                ref = MIRDEEP_RESULTS + "BowtieRefs/sample_{sample}_precursors.txt"
    output:     MIRDEEP_RESULTS + "sample_{sample}_precursors.aln"
    message:    "Aligning reads to predicted hairpin sequences (putative pre-miRNAs)."
    params:     MIRDEEP_RESULTS + "BowtieRefs/sample_{sample}_precursors"
    shell:
        """
        bowtie -a -v 0 -p 8 {params} -f {input.fasta} > {output}
        """

rule miRDeep_imitate_filtered_fasta:
    input:      script = MIRDEEP_SCRIPTS + "filter_alignments.pl",
                blast = MIRDEEP_RESULTS + "sample_{sample}_miRDeep15.bst",
                fasta = FILTERED_FASTQ + "sample_{sample}_miRDeep_ready.fasta"
    output:     MIRDEEP_RESULTS + "sample_{sample}_filtered.fa"
    message:    "Creating fasta file without rRNAs and tRNAs."
    shell:
        """
        {input.script} {input.blast} -b {input.fasta} > {output}
        """

rule miRDeep_hairpins_to_index:
    input:      MIRDEEP_RESULTS + "sample_{sample}_precursors.fa"
    output:     MIRDEEP_RESULTS + "BowtieRefs/sample_{sample}_precursors.txt"
    message:    "Building bowtie index for all putative pre-miRNAs"
    shell:
        """
        python Scripts/build_database.py {input} {output}
        """

rule miRDeep_RNAfold:
    input:      MIRDEEP_RESULTS + "sample_{sample}_precursors.fa"
    output:     MIRDEEP_RESULTS + "sample_{sample}_structures.txt"
    message:    "Folding putative pre-miRNAs to predict their structure."
    shell:
        """
        cat {input} | RNAfold --noPS > {output}
        """

rule miRDeep_excise_hairpins:
    input:      MIRDEEP_RESULTS + "sample_{sample}_miRDeep15.bst"
    output:     MIRDEEP_RESULTS + "sample_{sample}_precursors.fa"
    message:    "Excising putative pre-miRNAs from the genome."
    shell:
        """
        Scripts/miRDP/excise_candidate.pl {CARROT_GENOME} {input} 250 >{output}
        """

rule miRDeep_sort_cutoff:
    input:      MIRDEEP_RESULTS + "sample_{sample}_miRDeep.bst"
    output:     MIRDEEP_RESULTS + "sample_{sample}_miRDeep15.bst"
    message:    "Filtering alignment to set the maximum miRNA family size to 15 members."
    shell:
        """
        Scripts/miRDP/filter_alignments.pl {input} -c 15 > {output}
        """

rule miRDeep_convert_to_blast1:
    input:      bowtie = BOWTIE_OUTPUT + "sample_{sample}_miRDeep.aln",
                fasta_mirdeep = FILTERED_FASTQ + "sample_{sample}_miRDeep_ready.fasta"
    output:     MIRDEEP_RESULTS + "sample_{sample}_miRDeep.bst"
    message:    "Converting bowtie-output to BLAST."
    shell:
        """
        Scripts/miRDP/convert_bowtie_to_blast.pl {input.bowtie} {input.fasta_mirdeep} {CARROT_GENOME} > {output}
        """

rule miRDeep_align_to_genome:
    input:      FILTERED_FASTQ + "sample_{sample}_miRDeep_ready.fasta"
    output:     BOWTIE_OUTPUT + "sample_{sample}_miRDeep.aln"
    message:    "Aligning the reads to the genome."
    shell:
        """
        bowtie -a -v 0 -p 16 {CARROT_GENOME_BOWTIE} -f {input} > {output}
        """

rule make_mirdeep_fasta:
    input:      FILTERED_FASTQ + "sample_{sample}_filtered.fastq"
    output:     FILTERED_FASTQ + "sample_{sample}_miRDeep_ready.fasta"
    message:    "make the fasta files needed for mirdeep."
    params:     "filtered_read"
    shell:
        """
        python Scripts/make_mirdeep_fasa.py -s -i {input} -o {output} -n {params}
        """

#####
# compare the identified miRNAs with pakibase
        
rule miRNAs_miRBase_pakibase_report:
    input:  BOWTIE_OUTPUT + "sample_{sample}_miRBase-pakibase.sam"
    output: REPORTS + "miRNAs-to-pakibase_with-miRBase.txt"
    message:    "Summarising miRNAs to Barozai check."
    shell:
        """
        python Scripts/check_barozai_reads-miRBase.py
        """

rule map_miRBase_pakibase:
    input:  fastq = FILTERED_FASTQ + "sample_{sample}_miRBase_miRNAs.fasta"
    output: BOWTIE_OUTPUT + "sample_{sample}_miRBase-pakibase.sam"
    message:    "Aligning identified miRNAs to Barozai's miRNAs."
    params: MIRNA_BOWTIE_PARAMS
    shell:
        """
        bowtie2 -f {params} -x {PAKIBASE_M} -U {input} -S {output}
        """

#####
# make a list of miRNAs found with miRBase
rule table_mirbase_miRNA:
    input:      FILTERED_FASTQ + "sample_{sample}_miRBase_miRNAs.fasta"
    output:     RESULTS + "miRBase_found_{sample}.csv"
    message:    "creates a table of found miRNA for easy further study." 
    shell:
        """
        python {SCRIPTS}counts_to_table.py {input} {output}
        """

rule export_miRNAs:
    input:      fastq = FILTERED_FASTQ + "sample_{sample}_filtered.fastq",
                sam = BOWTIE_OUTPUT + "{sample}_mirbase.txt"
    output:     FILTERED_FASTQ + "sample_{sample}_miRBase_miRNAs.fasta"
    message:    "Making a fasta file (collapsed & counted) of all the miRNAs found by aligning sequencing reads to miRBase."
    shell:
        """
        python {SCRIPTS}make_miRBase_miRNAs.py {input.fastq} {input.sam} {output}
        """

rule export_bam_to_text_for_miRBase_miRNAs:
    input:      BOWTIE_OUTPUT + "{sample}_mirbase_sorted.bam"
    output:     BOWTIE_OUTPUT + "{sample}_mirbase.txt"
    message:    "Converting .bam file into .txt."
    shell:
        """
        samtools view {input} > {output}
        """
    
rule sort_bam_for_miRBase_miRNAs:
    input:      BOWTIE_OUTPUT + "{sample}_mirbase.bam"
    output:     BOWTIE_OUTPUT + "{sample}_mirbase_sorted.bam"
    params:     BOWTIE_OUTPUT + "{sample}_mirbase_sorted"
    message:    "Sorting .bam file."
    shell:
        """
        samtools sort {input} {params}
        """
    
rule make_bam_for_miRBase_miRNAs:
    input:      mirbase = DATABASES + "miRBase-cca-ath_mature.fa",
                sam = BOWTIE_OUTPUT + "sample_{sample}_miRbase.sam"
    output:     BOWTIE_OUTPUT + "{sample}_mirbase.bam"
    message:    "Converting .sam file into .bam."
    shell:
        """
        samtools view -bT {input.mirbase} {input.sam} > {output}
        """
#####
# which sequences by Barozai et al. (2013) do we have in our
# sequences bewteen 15 and 30 nt long?
rule pakibase_report:
    input:  BOWTIE_OUTPUT + "sample_{sample}_pakibase.sam"
    output: REPORTS + "miRNAs_from_barozai.txt"
    message:    "Summarising sequences to Barozai check."
    shell:
        """
        python Scripts/check_barozai_reads.py
        """

rule align_sequences_to_pakibase:
    input:  fastq = FILTERED_FASTQ + "sample_{sample}_filtered.fastq"
    output: BOWTIE_OUTPUT + "sample_{sample}_pakibase.sam"
    message:    "Aligning remaining sequences to Barozai's miRNAs."
    params: MIRNA_BOWTIE_PARAMS
    shell:
        """
        bowtie2 {params} -x {PAKIBASE_M} -U {input} -S {output}
        """

#####
'''
rule second_filter:
    input:      FILTERED_FASTQ + "sample_{sample}_no_tRNA_no_rRNA.fastq"
    output:     FILTERED_FASTQ + "sample_{sample}_filtered.fastq"
    params:     cutoff = "30",
                log_file = REPORTS + "removed_reads.csv"
    message:    "Second filtering step : <= 30 size cutoff"
    shell:
        """
        python Scripts/top_size_filter.py {input} {output} {params.cutoff} {params.log_file}
        """
'''

#label fastq
rule remove_tRNA_rRNA:
    input:      tRNA = BOWTIE_OUTPUT + "sample_{sample}_tRNA.sam",
                rRNA = BOWTIE_OUTPUT + "sample_{sample}_rRNA.sam",
                miRNA = BOWTIE_OUTPUT + "sample_{sample}_miRbase.sam",
                fastq = FILTERED_FASTQ + "sample_{sample}_no_spikein.fastq"
    output:     FILTERED_FASTQ + "sample_{sample}_filtered.fastq" #_no_tRNA_no_rRNA
    params:     labels = "tRNA,rRNA,miRNA",
                filters = "tRNA,rRNA",
                log_file = REPORTS + "removed_reads.csv"
    message:    "Labeling and removing the rRNA and tRNA of fastq."
    shell:
        """
        python Scripts/reads_filter.py -p -s {input.tRNA},{input.rRNA},{input.miRNA} -l {params.labels} -i {input.fastq} -n {output} -f {params.filters} -r {params.log_file}
        """

#map to rRNA
rule map_rRNAs:
    input:      FILTERED_FASTQ + "sample_{sample}_no_spikein.fastq"
    output:     BOWTIE_OUTPUT + "sample_{sample}_rRNA.sam"
    message:    "Aligning rRNA."
    shell:
        """
        bowtie2 {OTHER_BOWTIE_PARAMS} -x {RRNADB} -U {input} -S {output}
        """

#map to tRNA
rule map_tRNAs:
    input:      FILTERED_FASTQ + "sample_{sample}_no_spikein.fastq"
    output:     BOWTIE_OUTPUT + "sample_{sample}_tRNA.sam"
    message:    "Aligning tRNA."
    shell:
        """
        bowtie2 {OTHER_BOWTIE_PARAMS} -x {TRNADB} -U {input} -S {output}
        """

#####
#map to miRBase
rule map_miRNAs:
    input:      FILTERED_FASTQ + "sample_{sample}_no_spikein.fastq"
    output:     BOWTIE_OUTPUT + "sample_{sample}_miRbase.sam"
    message:    "Aligning miRNA."
    shell:
        """
        bowtie2 {MIRNA_BOWTIE_PARAMS} -x {MIRBASEDB} -U {input} -S {output}
        """

#####
#remove spikein sequences & select reads for length (maybe collapse similar sequences)
## (how do you select? do we put sizes at the top, for easy manipulation?)
rule first_filter:
    input:      alignment = BOWTIE_OUTPUT + "sample_{sample}_spikeins.sam",
                reads = FASTQ_FOLDER + "{sample}.fastq"
    output:     FILTERED_FASTQ + "sample_{sample}_no_spikein.fastq"
    params:     cutoff = "15,30",
                report = REPORTS + "spikein_report.csv",
                removed = REPORTS + "removed_reads.csv"
    message:    "First filtering step : <30, >15 size cutoff and no spike-ins."
    shell:
        """
        python Scripts/spikin_reads_filter.py -s -i {input.alignment} -q {input.reads} -o {output} -l {params.cutoff} -r {params.report} -p {params.removed}
        """

#####
#spikein analysis
rule spike_analysis:
    input:      BOWTIE_OUTPUT + "sample_{sample}_spikeins.sam"
    output:     REPORTS + "spike_in_analysis{sample}.csv"
    message:    "Performing spike-in analysis."
    shell:
        """
        python Scripts/SpikeAnalysis.py {input} {output}
        """

#map to spikeins
rule map_spikeins:
    input:      FASTQ_FOLDER +"{sample}.fastq"
    output:     BOWTIE_OUTPUT + "sample_{sample}_spikeins.sam"
    message:    "Aligning spike-ins."
    shell:
        """
        bowtie2 {SPIKE_BOWTIE_PARAMS} -x {SPIKEDB} -U {input} -S {output}
        """

#####
