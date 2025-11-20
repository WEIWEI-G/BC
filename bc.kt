#bash
REF_GENOME="/mnt/d/reference/human_grch38/Homo_sapiens.GRCh38.dna.prima
ry_assembly.fa"
KNOWN_SITES_MILLS="/mnt/d/reference/human_grch38/Mills_and_1000G_gold_s
tandard.indels.GRCh38.sites.vcf.gz"
KNOWN_SITES_1000G="/mnt/d/reference/human_grch38/1000G_phase1.indels.hg
19.sites.vcf.gz"
FASTQ_R1="01_raw_data/your_sample_R1.fq.gz"
FASTQ_R2="01_raw_data/your_sample_R2.fq.gz"
SAMPLE_NAME="my_sample"
RG_ID="unique_run_id_001" 
RG_PL="ILLUMINA" 
THREADS=8 
PREFIX="03_alignment_results/${SAMPLE_NAME}"
SORTED_BAM="${PREFIX}.sorted.bam"
DEDUP_BAM="${PREFIX}.sorted.markdup.bam"
DEDUP_METRICS="${PREFIX}.markdup_metrics.txt"
RECAL_TABLE="03_alignment_results/${SAMPLE_NAME}.recal_data.table"
FINAL_BAM="04_final_bam/${SAMPLE_NAME}.analysis_ready.bam"
RG_STRING="@RG\\tID:${RG_ID}\\tSM:${SAMPLE_NAME}\\tPL:${RG_PL}"

fastqc ${FASTQ_R1} ${FASTQ_R2} -o 02_qc_reports
bwa mem -t ${THREADS} -R "${RG_STRING}" ${REF_GENOME} ${FASTQ_R1}
${FASTQ_R2} | samtools sort -@ ${THREADS} -o ${SORTED_BAM} -
gatk MarkDuplicates \
 -I ${SORTED_BAM} \
 -O ${DEDUP_BAM} \
 -M ${DEDUP_METRICS}
samtools index ${DEDUP_BAM}
gatk BaseRecalibrator \
 -R ${REF_GENOME} \
 -I ${DEDUP_BAM} \
 --known-sites ${KNOWN_SITES_MILLS} \
 --known-sites ${KNOWN_SITES_1000G} \
 -O ${RECAL_TABLE}

echo "--- Step 5b: Applying BQSR model ---"
gatk ApplyBQSR \
-R ${REF_GENOME} \
 -I ${DEDUP_BAM} \
 --bqsr-recal-file ${RECAL_TABLE} \
 -O ${FINAL_BAM}

BAM_ANALYSIS_READY"${SAMPLE_ID}.analysis_ready.bam
GVCF_FILE"${SAMPLE_ID}.g.vcf.gz"
gatk HaplotypeCaller \
 -R ${REF_GENOME} \
 -I ${BAM_ANALYSIS_READY} \
 -O ${GVCF_FILE} \
 -ERC GVCF

RAW_VCF_FILE"${SAMPLE_ID}.raw.vcf.gz"
gatk GenotypeGVCFs \
 -R ${REF_GENOME} \
 -V ${GVCF_FILE} \
 -O ${RAW_VCF_FILE}
gatk VariantFiltration \
 -R ${REF_GENOME} \
 -V ${RAW_INDELS_FILE} \
 -O ${FILTERED_INDELS_FILE} \
 --filter-expression "QD < 2.0" --filter-name "QD2" \
 --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
 --filter-expression "FS > 200.0" --filter-name "FS200" \
 --filter-expression "ReadPosRankSum < -20.0" --filter-name 
"ReadPosRankSum-20

#R
packages_to_install <- c(
  "tidyverse",    # For data manipulation and visualization (ggplot2, dplyr, etc.)
  "maftools",     # For mutation analysis and visualization (oncoplots, lollipop plots)
  "pheatmap",     # For creating pretty heatmaps
  "ComplexHeatmap", # For advanced heatmaps and co-occurrence plots
  "fmsb",         # For creating radar charts
  "broom",        # For tidying model output
  "ragg"          # For better plot rendering
)

for (pkg in packages_to_install) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    if (pkg == "maftools") {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("maftools")
    }
  }
}

# Load the packages
library(tidyverse)
library(maftools)
library(pheatmap)
library(ComplexHeatmap)
library(fmsb)
library(broom)
library(ragg)

# Set a seed for reproducibility of random data generation
set.seed(42)
# We'll create a synthetic dataset of 1,262 patients
# that mirrors the characteristics described in the methods section.

n_patients <- 1262
patient_ids <- paste0("P", 1:n_patients)

clinical_data <- tibble(
  Tumor_Sample_Barcode = patient_ids,
  age_at_diagnosis = round(rnorm(n_patients, mean = 48, sd = 10)),
  menopause_status = sample(c("Pre", "Post"), n_patients, replace = TRUE, prob = c(0.55, 0.45)),
  ER_status = sample(c("Positive", "Negative"), n_patients, replace = TRUE, prob = c(0.7, 0.3)),
  PR_status = sample(c("Positive", "Negative"), n_patients, replace = TRUE, prob = c(0.6, 0.4)),
  HER2_status = sample(c("Positive", "Negative"), n_patients, replace = TRUE, prob = c(0.2, 0.8)),
  tumor_size_cm = round(rnorm(n_patients, mean = 2.5, sd = 1), 1),
  lymph_node_status = sample(c("Positive", "Negative"), n_patients, replace = TRUE, prob = c(0.4, 0.6)),
  molecular_subtype = sample(c("Luminal A", "Luminal B", "HER2-enriched", "Basal-like"), n_patients, replace = TRUE, prob = c(0.4, 0.3, 0.15, 0.15))
)

genes <- c("BRCA1", "BRCA2", "PALB2", "TP53", "ATM", "CHEK2", "PTEN")
gene_probs <- c(0.35, 0.30, 0.12, 0.08, 0.07, 0.06, 0.02)
variant_class <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site")
variant_probs <- c(0.4, 0.25, 0.15, 0.1, 0.1)

carrier_patients <- sample(patient_ids, size = round(n_patients * 0.15))

maf_data <- map_dfr(carrier_patients, ~{
  tibble(
    Hugo_Symbol = sample(genes, 1, prob = gene_probs),
    Tumor_Sample_Barcode = .x,
    Variant_Classification = sample(variant_class, 1, prob = variant_probs),
    Chromosome = sample(1:22, 1),
    Start_Position = sample(100000:5000000, 1)
  )
}) %>%
  mutate(
    End_Position = Start_Position,
    Reference_Allele = "A",
    Tumor_Seq_Allele2 = sample(c("T", "C", "G"), n(), replace = TRUE),
    Variant_Type = "SNP"
  )


clinical_data <- clinical_data %>%
  mutate(
    pathogenic_variant_carrier = ifelse(Tumor_Sample_Barcode %in% carrier_patients, "Carrier", "Non-carrier")
  )

cat("### Simulated Clinical Data (First 6 Rows):\n")
print(head(clinical_data))

cat("\n### Simulated MAF Data (First 6 Rows):\n")
print(head(maf_data))

cat("\n\n--- Starting Statistical Analysis ---\n")

cat("\n### Pathogenic Variant Status vs. Molecular Subtype\n")
contingency_table <- table(clinical_data$pathogenic_variant_carrier, clinical_data$molecular_subtype)
print(contingency_table)

chi_sq_result <- chisq.test(contingency_table)
cat("\nChi-squared test result:\n")
print(chi_sq_result)

cat("\n### Age at Diagnosis vs. Pathogenic Variant Status\n")
t_test_result <- t.test(age_at_diagnosis ~ pathogenic_variant_carrier, data = clinical_data)
cat("T-test result for age:\n")
print(t_test_result)

cat("\n### Tumor Size vs. Molecular Subtype\n")
anova_result <- aov(tumor_size_cm ~ molecular_subtype, data = clinical_data)
cat("ANOVA result for tumor size:\n")
print(summary(anova_result))











