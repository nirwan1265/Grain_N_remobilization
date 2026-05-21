################################################################################
### JOIN GO ANNOTATION TO SUPPLEMENTARY GWAS TABLE
################################################################################

library(readr)
library(dplyr)

################################################################################
### CONFIGURATION
################################################################################

SUPP_TABLE <- "tables/supplementary/SuppTable_GWAS_annotation_soilN_amino.csv"
GO_TABLE <- "tables/supplementary/GO_grain_N.txt"
OUTPUT_FILE <- "tables/supplementary/SuppTable_GWAS_annotation_soilN_amino_with_GO.csv"

################################################################################
### LOAD INPUTS
################################################################################

cat("\n=== Loading tables ===\n")
supp_df <- readr::read_csv(SUPP_TABLE, show_col_types = FALSE)
go_df <- readr::read_tsv(GO_TABLE, show_col_types = FALSE)

required_supp_cols <- c("GeneName", "Phenotypes", "Phenotypes_count", "pvalues")
required_go_cols <- c("GeneID", "GeneName", "Family_Subfamily", "Protein_Class", "GO_MF", "GO_BO", "GO_CC")

missing_supp <- setdiff(required_supp_cols, colnames(supp_df))
missing_go <- setdiff(required_go_cols, colnames(go_df))

if (length(missing_supp) > 0) {
  stop("Missing supplementary-table columns: ", paste(missing_supp, collapse = ", "))
}

if (length(missing_go) > 0) {
  stop("Missing GO-table columns: ", paste(missing_go, collapse = ", "))
}

################################################################################
### JOIN
################################################################################

cat("Supplementary rows:", nrow(supp_df), "\n")
cat("GO rows:", nrow(go_df), "\n")

go_clean <- go_df %>%
  dplyr::distinct(GeneID, .keep_all = TRUE)

joined_df <- supp_df %>%
  dplyr::rename(GeneID = GeneName) %>%
  dplyr::left_join(
    go_clean %>%
      dplyr::rename(
        GeneSymbol = GeneName
      ) %>%
      dplyr::select(
        GeneID,
        GeneSymbol,
        Family_Subfamily,
        Protein_Class,
        GO_MF,
        GO_BO,
        GO_CC
      ),
    by = "GeneID"
  ) %>%
  dplyr::select(
    GeneID,
    GeneSymbol,
    Phenotypes,
    Phenotypes_count,
    pvalues,
    Family_Subfamily,
    Protein_Class,
    GO_MF,
    GO_BO,
    GO_CC
  )

################################################################################
### SAVE
################################################################################

readr::write_csv(joined_df, OUTPUT_FILE)

matched_rows <- sum(!is.na(joined_df$GeneSymbol))
unmatched_rows <- sum(is.na(joined_df$GeneSymbol))

cat("\n=== Done ===\n")
cat("Saved:", OUTPUT_FILE, "\n")
cat("Rows:", nrow(joined_df), "\n")
cat("Matched rows:", matched_rows, "\n")
cat("Unmatched rows:", unmatched_rows, "\n")
cat("\nPreview:\n")
print(utils::head(joined_df, 10))
