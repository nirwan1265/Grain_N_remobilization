################################################################################
### AMINO-ACID GWAS EXPLORER
################################################################################

library(shiny)
library(DT)
library(dplyr)
library(tibble)
library(ggplot2)

################################################################################
### PATHS AND CONFIGURATION
################################################################################

get_app_home <- function() {
  frame_files <- vapply(sys.frames(), function(x) {
    if (!is.null(x$ofile)) {
      return(x$ofile)
    }
    ""
  }, character(1))

  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) {
    return(dirname(normalizePath(frame_files[[1]], mustWork = TRUE)))
  }

  normalizePath(getwd(), mustWork = TRUE)
}

APP_HOME <- get_app_home()
REPO_ROOT <- normalizePath(file.path(APP_HOME, "..", ".."), mustWork = TRUE)

resolve_required_path <- function(label, candidates) {
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      paste0(
        "Missing required file for ", label, ".\n",
        "Expected one of:\n- ", paste(candidates, collapse = "\n- "), "\n",
        "For shinyapps.io, include the file inside the app bundle (recommended: data/supplementary/)."
      ),
      call. = FALSE
    )
  }
  normalizePath(existing[[1]], mustWork = TRUE)
}

MANHATTAN_DIR <- normalizePath(file.path(APP_HOME, "Figs", "manhattan"), mustWork = TRUE)
QQ_DIR <- normalizePath(file.path(APP_HOME, "Figs", "qq"), mustWork = TRUE)

# Prefer app-local files (works for deployment), keep repo-root fallback for local runs.
BEST_HITS_PATH <- resolve_required_path(
  label = "best-hit GWAS table",
  candidates = c(
    file.path(APP_HOME, "data", "supplementary", "SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv"),
    file.path(APP_HOME, "..", "..", "tables", "supplementary", "SuppTable_amino_gwas_gene_best_by_phenotype_25kb.csv")
  )
)
ANNOTATION_PATH <- resolve_required_path(
  label = "GWAS annotation table",
  candidates = c(
    file.path(APP_HOME, "data", "supplementary", "SuppTable1_GWAS_annotation_soilN_amino_with_GO.csv"),
    file.path(APP_HOME, "..", "..", "tables", "supplementary", "SuppTable1_GWAS_annotation_soilN_amino_with_GO.csv")
  )
)
FAA_PATH <- resolve_required_path(
  label = "FAA phenotype table",
  candidates = c(
    file.path(APP_HOME, "data", "supplementary", "SuppTable_FAA.csv"),
    file.path(REPO_ROOT, "data", "SuppTable_FAA_BLUPs.csv")
  )
)

DEFAULT_INDIVIDUAL_TRAIT <- "P"
DEFAULT_COMBINED_SCOPE <- "__all__"

TRAIT_LABELS <- c(
  A = "Alanine",
  C = "Cysteine",
  D = "Aspartate",
  E = "Glutamate",
  F = "Phenylalanine",
  G = "Glycine",
  H = "Histidine",
  I = "Isoleucine",
  K = "Lysine",
  L = "Leucine",
  M = "Methionine",
  N = "Asparagine",
  P = "Proline",
  Q = "Glutamine",
  R = "Arginine",
  S = "Serine",
  T = "Threonine",
  V = "Valine",
  W = "Tryptophan",
  Y = "Tyrosine",
  Total_N = "Total nitrogen",
  Total_PBAA = "Total protein-bound amino acids"
)

################################################################################
### HELPERS
################################################################################

display_trait_name <- function(trait_name) {
  clean_name <- gsub("_+$", "", trait_name)

  if (trait_name %in% names(TRAIT_LABELS)) {
    paste0(clean_name, " - ", TRAIT_LABELS[[trait_name]])
  } else {
    clean_name
  }
}

trait_group <- function(trait_name) {
  sub("[.].*$", "", trait_name)
}

display_group_name <- function(group_name) {
  clean_name <- gsub("_+$", "", group_name)

  if (group_name %in% names(TRAIT_LABELS)) {
    paste0(clean_name, " - ", TRAIT_LABELS[[group_name]])
  } else {
    clean_name
  }
}

list_plot_traits <- function(dir_path, suffix) {
  files <- list.files(dir_path, full.names = FALSE)
  matches <- files[grepl(paste0(suffix, "[.]png$"), files)]
  sort(sub(paste0(suffix, "[.]png$"), "", matches))
}

empty_individual_table <- function() {
  tibble(
    GeneID = character(),
    GeneSymbol = character(),
    Phenotype = character(),
    Model = character(),
    SNP = character(),
    Chr = integer(),
    SNP_Pos = numeric(),
    P.value = numeric(),
    log10_P = numeric(),
    Relation = character(),
    Distance_to_Gene_bp = numeric(),
    Family_Subfamily = character(),
    Protein_Class = character(),
    GO_MF = character(),
    GO_BO = character(),
    GO_CC = character()
  )
}

empty_combined_gene_table <- function() {
  tibble(
    GeneID = character(),
    GeneSymbol = character(),
    GeneChr = character(),
    GeneStart = numeric(),
    GeneEnd = numeric(),
    Phenotypes = character(),
    Occurrences = integer(),
    Highest_log10P = numeric(),
    Top_Pvalue = numeric(),
    Top_Model = character(),
    Top_SNP = character(),
    Top_Relation = character(),
    Top_Distance_bp = numeric(),
    Family_Subfamily = character(),
    Protein_Class = character(),
    GO_MF = character(),
    GO_BO = character(),
    GO_CC = character()
  )
}

empty_phenotype_table <- function() {
  tibble(
    Phenotype = character(),
    Gene_count = integer(),
    Hit_rows = integer(),
    Highest_log10P = numeric(),
    Top_Gene = character(),
    Models = character()
  )
}

empty_faa_distribution <- function() {
  tibble(
    taxa = character(),
    value = numeric()
  )
}

build_plot_src <- function(plot_type, trait_name) {
  suffix <- if (identical(plot_type, "manhattan")) {
    "_manhattan.png"
  } else {
    "_qq.png"
  }

  dir_path <- if (identical(plot_type, "manhattan")) {
    MANHATTAN_DIR
  } else {
    QQ_DIR
  }

  file_name <- paste0(trait_name, suffix)
  file_path <- file.path(dir_path, file_name)

  if (!file.exists(file_path)) {
    return(NULL)
  }

  paste0("amino-figs/", plot_type, "/", utils::URLencode(file_name, reserved = TRUE))
}

build_faa_distribution <- function(faa_df, trait_name) {
  if (!trait_name %in% names(faa_df)) {
    return(empty_faa_distribution())
  }

  tibble(
    taxa = faa_df$taxa,
    value = as.numeric(faa_df[[trait_name]])
  ) %>%
    dplyr::filter(!is.na(value), is.finite(value))
}

make_faa_distribution_plot <- function(distribution_df, trait_name) {
  median_value <- stats::median(distribution_df$value, na.rm = TRUE)
  mean_value <- mean(distribution_df$value, na.rm = TRUE)

  ggplot(distribution_df, aes(x = value)) +
    geom_histogram(
      bins = 28,
      fill = "#6f8f72",
      color = "white",
      linewidth = 0.2
    ) +
    geom_vline(
      xintercept = median_value,
      color = "#2e4d37",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    labs(
      title = display_trait_name(trait_name),
      subtitle = paste0(
        "FAA distribution across ", nrow(distribution_df),
        " lines | median = ", formatC(median_value, digits = 2, format = "f"),
        " | mean = ", formatC(mean_value, digits = 2, format = "f")
      ),
      x = "Phenotype value",
      y = "Line count"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 10, color = "grey35"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10, color = "black"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

join_annotation <- function(df, annotation_lookup) {
  df %>%
    dplyr::left_join(annotation_lookup, by = "GeneID") %>%
    dplyr::mutate(
      GeneSymbol = dplyr::coalesce(AnnotationGeneSymbol, GeneSymbol)
    ) %>%
    dplyr::select(-AnnotationGeneSymbol)
}

build_individual_hits <- function(best_hits_df, trait_name, cutoff, annotation_lookup) {
  filtered_hits <- best_hits_df %>%
    dplyr::filter(
      Phenotype == trait_name,
      log10_P >= cutoff
    ) %>%
    dplyr::arrange(dplyr::desc(log10_P), GeneSymbol, SNP)

  if (nrow(filtered_hits) == 0) {
    return(empty_individual_table())
  }

  join_annotation(filtered_hits, annotation_lookup)
}

build_combined_gene_hits <- function(best_hits_df, selected_traits, cutoff, annotation_lookup) {
  filtered_hits <- best_hits_df %>%
    dplyr::filter(
      Phenotype %in% selected_traits,
      log10_P >= cutoff
    ) %>%
    dplyr::mutate(Phenotype = factor(Phenotype, levels = selected_traits))

  if (nrow(filtered_hits) == 0) {
    return(empty_combined_gene_table())
  }

  counts_df <- filtered_hits %>%
    dplyr::arrange(Phenotype, dplyr::desc(log10_P)) %>%
    dplyr::group_by(GeneID, GeneSymbol, GeneChr, GeneStart, GeneEnd) %>%
    dplyr::summarise(
      Phenotypes = paste(unique(as.character(Phenotype)), collapse = ";"),
      Occurrences = dplyr::n_distinct(Phenotype),
      .groups = "drop"
    )

  top_rows <- filtered_hits %>%
    dplyr::group_by(GeneID, GeneSymbol, GeneChr, GeneStart, GeneEnd) %>%
    dplyr::slice_max(log10_P, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      GeneID,
      GeneSymbol,
      GeneChr,
      GeneStart,
      GeneEnd,
      Highest_log10P = round(log10_P, 2),
      Top_Pvalue = signif(P.value, 3),
      Top_Model = Model,
      Top_SNP = SNP,
      Top_Relation = Relation,
      Top_Distance_bp = Distance_to_Gene_bp
    )

  counts_df %>%
    dplyr::select(-GeneSymbol) %>%
    dplyr::left_join(top_rows, by = c("GeneID", "GeneChr", "GeneStart", "GeneEnd")) %>%
    join_annotation(annotation_lookup = annotation_lookup) %>%
    dplyr::arrange(dplyr::desc(Occurrences), dplyr::desc(Highest_log10P), GeneSymbol)
}

build_phenotype_summary <- function(best_hits_df, selected_traits, cutoff, annotation_lookup) {
  filtered_hits <- best_hits_df %>%
    dplyr::filter(
      Phenotype %in% selected_traits,
      log10_P >= cutoff
    )

  if (nrow(filtered_hits) == 0) {
    return(empty_phenotype_table())
  }

  top_rows <- filtered_hits %>%
    dplyr::group_by(Phenotype) %>%
    dplyr::slice_max(log10_P, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    join_annotation(annotation_lookup = annotation_lookup) %>%
    dplyr::transmute(
      Phenotype,
      Top_Gene = GeneSymbol,
      Highest_log10P = round(log10_P, 2)
    )

  filtered_hits %>%
    dplyr::group_by(Phenotype) %>%
    dplyr::summarise(
      Gene_count = dplyr::n_distinct(GeneID),
      Hit_rows = dplyr::n(),
      Models = paste(sort(unique(Model)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::left_join(top_rows, by = "Phenotype") %>%
    dplyr::mutate(
      Phenotype = as.character(Phenotype)
    ) %>%
    dplyr::arrange(dplyr::desc(Gene_count), dplyr::desc(Highest_log10P), Phenotype)
}

################################################################################
### STATIC DATA
################################################################################

best_hits_raw <- read.csv(BEST_HITS_PATH, stringsAsFactors = FALSE, check.names = FALSE)
annotation_lookup <- read.csv(ANNOTATION_PATH, stringsAsFactors = FALSE, check.names = FALSE)
names(annotation_lookup)[1] <- "GeneID"
annotation_lookup <- annotation_lookup %>%
  dplyr::transmute(
    GeneID,
    AnnotationGeneSymbol = GeneSymbol,
    Family_Subfamily,
    Protein_Class,
    GO_MF,
    GO_BO,
    GO_CC
  ) %>%
  dplyr::distinct(GeneID, .keep_all = TRUE)
faa_raw <- read.csv(FAA_PATH, stringsAsFactors = FALSE, check.names = FALSE)
names(faa_raw)[1] <- "taxa"

AVAILABLE_TRAITS <- sort(Reduce(
  intersect,
  list(
    unique(best_hits_raw$Phenotype),
    list_plot_traits(MANHATTAN_DIR, "_manhattan"),
    list_plot_traits(QQ_DIR, "_qq")
  )
))

TRAIT_GROUPS <- tibble(
  Phenotype = AVAILABLE_TRAITS,
  Group = vapply(AVAILABLE_TRAITS, trait_group, character(1))
)

AVAILABLE_GROUPS <- sort(unique(TRAIT_GROUPS$Group))
GROUP_CHOICES <- stats::setNames(
  AVAILABLE_GROUPS,
  vapply(AVAILABLE_GROUPS, display_group_name, character(1))
)
TRAIT_CHOICES <- stats::setNames(
  AVAILABLE_TRAITS,
  vapply(AVAILABLE_TRAITS, display_trait_name, character(1))
)

DEFAULT_INDIVIDUAL_TRAIT <- if (DEFAULT_INDIVIDUAL_TRAIT %in% AVAILABLE_TRAITS) {
  DEFAULT_INDIVIDUAL_TRAIT
} else {
  AVAILABLE_TRAITS[[1]]
}
DEFAULT_INDIVIDUAL_GROUP <- trait_group(DEFAULT_INDIVIDUAL_TRAIT)

traits_for_group <- function(group_name) {
  if (identical(group_name, "__all__")) {
    return(AVAILABLE_TRAITS)
  }

  TRAIT_GROUPS$Phenotype[TRAIT_GROUPS$Group == group_name]
}

resource_paths <- shiny::resourcePaths()
current_fig_path <- unname(resource_paths["amino-figs"])

if (length(current_fig_path) == 0 || !identical(current_fig_path[[1]], file.path(APP_HOME, "Figs"))) {
  shiny::addResourcePath("amino-figs", file.path(APP_HOME, "Figs"))
}

################################################################################
### UI
################################################################################

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        background: #f5f6f8;
        color: #1f1f1f;
      }
      .app-title {
        font-size: 28px;
        font-weight: 700;
        margin-bottom: 6px;
      }
      .app-subtitle {
        font-size: 14px;
        color: #4d4d4d;
        margin-bottom: 20px;
      }
      .panel-card {
        background: white;
        border: 1px solid #d9d9d9;
        border-radius: 12px;
        padding: 16px 18px;
        box-shadow: 0 4px 14px rgba(0, 0, 0, 0.04);
      }
      .top-panel-card {
        height: 100%;
      }
      .card-title {
        font-size: 18px;
        font-weight: 700;
        margin-bottom: 12px;
      }
      .summary-strip {
        font-size: 13px;
        color: #444444;
        margin-bottom: 10px;
      }
      .left-controls {
        position: sticky;
        top: 16px;
      }
      .control-note {
        font-size: 12px;
        color: #555555;
        line-height: 1.45;
        margin-top: 10px;
      }
      .plot-image {
        width: 100%;
        height: 320px;
        max-height: 320px;
        object-fit: contain;
        display: block;
        background: white;
        border: 1px solid #ebebeb;
        border-radius: 10px;
      }
      .tab-content {
        padding-top: 18px;
      }
      .nav-tabs {
        border-bottom: 1px solid #d9d9d9;
      }
      .nav-tabs > li > a {
        font-weight: 600;
        color: #38506a;
      }
      .shiny-output-error-validation {
        color: #7a2d2d;
        font-weight: 600;
      }
      @media (max-width: 991px) {
        .left-controls {
          position: static;
        }
      }
    "))
  ),
  fluidRow(
    column(
      width = 12,
      div(class = "app-title", "Amino-Acid GWAS Explorer"),
      div(
        class = "app-subtitle",
        "Use the individual tab for one GWAS result at a time, or switch to the combined tab to count recurring gene hits across all results."
      )
    )
  ),
  tabsetPanel(
    id = "main_tabs",
    tabPanel(
      "Individual GWAS",
      fluidRow(
        column(
          width = 3,
          div(
            class = "panel-card left-controls",
            div(class = "card-title", "GWAS Controls"),
            selectInput(
              inputId = "individual_group",
              label = "Select amino acid",
              choices = GROUP_CHOICES,
              selected = DEFAULT_INDIVIDUAL_GROUP
            ),
            selectInput(
              inputId = "individual_trait",
              label = "Select GWAS result",
              choices = TRAIT_CHOICES,
              selected = DEFAULT_INDIVIDUAL_TRAIT
            ),
            selectInput(
              inputId = "individual_cutoff",
              label = HTML("Select <code>-log10(p)</code> threshold"),
              choices = c("5" = 5, "6" = 6, "7" = 7),
              selected = 7
            ),
            div(
              class = "control-note",
              "Pick one amino-acid family, then one specific GWAS result. The figures on the right update for that single result, and the table below shows the matching gene annotations."
            )
          )
        ),
        column(
          width = 9,
          fluidRow(
            column(
              width = 5,
              div(
                class = "panel-card top-panel-card",
                div(class = "card-title", "Manhattan Plot"),
                uiOutput("individual_manhattan_plot")
              )
            ),
            column(
              width = 3,
              div(
                class = "panel-card top-panel-card",
                div(class = "card-title", "QQ Plot"),
                uiOutput("individual_qq_plot")
              )
            ),
            column(
              width = 4,
              div(
                class = "panel-card top-panel-card",
                div(class = "card-title", "FAA Distribution"),
                plotOutput("individual_faa_distribution", height = "320px")
              )
            )
          ),
          tags$div(style = "height: 16px;"),
          div(
            class = "panel-card",
            div(class = "card-title", "GWAS Annotation Table"),
            uiOutput("individual_summary_strip"),
            DTOutput("individual_table")
          )
        )
      )
    ),
    tabPanel(
      "Combined Gene Hits",
      fluidRow(
        column(
          width = 3,
          div(
            class = "panel-card left-controls",
            div(class = "card-title", "Gene Hit Filters"),
            selectInput(
              inputId = "combined_scope",
              label = "Combine results from",
              choices = c("All GWAS results" = "__all__", GROUP_CHOICES),
              selected = DEFAULT_COMBINED_SCOPE
            ),
            selectInput(
              inputId = "combined_cutoff",
              label = HTML("Select <code>-log10(p)</code> threshold"),
              choices = c("5" = 5, "6" = 6, "7" = 7),
              selected = 7
            ),
            div(
              class = "control-note",
              "This tab counts how often genes recur across the chosen GWAS set. Use 'All GWAS results' for the full combined view, or switch to one amino-acid group to narrow it down."
            )
          )
        ),
        column(
          width = 9,
          div(
            class = "panel-card",
            div(class = "card-title", "Gene and Phenotype Summaries"),
            uiOutput("combined_summary_strip"),
            tabsetPanel(
              id = "combined_summary_tabs",
              tabPanel("Top Genes", DTOutput("combined_gene_table")),
              tabPanel("Phenotypes", DTOutput("combined_phenotype_table"))
            )
          )
        )
      )
    )
  )
)

################################################################################
### SERVER
################################################################################

server <- function(input, output, session) {
  observe({
    group_traits <- traits_for_group(input$individual_group)
    selected_trait <- isolate(input$individual_trait)

    if (!selected_trait %in% group_traits) {
      selected_trait <- group_traits[[1]]
    }

    updateSelectInput(
      session,
      "individual_trait",
      choices = stats::setNames(group_traits, vapply(group_traits, display_trait_name, character(1))),
      selected = selected_trait
    )
  })

  selected_combined_traits <- reactive({
    req(input$combined_scope)
    traits_for_group(input$combined_scope)
  })

  individual_hits <- reactive({
    req(input$individual_trait)
    build_individual_hits(
      best_hits_df = best_hits_raw,
      trait_name = input$individual_trait,
      cutoff = as.numeric(input$individual_cutoff),
      annotation_lookup = annotation_lookup
    )
  })

  individual_faa_distribution <- reactive({
    req(input$individual_trait)
    build_faa_distribution(
      faa_df = faa_raw,
      trait_name = input$individual_trait
    )
  })

  combined_gene_hits <- reactive({
    build_combined_gene_hits(
      best_hits_df = best_hits_raw,
      selected_traits = selected_combined_traits(),
      cutoff = as.numeric(input$combined_cutoff),
      annotation_lookup = annotation_lookup
    )
  })

  combined_phenotype_hits <- reactive({
    build_phenotype_summary(
      best_hits_df = best_hits_raw,
      selected_traits = selected_combined_traits(),
      cutoff = as.numeric(input$combined_cutoff),
      annotation_lookup = annotation_lookup
    )
  })

  output$individual_manhattan_plot <- renderUI({
    req(input$individual_trait)

    plot_src <- build_plot_src("manhattan", input$individual_trait)
    validate(
      need(!is.null(plot_src), "No Manhattan plot file was found for the selected GWAS result.")
    )

    tags$img(
      src = plot_src,
      alt = paste("Manhattan plot for", display_trait_name(input$individual_trait)),
      class = "plot-image"
    )
  })

  output$individual_qq_plot <- renderUI({
    req(input$individual_trait)

    plot_src <- build_plot_src("qq", input$individual_trait)
    validate(
      need(!is.null(plot_src), "No QQ plot file was found for the selected GWAS result.")
    )

    tags$img(
      src = plot_src,
      alt = paste("QQ plot for", display_trait_name(input$individual_trait)),
      class = "plot-image"
    )
  })

  output$individual_faa_distribution <- renderPlot({
    req(input$individual_trait)

    distribution_df <- individual_faa_distribution()
    validate(
      need(
        nrow(distribution_df) > 0,
        paste0(
          "No FAA distribution is available for ",
          display_trait_name(input$individual_trait),
          " in SuppTable_FAA.csv."
        )
      )
    )

    make_faa_distribution_plot(
      distribution_df = distribution_df,
      trait_name = input$individual_trait
    )
  }, res = 120)

  output$individual_summary_strip <- renderUI({
    df <- individual_hits()

    tags$div(
      class = "summary-strip",
      tags$strong(display_trait_name(input$individual_trait)),
      " | ",
      tags$strong(dplyr::n_distinct(df$GeneID)),
      " genes meeting the current cutoff | ",
      tags$strong(nrow(df)),
      " annotation rows"
    )
  })

  output$combined_summary_strip <- renderUI({
    selected_traits <- selected_combined_traits()
    gene_df <- combined_gene_hits()
    phenotype_df <- combined_phenotype_hits()

    tags$div(
      class = "summary-strip",
      tags$strong(length(selected_traits)),
      " GWAS results combined | ",
      tags$strong(nrow(gene_df)),
      " recurring genes | ",
      tags$strong(nrow(phenotype_df)),
      " phenotypes represented"
    )
  })

  output$individual_table <- renderDT({
    summary_df <- individual_hits() %>%
      dplyr::rename(
        `Gene ID` = GeneID,
        `Gene symbol` = GeneSymbol,
        Phenotype = Phenotype,
        Model = Model,
        SNP = SNP,
        Chromosome = Chr,
        `SNP position` = SNP_Pos,
        `p-value` = P.value,
        `-log10(p)` = log10_P,
        Relation = Relation,
        `Distance to gene (bp)` = Distance_to_Gene_bp,
        `Family / subfamily` = Family_Subfamily,
        `Protein class` = Protein_Class,
        `GO MF` = GO_MF,
        `GO BP` = GO_BO,
        `GO CC` = GO_CC
      )

    DT::datatable(
      summary_df,
      rownames = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 12,
        lengthMenu = c(12, 25, 50, 100),
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv"),
        autoWidth = TRUE
      )
    )
  })

  output$combined_gene_table <- renderDT({
    summary_df <- combined_gene_hits() %>%
      dplyr::rename(
        `Gene ID` = GeneID,
        `Gene symbol` = GeneSymbol,
        Chromosome = GeneChr,
        `Gene start` = GeneStart,
        `Gene end` = GeneEnd,
        Phenotypes = Phenotypes,
        Occurrences = Occurrences,
        `Highest -log10(p)` = Highest_log10P,
        `Top p-value` = Top_Pvalue,
        `Top model` = Top_Model,
        `Top SNP` = Top_SNP,
        `Top relation` = Top_Relation,
        `Top distance (bp)` = Top_Distance_bp,
        `Family / subfamily` = Family_Subfamily,
        `Protein class` = Protein_Class,
        `GO MF` = GO_MF,
        `GO BP` = GO_BO,
        `GO CC` = GO_CC
      )

    DT::datatable(
      summary_df,
      rownames = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 12,
        lengthMenu = c(12, 25, 50, 100),
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv"),
        autoWidth = TRUE
      )
    )
  })

  output$combined_phenotype_table <- renderDT({
    summary_df <- combined_phenotype_hits() %>%
      dplyr::mutate(
        Phenotype = vapply(Phenotype, display_trait_name, character(1))
      ) %>%
      dplyr::rename(
        `GWAS result` = Phenotype,
        `Gene count` = Gene_count,
        `Annotation rows` = Hit_rows,
        `Highest -log10(p)` = Highest_log10P,
        `Top gene` = Top_Gene,
        Models = Models
      )

    DT::datatable(
      summary_df,
      rownames = FALSE,
      filter = "top",
      extensions = "Buttons",
      options = list(
        pageLength = 12,
        lengthMenu = c(12, 25, 50, 100),
        scrollX = TRUE,
        dom = "Bfrtip",
        buttons = c("copy", "csv"),
        autoWidth = TRUE
      )
    )
  })
}

################################################################################
### APP OBJECT
################################################################################

app <- shinyApp(ui = ui, server = server)
app
