################################################################################
### AMINO-ACID GWAS EXPLORER
################################################################################

library(shiny)
library(DT)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(IRanges)

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

AMINO_DIR <- "/Users/nirwantandukar/Documents/Research/results/GWAS/Sarah_amino_acid/N_grain/Phenotypes_GWAS_Grain"
REFERENCE_GFF <- "/Users/nirwantandukar/Library/Mobile Documents/com~apple~CloudDocs/Research/Data/Maize/Maize.annotation/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
AVAILABLE_TRAITS <- c("D", "E", "N", "P", "Q", "Total_N", "Total_PBAA")
DEFAULT_TRAITS <- AVAILABLE_TRAITS
DEFAULT_PLOT_TRAIT <- "P"

TRAIT_LABELS <- c(
  D = "Aspartate",
  E = "Glutamate",
  N = "Asparagine",
  P = "Proline",
  Q = "Glutamine",
  Total_N = "Total nitrogen",
  Total_PBAA = "Total protein-bound amino acids"
)

WINDOW_BP <- 25000
MAX_BG_POINTS <- 150000
MODEL_COLORS <- c(
  MLM = "#0072B2",
  MLMM = "#E69F00",
  BLINK = "#009E73",
  FarmCPU = "#D55E00"
)

EMPTY_GENE_TABLE <- tibble(
  GeneID = character(),
  GeneSymbol = character(),
  Phenotypes = character(),
  Phenotypes_count = integer(),
  pvalues = character()
)

################################################################################
### THEME AND HELPERS
################################################################################

plot_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey35"),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.margin = margin(12, 12, 12, 12)
  )

trait_choice_labels <- stats::setNames(
  AVAILABLE_TRAITS,
  paste0(AVAILABLE_TRAITS, " \u2014 ", TRAIT_LABELS[AVAILABLE_TRAITS])
)

pretty_trait_name <- function(trait_name) {
  if (trait_name %in% names(TRAIT_LABELS)) {
    paste0(trait_name, " \u2014 ", TRAIT_LABELS[[trait_name]])
  } else {
    trait_name
  }
}

thin_rows <- function(df, max_n) {
  if (is.null(df) || nrow(df) <= max_n) {
    return(df)
  }

  keep_idx <- unique(round(seq(1, nrow(df), length.out = max_n)))
  df[keep_idx, , drop = FALSE]
}

get_gene_label <- function(genes_gr) {
  gene_id <- as.character(mcols(genes_gr)$ID)

  if ("Name" %in% colnames(mcols(genes_gr))) {
    gene_name <- as.character(mcols(genes_gr)$Name)
    ifelse(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  } else {
    gene_id
  }
}

annotate_hits <- function(df, phenotype_name, model_name, genes_gr, window_bp = 25000) {
  if (is.null(df) || nrow(df) == 0) {
    return(NULL)
  }

  snps <- GRanges(
    seqnames = Rle(paste0("chr", df$Chr)),
    ranges = IRanges(df$Pos, df$Pos)
  )
  mcols(snps)$SNP <- df$SNP
  mcols(snps)$P.value <- df$P.value
  mcols(snps)$SNP_Pos <- df$Pos

  extended <- snps
  start(extended) <- pmax(1L, start(snps) - window_bp)
  end(extended) <- end(snps) + window_bp

  overlaps <- findOverlaps(genes_gr, extended, ignore.strand = TRUE)
  if (length(overlaps) == 0) {
    return(NULL)
  }

  gene_idx <- queryHits(overlaps)
  snp_idx <- subjectHits(overlaps)
  gene_label <- get_gene_label(genes_gr)
  snp_pos <- mcols(extended)$SNP_Pos[snp_idx]

  tibble(
    Phenotype = phenotype_name,
    Model = model_name,
    GeneID = as.character(mcols(genes_gr)$ID[gene_idx]),
    GeneSymbol = gene_label[gene_idx],
    SNP = mcols(extended)$SNP[snp_idx],
    SNP_Pos = snp_pos,
    P.value = mcols(extended)$P.value[snp_idx],
    log10_P = -log10(mcols(extended)$P.value[snp_idx])
  )
}

collapse_best_per_phenotype <- function(annotation_df) {
  if (is.null(annotation_df) || nrow(annotation_df) == 0) {
    return(NULL)
  }

  annotation_df %>%
    dplyr::group_by(GeneID, GeneSymbol, Phenotype) %>%
    dplyr::slice_min(P.value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(GeneID, GeneSymbol, Phenotype, log10_P)
}

make_gene_summary <- function(best_hits_df, selected_traits) {
  if (is.null(best_hits_df) || nrow(best_hits_df) == 0) {
    return(EMPTY_GENE_TABLE)
  }

  best_hits_df %>%
    dplyr::mutate(Phenotype = factor(Phenotype, levels = selected_traits)) %>%
    dplyr::arrange(GeneSymbol, Phenotype) %>%
    dplyr::group_by(GeneID, GeneSymbol) %>%
    dplyr::summarise(
      Phenotypes = paste(as.character(Phenotype), collapse = ";"),
      Phenotypes_count = dplyr::n(),
      pvalues = paste(formatC(log10_P, digits = 2, format = "f"), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::select(
      GeneID,
      GeneSymbol,
      Phenotypes,
      Phenotypes_count,
      pvalues
    ) %>%
    dplyr::arrange(dplyr::desc(Phenotypes_count), GeneSymbol)
}

load_trait_raw <- function(trait_name, cache_env) {
  cache_key <- paste0("raw__", trait_name)
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  file_path <- file.path(AMINO_DIR, paste0(trait_name, ".csv"))
  raw_dt <- data.table::fread(
    file_path,
    select = c("SNP", "Chr", "Pos", "P.value", "model"),
    showProgress = FALSE
  )

  raw_df <- raw_dt %>%
    as_tibble() %>%
    dplyr::mutate(
      Chr = as.integer(Chr),
      Pos = as.numeric(Pos),
      P.value = as.numeric(P.value),
      model = as.character(model),
      log10_P = -log10(P.value)
    ) %>%
    dplyr::filter(
      !is.na(Chr),
      !is.na(Pos),
      !is.na(P.value),
      is.finite(P.value),
      P.value > 0,
      is.finite(log10_P)
    )

  assign(cache_key, raw_df, envir = cache_env)
  raw_df
}

build_plot_data <- function(trait_name, cache_env) {
  cache_key <- paste0("plot__", trait_name)
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  raw_df <- load_trait_raw(trait_name, cache_env)
  best_df <- raw_df %>%
    dplyr::group_by(SNP, Chr, Pos) %>%
    dplyr::slice_min(P.value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::rename(best_model = model)

  chr_tbl <- best_df %>%
    dplyr::group_by(Chr) %>%
    dplyr::summarise(chr_len = max(Pos, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Chr) %>%
    dplyr::mutate(
      offset = dplyr::lag(cumsum(chr_len), default = 0),
      center = offset + chr_len / 2
    )

  plot_df <- best_df %>%
    dplyr::left_join(chr_tbl %>% dplyr::select(Chr, offset), by = "Chr") %>%
    dplyr::mutate(pos_cum = Pos + offset) %>%
    dplyr::arrange(Chr, Pos)

  bg_df <- plot_df %>%
    dplyr::filter(log10_P < 5)
  sig_df <- plot_df %>%
    dplyr::filter(log10_P >= 5)

  plot_data <- list(
    bg_df = thin_rows(bg_df, MAX_BG_POINTS),
    sig_df = sig_df,
    chr_tbl = chr_tbl
  )

  assign(cache_key, plot_data, envir = cache_env)
  plot_data
}

build_trait_gene_hits <- function(trait_name, cutoff, genes_gr, cache_env) {
  cache_key <- paste0("hits__", trait_name, "__", cutoff)
  if (exists(cache_key, envir = cache_env, inherits = FALSE)) {
    return(get(cache_key, envir = cache_env, inherits = FALSE))
  }

  raw_df <- load_trait_raw(trait_name, cache_env) %>%
    dplyr::filter(log10_P >= cutoff)

  if (nrow(raw_df) == 0) {
    assign(cache_key, NULL, envir = cache_env)
    return(NULL)
  }

  annotation_list <- lapply(unique(raw_df$model), function(model_name) {
    model_df <- raw_df %>%
      dplyr::filter(model == model_name) %>%
      dplyr::select(SNP, Chr, Pos, P.value)

    annotate_hits(
      df = model_df,
      phenotype_name = trait_name,
      model_name = model_name,
      genes_gr = genes_gr,
      window_bp = WINDOW_BP
    )
  })

  best_hits <- collapse_best_per_phenotype(dplyr::bind_rows(annotation_list))
  assign(cache_key, best_hits, envir = cache_env)
  best_hits
}

make_manhattan_plot <- function(trait_name, cutoff, plot_data) {
  trait_title <- pretty_trait_name(trait_name)

  ggplot() +
    geom_point(
      data = plot_data$bg_df,
      aes(x = pos_cum, y = log10_P, color = factor(Chr %% 2)),
      alpha = 0.4,
      size = 0.45,
      show.legend = FALSE
    ) +
    scale_color_manual(values = c("0" = "grey65", "1" = "grey40")) +
    geom_point(
      data = plot_data$sig_df,
      aes(x = pos_cum, y = log10_P, color = best_model),
      alpha = 0.85,
      size = 1.15
    ) +
    scale_x_continuous(
      breaks = plot_data$chr_tbl$center,
      labels = plot_data$chr_tbl$Chr,
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_color_manual(values = MODEL_COLORS, breaks = names(MODEL_COLORS)) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "black", linewidth = 0.45) +
    geom_hline(yintercept = 6, linetype = "dashed", color = "grey30", linewidth = 0.45) +
    geom_hline(yintercept = 7, linetype = "solid", color = "black", linewidth = 0.55) +
    annotate(
      "label",
      x = Inf,
      y = cutoff,
      label = paste0("-log10(p) = ", cutoff),
      hjust = 1.05,
      vjust = -0.3,
      size = 4,
      label.size = NA,
      fill = "white"
    ) +
    labs(
      title = trait_title,
      subtitle = "Use the phenotype selector on the right to flip through selected traits",
      x = "Chromosome",
      y = expression(-log[10](italic(p))),
      color = NULL
    ) +
    plot_theme
}

################################################################################
### STATIC DATA
################################################################################

genes_gr <- rtracklayer::import(REFERENCE_GFF)
genes_gr <- genes_gr[mcols(genes_gr)$type == "gene"]
trait_cache <- new.env(parent = emptyenv())

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
      .right-controls {
        position: sticky;
        top: 16px;
      }
      .btn-blockish {
        width: 100%;
        margin-bottom: 8px;
      }
      .control-note {
        font-size: 12px;
        color: #555555;
        line-height: 1.45;
        margin-top: 10px;
      }
      .shiny-output-error-validation {
        color: #7a2d2d;
        font-weight: 600;
      }
    "))
  ),
  fluidRow(
    column(
      width = 12,
      div(class = "app-title", "Amino-Acid GWAS Explorer"),
      div(
        class = "app-subtitle",
        "Browse current amino-acid GWAS phenotypes, flip through Manhattan plots, and collapse candidate genes across one or more selected traits."
      )
    )
  ),
  fluidRow(
    column(
      width = 8,
      div(
        class = "panel-card",
        div(class = "card-title", "Manhattan Plot"),
        plotOutput("manhattan_plot", height = "460px")
      ),
      tags$div(style = "height: 16px;"),
      div(
        class = "panel-card",
        div(class = "card-title", "Gene Summary"),
        uiOutput("summary_strip"),
        DTOutput("gene_table")
      )
    ),
    column(
      width = 4,
      div(
        class = "panel-card right-controls",
        div(class = "card-title", "GWAS Options"),
        selectizeInput(
          inputId = "traits",
          label = "Select one or more phenotypes",
          choices = trait_choice_labels,
          selected = DEFAULT_TRAITS,
          multiple = TRUE,
          options = list(placeholder = "Choose amino-acid traits")
        ),
        fluidRow(
          column(
            width = 6,
            actionButton("select_all", "All current", class = "btn-blockish")
          ),
          column(
            width = 6,
            actionButton("clear_traits", "Clear", class = "btn-blockish")
          )
        ),
        radioButtons(
          inputId = "summary_scope",
          label = "Lower-table summary",
          choices = c(
            "Selected phenotypes (sum some)" = "selected",
            "All current phenotypes (sum all)" = "all"
          ),
          selected = "selected"
        ),
        selectInput(
          inputId = "cutoff",
          label = HTML("Minimum <code>-log10(p)</code>"),
          choices = c("5" = 5, "6" = 6, "7" = 7),
          selected = 7
        ),
        selectInput(
          inputId = "plot_trait",
          label = "Trait shown in the Manhattan plot",
          choices = trait_choice_labels,
          selected = DEFAULT_PLOT_TRAIT
        ),
        div(
          class = "control-note",
          "If you select multiple phenotypes, the top panel stays on one phenotype at a time while the lower table collapses genes across your chosen set. We can add the ratio phenotypes later once you decide which ones you want in here."
        )
      )
    )
  )
)

################################################################################
### SERVER
################################################################################

server <- function(input, output, session) {
  observeEvent(input$select_all, {
    updateSelectizeInput(session, "traits", selected = AVAILABLE_TRAITS)
  })

  observeEvent(input$clear_traits, {
    updateSelectizeInput(session, "traits", selected = character(0))
  })

  selected_traits <- reactive({
    AVAILABLE_TRAITS[AVAILABLE_TRAITS %in% input$traits]
  })

  table_traits <- reactive({
    if (identical(input$summary_scope, "all")) {
      AVAILABLE_TRAITS
    } else {
      selected_traits()
    }
  })

  observe({
    plot_choices <- selected_traits()
    if (length(plot_choices) == 0) {
      plot_choices <- AVAILABLE_TRAITS
    }

    selected_plot <- isolate(input$plot_trait)
    if (!selected_plot %in% plot_choices) {
      if (DEFAULT_PLOT_TRAIT %in% plot_choices) {
        selected_plot <- DEFAULT_PLOT_TRAIT
      } else {
        selected_plot <- plot_choices[[1]]
      }
    }

    updateSelectInput(
      session,
      "plot_trait",
      choices = stats::setNames(plot_choices, paste0(plot_choices, " \u2014 ", TRAIT_LABELS[plot_choices])),
      selected = selected_plot
    )
  })

  plot_data <- reactive({
    req(input$plot_trait)
    build_plot_data(input$plot_trait, trait_cache)
  })

  gene_summary <- reactive({
    traits_to_use <- table_traits()

    validate(
      need(
        length(traits_to_use) > 0,
        "Select at least one phenotype on the right, or switch the summary to 'All current phenotypes'."
      )
    )

    cutoff_value <- as.numeric(input$cutoff)

    withProgress(message = "Building gene summary", value = 0, {
      hit_list <- vector("list", length(traits_to_use))

      for (i in seq_along(traits_to_use)) {
        incProgress(1 / max(length(traits_to_use), 1), detail = traits_to_use[[i]])
        hit_list[[i]] <- build_trait_gene_hits(
          trait_name = traits_to_use[[i]],
          cutoff = cutoff_value,
          genes_gr = genes_gr,
          cache_env = trait_cache
        )
      }

      combined_hits <- dplyr::bind_rows(hit_list)
      make_gene_summary(combined_hits, traits_to_use)
    })
  })

  output$manhattan_plot <- renderPlot({
    req(input$plot_trait)
    make_manhattan_plot(
      trait_name = input$plot_trait,
      cutoff = as.numeric(input$cutoff),
      plot_data = plot_data()
    )
  }, res = 120)

  output$summary_strip <- renderUI({
    summary_df <- gene_summary()
    traits_to_use <- table_traits()

    tags$div(
      class = "summary-strip",
      tags$strong(length(traits_to_use)),
      " phenotypes in summary | ",
      tags$strong(nrow(summary_df)),
      " genes meeting the current cutoff | annotation window = ",
      WINDOW_BP / 1000,
      " kb each side"
    )
  })

  output$gene_table <- renderDT({
    summary_df <- gene_summary() %>%
      dplyr::rename(
        `Gene ID` = GeneID,
        `Gene symbol` = GeneSymbol,
        `Phenotypes` = Phenotypes,
        `Phenotype count` = Phenotypes_count,
        `-log10(p)` = pvalues
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
