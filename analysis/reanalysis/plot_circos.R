# =============================================================================
# Circos Plot Visualisation of DGSA Results
# =============================================================================
# Produces circos plots comparing QuSAGE and dearseq results across timepoints,
# for both positively and negatively regulated gene sets.
# =============================================================================

# ── Packages ──────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(fs)
  library(tidyverse)
  library(circlize)
})

# ── Paths ─────────────────────────────────────────────────────────────────────

p_results_dearseq <- fs::path("output", "results", "reanalysis", "dearseq_dgsa_results_processed.rds")
p_results_qusage  <- fs::path("output", "results", "reanalysis",  "qusage_dgsa_results_processed.rds")
figures_folder    <- fs::path("output", "figures", "reanalysis")

# ── Load and combine results ───────────────────────────────────────────────────

results_df <- bind_rows(
  readRDS(p_results_dearseq),
  readRDS(p_results_qusage)
)

# =============================================================================
# CONFIGURATION
# =============================================================================

# Gene set aggregates to include in all plots
btm_aggregates <- c(
  "Antigen Presentation",
  "Inflammatory/TLR/Chemokines",
  "Interferon/Antiviral Sensing",
  "Monocytes",
  "DC Activation",
  "Neutrophils",
  "NK Cells",
  "Signal Transduction",
  "ECM And Migration",
  "Energy Metabolism",
  "Cell Cycle",
  "Platelets",
  "T Cells",
  "B Cells",
  "Plasma Cells"
)

# Arguments shared across every plot_circos call
circos_defaults <- list(
  conditions         = levels(results_df$condition),
  aggregates_name    = btm_aggregates,
  p_correction       = "BH",
  p_approach         = "global",
  p_threshold        = 0.05,
  filter_variable    = "none",
  filter_mode        = "user",
  user_threshold     = 0.5,
  quantile_threshold = 0.5,
  scores             = "fc.score",
  correlation        = "corr.mean",
  ring               = "expression",
  order              = "set_all",
  quantile_scoreclip = 0.995,
  legend             = FALSE,
  title_size         = 6,
  title_line         = -4
)

# =============================================================================
# PLOT FUNCTIONS
# =============================================================================

# Render a circos plot for one method x timepoint x arc combination.
# All arguments not listed here are drawn from circos_defaults.
plot_circos <- function(method_name,
                        conditions,
                        timepoint,
                        aggregates_name,
                        p_correction       = c("BH", "bonferroni", "holm", "hommel", "hochberg", "BY"),
                        p_approach         = c("global", "withinTime", "withinComparison"),
                        p_threshold        = 0.05,
                        filter_variable    = c("none", "fc.score", "activation.score"),
                        filter_mode        = c("user", "data"),
                        user_threshold     = 0.5,
                        quantile_threshold = 0.5,
                        scores             = c("fc.score", "activation.score"),
                        correlation        = c("corr.mean", "mean.corr"),
                        arc                = c("any", "direction", "positive", "negative"),
                        ring               = c("all", "expression", "none"),
                        order              = c("set_all", "set_available", "cluster"),
                        quantile_scoreclip = 0.995,
                        legend             = TRUE,
                        plot_title         = NULL,
                        title_size         = 2.5,
                        title_line         = 0) {
  
  # Subset to selected method, conditions, and aggregates
  results_df_circos <- results_df %>%
    filter(
      gs.aggregate %in% aggregates_name,
      condition    %in% conditions,
      method       == method_name
    ) %>%
    mutate(gs.aggregate = droplevels(gs.aggregate))
  
  # Global gene set ordering within aggregates
  global_order <- results_df_circos %>%
    select(gs.name, gs.aggregate) %>%
    distinct() %>%
    arrange(gs.aggregate) %>%
    pull(gs.name)
  
  results_df_circos <- results_df_circos %>%
    mutate(gs.name = factor(gs.name, levels = global_order))
  
  # Determine which conditions to show
  if (order %in% c("set_available", "cluster")) {
    condition_names <- results_df_circos %>%
      filter(time == timepoint) %>%
      pull(condition) %>%
      unique()
  } else {
    condition_names <- unique(results_df$condition)
  }
  
  condition_order <- levels(results_df$condition)
  condition_order <- condition_order[condition_order %in% condition_names]
  
  # Build condition metadata (ordering, colour, availability)
  if (order %in% c("set_all", "set_available")) {
    condition_circos_metadata <- data.frame(
      order     = match(condition_names, condition_order),
      condition = factor(condition_names, levels = condition_order),
      stringsAsFactors = FALSE
    ) %>%
      arrange(condition)
    
  } else if (order == "cluster") {
    mat_df <- results_df_circos %>%
      filter(time == timepoint) %>%
      select(condition, gs.name, score = .data[[scores]]) %>%
      pivot_wider(names_from = gs.name, values_from = score) %>%
      mutate(across(-condition, ~ tidyr::replace_na(.x, 0)))
    
    mat <- mat_df %>%
      column_to_rownames("condition") %>%
      as.matrix()
    
    hc                <- hclust(dist(mat, method = "euclidean"), method = "ward.D2")
    ordered_conditions <- hc$labels[hc$order]
    
    condition_circos_metadata <- data.frame(
      order     = match(condition_names, ordered_conditions),
      condition = factor(condition_names, levels = ordered_conditions),
      stringsAsFactors = FALSE
    ) %>%
      arrange(condition)
  }
  
  # Per-timepoint availability columns (expression + immune response)
  unique_times    <- unique(results_df_circos$time)
  condition_colors <- results_df %>%
    select(condition, condition.colour) %>%
    distinct() %>%
    deframe()
  
  for (t in unique_times) {
    conditions_expr_at_time <- results_df_circos %>%
      filter(time == t) %>%
      pull(condition) %>%
      unique()
    
    conditions_resp_at_time <- results_df %>%
      filter(time == t, !is.na(corr.mean)) %>%
      pull(condition) %>%
      unique()
    
    expr_col  <- paste0("condition_expr_available_", t)
    resp_col  <- paste0("condition_response_available_", t)
    color_col <- paste0("condition_color_", t)
    text_col  <- paste0("text_color_", t)
    
    condition_circos_metadata <- condition_circos_metadata %>%
      mutate(
        !!expr_col  := condition %in% conditions_expr_at_time,
        !!resp_col  := condition %in% conditions_resp_at_time,
        !!color_col := if_else(.data[[expr_col]], condition_colors[as.character(condition)], "#D3D3D3"),
        !!text_col  := if_else(.data[[expr_col]], "black", "#9999a1")
      )
  }
  
  condition_circos_metadata$xmin <- 0
  condition_circos_metadata$xmax <- length(global_order)
  
  # Fixed gene set x-positions
  gene_set_positions <- setNames(
    seq(0.5, length(global_order) - 0.5, length.out = length(global_order)),
    global_order
  )
  
  # Aggregate colour and proportion metadata (for the aggregate ring)
  aggregate_colors <- results_df %>%
    select(gs.aggregate, gs.colour) %>%
    distinct() %>%
    { setNames(.$gs.colour, .$gs.aggregate) }
  
  aggregate_distinct <- results_df_circos %>%
    select(gs.name, gs.aggregate) %>%
    distinct()
  
  pt <- prop.table(table(aggregate_distinct$gs.aggregate))
  
  aggregate_proportions <- data.frame(
    gs.aggregate = names(pt),
    prop         = as.vector(pt),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      colour    = aggregate_colors[gs.aggregate],
      cum_start = c(0, head(cumsum(prop), -1)),
      cum_end   = cumsum(prop)
    )
  
  # ── Significance filtering ────────────────────────────────────────────────
  
  results_df_significant <- if (p_correction == "none") {
    results_df_circos %>% filter(rawPval < p_threshold, time == timepoint)
  } else {
    results_df_circos %>%
      filter(
        .data[[paste0(p_approach, ".adjPval_", p_correction)]] < p_threshold,
        time == timepoint
      )
  }
  
  if (filter_variable != "none" && filter_mode == "data") {
    if (p_approach == "global") {
      threshold <- results_df_circos %>%
        pull(filter_variable) %>% abs() %>% quantile(quantile_threshold)
      results_df_significant <- results_df_significant %>%
        filter(abs(.data[[filter_variable]]) > threshold)
      
    } else if (p_approach == "withinTime") {
      threshold <- results_df_circos %>%
        filter(time == timepoint) %>%
        pull(filter_variable) %>% abs() %>% quantile(quantile_threshold)
      results_df_significant <- results_df_significant %>%
        filter(abs(.data[[filter_variable]]) > threshold)
      
    } else if (p_approach == "withinCondition") {
      results_df_significant <- results_df_circos %>%
        filter(time == timepoint) %>%
        group_by(condition) %>%
        mutate(abs_score = abs(.data[[filter_variable]])) %>%
        filter(abs_score >= quantile(abs_score, probs = quantile_threshold, na.rm = TRUE)) %>%
        ungroup()
    }
  } else if (filter_variable != "none" && filter_mode == "user") {
    results_df_significant <- results_df_significant %>%
      filter(abs(.data[[filter_variable]]) > user_threshold)
  }
  
  # ── Arc (link) computation ────────────────────────────────────────────────
  
  # Helper: build shared_links for arcs based on a pre-filtered significant set
  build_links_from_modules <- function(df_sig) {
    modules <- df_sig %>%
      group_by(gs.name) %>%
      filter(n() > 1) %>%
      summarise(condition_names = list(unique(condition)), .groups = "drop")
    
    if (nrow(modules) == 0) return(NULL)
    
    modules %>%
      mutate(pairs = purrr::map(condition_names, ~ combn(.x, 2, simplify = FALSE))) %>%
      select(gs.name, pairs) %>%
      unnest(pairs) %>%
      unnest_wider(pairs, names_sep = "_") %>%
      rename(condition1 = pairs_1, condition2 = pairs_2) %>%
      mutate(position = gene_set_positions[gs.name])
  }
  
  shared_links <- if (arc == "any") {
    build_links_from_modules(results_df_significant)
    
  } else if (arc == "positive") {
    build_links_from_modules(filter(results_df_significant, .data[[scores]] >= 0))
    
  } else if (arc == "negative") {
    build_links_from_modules(filter(results_df_significant, .data[[scores]] <= 0))
    
  } else if (arc == "direction") {
    col1 <- paste0(scores, 1)
    col2 <- paste0(scores, 2)
    
    results_df_significant %>%
      select(gs.name, condition, .data[[scores]]) %>%
      inner_join(
        results_df_significant %>% select(gs.name, condition, .data[[scores]]),
        by = "gs.name", suffix = c("1", "2"), relationship = "many-to-many"
      ) %>%
      mutate(!!col1 := as.numeric(.data[[col1]]),
             !!col2 := as.numeric(.data[[col2]])) %>%
      filter(condition1 != condition2,
             .data[[col1]] * .data[[col2]] > 0) %>%
      transmute(gs.name, condition1, condition2,
                position = gene_set_positions[gs.name]) %>%
      distinct() %>%
      { if (nrow(.) == 0) NULL else . }
  }
  
  any_links <- !is.null(shared_links) && nrow(shared_links) > 0
  
  # ── Initialise circos ─────────────────────────────────────────────────────
  
  circos.clear()
  circos.par(
    cell.padding  = c(0, 0, 0, 0),
    track.margin  = c(0, 0.01),
    start.degree  = 81,
    gap.degree    = 2
  )
  circos.par("canvas.xlim" = c(-1.2, 1.2), "canvas.ylim" = c(-1.2, 1.2))
  
  circos.initialize(
    factors = condition_circos_metadata$condition,
    xlim    = cbind(condition_circos_metadata$xmin, condition_circos_metadata$xmax)
  )
  
  colour_var <- paste0("condition_color_", timepoint)
  text_var   <- paste0("text_color_",      timepoint)
  expr_var   <- paste0("condition_expr_available_", timepoint)
  
  if (!is.null(plot_title))
    title(plot_title, line = title_line, cex.main = title_size)
  
  # ── Layer 1: condition segments ───────────────────────────────────────────
  
  suppressMessages({
    circos.trackPlotRegion(
      ylim         = c(0, 1),
      factors      = condition_circos_metadata$condition,
      track.height = 0.1,
      bg.lwd       = 1,
      panel.fun    = function(x, y) {
        name  <- get.cell.meta.data("sector.index")
        i     <- get.cell.meta.data("sector.numeric.index")
        xlim  <- get.cell.meta.data("xlim")
        ylim  <- get.cell.meta.data("ylim")
        theta <- circlize(mean(xlim), 1.3)[1, 1] %% 360
        dd    <- ifelse(theta > 200 && theta < 340, "outside", "inside")
        
        # In the panel.fun of Layer 1, replace the circos.text call with:
        circos.text(
          x          = mean(xlim),
          y          = ylim[2] + 0.4,
          labels     = name,
          facing     = "bending.outside",   # curves letters along the arc
          cex        = 1.5,
          adj        = c(0.5, 0.5),
          niceFacing = TRUE,                 # flips text on the bottom half so it reads naturally
          col        = condition_circos_metadata[[text_var]][i]
        )
        circos.rect(
          xleft   = xlim[1], ybottom = ylim[1],
          xright  = xlim[2], ytop    = ylim[2],
          col     = condition_circos_metadata[[colour_var]][i],
          border  = "black"
        )
        circos.axis(labels = FALSE, major.tick = FALSE)
      }
    )
  })
  
  # ── Layer 2: correlation with immune response (ring == "all" only) ─────────
  
  if (ring == "all") {
    resp_var <- paste0("condition_response_available_", timepoint)
    
    suppressMessages({
      cor_data <- results_df_circos %>%
        filter(time == timepoint, !is.na(.data[[correlation]])) %>%
        group_by(condition) %>%
        summarise(
          cor_scores = list(.data[[correlation]]),
          gs.names   = list(gs.name),
          .groups    = "drop"
        ) %>%
        split(.$condition)
      
      cor_data <- cor_data[sapply(cor_data, nrow) > 0]
      
      circos.trackPlotRegion(
        ylim         = c(-1, 1),
        factors      = condition_circos_metadata$condition,
        track.height = 0.1,
        bg.lwd       = 1,
        bg.border    = ifelse(condition_circos_metadata[[resp_var]], "black", "white"),
        panel.fun    = function(x, y) {
          name <- get.cell.meta.data("sector.index")
          if (!name %in% names(cor_data)) return()
          
          sec_data          <- cor_data[[name]]
          correlation_scores <- sec_data$cor_scores[[1]]
          gs.names          <- sec_data$gs.names[[1]]
          positions         <- gene_set_positions[as.character(gs.names)]
          
          circos.rect(
            xleft   = positions - 0.5, ybottom = 0,
            xright  = positions + 0.5, ytop    = correlation_scores,
            col     = ifelse(correlation_scores > 0, "purple", "orange"),
            border  = NA
          )
        }
      )
    })
  }
  
  # ── Layer 3: activation/FC scores ─────────────────────────────────────────
  
  if (ring %in% c("all", "expression")) {
    all_scores <- results_df_circos %>%
      filter(time == timepoint) %>%
      pull(!!sym(scores)) %>%
      as.numeric()
    
    clip_threshold    <- quantile(abs(all_scores), quantile_scoreclip)
    max_clipped_score <- max(abs(pmin(pmax(all_scores, -clip_threshold), clip_threshold)))
    
    act_data <- results_df_circos %>%
      filter(time == timepoint) %>%
      group_by(condition, gs.name) %>%
      summarise(avg_score = mean(!!sym(scores), na.rm = TRUE), .groups = "drop") %>%
      group_by(condition) %>%
      summarise(
        gs.names       = list(gs.name),
        raw_scores     = list(avg_score),
        .groups        = "drop"
      ) %>%
      mutate(clipped_scores = purrr::map(
        raw_scores, ~ pmin(pmax(.x, -clip_threshold), clip_threshold)
      )) %>%
      split(.$condition)
    
    act_data <- act_data[sapply(act_data, nrow) > 0]
    
    suppressMessages({
      circos.trackPlotRegion(
        ylim         = c(-1, 1),
        factors      = condition_circos_metadata$condition,
        track.height = 0.1,
        bg.lwd       = 1,
        bg.border    = ifelse(condition_circos_metadata[[expr_var]], "black", "white"),
        panel.fun    = function(x, y) {
          sector_name <- get.cell.meta.data("sector.index")
          ylim        <- get.cell.meta.data("ylim")
          if (!sector_name %in% names(act_data)) return()
          
          sector_act    <- act_data[[sector_name]]
          clipped       <- sector_act$clipped_scores[[1]]
          gs.names      <- sector_act$gs.names[[1]]
          positions     <- gene_set_positions[match(gs.names, names(gene_set_positions))]
          keep          <- !is.na(positions)
          positions     <- positions[keep]
          clipped       <- clipped[keep]
          
          bar_lengths <- ifelse(
            clipped < 0,
            clipped * -(ylim[1] / max_clipped_score),
            clipped *  (ylim[2] / max_clipped_score)
          )
          
          circos.rect(
            xleft   = positions - 0.5, ybottom = 0,
            xright  = positions + 0.5, ytop    = bar_lengths,
            col     = ifelse(clipped > 0, "red", "blue"),
            border  = NA
          )
        }
      )
    })
  }
  
  # ── Layer 4: gene set aggregate colour ring ────────────────────────────────
  
  suppressMessages({
    circos.trackPlotRegion(
      ylim         = c(0, 1),
      factors      = condition_circos_metadata$condition,
      track.height = 0.05,
      bg.lwd       = 1,
      bg.border    = ifelse(condition_circos_metadata[[expr_var]], "black", "white"),
      panel.fun    = function(x, y) {
        sector_index       <- get.cell.meta.data("sector.index")
        xlim               <- get.cell.meta.data("xlim")
        condition_available <- condition_circos_metadata[
          condition_circos_metadata$condition == sector_index, expr_var
        ]
        
        x_start <- xlim[1] + (xlim[2] - xlim[1]) * aggregate_proportions$cum_start
        x_end   <- xlim[1] + (xlim[2] - xlim[1]) * aggregate_proportions$cum_end
        
        if (condition_available) {
          col_fill   <- aggregate_proportions$colour
          border_col <- rep(NA, length(col_fill))
        } else {
          col_fill   <- rep("white", nrow(aggregate_proportions))
          border_col <- rep("white", nrow(aggregate_proportions))
        }
        
        circos.rect(
          xleft = x_start, ybottom = 0,
          xright = x_end,  ytop    = 1,
          col    = col_fill, border = border_col
        )
      }
    )
  })
  
  # ── Arcs ──────────────────────────────────────────────────────────────────
  
  if (any_links) {
    suppressMessages({
      shared_links <- shared_links %>%
        left_join(
          distinct(results_df_circos, gs.name, gs.aggregate),
          by = "gs.name"
        ) %>%
        mutate(
          gs.aggregate = as.character(gs.aggregate),
          link_colour  = aggregate_colors[gs.aggregate]
        )
      
      mapply(
        function(c1, c2, pos, col) {
          circos.link(
            sector.index1 = c1, point1 = c(pos, pos),
            sector.index2 = c2, point2 = c(pos, pos),
            col = adjustcolor(col, alpha.f = 0.3)
          )
        },
        shared_links$condition1,
        shared_links$condition2,
        shared_links$position,
        shared_links$link_colour
      )
    })
  }
  
  # ── Legend (inline, when requested) ───────────────────────────────────────
  
  if (legend) {
    draw_circos_legend(aggregates_name = aggregates_name, ring = ring)
  }
}

# -----------------------------------------------------------------------------

# Draw the circos legend, either populated or as a blank placeholder to
# reserve space in a multi-panel layout.
draw_circos_legend <- function(aggregates_name,
                               ring        = c("all", "expression", "none"),
                               placeholder = FALSE,
                               scale       = 1,
                               font_scale  = 1) {
  ring <- match.arg(ring)
  
  aggregate_colors <- results_df %>%
    select(gs.aggregate, gs.colour) %>%
    distinct() %>%
    { setNames(.$gs.colour, .$gs.aggregate) }
  
  # Build label / colour / spacing vectors depending on ring type
  if (ring == "all") {
    legend_labels      <- c(
      expression(bold("Gene Set Aggregate")), aggregates_name,
      expression(bold("Direction of Regulation")), "Upregulated", "Downregulated",
      expression(bold("Correlation with Ab Response")), "Positive", "Negative"
    )
    legend_colours     <- c(NA, aggregate_colors[aggregates_name], NA, "red", "blue", NA, "purple", "orange")
    legend_cex_base    <- c(1.5, rep(1.2, length(aggregates_name)), 1.5, 1.2, 1.2, 1.5, 1.2, 1.2)
    legend_yisp_base   <- c(1,   rep(0.8, length(aggregates_name)), 1.5, 0.8, 0.8, 1.5, 0.8, 0.8)
    
  } else if (ring == "expression") {
    legend_labels      <- c(
      expression(bold("Gene Set Aggregate")), aggregates_name,
      expression(bold("Direction of Regulation")), "Upregulated", "Downregulated"
    )
    legend_colours     <- c(NA, aggregate_colors[aggregates_name], NA, "red", "blue")
    legend_cex_base    <- c(1.5, rep(1.2, length(aggregates_name)), 1.5, 1.2, 1.2)
    legend_yisp_base   <- c(1,   rep(0.8, length(aggregates_name)), 1.5, 0.8, 0.8)
    
  } else {
    legend_labels      <- c(expression(bold("Gene Set Aggregate")), aggregates_name)
    legend_colours     <- c(NA, aggregate_colors[aggregates_name])
    legend_cex_base    <- c(1.5, rep(1.2, length(aggregates_name)))
    legend_yisp_base   <- c(1,   rep(0.8, length(aggregates_name)))
  }
  
  legend_cex  <- legend_cex_base  * scale * font_scale
  legend_yisp <- legend_yisp_base * scale * font_scale
  text_width  <- 1.2 * max(strwidth(legend_labels, cex = max(legend_cex)))
  
  if (placeholder) {
    legend(
      "center",
      legend     = rep(" ", length(legend_labels)),
      fill       = rep(NA,  length(legend_labels)),
      border     = rep(NA,  length(legend_labels)),
      cex        = legend_cex,
      y.intersp  = legend_yisp,
      text.width = text_width,
      ncol       = 1, bty = "n", xjust = 0,
      text.col   = rep(NA, length(legend_labels)),
      title      = NULL
    )
  } else {
    legend(
      "center",
      legend     = legend_labels,
      fill       = legend_colours,
      border     = "white",
      cex        = legend_cex,
      y.intersp  = legend_yisp,
      text.width = text_width,
      ncol       = 1, bty = "n", xjust = 0,
      text.col   = rep("black", length(legend_labels)),
      title      = NULL
    )
  }
}

# -----------------------------------------------------------------------------

# Draw a bold row annotation label in the centre of an empty plot region,
# or reserve blank space when placeholder = TRUE.
plot_row_annotation <- function(label       = NULL,
                                size        = 5,
                                placeholder = FALSE) {
  par(mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
       xlim = c(0, 1), ylim = c(0, 1))
  if (!placeholder && !is.null(label))
    text(0.5, 0.5, label, cex = size, font = 2)
}

# =============================================================================
# RENDERING HELPERS
# =============================================================================

# Thin wrapper: call plot_circos merging shared defaults with per-call overrides
plot_circos_default <- function(...) {
  do.call(plot_circos, c(circos_defaults, list(...)))
}

# Render one four-panel row: annotation | QuSAGE | dearseq | legend
render_row <- function(day,
                       arc,
                       legend_placeholder,
                       legend_scale  = 1.3,
                       qusage_title  = NULL,
                       dearseq_title = NULL) {
  plot_row_annotation(paste0("Day ", day), size = 6)
  
  circos.clear()
  plot_circos_default(
    method_name = "qusage",
    timepoint   = day,
    arc         = arc,
    plot_title  = qusage_title
  )
  
  circos.clear()
  plot_circos_default(
    method_name = "dearseq",
    timepoint   = day,
    arc         = arc,
    plot_title  = dearseq_title
  )
  
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  draw_circos_legend(
    aggregates_name = btm_aggregates,
    ring            = "expression",
    placeholder     = legend_placeholder,
    scale           = legend_scale
  )
}

# Open a single-day PDF (1 row × 4 columns) and render one day
render_single_day_pdf <- function(filename, day, arc) {
  pdf(fs::path(figures_folder, filename), width = 24, height = 9.5)
  on.exit(dev.off())
  
  layout(matrix(1:4, nrow = 1), widths = c(0.1, 0.35, 0.35, 0.25))
  par(mar = rep(0, 4))
  
  render_row(
    day                = day,
    arc                = arc,
    legend_placeholder = FALSE,
    legend_scale       = 1.3,
    qusage_title       = "QuSAGE",
    dearseq_title      = "dearseq"
  )
}

# Open a comparison PDF (3 rows × 4 columns) and render days 1, 3, 7.
# The legend is shown on day 3 and blank on days 1 and 7.
render_comparison_pdf <- function(filename, arc) {
  pdf(fs::path(figures_folder, filename), width = 27, height = 25)
  on.exit(dev.off())
  
  layout(
    matrix(1:12, nrow = 3, byrow = TRUE),
    widths  = c(0.1, 0.35, 0.35, 0.35),
    heights = c(0.33, 0.33, 0.33)
  )
  par(mar = rep(0, 4))
  
  render_row(1, arc, legend_placeholder = TRUE,  legend_scale = 1.3,
             qusage_title = "QuSAGE", dearseq_title = "dearseq")
  render_row(3, arc, legend_placeholder = FALSE, legend_scale = 1.5)
  render_row(7, arc, legend_placeholder = TRUE,  legend_scale = 1.3)
}

# =============================================================================
# GENERATE FIGURES
# =============================================================================

# Multi-day comparison figures
render_comparison_pdf("circos_comparison.pdf",          arc = "positive")
render_comparison_pdf("circos_comparison_negative.pdf", arc = "negative")

# Single-day figures
single_day_configs <- tribble(
  ~filename,                   ~day, ~arc,
  "circos_day1.pdf",           1,    "positive",
  "circos_day3.pdf",           3,    "positive",
  "circos_day7.pdf",           7,    "positive",
  "circos_day1_negative.pdf",  1,    "negative",
  "circos_day3_negative.pdf",  3,    "negative",
  "circos_day7_negative.pdf",  7,    "negative"
)

purrr::pwalk(single_day_configs, render_single_day_pdf)