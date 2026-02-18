# ============================================================
# CLUSTERING + HEATMAP PIPELINE (refactored)
# - single place to set analysis arguments
# - removes duplicated code blocks
# - saves every figure with parameter-tagged filenames
# ============================================================

# -----------------------------
# 0) USER CONFIG (EDIT HERE)
# -----------------------------
CFG <- list(
  # Paths
  workdir = "C:/Users/yvandenboer/OneDrive - Institute of Natural Sciences/1_KBIN/99_General/6_Colleagues/2026.02.11_Katja_R",
  excel_file = "PhD_DLS_SLS_countings.xlsx",
  sheets = list(own = 1, all = 2),

  # Output folders
  out_root = "Plots",
  out_dir_cluster = file.path("Plots", "Cluster"),
  out_dir_heatmap = file.path("Plots", "Heatmap"),

  # General data cleaning
  drop_taxa = c("Sterile", "No ostracods", "fragmented"),
  drop_localities = c("SLS"),                 # the generic locality label to remove
  drop_bad_lit_localities = c(),
  drop_family_level_suffix = "ae$",           # remove entries ending with '...ae' at genus level

  # Default analysis settings (used unless overridden in runs below)
  defaults = list(
    level = "species",                        # "species" or "genus"
    use_all_data = TRUE,                      # TRUE = union PA across own+published; FALSE = own abundances
    plot_points = "Locality",                 # "Locality" or "Basin"
    include_groups = TRUE,                    # keep Genus/Taxon entries containing "Group"
    hull_group_var = "LakeSystem",            # "LakeSystem", "Basin", or NULL
    min_taxa_present = 3L,                    # filter localities with fewer present taxa
    exclude_locality = c(),
    k_override = NULL,                        # force k; otherwise best silhouette k
    seed = 123                                # for NMDS/metaMDS reproducibility
  ),

  # Which plots to make
  make = list(
    dendrogram = TRUE,
    silhouette = TRUE,
    pam = TRUE,
    pcoa = TRUE,
    nmds = TRUE
  ),

  # Runs: you asked for:
  # - by Locality
  # - species + genus
  # - all_data + own_data
  # - groups included + excluded
  # - min 3 taxa per locality (min_taxa_present = 3)
  # - no min occurrence threshold for heatmap
  runs = list(
    clustering = list(
      # species
      list(name = "species_all_groups",    level = "species", use_all_data = TRUE,  plot_points = "Locality", include_groups = TRUE,  hull_group_var = "LakeSystem", min_taxa_present = 3L)
      # list(name = "species_all_nogroups",  level = "species", use_all_data = TRUE,  plot_points = "Locality", include_groups = FALSE, hull_group_var = "LakeSystem", min_taxa_present = 3L),
      # list(name = "species_own_groups",    level = "species", use_all_data = FALSE, plot_points = "Locality", include_groups = TRUE,  hull_group_var = "LakeSystem", min_taxa_present = 3L),
      # list(name = "species_own_nogroups",  level = "species", use_all_data = FALSE, plot_points = "Locality", include_groups = FALSE, hull_group_var = "LakeSystem", min_taxa_present = 3L),
      
      # genus
      # list(name = "genus_all_groups",      level = "genus",   use_all_data = TRUE,  plot_points = "Locality", include_groups = TRUE,  hull_group_var = "LakeSystem", min_taxa_present = 3L),
      # list(name = "genus_all_nogroups",    level = "genus",   use_all_data = TRUE,  plot_points = "Locality", include_groups = FALSE, hull_group_var = "LakeSystem", min_taxa_present = 3L),
      # list(name = "genus_own_groups",      level = "genus",   use_all_data = FALSE, plot_points = "Locality", include_groups = TRUE,  hull_group_var = "LakeSystem", min_taxa_present = 3L),
      # list(name = "genus_own_nogroups",    level = "genus",   use_all_data = FALSE, plot_points = "Locality", include_groups = FALSE, hull_group_var = "LakeSystem", min_taxa_present = 3L)
    ),

    heatmap = list(
      # Source-coded heatmaps (0/1/2/3), by Locality:
      list(name = "species_groups",   level = "species", use_all_data = TRUE, plot_points = "Locality", include_groups = TRUE, min_taxa_present = 3L, lake_system_bar_height_mm = 3),
      list(name = "species_nogroups", level = "species", use_all_data = TRUE, plot_points = "Locality", include_groups = FALSE, min_taxa_present = 3L, lake_system_bar_height_mm = 3),
      list(name = "genus_groups",     level = "genus", use_all_data = TRUE,   plot_points = "Locality", include_groups = TRUE, min_taxa_present = 3L, lake_system_bar_height_mm = 3),
      list(name = "genus_nogroups",   level = "genus", use_all_data = TRUE,   plot_points = "Locality", include_groups = FALSE, min_taxa_present = 3L, lake_system_bar_height_mm = 3),
      # Own-data heatmaps (exclude literature; produces raw/rel/hel/PA heatmaps in one go):
      list(name = "species_groups_own",   level = "species", use_all_data = FALSE, plot_points = "Locality", include_groups = TRUE,  min_taxa_present = 3L, lake_system_bar_height_mm = 3),
      list(name = "genus_groups_own",     level = "genus",   use_all_data = FALSE, plot_points = "Locality", include_groups = TRUE,  min_taxa_present = 3L, lake_system_bar_height_mm = 3)

    )
  )
)

# -----------------------------
# 1) SETUP
# -----------------------------
setwd(CFG$workdir)

if (!requireNamespace("ragg", quietly = TRUE)) install.packages("ragg")

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(vegan)
  library(cluster)
  library(ggdendro)
  library(dendextend)
  library(ggnewscale)
  library(ComplexHeatmap)
  library(circlize)
})

dir.create(CFG$out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_dir_cluster, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_dir_heatmap, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 2) HELPERS

metadata_mapping <- NULL
# -----------------------------

slugify <- function(x) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+", "-", x)
  x <- gsub("[^A-Za-z0-9._-]+", "", x)
  x <- gsub("-{2,}", "-", x)
  x <- gsub("(^-|-$)", "", x)
  x
}


format_level_label <- function(level) {
  level <- match.arg(level, c("genus", "species"))
  paste0(level, "-level")
}

format_transform_label <- function(nm) {
  switch(nm,
         pa_jc = "P/A (Jaccard)",
         raw_bc = "Abundance (Bray-Curtis)",
         rel_bc = "Relative abundance (Bray-Curtis)",
         hel_bc = "Abundance (Hellinger)",
         nm)
}

save_tiff <- function(filename, plot, width = 2500, height = 2000, res = 300) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ragg::agg_tiff(filename, width = width, height = height, units = "px", res = res, compression = "lzw")
  print(plot)
  grDevices::dev.off()
}

drop_groups_long <- function(df, taxon_col, include_groups = TRUE) {
  if (isTRUE(include_groups)) return(df)
  df %>% dplyr::filter(!grepl("Group", .data[[taxon_col]], ignore.case = TRUE))
}

presence_absence <- function(mat) (mat > 0) * 1

filter_min_taxa_present <- function(mat, min_taxa_present = 1L) {
  min_taxa_present <- as.integer(min_taxa_present)
  if (is.na(min_taxa_present) || min_taxa_present < 1L) return(mat)
  keep <- rowSums(presence_absence(mat)) >= min_taxa_present
  mat2 <- mat[keep, , drop = FALSE]
  if (nrow(mat2) == 0L) stop("After filtering by min_taxa_present, no rows remain.")
  mat2
}


drop_zero_taxa <- function(mat) {
  # Remove taxa/columns that are absent in all retained localities.
  # This matters after you filter out low-richness localities: taxa that only occurred
  # in the removed localities become all-zeros and should not appear in heatmaps/plots.
  keep <- colSums(presence_absence(mat)) > 0
  mat2 <- mat[, keep, drop = FALSE]
  if (ncol(mat2) == 0L) stop("After dropping zero-occurrence taxa, no columns remain.")
  mat2
}

exclude_matrix_rows_by_id <- function(mat, exclude = NULL) {
  if (is.null(exclude) || length(exclude) == 0) return(mat)
  ex <- as.character(exclude)
  mat2 <- mat[!rownames(mat) %in% ex, drop = FALSE]
  if (nrow(mat2) == 0L) stop("After excluding ids, no rows remain.")
  mat2
}

rel_abund <- function(mat) {
  rs <- rowSums(mat)
  mat / rs
}

hellinger <- function(mat) {
  rs <- rowSums(mat)
  sqrt(mat / rs)
}

# -----------------------------
# 3) DATA IMPORT + CLEANING (single source of truth)
# -----------------------------

clean_taxa_table <- function(df) {
  df %>%
    mutate(Locality = str_trim(str_remove(Locality, "\\s*\\(.*\\)$"))) %>%
    filter(!Taxon %in% CFG$drop_taxa,
           !Locality %in% CFG$drop_localities) %>%
    mutate(
      Basin = str_remove(Basin, " Basin"),
      Taxon = case_match(
        Taxon,
        "Candona basic ???" ~ "Candona sp. (Trepča)",
        "Candona Group T" ~ "Candona sp. T",
        .default = str_remove(Taxon, "\\s*\\?$")
      ),
      Genus = ifelse(
        str_detect(Taxon, "Group"),
        Taxon,
        str_remove(Taxon, "[ ].*$")
      ) %>% str_remove("[ ]*sp\\.") %>% str_trim()
    ) %>%
    filter(!str_detect(Genus, CFG$drop_family_level_suffix))
}

own_raw <- read_excel(CFG$excel_file, sheet = CFG$sheets$own) %>%
  select(!c(Sample, Fraction_um, `Adult/juvenile`, Notes)) %>%
  rename("Taxon_rank" = "Taxon rank") %>%
  select(!Taxon_rank) %>%
  clean_taxa_table()

# own observations table for source coding in sheet 2
own_obs <- own_raw %>%
  distinct(Locality, Taxon) %>%
  mutate(own_data_observation = TRUE)

all_raw <- read_excel(CFG$excel_file, sheet = CFG$sheets$all) %>%
  clean_taxa_table() %>%
  filter(!Locality %in% CFG$drop_bad_lit_localities) %>%
  distinct(Locality, Taxon, .keep_all = TRUE) %>%
  left_join(own_obs, by = c("Locality", "Taxon")) %>%
  mutate(own_data_observation = coalesce(own_data_observation, FALSE))

# -----------------------------
# 4) MATRIX + METADATA BUILDERS
# -----------------------------

build_metadata_mapping <- function(level = c("genus", "species"),
                                   use_all_data = TRUE,
                                   plot_points = c("Locality", "Basin"),
                                   include_groups = TRUE) {
  level <- match.arg(level)
  plot_points <- match.arg(plot_points)
  use_all_data <- isTRUE(use_all_data)
  id_col <- if (plot_points == "Basin") "Basin" else "Locality"

  if (use_all_data) {
    long_tbl <- if (level == "genus") {
      all_raw %>% select(Locality, Basin, LakeSystem, Genus, own_data_observation) %>% distinct()
    } else {
      all_raw %>% select(Locality, Basin, LakeSystem, Taxon, own_data_observation) %>% distinct()
    }

    long_tbl <- if (level == "genus") drop_groups_long(long_tbl, "Genus", include_groups) else drop_groups_long(long_tbl, "Taxon", include_groups)

    long_tbl %>%
      mutate(
        ID = as.character(.data[[id_col]]),
        Basin = as.character(Basin),
        LakeSystem = as.character(LakeSystem),
        own_data_observation = coalesce(own_data_observation, FALSE)
      ) %>%
      group_by(ID) %>%
      summarise(
        Basin = first(stats::na.omit(Basin)),
        LakeSystem = first(stats::na.omit(LakeSystem)),
        own_data_observation = any(own_data_observation),
        .groups = "drop"
      ) %>%
      mutate(
        Basin = if (plot_points == "Basin") ID else Basin,
        Source = factor(if_else(own_data_observation, "Own observations", "Literature"),
                        levels = c("Own observations", "Literature"))
      ) %>%
      rename(Locality = ID)
  } else {
    long_tbl <- if (level == "genus") {
      own_raw %>% select(Locality, Basin, LakeSystem, Genus, Count) %>% distinct()
    } else {
      own_raw %>% select(Locality, Basin, LakeSystem, Taxon, Count) %>% distinct()
    }
    long_tbl <- if (level == "genus") drop_groups_long(long_tbl, "Genus", include_groups) else drop_groups_long(long_tbl, "Taxon", include_groups)

    long_tbl %>%
      mutate(ID = as.character(.data[[id_col]])) %>%
      group_by(ID) %>%
      summarise(
        Basin = first(stats::na.omit(as.character(Basin))),
        LakeSystem = first(stats::na.omit(as.character(LakeSystem))),
        own_data_observation = TRUE,
        .groups = "drop"
      ) %>%
      mutate(
        Basin = if (plot_points == "Basin") ID else Basin,
        Source = factor("Own observations", levels = c("Own observations", "Literature"))
      ) %>%
      rename(Locality = ID)
  }
}

get_comm_matrix <- function(level = c("genus", "species"),
                            use_all_data = TRUE,
                            plot_points = c("Locality", "Basin"),
                            include_groups = TRUE) {
  level <- match.arg(level)
  plot_points <- match.arg(plot_points)
  use_all_data <- isTRUE(use_all_data)
  id_col <- if (plot_points == "Basin") "Basin" else "Locality"

  if (use_all_data) {
    # union presence across sources: presence = 1 if in own OR published
    if (level == "genus") {
      own_pa <- own_raw %>%
        select(!!sym(id_col), Genus) %>% rename(ID = !!sym(id_col), Tax = Genus) %>%
        distinct() %>% mutate(Count = 1L) %>%
        drop_groups_long(., "Tax", include_groups)

      pub_pa <- all_raw %>%
        select(!!sym(id_col), Genus) %>% rename(ID = !!sym(id_col), Tax = Genus) %>%
        distinct() %>% mutate(Count = 1L) %>%
        drop_groups_long(., "Tax", include_groups)
    } else {
      own_pa <- own_raw %>%
        select(!!sym(id_col), Taxon) %>% rename(ID = !!sym(id_col), Tax = Taxon) %>%
        dplyr::mutate(Tax = trimws(as.character(Tax))) %>%
        dplyr::filter(grepl("\\s+", Tax)) %>%
        distinct() %>% mutate(Count = 1L) %>%
        drop_groups_long(., "Tax", include_groups)

      pub_pa <- all_raw %>%
        select(!!sym(id_col), Taxon) %>% rename(ID = !!sym(id_col), Tax = Taxon) %>%
        dplyr::mutate(Tax = trimws(as.character(Tax))) %>%
        dplyr::filter(grepl("\\s+", Tax)) %>%
        distinct() %>% mutate(Count = 1L) %>%
        drop_groups_long(., "Tax", include_groups)
    }

    mat_df <- bind_rows(own_pa, pub_pa) %>%
      group_by(ID, Tax) %>%
      summarise(Count = max(Count), .groups = "drop") %>%
      pivot_wider(names_from = Tax, values_from = Count, values_fill = 0) %>%
      as.data.frame()

    rownames(mat_df) <- mat_df$ID
    return(as.matrix(select(mat_df, -ID)))
  }

  # own only: abundance matrix
  if (level == "genus") {
    long_tbl <- own_raw %>%
      select(!!sym(id_col), Genus, Count) %>% rename(ID = !!sym(id_col), Tax = Genus) %>%
      drop_groups_long(., "Tax", include_groups) %>%
      group_by(ID, Tax) %>% summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")
  } else {
    long_tbl <- own_raw %>%
      select(!!sym(id_col), Taxon, Count) %>% rename(ID = !!sym(id_col), Tax = Taxon) %>%
      dplyr::mutate(Tax = trimws(as.character(Tax))) %>%
      dplyr::filter(grepl("\\s+", Tax)) %>%
      drop_groups_long(., "Tax", include_groups) %>%
      group_by(ID, Tax) %>% summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop")
  }

  mat_df <- long_tbl %>%
    pivot_wider(names_from = Tax, values_from = Count, values_fill = 0) %>%
    as.data.frame()
  rownames(mat_df) <- mat_df$ID
  as.matrix(select(mat_df, -ID))
}

# -----------------------------
# 5) PLOTS
# -----------------------------

hull_palette <- function(levels_vec) {
  lv <- as.character(levels_vec)
  base <- c(
    DLS = "#1b9e77",
    SLS = "#7570b3",
    PBS = "#d95f02",
    PLS = "#e7298a"
  )
  if (all(lv %in% names(base))) return(base[lv])
  setNames(scales::hue_pal(h = c(110, 470))(length(lv)), lv)
}

cluster_palette <- function(n) scales::hue_pal(h = c(0, 260))(n)

# Lighter/pastel variant of a rainbow-like qualitative palette.
# Uses HCL (Hue-Chroma-Luminance) to keep colours readable but soft.
# Default luminance is slightly lowered (darker) to improve contrast.
pastel_rainbow <- function(n, l = 68, c = 58, start = 10, end = 350, alpha = 1) {
  n <- as.integer(n)
  if (is.na(n) || n < 1L) return(character())
  hues <- seq(from = start, to = end, length.out = n + 1L)[seq_len(n)]
  grDevices::hcl(h = hues, c = c, l = l, alpha = alpha)
}

compact <- function(x) Filter(Negate(is.null), x)

# Cluster palette used consistently across:
# - taxa dendrograms
# - locality/basin dendrograms (standalone plots)
# - heatmap column dendrogram
#
# The mapping is purely by cluster id ("1".."k"), so the same parameter set and
# same k yields the same colours across all dendrogram contexts.
cluster_palette_named <- function(k, palette_fun = pastel_rainbow) {
  k <- as.integer(k)
  if (is.na(k) || k < 1L) stop("k must be a positive integer")
  stats::setNames(palette_fun(k), as.character(seq_len(k)))
}

hull_legend_title <- function(hull_group_var) {
  if (is.null(hull_group_var)) return(NULL)
  if (identical(hull_group_var, "LakeSystem")) return("Lake System")
  if (identical(hull_group_var, "Basin")) return("Basin")
  hull_group_var
}

make_hulls_poly_and_pairs <- function(df, x = "Axis1", y = "Axis2", group_var) {
  if (is.null(group_var)) return(list(poly = NULL, pairs = NULL))
  if (!group_var %in% names(df)) stop("hull_group_var not found in plotting data: ", group_var)
  
  df2 <- df %>% filter(!is.na(.data[[group_var]]), !is.na(.data[[x]]), !is.na(.data[[y]]))
  if (nrow(df2) == 0) return(list(poly = NULL, pairs = NULL))
  
  groups_n <- df2 %>%
    count(.data[[group_var]], name = "n") %>%
    rename(HullGroup = !!group_var)
  
  # Polygons: n >= 3, closed
  poly <- df2 %>%
    inner_join(groups_n %>% filter(n >= 3), by = setNames("HullGroup", group_var)) %>%
    group_by(.data[[group_var]]) %>%
    group_modify(\(d, key) {
      stopifnot("Axis1" %in% names(d), "Axis2" %in% names(d))
      
      d_xy <- dplyr::filter(d, !is.na(.data[["Axis1"]]), !is.na(.data[["Axis2"]]))
      if (nrow(d_xy) < 3) return(dplyr::slice(d_xy, 0))
      
      idx <- grDevices::chull(d_xy[["Axis1"]], d_xy[["Axis2"]])
      
      h <- d_xy[idx, , drop = FALSE]          # NOTE: comma => rows
      h <- dplyr::bind_rows(h, h[1, , drop = FALSE])  # close polygon
      h
    }) %>%
    ungroup() %>%
    rename(HullGroup = !!group_var)
  
  if (nrow(poly) == 0) poly <- NULL
  
  # Pairs: n == 2 -> segment between points
  pairs <- df2 %>%
    inner_join(groups_n %>% filter(n == 2), by = setNames("HullGroup", group_var)) %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      x = first(.data[[x]]),
      y = first(.data[[y]]),
      xend = dplyr::nth(.data[[x]], 2),
      yend = dplyr::nth(.data[[y]], 2),
      .groups = "drop"
    ) %>%
    rename(HullGroup = !!group_var)
  
  if (nrow(pairs) == 0) pairs <- NULL
  
  list(poly = poly, pairs = pairs)
}

color_dendrogram_by_group <- function(hc_obj, group_var) {
  if (is.null(group_var) || is.null(metadata_mapping)) return(as.dendrogram(hc_obj))
  
  group_var <- as.character(group_var)
  if (length(group_var) != 1L) stop("hull_group_var/group_var must be a single column name (length 1).")
  if (!group_var %in% names(metadata_mapping)) stop("group_var not in metadata_mapping: ", group_var)
  
  lab <- hc_obj$labels
  
  group_map <- metadata_mapping %>%
    dplyr::select(Locality, Group = dplyr::all_of(group_var)) %>%
    dplyr::distinct(Locality, .keep_all = TRUE)
  
  group_by_lab <- tibble::tibble(Locality = lab) %>%
    dplyr::left_join(group_map, by = "Locality") %>%
    dplyr::mutate(Group = as.character(Group)) %>%
    dplyr::pull(Group)
  
  groups <- sort(unique(stats::na.omit(group_by_lab)))
  pal <- if (length(groups) > 0) stats::setNames(scales::hue_pal()(length(groups)), groups) else character()
  na_col <- "grey60"
  
  dend <- as.dendrogram(hc_obj)
  
  set_branch_cols_inplace <- function(node) {
    if (stats::is.leaf(node)) {
      lbl <- attr(node, "label")
      g <- group_by_lab[match(lbl, lab)]
      col <- if (is.na(g) || !nzchar(g) || !(g %in% names(pal))) na_col else pal[[g]]
      attr(node, "edgePar") <- c(attr(node, "edgePar"), list(col = col, lwd = 0.7))
      return(list(node = node, groups = g))
    }
    
    gs_all <- character(0)
    for (i in seq_along(node)) {
      child <- set_branch_cols_inplace(node[[i]])
      node[[i]] <- child$node
      gs_all <- c(gs_all, child$groups)
    }
    
    gs_all_non_na <- gs_all[!is.na(gs_all)]
    g_unique <- unique(gs_all_non_na)
    col <- if (length(g_unique) == 1L) {
      g <- g_unique[[1]]
      if (!nzchar(g) || !(g %in% names(pal))) na_col else pal[[g]]
    } else {
      na_col
    }
    
    attr(node, "edgePar") <- c(attr(node, "edgePar"), list(col = col, lwd = 0.7))
    list(node = node, groups = gs_all)
  }
  
  dend2 <- set_branch_cols_inplace(dend)$node
  attr(dend2, "group_palette") <- pal
  attr(dend2, "na_col") <- na_col
  dend2
}

make_hulls_poly <- function(df, x = "Axis1", y = "Axis2", group_var) {
  if (is.null(group_var) || !group_var %in% names(df)) return(NULL)
  
  df2 <- df |>
    dplyr::filter(
      !is.na(.data[[group_var]]),
      !is.na(.data[[x]]),
      !is.na(.data[[y]])
    )
  
  if (nrow(df2) == 0) return(NULL)
  
  # Keep only groups with >= 3 points (valid polygon hull)
  df2 <- df2 |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::filter(dplyr::n() >= 3) |>
    dplyr::ungroup()
  
  if (nrow(df2) == 0) return(NULL)
  
  df2 |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::slice(chull(.data[[x]], .data[[y]])) |>
    dplyr::ungroup()
}

make_hulls_line <- function(df, x = "Axis1", y = "Axis2", group_var) {
  if (is.null(group_var) || !group_var %in% names(df)) return(NULL)
  
  df2 <- df |>
    dplyr::filter(
      !is.na(.data[[group_var]]),
      !is.na(.data[[x]]),
      !is.na(.data[[y]])
    )
  
  if (nrow(df2) == 0) return(NULL)
  
  # Keep only groups with exactly 2 points (draw a line)
  df2 <- df2 |>
    dplyr::group_by(.data[[group_var]]) |>
    dplyr::filter(dplyr::n() == 2) |>
    dplyr::ungroup()
  
  if (nrow(df2) == 0) return(NULL)
  
  df2
}


# -----------------------------
# DENDROGRAM COLOURING HELPERS
# (must be defined before plot_dendrogram())
# -----------------------------

# Colour dendrogram edges by cluster membership (branches only; labels handled separately).
color_dend_by_membership <- function(dend, membership, palette, mixed_col = "grey40", lwd_col = 2.5, lwd_default = 1) {
  stopifnot(inherits(dend, "dendrogram"))
  if (is.null(names(membership))) stop("membership must be a named vector (names are leaf labels).")
  if (is.null(names(palette))) stop("palette must be a named vector (names are cluster IDs as character).")

  rec <- function(node) {
    if (is.leaf(node)) {
      lab <- labels(node)
      if (!lab %in% names(membership)) stop("membership missing label: ", lab)
      id <- as.character(membership[[lab]])
      if (!id %in% names(palette)) stop("palette missing cluster id: ", id)
      attr(node, "leaf_ids") <- id
      # edgePar on leaves controls the terminal branch (to the leaf)
      ep <- attr(node, "edgePar")
      if (is.null(ep) || !is.list(ep)) ep <- list()
      ep <- modifyList(ep, list(col = unname(palette[[id]]), lwd = lwd_col))
      attr(node, "edgePar") <- ep
      return(node)
    }

    node[] <- lapply(node, rec)
    ids <- unique(unlist(lapply(node, function(x) attr(x, "leaf_ids"))))
    attr(node, "leaf_ids") <- ids

    col_use <- if (length(ids) == 1L) unname(palette[[ids]]) else mixed_col
    ep <- attr(node, "edgePar")
    if (is.null(ep) || !is.list(ep)) ep <- list()
    lwd_use <- if (length(ids) == 1L) lwd_col else lwd_default
    ep <- modifyList(ep, list(col = col_use, lwd = lwd_use))
    attr(node, "edgePar") <- ep
    node
  }

  dend2 <- rec(dend)

  # Clean up helper attribute
  clean <- function(node) {
    attr(node, "leaf_ids") <- NULL
    if (!is.leaf(node)) node[] <- lapply(node, clean)
    node
  }
  clean(dend2)
}

# Ensure internal node heights are strictly increasing (avoids drawing artefacts with ties).
bump_zero_heights <- function(dend, eps = 1e-6) {
  rec <- function(node) {
    if (is.leaf(node)) return(node)
    node[] <- lapply(node, rec)
    ch <- vapply(node, function(x) attr(x, "height"), numeric(1))
    h  <- attr(node, "height")
    if (is.null(h)) h <- max(ch, na.rm = TRUE)
    min_ok <- max(ch, na.rm = TRUE)
    if (!is.finite(h) || h <= min_ok) h <- min_ok + eps
    attr(node, "height") <- h
    node
  }
  rec(dend)
}

# Build dendrogram from hclust and colour branches according to a leaf->cluster mapping (or cutree()).
color_cluster_leaf_branches <- function(hc_obj,
                                       k,
                                       leaf_group_ids = NULL,
                                       palette_k = NULL,
                                       mixed_col = "grey70",
                                       bump_zero = FALSE,
                                       eps = 1e-6) {
  stopifnot(inherits(hc_obj, "hclust"))
  dend <- as.dendrogram(hc_obj)
  if (isTRUE(bump_zero)) dend <- bump_zero_heights(dend, eps = eps)

  if (is.null(leaf_group_ids)) {
    leaf_group_ids <- cutree(hc_obj, k = k)
  }
  if (is.null(names(leaf_group_ids))) stop("leaf_group_ids must be a named vector (names are leaf labels).")

  labs <- labels(dend)
  leaf_group_ids <- leaf_group_ids[labs]
  if (anyNA(leaf_group_ids)) {
    missing <- labs[is.na(leaf_group_ids)]
    stop("leaf_group_ids missing these labels: ", paste(missing[1:min(20, length(missing))], collapse = ", "))
  }

  ids <- as.integer(leaf_group_ids)
  id_levels <- sort(unique(ids))

  if (is.null(palette_k)) {
    palette_k <- setNames(grDevices::hcl.colors(length(id_levels), palette = "Dynamic"),
                          as.character(id_levels))
  }
  if (!all(as.character(id_levels) %in% names(palette_k))) {
    miss <- setdiff(as.character(id_levels), names(palette_k))
    stop("palette_k missing IDs: ", paste(miss, collapse = ", "))
  }

  color_dend_by_membership(dend, stats::setNames(ids, labs), palette_k, mixed_col = mixed_col)
}

plot_dendrogram <- function(hc_obj, k, title,
                            hull_group_var = NULL,
                            exclude_dendrogram_loc = NULL) {
  if (length(hc_obj$labels) < 2L) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title = title, " (dendrogram unavailable: <2 localities)"))
  }

  cl_pal <- cluster_palette_named(k)
  membership <- stats::cutree(hc_obj, k = k)
  dend <- color_cluster_leaf_branches(
    hc_obj,
    k = k,
    leaf_group_ids = membership,
    palette_k = cl_pal,
    mixed_col = "grey70",
    bump_zero = any(hc_obj$height <= 0)
  )

  gd <- dendextend::as.ggdend(dend)
  seg <- gd$segments
  if (is.null(seg) || nrow(seg) == 0L) {
    gd <- dendextend::as.ggdend(as.dendrogram(hc_obj))
    seg <- gd$segments
    seg$col <- "grey30"
  }
  if (!"col" %in% names(seg)) seg$col <- "grey30"

  y_max_raw <- suppressWarnings(max(c(seg$y, seg$yend), na.rm = TRUE))
  if (!is.finite(y_max_raw) || y_max_raw <= 0) y_max_raw <- 1
  if (y_max_raw < 1) {
    sf <- 1 / y_max_raw
    seg$y    <- seg$y * sf
    seg$yend <- seg$yend * sf
  }

  seg$y    <- pmax(seg$y, 0)
  seg$yend <- pmax(seg$yend, 0)

  labs_df <- tibble::tibble(
    Locality = hc_obj$labels[hc_obj$order],
    x = seq_along(hc_obj$order),
    Cluster = as.character(unname(stats::cutree(hc_obj, k = k)[hc_obj$labels[hc_obj$order]]))
  ) %>%
    dplyr::left_join(metadata_mapping, by = "Locality") %>%
    dplyr::mutate(
      LakeSystem = as.character(LakeSystem),
      LakeSystem_plot = dplyr::if_else(
        is.na(LakeSystem) | LakeSystem == "Unknown",
        NA_character_, LakeSystem
      )
    )

  if (!is.null(exclude_dendrogram_loc) && length(exclude_dendrogram_loc) > 0) {
    labs_df <- labs_df %>% dplyr::filter(!Locality %in% as.character(exclude_dendrogram_loc))
  }

  y_max <- max(c(seg$y, seg$yend), na.rm = TRUE)
  y_scale <- max(y_max, 1)

  ## ---- VISUAL TIP GAP (key edit) ----
  tip_gap_frac <- 0.05
  tip_shift <- tip_gap_frac * y_scale
  seg$y    <- seg$y + tip_shift
  seg$yend <- seg$yend + tip_shift
  ## ----------------------------------

  node_gap_frac  <- 0.09
  label_gap_frac <- 0.06

  tip_y   <- 0
  node_y  <- -node_gap_frac * y_scale
  label_y <- node_y - label_gap_frac * y_scale
  y_min   <- label_y - 0.30 * y_scale

  tip_pos <- if ("label" %in% names(seg)) {
    seg %>%
      dplyr::filter(!is.na(.data$label), .data$x == .data$xend, .data$yend == tip_shift) %>%
      dplyr::distinct(.data$label, .data$xend)
  } else {
    gd$labels %>%
      dplyr::transmute(label = .data$label, xend = .data$x)
  }

  leaf_ext <- tip_pos %>%
    dplyr::left_join(
      labs_df %>% dplyr::select(Locality, Cluster),
      by = c("label" = "Locality")
    ) %>%
    dplyr::mutate(
      x = xend,
      y = tip_y,
      yend = node_y,
      col = unname(cl_pal[as.character(Cluster)])
    )

  labs_df$Cluster <- factor(labs_df$Cluster, levels = names(cl_pal))

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = seg,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = I(col)),
      linewidth = 1
    ) +
    ggplot2::geom_segment(
      data = leaf_ext,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = I(col)),
      linewidth = 1
    ) +
    ggplot2::scale_color_identity(guide = "none") +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(
      data = labs_df,
      ggplot2::aes(x = x, y = node_y, shape = Source, color = LakeSystem_plot),
      size = 2.1
    ) +
    ggplot2::scale_color_manual(
      values = hull_palette(sort(unique(stats::na.omit(as.character(labs_df$LakeSystem_plot))))),
      na.value = "grey60",
      guide = ggplot2::guide_legend(title = "Lake system"),
      na.translate = FALSE
    ) +
    ggnewscale::new_scale_color() +
    ggplot2::geom_text(
      data = labs_df,
      ggplot2::aes(x = x, y = label_y, label = Locality),
      angle = 90, hjust = 1, vjust = 0.5, size = 2.5, color = "grey20"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Own observations" = 16, "Literature" = 17),
      drop = TRUE, na.translate = FALSE
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 5.5, 0, 5.5, unit = "pt"),
      panel.border = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_min, max(y_max + tip_shift, 1)),
      clip = "off"
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(color = "grey30"))
    )
}

pcoa_df <- function(dist_obj, mat, hc_obj, k) {
  fit <- stats::cmdscale(dist_obj, eig = TRUE, k = 2)

  tibble::tibble(
    Locality = as.character(rownames(mat)),
    Axis1    = fit$points[, 1],
    Axis2    = fit$points[, 2],
    Cluster  = factor(stats::cutree(hc_obj, k = k)[rownames(mat)])
  ) %>%
    dplyr::mutate(Locality = stringr::str_trim(Locality)) %>%
    dplyr::left_join(
      metadata_mapping %>%
        dplyr::mutate(Locality = stringr::str_trim(as.character(Locality))) %>%
        dplyr::select(Locality, Basin, LakeSystem, own_data_observation, Source),
      by = "Locality"
    ) %>%
    dplyr::mutate(
      own_data_observation = dplyr::coalesce(own_data_observation, FALSE),
      Source = dplyr::coalesce(Source, factor(
        dplyr::if_else(own_data_observation, "Own observations", "Literature"),
        levels = c("Own observations", "Literature")
      ))
    )
}


plot_pcoa <- function(dist_obj, mat, hc_obj, k,
                      title,
                      hull_group_var = NULL,
                      hull_alpha = 0.15,
                      xlim = NULL, ylim = NULL,
                      exclude_pcoa_loc = NULL) {

  df <- pcoa_df(dist_obj, mat, hc_obj, k)

  if (!is.null(exclude_pcoa_loc) && length(exclude_pcoa_loc) > 0) {
    df <- df %>% dplyr::filter(!Locality %in% as.character(exclude_pcoa_loc))
  }

  df$Cluster <- droplevels(df$Cluster)

  hulls <- make_hulls_poly_and_pairs(df, x = "Axis1", y = "Axis2", group_var = hull_group_var)
  hull_name <- hull_legend_title(hull_group_var)

  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "PCoA axis 1", y = "PCoA axis 2") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(hulls$poly) || !is.null(hulls$pairs)) {
    hull_levels <- sort(unique(c(
      if (!is.null(hulls$poly))  as.character(hulls$poly$HullGroup)  else character(0),
      if (!is.null(hulls$pairs)) as.character(hulls$pairs$HullGroup) else character(0)
    )))
    hp <- hull_palette(hull_levels)

    if (!is.null(hulls$poly)) {
      p <- p +
        ggplot2::geom_polygon(
          data = hulls$poly,
          ggplot2::aes(x = Axis1, y = Axis2, group = HullGroup, fill = HullGroup),
          inherit.aes = FALSE,
          alpha = hull_alpha,
          color = NA
        ) +
        ggplot2::scale_fill_manual(values = hp, guide = "none")
    }

    p <- p + ggnewscale::new_scale_color()

    if (!is.null(hulls$poly)) {
      p <- p +
        ggplot2::geom_path(
          data = hulls$poly,
          ggplot2::aes(x = Axis1, y = Axis2, group = HullGroup, color = HullGroup),
          inherit.aes = FALSE,
          linewidth = 0.35
        )
    }

    if (!is.null(hulls$pairs)) {
      p <- p +
        ggplot2::geom_segment(
          data = hulls$pairs,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = HullGroup),
          inherit.aes = FALSE,
          linewidth = 0.5
        )
    }

    p <- p + ggplot2::scale_color_manual(values = hp, name = hull_name)
  }

  p <- p +
    ggnewscale::new_scale_color() +
    ggplot2::geom_text(
      data = df,
      ggplot2::aes(Axis1, Axis2, label = Locality, color = Cluster),
      vjust = -1.25, size = 3, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(Axis1, Axis2, color = Cluster, shape = Source),
      size = 3,
      show.legend = c(color = FALSE, shape = TRUE)
    ) +
    ggplot2::scale_color_manual(values = cluster_palette(length(levels(df$Cluster))), guide = "none") +
    ggplot2::scale_shape_manual(values = c("Own observations" = 16, "Literature" = 17), drop = FALSE) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(color = "grey30")))

  if (!is.null(xlim) || !is.null(ylim)) p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  p
}


plot_nmds <- function(mat, dist_method_for_nmds, hc_obj, k,
                      title,
                      hull_group_var = NULL,
                      hull_alpha = 0.15,
                      seed = 123, nmds_k = 2, trymax = 100,
                      exclude_nmds_loc = NULL) {

  mat_fit <- exclude_matrix_rows_by_id(mat, exclude = exclude_nmds_loc)
  if (nrow(mat_fit) < 2L) stop("NMDS requires at least 2 points after exclude_nmds_loc filtering.")

  set.seed(seed)
  fit <- vegan::metaMDS(mat_fit, distance = dist_method_for_nmds, k = nmds_k, trymax = trymax)

  df <- tibble::tibble(
    Locality = as.character(rownames(mat_fit)),
    Axis1    = fit$points[, 1],
    Axis2    = fit$points[, 2],
    Cluster  = factor(stats::cutree(hc_obj, k = k)[rownames(mat_fit)])
  ) %>%
    dplyr::mutate(Locality = stringr::str_trim(Locality)) %>%
    dplyr::left_join(
      metadata_mapping %>%
        dplyr::mutate(Locality = stringr::str_trim(as.character(Locality))) %>%
        dplyr::select(Locality, Basin, LakeSystem, own_data_observation, Source),
      by = "Locality"
    ) %>%
    dplyr::mutate(
      own_data_observation = dplyr::coalesce(own_data_observation, FALSE),
      Source = dplyr::coalesce(Source, factor(
        dplyr::if_else(own_data_observation, "Own observations", "Literature"),
        levels = c("Own observations", "Literature")
      ))
    )

  df$Cluster <- droplevels(df$Cluster)

  hulls <- make_hulls_poly_and_pairs(df, x = "Axis1", y = "Axis2", group_var = hull_group_var)
  hull_name <- hull_legend_title(hull_group_var)

  stress_txt <- paste0("stress: ", formatC(fit$stress, format = "f", digits = 3))

  p <- ggplot2::ggplot() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = title, x = "NMDS axis 1", y = "NMDS axis 2") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(hulls$poly) || !is.null(hulls$pairs)) {
    hull_levels <- sort(unique(c(
      if (!is.null(hulls$poly))  as.character(hulls$poly$HullGroup)  else character(0),
      if (!is.null(hulls$pairs)) as.character(hulls$pairs$HullGroup) else character(0)
    )))
    hp <- hull_palette(hull_levels)

    if (!is.null(hulls$poly)) {
      p <- p +
        ggplot2::geom_polygon(
          data = hulls$poly,
          ggplot2::aes(x = Axis1, y = Axis2, group = HullGroup, fill = HullGroup),
          inherit.aes = FALSE,
          alpha = hull_alpha,
          color = NA,
          show.legend = FALSE
        ) +
        ggplot2::scale_fill_manual(values = hp, guide = "none")
    }

    p <- p + ggnewscale::new_scale_color()

    if (!is.null(hulls$poly)) {
      p <- p +
        ggplot2::geom_path(
          data = hulls$poly,
          ggplot2::aes(x = Axis1, y = Axis2, group = HullGroup, color = HullGroup),
          inherit.aes = FALSE,
          linewidth = 0.35
        )
    }

    if (!is.null(hulls$pairs)) {
      p <- p +
        ggplot2::geom_segment(
          data = hulls$pairs,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = HullGroup),
          inherit.aes = FALSE,
          linewidth = 0.5
        )
    }

    p <- p + ggplot2::scale_color_manual(values = hp, name = hull_name)
  }

  p +
    ggnewscale::new_scale_color() +
    ggplot2::geom_text(
      data = df,
      ggplot2::aes(Axis1, Axis2, label = Locality, color = Cluster),
      vjust = -0.8, size = 3, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(Axis1, Axis2, color = Cluster, shape = Source),
      size = 3,
      show.legend = c(color = FALSE, shape = TRUE)
    ) +
    ggplot2::scale_color_manual(values = cluster_palette(length(levels(df$Cluster))), guide = "none") +
    ggplot2::scale_shape_manual(values = c("Own observations" = 16, "Literature" = 17), drop = FALSE) +
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(color = "grey30"))) +
    ggplot2::annotate("text", x = -Inf, y = -Inf, label = stress_txt, hjust = -0.05, vjust = -0.6, size = 3)
}


silhouette_curve_hclust <- function(dist_obj, hc_obj) {
  n <- attr(dist_obj, "Size")
  if (is.null(n) || n < 3) return(tibble(k = integer(), mean_silhouette = numeric()))
  ks <- 2:(n - 1)
  ms <- sapply(ks, function(k) {
    cl <- cutree(hc_obj, k = k)
    if (length(unique(cl)) < 2) return(NA_real_)
    sil <- silhouette(cl, dist_obj)
    mean(sil[, 3], na.rm = TRUE)
  })
  tibble(k = ks, mean_silhouette = ms)
}

plot_silhouette <- function(sil_df, title) {
  ggplot(sil_df, aes(k, mean_silhouette)) +
    geom_line() + geom_point() +
    labs(title = title, x = "k", y = "Mean silhouette") +
    theme_bw()
}

# -----------------------------
# 6) CLUSTERING PIPELINE
# -----------------------------

run_clustering <- function(run_params) {
  # merge defaults with run-specific overrides
  p <- modifyList(CFG$defaults, run_params)

  level <- match.arg(p$level, c("genus","species"))
  plot_points <- match.arg(p$plot_points, c("Locality","Basin"))
  use_all_data <- isTRUE(p$use_all_data)

  mat <- get_comm_matrix(level = level, use_all_data = use_all_data, plot_points = plot_points, include_groups = isTRUE(p$include_groups))
  mat <- filter_min_taxa_present(mat, p$min_taxa_present)
  mat <- exclude_matrix_rows_by_id(mat, p$exclude_locality)
  mat <- drop_zero_taxa(mat)

  meta_map <- build_metadata_mapping(level = level, use_all_data = use_all_data, plot_points = plot_points, include_groups = isTRUE(p$include_groups))
  metadata_mapping <<- meta_map

  # distances:
  mats <- list()
  if (use_all_data) {
    # PA-only rule for all_data
    mats$pa_jc <- presence_absence(mat)
  } else {
    mats$raw_bc <- mat
    mats$rel_bc <- rel_abund(mat)
    mats$hel_bc <- hellinger(mat)
    mats$pa_jc  <- presence_absence(mat)
  }

  file_prefix <- file.path(CFG$out_dir_cluster, run_params$name)

  out <- list(params = p, mats = list(), plots = list())

  for (nm in names(mats)) {
    m <- mats[[nm]]
    dist_obj <- vegdist(m, method = if (nm == "pa_jc") "jaccard" else "bray", binary = (nm == "pa_jc"))
    hc_obj <- hclust(dist_obj, method = "ward.D2")

    sil_df <- if (isTRUE(CFG$make$silhouette)) silhouette_curve_hclust(dist_obj, hc_obj) else NULL
    k_best <- if (!is.null(p$k_override)) as.integer(p$k_override) else {
      if (!is.null(sil_df) && nrow(sil_df) > 0) sil_df$k[which.max(sil_df$mean_silhouette)] else 2L
    }

    data_label <- if (use_all_data) "All data" else "Own data"
    level_label <- format_level_label(level)
    transform_label <- format_transform_label(nm)
    title_core <- paste0(data_label, ", ", level_label, ", ", transform_label, " (k = ", k_best, ")")

    if (isTRUE(CFG$make$silhouette) && !is.null(sil_df) && nrow(sil_df) > 0) {
      ps <- plot_silhouette(sil_df, title = paste0("Silhouette: ", title_core))
      save_tiff(paste0(file_prefix, "_silhouette_", nm, "_minTaxa", p$min_taxa_present, ".tiff"), ps, width = 2400, height = 1800)
      out$plots[[paste0(nm,"_sil")]] <- ps
    }

    if (isTRUE(CFG$make$pcoa)) {
      pp <- plot_pcoa(dist_obj, m, hc_obj, k = k_best, title = paste0("PCoA: ", title_core), hull_group_var = p$hull_group_var)
      save_tiff(paste0(file_prefix, "_PCoA_", nm, "_minTaxa", p$min_taxa_present, ".tiff"), pp)
      out$plots[[paste0(nm,"_pcoa")]] <- pp
    }

    if (isTRUE(CFG$make$nmds)) {
      pn <- plot_nmds(m, dist_method_for_nmds = "bray", hc_obj, k = k_best, title = paste0("NMDS: ", title_core),
                      hull_group_var = p$hull_group_var, seed = p$seed)
      save_tiff(paste0(file_prefix, "_NMDS_", nm, "_minTaxa", p$min_taxa_present, ".tiff"), pn)
      out$plots[[paste0(nm,"_nmds")]] <- pn
    }

    if (isTRUE(CFG$make$dendrogram)) {
      pd <- plot_dendrogram(
        hc_obj,
        k = k_best,
        title = paste0("Dendrogram: ", title_core),
        hull_group_var = p$hull_group_var,
        exclude_dendrogram_loc = p$exclude_locality
      )
      save_tiff(paste0(file_prefix, "_dendrogram_", nm, "_minTaxa", p$min_taxa_present, ".tiff"), pd, width = 2600, height = 1800)
      out$plots[[paste0(nm,"_dend")]] <- pd
    }

    out$mats[[nm]] <- list(mat = m, dist = dist_obj, hc = hc_obj, k = k_best)
  }

  out
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# -----------------------------
# 7) HEATMAPS (taxa x locality with clustering-aligned dendrogram styling)
# Rules:
# - species-level heatmap EXCLUDES genus-only / higher-taxon records
# - genus-level heatmap combines genus-only observations with any species in that genus
#   (i.e., collapse to Genus, same as clustering)
# - cells are Present / Absent (sage green / almost white)
# - two locality annotation bars: LakeSystem + Source (Own observations vs Literature)
# -----------------------------

is_species_taxon <- function(x) {
  # Your Taxon format is typically either:
  # - species-level: "bihacensis, Amplocypris"  (contains a comma)
  # - genus/higher: "Amplocypris" / "Candonidae" / "Candona sp. T" (no comma)
  # For species heatmaps, we only keep the comma-formatted records.
  (str_detect(x, " ") | (str_detect(x, " sp.")) & !str_detect(x, "Group"))
}

make_presence_matrix <- function(level = c("genus","species"),
                                 plot_points = c("Locality","Basin"),
                                 include_groups = TRUE) {
  level <- match.arg(level)
  plot_points <- match.arg(plot_points)

  id_col <- if (plot_points == "Basin") "Basin" else "Locality"

  if (level == "genus") {
    long_tbl <- bind_rows(
      own_raw %>% mutate(.src = "own"),
      all_raw %>% mutate(.src = if_else(own_data_observation, "own", "pub"))
    ) %>%
      select(!!sym(id_col), Genus) %>%
      rename(ID = !!sym(id_col), Tax = Genus) %>%
      distinct() %>%
      mutate(PA = 1L) %>%
      drop_groups_long(., "Tax", include_groups)
  } else {
    # species-only: exclude genus-level and higher taxa
    long_tbl <- bind_rows(
      own_raw %>% mutate(.src = "own"),
      all_raw %>% mutate(.src = if_else(own_data_observation, "own", "pub"))
    ) %>%
      filter(is_species_taxon(Taxon)) %>%
      select(!!sym(id_col), Taxon) %>%
      rename(ID = !!sym(id_col), Tax = Taxon) %>%
      distinct() %>%
      mutate(PA = 1L) %>%
      drop_groups_long(., "Tax", include_groups)
  }

  mat_df <- long_tbl %>%
    group_by(ID, Tax) %>%
    summarise(PA = max(PA), .groups = "drop") %>%
    pivot_wider(names_from = Tax, values_from = PA, values_fill = 0) %>%
    as.data.frame()

  rownames(mat_df) <- mat_df$ID
  as.matrix(select(mat_df, -ID))
}


save_heatmaps <- function(run_params) {
  p <- modifyList(
    list(
      level = "species",
      use_all_data = TRUE,                # TRUE = own + literature (PA only); FALSE = own-only (raw/rel/hel/PA)
      plot_points = "Locality",
      include_groups = TRUE,
      annotation_type = "shape",          # "shape" (default) or "bar"
      min_taxa_present = 1L,
      lake_system_bar_height_mm = 3,
      source_bar_height_mm = 3,
      # taxa-present count annotation (top)
      taxa_count_bar_height_mm = 6,
      taxa_count_text_height_mm = 3,
      # control plot height scaling (inches per row + margins)
      row_height_in = 0.12,   # ~3.0 mm per row
      base_height_in = 2.0,   # fixed margins/title/legends
      min_height_in = 4.0,
      max_height_in = 30.0
    ),
    run_params
  )

  level <- match.arg(p$level, c("genus", "species"))
  plot_points <- match.arg(p$plot_points, c("Locality", "Basin"))
  use_all_data <- isTRUE(p$use_all_data)

  # community matrix: rows = localities/basins, cols = taxa
  mat <- get_comm_matrix(level = level, use_all_data = use_all_data, plot_points = plot_points, include_groups = isTRUE(p$include_groups))
  mat <- filter_min_taxa_present(mat, p$min_taxa_present)
  mat <- exclude_matrix_rows_by_id(mat, p$exclude_locality %||% c())
  mat <- drop_zero_taxa(mat)

  if (nrow(mat) < 2L || ncol(mat) < 2L) {
    warning(sprintf("Heatmap skipped: matrix is %d sites x %d taxa after filtering (need >=2 each).", nrow(mat), ncol(mat)))
    return(invisible(NULL))
  }

  # build the same set of transformed matrices as clustering
  mats <- list()
  if (use_all_data) {
    mats$pa_jc <- presence_absence(mat)
  } else {
    mats$raw_bc <- mat
    mats$rel_bc <- rel_abund(mat)
    mats$hel_bc <- hellinger(mat)
    mats$pa_jc  <- presence_absence(mat)
  }

  # helper: choose best k via silhouette, identical to clustering logic
  best_k <- function(dist_obj, hc_obj) {
    sil_df <- silhouette_curve_hclust(dist_obj, hc_obj)
    if (!is.null(p$k_override)) return(as.integer(p$k_override))
    if (!is.null(sil_df) && nrow(sil_df) > 0) return(as.integer(sil_df$k[which.max(sil_df$mean_silhouette)]))
    2L
  }

  # helper: build a wide qualitative palette for branch colouring.
  # No recycling: request exactly n colours.
  make_cluster_cols <- function(n) {
    n <- as.integer(n)
    if (is.na(n) || n < 1L) stop("n must be a positive integer.")
    pastel_rainbow(n)
  }

  # helper: for presence/absence taxa matrices, cluster unique 0/1 profiles first,
  # then map profile-cluster IDs back to taxa so identical profiles always share a group.
  pa_profile_cluster_ids <- function(m_tax_x_loc_bin, k, hclust_method = "ward.D2") {
    if (!is.matrix(m_tax_x_loc_bin)) m_tax_x_loc_bin <- as.matrix(m_tax_x_loc_bin)
    m_bin <- ifelse(m_tax_x_loc_bin > 0, 1L, 0L)

    keys <- apply(m_bin, 1, function(x) paste0(as.integer(x), collapse = "|"))
    key_to_first <- !duplicated(keys)
    m_unique <- m_bin[key_to_first, , drop = FALSE]
    key_unique <- keys[key_to_first]

    rownames(m_unique) <- key_unique

    if (nrow(m_unique) == 1L) {
      cl_unique <- stats::setNames(1L, key_unique)
    } else {
      d_unique <- vegan::vegdist(m_unique, method = "jaccard", binary = TRUE)
      hc_unique <- stats::hclust(d_unique, method = hclust_method)
      k_use <- max(1L, min(as.integer(k), nrow(m_unique)))
      cl_unique <- stats::cutree(hc_unique, k = k_use)
      names(cl_unique) <- rownames(m_unique)
    }

    cl_taxa <- unname(cl_unique[keys])
    stats::setNames(as.integer(cl_taxa), rownames(m_bin))
  }

  # helper: colour dendrogram edges by leaf groups using a label-based mapping.
  # This avoids order-dependent mismatches in trees with many distance ties.
  # helper: colour dendrogram edges by leaf groups using a label-based mapping.
# This avoids order-dependent mismatches in trees with many distance ties and allows
# different linewidths for single-cluster vs mixed branches (used by ComplexHeatmap).
color_dend_by_membership <- function(dend, membership, palette, mixed_col = "grey40",
                                     lwd_col = 2.5, lwd_default = 1) {
  stopifnot(inherits(dend, "dendrogram"))
  if (is.null(names(membership))) stop("membership must be a named vector (names are leaf labels).")
  if (is.null(names(palette))) stop("palette must be a named vector (names are cluster IDs as character).")

  rec <- function(node) {
    if (is.leaf(node)) {
      lab <- labels(node)
      if (!lab %in% names(membership)) stop("membership missing label: ", lab)
      id <- as.character(membership[[lab]])
      if (!id %in% names(palette)) stop("palette missing cluster id: ", id)
      attr(node, "leaf_ids") <- id
      # edgePar on leaves controls the terminal branch (to the leaf)
      ep <- attr(node, "edgePar")
      if (is.null(ep) || !is.list(ep)) ep <- list()
      ep <- modifyList(ep, list(col = unname(palette[[id]]), lwd = lwd_col))
      attr(node, "edgePar") <- ep
      return(node)
    }

    node[] <- lapply(node, rec)
    ids <- unique(unlist(lapply(node, function(x) attr(x, "leaf_ids"))))
    attr(node, "leaf_ids") <- ids

    col_use <- if (length(ids) == 1L) unname(palette[[ids]]) else mixed_col
    lwd_use <- if (length(ids) == 1L) lwd_col else lwd_default

    ep <- attr(node, "edgePar")
    if (is.null(ep) || !is.list(ep)) ep <- list()
    # IMPORTANT: ComplexHeatmap uses 'lwd' (not 'linewidth') when converting to grid grobs.
    ep <- modifyList(ep, list(col = col_use, lwd = lwd_use))
    attr(node, "edgePar") <- ep
    node
  }

  dend2 <- rec(dend)

  # Clean helper attribute
  clean <- function(node) {
    attr(node, "leaf_ids") <- NULL
    if (!is.leaf(node)) node[] <- lapply(node, clean)
    node
  }
  clean(dend2)
}

  # helper: ensure dendrogram has drawable leaf-adjacent segments by enforcing
  # strictly increasing heights (useful for PA + Jaccard where many merges occur at height 0).
  bump_zero_heights <- function(dend, eps = 1e-6) {
    stopifnot(inherits(dend, "dendrogram"))
    rec <- function(node) {
      if (is.leaf(node)) {
        if (is.null(attr(node, "height"))) attr(node, "height") <- 0
        return(node)
      }
      for (k in seq_along(node)) node[[k]] <- rec(node[[k]])

      ch <- sapply(node, function(x) attr(x, "height"))
      h  <- attr(node, "height")
      if (is.null(h)) h <- max(ch)

      min_ok <- max(ch)
      if (is.na(h) || h <= min_ok) h <- min_ok + eps

      attr(node, "height") <- h
      node
    }
    rec(dend)
  }

  # wrapper used throughout the script: build a dendrogram from an hclust and colour
  # branches based on an explicit leaf->cluster mapping (or cutree if none provided).
  # Uses membership-driven colouring to avoid PA + Jaccard tie artefacts.
  color_cluster_leaf_branches <- function(hc_obj,
                                         k,
                                         leaf_group_ids = NULL,
                                         palette_k = NULL,
                                         mixed_col = "grey70",
                                         bump_zero = FALSE,
                                         eps = 1e-6) {
    stopifnot(inherits(hc_obj, "hclust"))
    dend <- as.dendrogram(hc_obj)
    if (isTRUE(bump_zero)) dend <- bump_zero_heights(dend, eps = eps)

    if (is.null(leaf_group_ids)) {
      leaf_group_ids <- cutree(hc_obj, k = k)
    }

    # align to dendrogram labels
    labs <- labels(dend)
    if (is.null(names(leaf_group_ids))) {
      stop("leaf_group_ids must be a named vector with names matching leaf labels")
    }
    leaf_group_ids <- leaf_group_ids[labs]
    if (anyNA(leaf_group_ids)) {
      missing <- labs[is.na(leaf_group_ids)]
      stop("leaf_group_ids missing these labels: ", paste(missing[1:min(20, length(missing))], collapse = ", "))
    }

    ids <- as.integer(leaf_group_ids)
    id_levels <- sort(unique(ids))

    if (is.null(palette_k)) {
      palette_k <- setNames(grDevices::hcl.colors(length(id_levels), palette = "Dynamic"), as.character(id_levels))
    }

    # ensure palette covers all ids
    if (!all(as.character(id_levels) %in% names(palette_k))) {
      miss <- setdiff(as.character(id_levels), names(palette_k))
      stop("palette_k missing IDs: ", paste(miss, collapse = ", "))
    }

    color_dend_by_membership(dend, leaf_group_ids, palette_k, mixed_col = mixed_col)
  }



# metadata for locality annotations (match clustering behaviour)
  meta_map <- build_metadata_mapping(level = level, use_all_data = use_all_data, plot_points = plot_points, include_groups = isTRUE(p$include_groups))
  metadata_mapping <<- meta_map

  # palettes
  # Sage-green only (no blue): light -> dark
  sage_palette <- c("#f2f5ee", "#dbe6d4", "#a7c2a1", "#5f7f57")
  ls_levels <- c("DLS", "SLS", "PBS")
  ls_cols <- setNames(hull_palette(ls_levels), ls_levels)

  src_levels <- c("Own observations", "Literature")
  src_cols <- setNames(c(sage_palette[3], sage_palette[4]), src_levels)

  out <- list(params = p, heatmaps = list())

  for (nm in names(mats)) {
    m_loc_x_tax <- mats[[nm]]              # rows = sites, cols = taxa
    m_tax_x_loc <- t(m_loc_x_tax)          # rows = taxa, cols = sites (heatmap orientation)

    # drop taxa/sites that are completely zero after transforms
    if (nm == "pa_jc") {
      keep_tax <- rowSums(m_tax_x_loc > 0) > 0
      keep_loc <- colSums(m_tax_x_loc > 0) > 0
    } else {
      keep_tax <- rowSums(m_tax_x_loc, na.rm = TRUE) > 0
      keep_loc <- colSums(m_tax_x_loc, na.rm = TRUE) > 0
    }
    m_tax_x_loc <- m_tax_x_loc[keep_tax, keep_loc, drop = FALSE]
    m_loc_x_tax <- t(m_tax_x_loc)

    # Ensure taxa labels are unique (prevents PA taxa dendrogram colouring issues with duplicated labels)
    taxa <- rownames(m_tax_x_loc)
    taxa_u <- make.unique(taxa, sep = " •dup")
    if (!identical(taxa, taxa_u)) {
      rownames(m_tax_x_loc) <- taxa_u
      m_loc_x_tax <- t(m_tax_x_loc)
    }

    if (nrow(m_tax_x_loc) < 2L || ncol(m_tax_x_loc) < 2L) {
      warning(sprintf("Heatmap skipped (%s): matrix is %d taxa x %d sites after filtering (need >=2 each).", nm, nrow(m_tax_x_loc), ncol(m_tax_x_loc)))
      next
    }

    # distances + clustering parameters aligned with clustering analysis
    dist_method <- if (nm == "pa_jc") "jaccard" else "bray"
    binary_flag <- (nm == "pa_jc")

    dist_tax <- vegan::vegdist(m_tax_x_loc, method = dist_method, binary = binary_flag)
    hc_tax <- stats::hclust(dist_tax, method = "ward.D2")
    k_tax <- best_k(dist_tax, hc_tax)

    # Build row dendrogram and colour edges by membership (label-based mapping).
    # For PA+Jaccard, bump zero-height merges so leaf-adjacent segments are drawable.
    if (nm == "pa_jc") {
      # identical 0/1 profiles receive the same cluster id
      membership_tax <- pa_profile_cluster_ids(m_tax_x_loc, k = k_tax, hclust_method = "ward.D2")

      ids <- sort(unique(as.integer(membership_tax)))
      palette_k <- setNames(make_cluster_cols(length(ids)), as.character(ids))

      row_dend <- as.dendrogram(hc_tax)
      row_dend <- bump_zero_heights(row_dend, eps = 1e-6)
      row_dend <- color_dend_by_membership(row_dend, membership_tax, palette_k, mixed_col = "grey70")
    } else {
      membership_tax <- stats::cutree(hc_tax, k = k_tax)

      ids <- sort(unique(as.integer(membership_tax)))
      palette_k <- setNames(make_cluster_cols(length(ids)), as.character(ids))

      row_dend <- as.dendrogram(hc_tax)
      row_dend <- color_dend_by_membership(row_dend, membership_tax, palette_k, mixed_col = "grey70")
    }


    dist_loc <- vegan::vegdist(m_loc_x_tax, method = dist_method, binary = binary_flag)
    hc_loc <- stats::hclust(dist_loc, method = "ward.D2")
    k_loc <- best_k(dist_loc, hc_loc)
    membership_loc <- stats::cutree(hc_loc, k = k_loc)
    pal_loc <- cluster_palette_named(k_loc)
    col_dend <- color_cluster_leaf_branches(
      hc_loc,
      k_loc,
      leaf_group_ids = membership_loc,
      palette_k = pal_loc,
      mixed_col = "grey70",
      bump_zero = (nm == "pa_jc")
    )

    # locality annotations
    loc_ids <- colnames(m_tax_x_loc)
    meta_sub <- meta_map %>%
      dplyr::filter((if (plot_points == "Basin") Basin else Locality) %in% loc_ids) %>%
      dplyr::distinct(Basin, Locality, LakeSystem, Source)

    key_col <- if (plot_points == "Basin") "Basin" else "Locality"
    meta_sub <- meta_sub %>% dplyr::mutate(Key = .data[[key_col]])

    loc_group  <- setNames(as.character(meta_sub$LakeSystem), meta_sub$Key)[loc_ids]
    loc_source <- setNames(as.character(meta_sub$Source),     meta_sub$Key)[loc_ids]

    taxa_present_n <- colSums(m_tax_x_loc > 0)

    top_ha <- ComplexHeatmap::HeatmapAnnotation(
      Richness = ComplexHeatmap::anno_text(
        taxa_present_n,
        rot = 60,
        gp = grid::gpar(fontsize = 8),
        location = 0.15,
        just = "left"
      ),
      show_annotation_name = FALSE
    )

    ann_type <- match.arg(as.character(p$annotation_type), c("bar", "shape"))

    if (ann_type == "bar") {
      # Bar annotations with per-locality dotted separators and a surrounding border.
      bottom_ha <- ComplexHeatmap::HeatmapAnnotation(
        LakeSystem = ComplexHeatmap::anno_simple(
          loc_group,
          col = ls_cols,
      border = (ann_type == "bar"),
          gp = grid::gpar(col = "black", lty = "dotted", lwd = 0.5)
        ),
        Source = ComplexHeatmap::anno_simple(
          loc_source,
          col = src_cols,
          border = (ann_type == "bar"),
          gp = grid::gpar(col = "black", lty = "dotted", lwd = 0.5)
        ),
        show_annotation_name = FALSE,
        annotation_height = grid::unit.c(
          grid::unit(as.numeric(p$lake_system_bar_height_mm), "mm"),
          grid::unit(as.numeric(p$source_bar_height_mm), "mm")
        ),
        show_legend = FALSE,
        border = FALSE
      )
    } else {
      # Shape annotation: point colour = LakeSystem, point shape = Source.
      src_chr <- as.character(loc_source)
      pch_vec <- ifelse(src_chr == "Own observations", 16, 17)
      col_vec <- unname(ls_cols[as.character(loc_group)])
      col_vec[is.na(col_vec)] <- "grey60"

      bottom_ha <- ComplexHeatmap::HeatmapAnnotation(
        LocalityInfo = ComplexHeatmap::anno_points(
          x = rep(1, length(loc_ids)),
          pch = pch_vec,
          gp = grid::gpar(col = col_vec),
          size = grid::unit(2.2, "mm"),
          axis = FALSE,
          border = FALSE
        ),
        show_annotation_name = FALSE,
        annotation_height = grid::unit(max(as.numeric(p$lake_system_bar_height_mm), as.numeric(p$source_bar_height_mm)), "mm"),
        show_legend = FALSE,
        border = FALSE,
        gp = grid::gpar(col = NA, fill = NA)
      )
    }

    # heatmap colors
    if (nm == "pa_jc") {
      mat_plot <- ifelse(m_tax_x_loc > 0, 1, 0)
      col_fun <- c("0" = sage_palette[1], "1" = sage_palette[4])
      lgd_main <- ComplexHeatmap::Legend(
        title = "Occurrence",
        at = c("0", "1"),
        labels = c("Absent", "Present"),
        legend_gp = grid::gpar(fill = unname(col_fun[c("0", "1")])),
        grid_height = grid::unit(4, "mm")
      )
      hm_name <- "Occurrence"
    } else {
      mat_plot <- m_tax_x_loc
      rng <- range(as.vector(mat_plot), finite = TRUE)
      if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) rng <- c(0, 1)
      mid <- mean(rng)
      col_fun <- circlize::colorRamp2(c(rng[1], mid, rng[2]), c(sage_palette[1], sage_palette[2], sage_palette[4]))
      at_vals <- unique(round(quantile(as.vector(mat_plot), probs = c(0, 0.5, 1), na.rm = TRUE), 3))
      lgd_main <- ComplexHeatmap::Legend(title = "Abundance", col_fun = col_fun, at = at_vals)
      hm_name <- "Abundance"
    }
    
    ht <- ComplexHeatmap::Heatmap(
      mat_plot,
      name = hm_name,
      col = col_fun,
      top_annotation = top_ha,
      bottom_annotation = bottom_ha,
      cluster_rows = row_dend,
      cluster_columns = col_dend,
      row_dend_reorder = FALSE,
      column_dend_reorder = FALSE,
      row_dend_width = grid::unit(35, "mm"),
      border = (ann_type == "bar"),
      layer_fun = function(j, i, x, y, w, h, fill) {
        grid::grid.rect(x, y, w, h, gp = grid::gpar(col = "#c9d3c0", fill = NA, lwd = 0.4, lty = "dotted"))
        # For quantitative heatmaps, print abundance values inside each cell.
        if (nm != "pa_jc") {
          v <- mat_plot[cbind(i, j)]  # vector of values aligned to i/j pairs
          
          ok <- !is.na(v) & (v > 0)
          if (any(ok)) {
            lab <- character(length(v))
            if (nm == "raw_bc") {
              lab[ok] <- as.character(as.integer(round(v[ok])))
            } else {
              lab[ok] <- formatC(v[ok], format = "f", digits = 2)
            }
            grid::grid.text(lab[ok], x[ok], y[ok], gp = grid::gpar(fontsize = 6))
          }
        }
      },
      row_dend_side = "left",
      row_names_side = "left",
      column_dend_side = "bottom",
      column_names_side = "bottom",
      show_column_names = TRUE,
      column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = 6),
      row_names_gp = grid::gpar(fontsize = 7, fontface = "italic"),
      show_heatmap_legend = FALSE
    )

    # legends for annotations
    ls_levels_present <- sort(unique(stats::na.omit(as.character(loc_group))))
    ls_levels_present <- intersect(ls_levels_present, names(ls_cols))
    src_levels_present <- sort(unique(stats::na.omit(as.character(loc_source))))
    src_levels_present <- intersect(src_levels_present, names(src_cols))

    # If the plot contains only one source level (e.g. own-only heatmaps),
    # omit the "Data source" legend entirely.
    show_source_legend <- length(src_levels_present) > 1L

    if (ann_type == "bar") {
      lgd_sys <- ComplexHeatmap::Legend(
        title = "Lake system",
        at = ls_levels_present,
        legend_gp = grid::gpar(fill = unname(ls_cols[ls_levels_present]))
      )
      lgd_src <- if (show_source_legend) {
        ComplexHeatmap::Legend(
          title = "Data source",
          at = src_levels_present,
          legend_gp = grid::gpar(fill = unname(src_cols[src_levels_present]))
        )
      } else {
        NULL
      }
    } else {
      lgd_sys <- ComplexHeatmap::Legend(
        title = "Lake system",
        at = ls_levels_present,
        type = "points",
        pch = rep(16, length(ls_levels_present)),
        legend_gp = grid::gpar(col = unname(ls_cols[ls_levels_present]), fill = NA))
      lgd_src <- if (show_source_legend) {
        ComplexHeatmap::Legend(
          title = "Data source",
          at = src_levels_present,
          type = "points",
          pch = unname(c("Own observations" = 16, "Literature" = 17)[src_levels_present]),
          legend_gp = grid::gpar(col = "grey30", fill = NA)
        )
      } else {
        NULL
      }
    }

    fn <- file.path(CFG$out_dir_heatmap, paste0(run_params$name, "_heatmap_", nm, "_minTaxa", p$min_taxa_present, ".tiff"))

    height_in <- max(p$min_height_in, min(p$max_height_in, p$base_height_in + p$row_height_in * nrow(m_tax_x_loc)))

    ragg::agg_tiff(
      filename = fn,
      width = 2400, height = max(1800, round(height_in * 300)),
      units = "px",
      res = 300,
      compression = "lzw"
    )
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_list = list(lgd_main),
      annotation_legend_list = compact(list(lgd_sys, lgd_src)),
      merge_legends = TRUE,
      padding = grid::unit(c(6, 6, 6, 6), "mm")
    )
    ComplexHeatmap::decorate_heatmap_body(hm_name, {
      grid::grid.rect(x = 0.5, y = 0.5, width = 1, height = 1, gp = grid::gpar(col = "grey30", fill = NA, lwd = 1))
      grid::grid.text(
        "Taxa richness",
        x = grid::unit(-22, "mm"),
        y = grid::unit(1, "npc") + grid::unit(2, "mm"),
        just = c("left", "bottom"),
        gp = grid::gpar(fontface = "bold", fontsize = 9)
      )
    })
    dev.off()

    out$heatmaps[[nm]] <- list(file = fn, heatmap = ht, k_tax = k_tax, k_loc = k_loc, height_in = height_in)
  }

  invisible(out)
}




# -----------------------------
# 8) RUN EVERYTHING
# -----------------------------

set.seed(CFG$defaults$seed)

results <- list(clustering = list(), heatmap = list())

for (r in CFG$runs$clustering) {
  results$clustering[[r$name %||% paste0("run_", length(results$clustering)+1)]] <- run_clustering(r)
}

for (r in CFG$runs$heatmap) {
  results$heatmap[[r$name %||% paste0("heatmap_", length(results$heatmap)+1)]] <- save_heatmaps(r)
}

results

