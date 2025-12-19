# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr, dplyr, tidyr, lubridate, purrr,
    ggplot2, ggtext, ggsankey, RColorBrewer,
    grid, gridExtra, Cairo, openxlsx, stringr
  )
})

# Load data
wd <- "./"; setwd(wd)
df_data     <- readRDS("data/clean data/anti_treat_index.RData")
df_baseline <- readRDS("data/clean data/baseline_outcomes_index.RData")

# Tidy baseline / remove antifungals
df_data <- df_data[!df_data$anti_group %in% c("Azole","Polyene","Echinocandin"),]

# anti_used (only CRA + sulbactam overrides; CRE Ceftazidime/avibactam; otherwise keep anti_group)
df_data <- df_data %>%
  mutate(
    anti_used = case_when(
      # CRA & sulbactam
      aci_car == 1 & str_detect(tolower(anti_names), "sulbactam") ~ "Sulbactam",
      
      # CRE & exact CAZ-AVI
      ent_car == 1 & str_detect(anti_names, regex("^Ceftazidime/avibactam$", ignore_case = TRUE)) ~ "Ceftazidime/avibactam",
      
      # CRE & other anti-pseudomonal PIP/TAZ group
      ent_car == 1 & 
        !str_detect(anti_names, regex("^Ceftazidime/avibactam$", ignore_case = TRUE)) &
        anti_group == "Anti-pseudomonal penicillin/Beta-lactamase inhibitor" ~ 
        "Other Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
      
      # Default
      TRUE ~ anti_group
    )
  )

# Per-pathogen priority for anti_used (full names, no abbreviations)
priority_map <- list(
  CRA   = c("Sulbactam","Polymyxin","Carbapenem",
            "Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
            "Other Beta-lactam/Beta-lactamase inhibitor",
            "Aminoglycoside","Third-generation cephalosporin",
            "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  `3GCRE` = c("Carbapenem","Polymyxin",
              "Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
              "Other Beta-lactam/Beta-lactamase inhibitor",
              "Aminoglycoside","Third-generation cephalosporin",
              "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  CRE   = c("Ceftazidime/avibactam","Polymyxin","Carbapenem",
            "Other Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
            "Other Beta-lactam/Beta-lactamase inhibitor",
            "Aminoglycoside","Third-generation cephalosporin",
            "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  CRP   = c("Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
            "Polymyxin","Carbapenem",
            "Other Beta-lactam/Beta-lactamase inhibitor","Aminoglycoside"),
  VRE   = c("Lipopeptide","Oxazolidinone","Penicillin",
            "Third-generation cephalosporin","Phosphonic"),
  MRSA  = c("Glycopeptide","Lipopeptide","Oxazolidinone",
            "Sulfonamide-trimethoprim-combination","Lincosamide")
)

# Build a complete priority vector for a subset (append unseen classes at the end)
priority_all_for <- function(key, pool_values) {
  base <- priority_map[[key]]
  others <- setdiff(sort(unique(pool_values)), base)
  c(base, others)
}

# Helpers to combine classes on the same day
combine_same_day <- function(v, priority=NULL){
  v <- unique(na.omit(v))
  if (!length(v)) return(NA_character_)
  if (is.null(priority) || !length(priority)) {
    return(paste(sort(v), collapse=", "))
  }
  ord <- match(v, priority)
  max_ord <- suppressWarnings(max(ord, na.rm = TRUE))
  if (!is.finite(max_ord)) max_ord <- 0
  ord[is.na(ord)] <- max_ord + seq_len(sum(is.na(ord)))
  v <- v[order(ord)]
  paste(v, collapse=", ") %>% gsub("\\s*-\\s*", "-", .) %>% trimws()
}

# Expand to patient-day, then select window days around ref dates
daily_antibiotics <- function(df,
                              id_col="recordid",
                              group_col="anti_used",
                              start_col="anti_start",
                              end_col  ="anti_end",
                              windows=list(
                                emp=list(ref_col="inf_onset", offsets=c(-1,0,1,2)),
                                def=list(ref_col="spec_date", offsets=c(3,4,5))
                              ),
                              same_day_priority=NULL){
  stopifnot(all(c(id_col,group_col,start_col,end_col) %in% names(df)))
  to_date <- function(x) as.Date(x)
  df <- df %>% mutate(
    !!start_col := to_date(.data[[start_col]]),
    !!end_col   := to_date(.data[[end_col]])
  )
  rx_base <- df %>% distinct(.data[[id_col]], .data[[group_col]], .data[[start_col]], .data[[end_col]])
  anti_by_day <- rx_base %>%
    filter(!is.na(.data[[start_col]]), !is.na(.data[[end_col]]), .data[[end_col]] >= .data[[start_col]]) %>%
    rowwise() %>%
    mutate(day_date = list(seq(.data[[start_col]], .data[[end_col]], by="day"))) %>%
    unnest(day_date) %>% ungroup() %>%
    group_by(.data[[id_col]], day_date) %>%
    summarise(anti = combine_same_day(.data[[group_col]], priority = same_day_priority), .groups="drop")
  
  mk_grid <- function(ref_col, offsets, prefix){
    ids <- df %>%
      select(all_of(c(id_col, ref_col))) %>% distinct() %>%
      mutate(!!ref_col := to_date(.data[[ref_col]])) %>% filter(!is.na(.data[[ref_col]]))
    ids %>%
      rowwise() %>%
      mutate(.dates = list(.data[[ref_col]] + lubridate::days(offsets)), .offs = list(offsets)) %>%
      unnest(c(.dates, .offs)) %>%
      ungroup() %>%
      mutate(win = prefix, rel = .offs,
             label = paste0(prefix, ifelse(.offs >= 0, paste0("+", .offs), .offs))) %>%
      rename(day_date = .dates, ref_date = !!ref_col)
  }
  grids <- purrr::imap(windows, ~ mk_grid(.x$ref_col, .x$offsets, .y)) %>% bind_rows()
  daily_long <- grids %>%
    left_join(anti_by_day, by = setNames(c(id_col,"day_date"), c(id_col,"day_date"))) %>%
    mutate(anti = tidyr::replace_na(anti, "None")) %>%
    select(all_of(c(id_col)), day_date, win, label, rel, anti)
  
  col_name <- function(win, rel){
    if (win=="emp"){
      if (rel<0) paste0("emp_m", abs(rel))
      else if (rel==0) "emp_0" else paste0("emp_p", rel)
    } else paste0("def_d", rel)
  }
  daily_wide <- daily_long %>%
    mutate(col = pmap_chr(list(win,rel), col_name)) %>%
    select(all_of(c(id_col)), col, anti) %>%
    distinct() %>%
    pivot_wider(names_from = col, values_from = anti)
  
  list(long=daily_long, wide=daily_wide)
}


# Outcomes (Def window: all None -> Outcome)
df_outcomes <- df_baseline %>%
  select(recordid, ho_discharge_date, mortality_date) %>%
  mutate(outcomes = dplyr::case_when(
    !is.na(mortality_date) ~ "Died",
    is.na(ho_discharge_date) & is.na(mortality_date) ~ "None",
    TRUE ~ "Discharged"
  )) %>%
  select(recordid, outcomes)

# Organism recordid subsets
aci_car_id  <- df_data$recordid[which(df_data$aci_car == 1)]
ent_thir_id <- df_data$recordid[which(df_data$ent_thir == 1)]
ent_car_id  <- df_data$recordid[which(df_data$ent_car == 1)]
pse_car_id  <- df_data$recordid[which(df_data$pse_car == 1)]
entc_van_id <- df_data$recordid[which(df_data$entc_van == 1)]
sa_meth_id  <- df_data$recordid[which(df_data$sa_meth == 1)]

subset_list <- list(
  list(key="CRA", ids = aci_car_id, title = "Carbapenem-resistant Acinetobacter spp.", space=30, height=20),
  list(key="3GCRE", ids = ent_thir_id, title = "Third-generation cephalosporin-resistant Enterobacterales", space=30, height=20),
  list(key="CRE",   ids = ent_car_id,  title = "Carbapenem-resistant Enterobacterales", space=30, height=18),
  list(key="CRP",   ids = pse_car_id,  title = "Carbapenem-resistant Pseudomonas spp.",  space=5, height=20),
  list(key="VRE",   ids = entc_van_id, title = "Vancomycin-resistant Enterococcus spp.",  space=2,  height=16),
  list(key="MRSA",  ids = sa_meth_id,  title = "Methicillin-resistant Staphylococcus aureus", space=5, height=16)
)

# Priority-based transformers
normalize_daily_local_factory <- function(priority_vec){
  force(priority_vec)
  function(s){
    s <- ifelse(is.na(s) | s=="", "None", s)
    if (s == "None") return("None")
    tokens <- strsplit(s, "\\s*,\\s*")[[1]]
    tokens <- gsub("\\s*-\\s*", "-", tokens)
    ord  <- match(tokens, priority_vec)
    maxo <- suppressWarnings(max(ord, na.rm=TRUE)); if(!is.finite(maxo)) maxo <- 0
    ord[is.na(ord)] <- maxo + seq_len(sum(is.na(ord)))
    paste(tokens[order(ord)], collapse = ", ")
  }
}

pick_base_by_priority_factory <- function(priority_vec){
  force(priority_vec)
  function(s){
    vapply(s, function(x){
      if (is.na(x) || x %in% c("None","Died","Discharged")) return(x)
      toks <- strsplit(x, "\\s*,\\s*")[[1]]
      toks <- gsub("\\s*-\\s*", "-", toks)
      ord  <- match(toks, priority_vec)
      maxo <- suppressWarnings(max(ord, na.rm=TRUE)); if(!is.finite(maxo)) maxo <- 0
      ord[is.na(ord)] <- maxo + seq_len(sum(is.na(ord)))
      toks[which.min(ord)]
    }, character(1))
  }
}

norm_combo_vec_factory <- function(priority_vec){
  force(priority_vec)
  function(s){
    vapply(s, function(x){
      if (is.na(x)) return(NA_character_)
      toks <- strsplit(x, "\\s*,\\s*")[[1]]
      toks <- gsub("\\s*-\\s*", "-", toks)
      ord  <- match(toks, priority_vec)
      maxo <- suppressWarnings(max(ord, na.rm=TRUE)); if(!is.finite(maxo)) maxo <- 0
      ord[is.na(ord)] <- maxo + seq_len(sum(is.na(ord)))
      paste(toks[order(ord)], collapse = ", ")
    }, character(1))
  }
}


# Build per-subset window sets
build_window_sets_for_ids <- function(id_vec, priority_vec){
  res_local <- daily_antibiotics(
    df = df_data %>% dplyr::filter(recordid %in% id_vec),
    id_col="recordid",
    group_col="anti_used",
    start_col="anti_start",
    end_col  ="anti_end",
    windows = list(
      emp=list(ref_col="inf_onset", offsets=c(-1,0,1,2)),
      def=list(ref_col="spec_date", offsets=c(3,4,5))
    ),
    same_day_priority = priority_vec
  )
  emp_def_wide_local <- res_local$wide
  emp_cols <- intersect(c("emp_m1","emp_0","emp_p1","emp_p2"), names(emp_def_wide_local))
  def_cols <- intersect(c("def_d3","def_d4","def_d5"), names(emp_def_wide_local))
  emp_def_wide_local <- emp_def_wide_local %>%
    mutate(across(all_of(c(emp_cols, def_cols)), as.character))
  normalize_daily_local <- normalize_daily_local_factory(priority_vec)
  
  # Empirical window set
  emp_final_local <- {
    long <- emp_def_wide_local %>%
      select(recordid, all_of(emp_cols)) %>%
      pivot_longer(cols = all_of(emp_cols), names_to="day", values_to="treatment") %>%
      mutate(treatment = vapply(treatment, normalize_daily_local, character(1)))
    all_none <- long %>% group_by(recordid) %>% summarise(all_none = all(treatment=="None"), .groups="drop")
    keep_non_none <- long %>% inner_join(all_none, by="recordid") %>%
      filter(!all_none, treatment!="None") %>%
      distinct(recordid, treatment)
    keep_none <- all_none %>% filter(all_none) %>%
      transmute(recordid, treatment="None")
    bind_rows(keep_non_none, keep_none) %>% rename(emp_treatment = treatment)
  }
  
  # Definitive window set (None -> outcome)
  def_final_local <- {
    long <- emp_def_wide_local %>%
      select(recordid, all_of(def_cols)) %>%
      pivot_longer(cols = all_of(def_cols), names_to="day", values_to="treatment") %>%
      mutate(treatment = vapply(treatment, normalize_daily_local, character(1)))
    all_none <- long %>% group_by(recordid) %>% summarise(all_none = all(treatment=="None"), .groups="drop")
    keep_non_none <- long %>% inner_join(all_none, by="recordid") %>%
      filter(!all_none, treatment!="None") %>%
      distinct(recordid, treatment) %>% rename(def_treatment = treatment)
    keep_outcome <- all_none %>% filter(all_none) %>%
      left_join(df_outcomes, by="recordid") %>%
      mutate(def_treatment = ifelse(is.na(outcomes), "None", outcomes)) %>%
      select(recordid, def_treatment) %>% distinct()
    bind_rows(keep_non_none, keep_outcome)
  }
  
  list(emp_final = emp_final_local,
       def_final = def_final_local,
       emp_def_wide = emp_def_wide_local)
}


# Pair-level exclusion
apply_pair_exclusion <- function(emp_df, def_df) {
  special <- c("None","Died","Discharged")
  joined <- emp_df %>%
    full_join(def_df, by = "recordid", relationship = "many-to-many") %>%
    mutate(
      emp_treatment = if_else(is.na(emp_treatment), "None", emp_treatment),
      def_treatment = if_else(is.na(def_treatment), "None", def_treatment)
    ) %>%
    # Drop Emp=None with Def in {None, Died, Discharged}
    filter(!(emp_treatment == "None" & def_treatment %in% special))
  
  list(
    joined = joined,
    emp = joined %>% distinct(recordid, emp_treatment),
    def = joined %>% distinct(recordid, def_treatment)
  )
}

# Global container for special nodes export
special_counts_global <- list()

# Sankey plotting with <5% -> Others + symmetric exception, and fixed totals
plot_sankey_based_for_ids <- function(id_vec, main_title, subset_key, space_value=40){
  special <- c("None","Died","Discharged")
  priority_vec <- priority_all_for(subset_key, pool_values = df_data$anti_used)
  
  ws <- build_window_sets_for_ids(id_vec, priority_vec)
  emp_final_local <- ws$emp_final
  def_final_local <- ws$def_final
  
  pick_base_by_priority_vec <- pick_base_by_priority_factory(priority_vec)
  norm_combo_vec  <- norm_combo_vec_factory(priority_vec)
  
  emp_pairs_raw <- emp_final_local %>% distinct(recordid, emp_treatment)
  def_pairs_raw <- def_final_local %>% distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  emp_pairs <- filt$emp
  def_pairs <- filt$def
  joined    <- filt$joined
  
  # Counts before merging
  mk_rx_counts_side <- function(recordids, trts){
    tibble(recordid = recordids, treatment = trts) %>%
      mutate(
        is_special = treatment %in% special,
        is_combo   = !is_special & grepl(",", treatment),
        node_name  = dplyr::case_when(
          # Special nodes (None / Died / Discharged)
          is_special ~ treatment,
          
          # CRE: any treatment containing Ceftazidime/avibactam
          # (monotherapy or in combination) is grouped into one node
          subset_key == "CRE" & grepl("Ceftazidime/avibactam", treatment, ignore.case = TRUE) ~
            "Ceftazidime/avibactam monotherapy or combinations",
          
          # Other combinations: base + "-based combinations"
          is_combo   ~ paste0(pick_base_by_priority_vec(treatment), "-based combinations"),
          
          # Single agents (non-special)
          TRUE       ~ treatment
        ),
        combo_sig  = ifelse(is_combo, norm_combo_vec(treatment), NA_character_)
      ) %>%
      distinct(recordid, node_name, combo_sig) %>%
      mutate(rx = 1L) %>%
      group_by(node_name) %>%
      summarise(rx = sum(rx), .groups = "drop")
  }
  
  emp_rx_tbl0 <- mk_rx_counts_side(emp_pairs$recordid, emp_pairs$emp_treatment)
  def_rx_tbl0 <- mk_rx_counts_side(def_pairs$recordid, def_pairs$def_treatment)
  
  # Log and stash special counts
  cat("Empirical special node counts:\n")
  print(emp_rx_tbl0 %>% dplyr::filter(node_name %in% special))
  cat("Definitive special node counts:\n")
  print(def_rx_tbl0 %>% dplyr::filter(node_name %in% special))
  special_counts_global[[paste0(subset_key, "_Empirical")]]  <- emp_rx_tbl0 %>%
    dplyr::filter(node_name %in% special) %>%
    dplyr::mutate(window = "Empirical", subset = subset_key)
  special_counts_global[[paste0(subset_key, "_Definitive")]] <- def_rx_tbl0 %>%
    dplyr::filter(node_name %in% special) %>%
    dplyr::mutate(window = "Definitive", subset = subset_key)
  
  # <5% -> Others rule (with symmetric exception)
  pct_threshold <- 0.05
  total_emp_rx <- sum(emp_rx_tbl0$rx[!emp_rx_tbl0$node_name %in% special], na.rm = TRUE)
  total_def_rx <- sum(def_rx_tbl0$rx[!def_rx_tbl0$node_name %in% special], na.rm = TRUE)
  
  emp_keep <- emp_rx_tbl0 %>%
    filter(!node_name %in% special) %>%
    mutate(pct = rx / ifelse(total_emp_rx > 0, total_emp_rx, 1)) %>%
    filter(pct > pct_threshold) %>%
    pull(node_name)
  
  def_keep_initial <- def_rx_tbl0 %>%
    filter(!node_name %in% special) %>%
    mutate(pct = rx / ifelse(total_def_rx > 0, total_def_rx, 1)) %>%
    filter(pct > pct_threshold) %>%
    pull(node_name)
  
  # Exceptions:
  # - Emp-kept must be kept on Def
  def_keep <- union(def_keep_initial, emp_keep)
  # - Def-kept must be kept on Emp (symmetric)
  emp_keep <- union(emp_keep, def_keep_initial)
  
  remap_emp <- function(nm) {
    if (nm %in% special) return(nm)
    
    # For CRE: always keep CAZ-AVI as its own node
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", nm, ignore.case = TRUE)) {
      return("Ceftazidime/avibactam monotherapy or combinations")
    }
    
    if (nm %in% emp_keep) nm else "Others"
  }
  
  remap_def <- function(nm) {
    if (nm %in% special) return(nm)
    
    # For CRE: always keep CAZ-AVI as its own node
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", nm, ignore.case = TRUE)) {
      return("Ceftazidime/avibactam monotherapy or combinations")
    }
    
    if (nm %in% def_keep) nm else "Others"
  }
  
  # Build raw node names from treatments (before adding counts)
  build_node_name <- function(t, side = c("emp","def")){
    side <- match.arg(side)
    is_special <- t %in% special
    is_combo   <- !is_special & grepl(",", t)
    if (is_special) return(t)
    if (is_combo) paste0(pick_base_by_priority_vec(t), "-based combinations") else t
  }
  joined_nodes <- joined %>%
    transmute(
      emp_node_raw = vapply(emp_treatment, build_node_name, character(1), side="emp"),
      def_node_raw = vapply(def_treatment, build_node_name, character(1), side="def")
    )
  
  # --- Build links (kept) ---
  joined_collapsed <- joined_nodes %>%
    mutate(
      emp_node = vapply(emp_node_raw, remap_emp, character(1)),
      def_node = vapply(def_node_raw, remap_def, character(1))
    ) %>%
    count(emp_node, def_node, name = "n")
  
  # --- SIDE COUNTS preserving combo signatures and Others granularity ---
  make_side_counts <- function(pairs_df, side = c("emp","def")) {
    side <- match.arg(side)
    raw_treat <- if (side == "emp") pairs_df$emp_treatment else pairs_df$def_treatment
    node_raw  <- vapply(raw_treat, build_node_name, character(1), side = side)
    is_combo  <- grepl("-based combinations$", node_raw)
    combo_sig <- ifelse(is_combo, norm_combo_vec(raw_treat), NA_character_)
    node_final <- if (side == "emp") {
      vapply(node_raw, remap_emp, character(1))
    } else {
      vapply(node_raw, remap_def, character(1))
    }
    # signature for dedup
    sig_final <- ifelse(
      node_final == "Others",
      ifelse(is_combo, paste0("C:", combo_sig), paste0("S:", node_raw)),
      ifelse(is_combo, paste0("C:", combo_sig), "S")
    )
    tibble(recordid = pairs_df$recordid,
           node_name = node_final,
           sig_final = sig_final) %>%
      distinct(recordid, node_name, sig_final) %>%
      count(node_name, name = "rx")
  }
  emp_rx_tbl <- make_side_counts(emp_pairs, "emp")
  def_rx_tbl <- make_side_counts(def_pairs, "def")
  
  # Build display labels
  emp_rx_map <- setNames(paste0(emp_rx_tbl$node_name, " (", emp_rx_tbl$rx, ")"), emp_rx_tbl$node_name)
  def_rx_map <- setNames(paste0(def_rx_tbl$node_name, " (", def_rx_tbl$rx, ")"), def_rx_tbl$node_name)
  
  get_base <- function(nm) sub("-based combinations$", "", nm)
  emp_based_only <- emp_rx_tbl %>% filter(grepl("-based combinations$", node_name))
  def_based_only <- def_rx_tbl %>% filter(grepl("-based combinations$", node_name))
  emp_based_map  <- setNames(paste0(emp_based_only$node_name, " (", emp_based_only$rx, ")"),
                             get_base(emp_based_only$node_name))
  def_based_map  <- setNames(paste0(def_based_only$node_name, " (", def_based_only$rx, ")"),
                             get_base(def_based_only$node_name))
  
  safe_pick <- function(map, key, fallback) {
    out <- unname(map[key]); if (is.null(out) || is.na(out)) fallback else out
  }
  
  label_emp <- function(t) {
    # CRE: force all CAZ-AVI (mono + combos) to use the merged label
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", t, ignore.case = TRUE)) {
      return(
        safe_pick(
          emp_rx_map,
          "Ceftazidime/avibactam monotherapy or combinations",
          "Ceftazidime/avibactam monotherapy or combinations (0)"
        )
      )
    }
    
    if (t %in% special) return(t)
    if (grepl(",", t))
      return(safe_pick(emp_based_map, pick_base_by_priority_vec(t),
                       paste0(pick_base_by_priority_vec(t), "-based combinations (0)")))
    safe_pick(emp_rx_map, t, paste0(t, " (0)"))
  }
  
  label_def <- function(t) {
    # CRE: same logic for definitive window
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", t, ignore.case = TRUE)) {
      return(
        safe_pick(
          def_rx_map,
          "Ceftazidime/avibactam monotherapy or combinations",
          "Ceftazidime/avibactam monotherapy or combinations (0)"
        )
      )
    }
    
    if (t %in% special) return(t)
    if (grepl(",", t))
      return(safe_pick(def_based_map, pick_base_by_priority_vec(t),
                       paste0(pick_base_by_priority_vec(t), "-based combinations (0)")))
    safe_pick(def_rx_map, t, paste0(t, " (0)"))
  }
  
  links_named <- joined_collapsed %>%
    transmute(
      emp_node = vapply(emp_node, label_emp, character(1)),
      def_node = vapply(def_node, label_def, character(1)),
      n
    )
  links_rep <- tidyr::uncount(links_named, weights = n) %>%
    mutate(emp_node = as.character(emp_node), def_node = as.character(def_node))
  
  df_long <- ggsankey::make_long(links_rep, emp_node, def_node)
  
  # Visible node bases
  strip_paren <- function(s) sub(" \\(.*\\)$", "", s)
  node_base <- function(s) { s <- sub(" \\(.*\\)$", "", s); sub("-based combinations$", "", s) }
  df_long$node_base <- node_base(as.character(df_long$node))
  
  # Totals for side headers (exclude special nodes)
  emp_node_rx_map <- setNames(emp_rx_tbl$rx, emp_rx_tbl$node_name)
  def_node_rx_map <- setNames(def_rx_tbl$rx, def_rx_tbl$node_name)
  visible_emp_nodes <- unique(strip_paren(as.character(df_long$node[df_long$x == "emp_node"])))
  visible_def_nodes <- unique(strip_paren(as.character(df_long$node[df_long$x == "def_node"])))
  N_emp_rx <- sum(emp_node_rx_map[setdiff(visible_emp_nodes, special)], na.rm = TRUE)
  N_def_rx <- sum(def_node_rx_map[setdiff(visible_def_nodes, special)], na.rm = TRUE)
  
  # Title (n patients)
  N_patient <- length(unique(id_vec))
  make_title_html <- function(main_title, n){
    title_core <- main_title
    if (grepl("Acinetobacter", main_title, ignore.case=TRUE)) {
      title_core <- "Carbapenem-resistant <i>Acinetobacter</i> spp."
    } else if (grepl("Pseudomonas", main_title, ignore.case=TRUE)) {
      title_core <- "Carbapenem-resistant <i>Pseudomonas</i> spp."
    } else if (grepl("Enterococcus", main_title, ignore.case=TRUE)) {
      title_core <- "Vancomycin-resistant <i>Enterococcus</i> spp."
    } else if (grepl("Staphylococcus aureus", main_title, ignore.case=TRUE)) {
      title_core <- "Methicillin-resistant <i>Staphylococcus aureus</i>"
    }
    sprintf("<b>%s (n = %s)</b>", title_core, format(n, big.mark = "", scientific = FALSE))
  }
  title_html <- make_title_html(main_title, N_patient)
  
  # Node ordering
  node_chr <- unique(as.character(df_long$node))
  
  # Special nodes
  special_nodes <- c(
    grep("^None( \\(|$)", node_chr, value = TRUE),
    grep("^Died( \\(|$)", node_chr, value = TRUE),
    grep("^Discharged( \\(|$)", node_chr, value = TRUE),
    grep("^Others( \\(|$)", node_chr, value = TRUE)
  )
  
  # Helpers
  strip_counts <- function(s) sub(" \\(.*\\)$", "", s)
  strip_combo  <- function(s) sub("-based combinations$", "", s)
  
  # Exact match by label while ignoring counts
  pick_nodes_exact <- function(label, pool) {
    pool_no_counts <- strip_counts(pool)
    pool[pool_no_counts == label]
  }
  
  # Custom display orders (no counts)
  cra_order <- rev(c(
    "Sulbactam-based combinations", "Polymyxin-based combinations",
    "Polymyxin", "Sulbactam",
    "Carbapenem-based combinations", "Carbapenem",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor"
  ))
  tgcre_order <- rev(c(
    "Carbapenem", "Carbapenem-based combinations",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor-based combinations",
    "Third-generation cephalosporin", "Third-generation cephalosporin-based combinations",
    "Aminoglycoside-based combinations"
  ))
  cre_order <- rev(c(
    "Ceftazidime/avibactam monotherapy or combinations",
    "Polymyxin-based combinations",
    "Carbapenem", "Carbapenem-based combinations",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor-based combinations",
    "Aminoglycoside-based combinations"
  ))
  crp_order <- rev(c(
    "Polymyxin-based combinations",
    "Carbapenem", "Carbapenem-based combinations",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor",
    "Anti-pseudomonal penicillin/Beta-lactamase inhibitor-based combinations",
    "Aminoglycoside-based combinations"
  ))
  custom_orders <- list(CRA = cra_order, `3GCRE` = tgcre_order, CRE = cre_order, CRP = crp_order)
  
  if (subset_key %in% names(custom_orders)) {
    # Custom nodes present in this plot (ignore counts)
    custom_nodes <- unique(unlist(lapply(custom_orders[[subset_key]], pick_nodes_exact, pool = node_chr), use.names = FALSE))
    
    # Original priority order (combo first, then single)
    priority_vec <- rev(priority_map[[subset_key]])
    ordered_nodes <- c()
    for (base in priority_vec) {
      single <- grep(paste0("^", base, "( \\(|$)"), node_chr, value = TRUE)
      combo  <- grep(paste0("^", base, "-based combinations"), node_chr, value = TRUE)
      ordered_nodes <- c(ordered_nodes, combo, single)
    }
    # Remove anything already placed by custom
    ordered_nodes <- setdiff(ordered_nodes, c(special_nodes, custom_nodes))
    
    # Remaining classes not in priority_vec: alphabetical by base, combos before singles
    remaining_nodes <- setdiff(node_chr, c(special_nodes, custom_nodes, ordered_nodes))
    remaining_bases <- sort(unique(strip_combo(strip_counts(remaining_nodes))))
    remaining_ordered <- c()
    for (base in remaining_bases) {
      combo  <- grep(paste0("^", base, "-based combinations"), remaining_nodes, value = TRUE)
      single <- grep(paste0("^", base, "( \\(|$)"),            remaining_nodes, value = TRUE)
      remaining_ordered <- c(remaining_ordered, combo, single)
    }
    
    final_levels <- c(special_nodes, remaining_ordered, ordered_nodes, custom_nodes)
    
  } else {
    # VRE/MRSA 
    priority_vec <- rev(priority_map[[subset_key]])
    ordered_nodes <- c()
    for (base in priority_vec) {
      single <- grep(paste0("^", base, "( \\(|$)"), node_chr, value = TRUE)
      combo  <- grep(paste0("^", base, "-based combinations"), node_chr, value = TRUE)
      ordered_nodes <- c(ordered_nodes, combo, single)
    }
    remaining_nodes <- setdiff(node_chr, c(special_nodes, ordered_nodes))
    remaining_nodes <- remaining_nodes[!strip_combo(strip_counts(remaining_nodes)) %in% priority_vec]
    remaining_bases <- sort(unique(strip_combo(strip_counts(remaining_nodes))))
    remaining_ordered <- c()
    for (base in remaining_bases) {
      combo  <- grep(paste0("^", base, "-based combinations"), remaining_nodes, value = TRUE)
      single <- grep(paste0("^", base, "( \\(|$)"), remaining_nodes, value = TRUE)
      remaining_ordered <- c(remaining_ordered, combo, single)
    }
    final_levels <- c(special_nodes, remaining_ordered, ordered_nodes)
  }
  
  df_long$node <- factor(df_long$node, levels = final_levels)
  
  # Colors: Others & special nodes = gray; combos = lighter shade
  col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  mycol <- unlist(mapply(RColorBrewer::brewer.pal, col_pals$maxcolors, rownames(col_pals)))
  vibrant_first <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                     "#FFD92F", "#A65628", "#F781BF", "#66C2A5", "#E7298A", 
                     "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D" )
  mycol <- c(vibrant_first, setdiff(mycol, vibrant_first))
  
  uniq_bases <- unique(df_long$node_base)
  color_bases <- setdiff(uniq_bases, c("None","Died","Discharged","Others"))
  if (length(color_bases) > length(mycol)) mycol <- rep(mycol, length.out = length(color_bases))
  base_colors <- setNames(mycol[seq_along(color_bases)], color_bases)
  base_colors[c("None","Died","Discharged")] <- "gray80"
  base_colors["Others"] <- "#35978f"
  
  lighten <- function(col, amount = 0.4) {
    col_rgb <- col2rgb(col, alpha = FALSE) / 255
    mixed   <- (1 - amount) * col_rgb + amount * 1
    rgb(mixed[1], mixed[2], mixed[3])
  }
  unique_nodes <- unique(as.character(df_long$node))
  node_colors <- vapply(seq_along(unique_nodes), function(i){
    b <- node_base(unique_nodes[i])
    col <- base_colors[b]
    if (grepl("-based combinations\\b", unique_nodes[i])) return(lighten(col, 0.4))
    unname(col)
  }, character(1))
  names(node_colors) <- unique_nodes
  
  # Plot
  ggplot(
    df_long,
    aes(x = x, next_x = next_x, node = node, next_node = next_node,
        fill = node, label = node)
  ) +
    geom_sankey(flow.alpha = 0.75, smooth = 9, width = 0.05, alpha = 1, space = space_value) +
    geom_sankey_label(
      aes(
        x = stage(x, after_stat = x + 0.1 * dplyr::case_when(x == 1 ~ -0.27, x == 2 ~ 0.27, .default = 0)),
        hjust = dplyr::case_when(x == "emp_node" ~ 1, x == "def_node" ~ 0, .default = 0.5)
      ),
      size = 6, color = "black", fill = NA, fontface = "bold",
      label.r = unit(0.3, "mm"), label.size = 0.2, space = space_value,
      family = "Times New Roman"
    ) +
    theme_void() +
    annotate("text", x = 0.95, y = Inf,
             label = paste0("Empirical prescriptions\n(n = ", N_emp_rx, ")"),
             hjust = 1, vjust = 1, size = 7, fontface = "bold", family = "Times New Roman") +
    annotate("text", x = 2.05, y = Inf,
             label = paste0("Definitive prescriptions\n(n = ", N_def_rx, ")"),
             hjust = 0, vjust = 1, size = 7, fontface = "bold", family = "Times New Roman") +
    labs(title = title_html) +
    theme(
      plot.title = ggtext::element_markdown(hjust = 0.5, size = 22, family = "Times New Roman"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = node_colors)
}

# Build plots, save, crop
plots <- lapply(subset_list, function(s) {
  plot_sankey_based_for_ids(
    id_vec = s$ids,
    main_title = s$title,
    subset_key = s$key,
    space_value = s$space
  )
})

heights <- vapply(subset_list, function(s) s$height, numeric(1))
dir.create("output/figure", recursive = TRUE, showWarnings = FALSE)

labels <- rep(letters[1:3], length.out = length(plots)) 
x_positions <- rep(c(0.128, 0.085, 0.072, 0.087, 0.087, 0.087), length.out = length(plots))

label_and_save <- function(p, label, x_pos, file, width = 47, height = 20) {
  lab_g <- grid::textGrob(
    paste0("(", label, ")"),
    x = grid::unit(x_pos, "npc"),
    y = grid::unit(0.9, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 20, fontface = "bold", family = "Times New Roman")
  )
  g <- gridExtra::arrangeGrob(p, top = lab_g)
  
  # Use grDevices::cairo_pdf
  grDevices::cairo_pdf(filename = file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.draw(g)
}

for (i in seq_along(plots)) {
  fn <- sprintf("output/figure/sankey_others_p%s.pdf", i)
  label_and_save(
    p = plots[[i]],
    label = labels[i],
    x_pos = x_positions[i],
    file = fn,
    width = 47,
    height = heights[i]
  )
}

pdfcrop_bin <- Sys.which("pdfcrop")
if (nzchar(pdfcrop_bin)) {
  pdfs <- list.files("output/figure", pattern = "^sankey_others_p\\d+\\.pdf$", full.names = TRUE)
  for (f in pdfs) system2(pdfcrop_bin, c(shQuote(f), shQuote(f)))
} else {
  message("pdfcrop not found; skipping cropping.")
}
