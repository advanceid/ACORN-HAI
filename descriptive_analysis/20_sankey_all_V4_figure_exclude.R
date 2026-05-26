# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr, dplyr, tidyr, stringr, lubridate, purrr,
    ggplot2, ggtext, ggsankey, RColorBrewer,
    grid, gridExtra, Cairo, openxlsx
  )
})

# Load data
wd <- "./"; setwd(wd)
df_data     <- readRDS("data/clean_data/anti_treat_index.RData")
df_baseline <- readRDS("data/clean_data/baseline_outcomes_index.RData")
episode_num <- readRDS("data/clean_data/episode_num_for_treatment.RData")

# Remove antifungals
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
        anti_group == "Anti-pseudomonal penicillin/beta-lactamase inhibitor" ~ 
        "Other anti-pseudomonal penicillin/beta-lactamase inhibitor",
      
      # Default
      TRUE ~ anti_group
    )
  )

# Per-pathogen priority for anti_used (full names, no abbreviations)
priority_map <- list(
  CRA   = c("Sulbactam","Polymyxin","Carbapenem",
            "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
            "Other beta-lactam/beta-lactamase inhibitor",
            "Aminoglycoside","Third-generation cephalosporin",
            "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  `3GCRE` = c("Carbapenem","Polymyxin",
              "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
              "Other beta-lactam/beta-lactamase inhibitor",
              "Aminoglycoside","Third-generation cephalosporin",
              "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  CRE   = c("Ceftazidime/avibactam","Polymyxin","Carbapenem",
            "Other anti-pseudomonal penicillin/beta-lactamase inhibitor",
            "Other beta-lactam/beta-lactamase inhibitor",
            "Aminoglycoside","Third-generation cephalosporin",
            "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  CRP   = c("Anti-pseudomonal penicillin/beta-lactamase inhibitor",
            "Polymyxin","Carbapenem",
            "Other beta-lactam/beta-lactamase inhibitor","Aminoglycoside"),
  VRE   = c("Lipopeptide","Oxazolidinone","Penicillin",
            "Third-generation cephalosporin","Phosphonic"),
  MRSA  = c("Glycopeptide","Lipopeptide","Oxazolidinone",
            "Sulfonamide-trimethoprim-combination","Lincosamide")
)

# Build complete priority vector for a subset (append unseen classes at the end)
priority_all_for <- function(key, pool_values) {
  base <- priority_map[[key]]
  others <- setdiff(sort(unique(pool_values)), base)
  c(base, others)
}

# Helper: combine multiple classes on the same day with a priority order
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

# Expand to patient-day, then pick window days around reference dates
daily_antibiotics <- function(df,
                              id_col="recordid",
                              group_col="anti_group",
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
      mutate(.dates = list(.data[[ref_col]] + days(offsets)), .offs = list(offsets)) %>%
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

# Outcomes (def all None -> Outcome)
df_outcomes <- df_baseline %>%
  select(recordid, ho_discharge_date, mortality_date) %>%
  mutate(outcomes = dplyr::case_when(
    !is.na(mortality_date) ~ "Died",
    is.na(ho_discharge_date) & is.na(mortality_date) ~ "None",
    TRUE ~ "Discharged"
  )) %>%
  select(recordid, outcomes)

# Organism subsets
aci_car_id  <- episode_num %>% filter(CRA == TRUE) %>% pull(recordid) %>% unique()
ent_thir_id <- episode_num %>% filter(`3GCRE` == TRUE) %>% pull(recordid) %>% unique()
ent_car_id  <- episode_num %>% filter(CRE == TRUE) %>% pull(recordid) %>% unique()
pse_car_id  <- episode_num %>% filter(CRP == TRUE) %>% pull(recordid) %>% unique()
entc_van_id <- episode_num %>% filter(VRE == TRUE) %>% pull(recordid) %>% unique()
sa_meth_id  <- episode_num %>% filter(MRSA == TRUE) %>% pull(recordid) %>% unique()

subset_list <- list(
  list(key="CRA",   ids = aci_car_id,  title = "Carbapenem-resistant Acinetobacter spp.", space=70, height=29),
  list(key="3GCRE", ids = ent_thir_id, title = "Third-generation cephalosporin-resistant Enterobacterales", space=60, height=29),
  list(key="CRE",   ids = ent_car_id,  title = "Carbapenem-resistant Enterobacterales", space=60, height=27),
  list(key="CRP",   ids = pse_car_id,  title = "Carbapenem-resistant Pseudomonas spp.",  space=10, height=27),
  list(key="VRE",   ids = entc_van_id, title = "Vancomycin-resistant Enterococcus spp.",  space=4,  height=21),
  list(key="MRSA",  ids = sa_meth_id,  title = "Methicillin-resistant Staphylococcus aureus", space=10, height=23)
)

# Factories using a local (subset) priority vector
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

# Build per-subset window sets based on anti_used and subset priority
build_window_sets_for_ids <- function(id_vec, priority_vec){
  all_ids <- tibble(recordid = unique(id_vec))
  
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
  
  emp_def_wide_local <- all_ids %>%
    left_join(res_local$wide, by = "recordid")
  
  emp_cols <- c("emp_m1","emp_0","emp_p1","emp_p2")
  def_cols <- c("def_d3","def_d4","def_d5")
  
  for (cc in c(emp_cols, def_cols)) {
    if (!cc %in% names(emp_def_wide_local)) {
      emp_def_wide_local[[cc]] <- NA_character_
    }
  }
  
  emp_def_wide_local <- emp_def_wide_local %>%
    mutate(across(all_of(c(emp_cols, def_cols)), as.character))
  
  normalize_daily_local <- normalize_daily_local_factory(priority_vec)
  
  emp_final_local <- {
    long <- emp_def_wide_local %>%
      select(recordid, all_of(emp_cols)) %>%
      pivot_longer(cols = all_of(emp_cols), names_to="day", values_to="treatment") %>%
      mutate(
        treatment = ifelse(is.na(treatment) | treatment == "", "None", treatment),
        treatment = vapply(treatment, normalize_daily_local, character(1))
      )
    
    all_none <- long %>%
      group_by(recordid) %>%
      summarise(all_none = all(treatment == "None"), .groups="drop")
    
    keep_non_none <- long %>%
      inner_join(all_none, by="recordid") %>%
      filter(!all_none, treatment != "None") %>%
      distinct(recordid, treatment)
    
    keep_none <- all_none %>%
      filter(all_none) %>%
      transmute(recordid, treatment = "None")
    
    bind_rows(keep_non_none, keep_none) %>%
      distinct(recordid, treatment) %>%
      rename(emp_treatment = treatment)
  }
  
  def_final_local <- {
    long <- emp_def_wide_local %>%
      select(recordid, all_of(def_cols)) %>%
      pivot_longer(cols = all_of(def_cols), names_to="day", values_to="treatment") %>%
      mutate(
        treatment = ifelse(is.na(treatment) | treatment == "", "None", treatment),
        treatment = vapply(treatment, normalize_daily_local, character(1))
      )
    
    all_none <- long %>%
      group_by(recordid) %>%
      summarise(all_none = all(treatment == "None"), .groups="drop")
    
    keep_non_none <- long %>%
      inner_join(all_none, by="recordid") %>%
      filter(!all_none, treatment != "None") %>%
      distinct(recordid, treatment) %>%
      rename(def_treatment = treatment)
    
    keep_outcome <- all_none %>%
      filter(all_none) %>%
      left_join(df_outcomes, by="recordid") %>%
      mutate(def_treatment = ifelse(is.na(outcomes), "None", outcomes)) %>%
      select(recordid, def_treatment) %>%
      distinct()
    
    bind_rows(keep_non_none, keep_outcome) %>%
      distinct(recordid, def_treatment)
  }
  
  list(
    emp_final = emp_final_local,
    def_final = def_final_local,
    emp_def_wide = emp_def_wide_local
  )
}

# Pair-level exclusion helper 
apply_pair_exclusion <- function(emp_df, def_df) {
  joined <- emp_df %>%
    full_join(def_df, by = "recordid", relationship = "many-to-many") %>%
    mutate(
      emp_treatment = if_else(is.na(emp_treatment), "None", emp_treatment),
      def_treatment = if_else(is.na(def_treatment), "None", def_treatment)
    )
  
  list(
    joined = joined,
    emp = joined %>% distinct(recordid, emp_treatment),
    def = joined %>% distinct(recordid, def_treatment)
  )
}


# Plot sankey for a subset (special nodes without counts)
special_counts_global <- list()

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
  
  mk_rx_counts_side <- function(recordids, trts){
    tibble(recordid = recordids, treatment = trts) %>%
      mutate(
        is_special = treatment %in% special,
        is_combo   = !is_special & grepl(",", treatment),
        node_name  = dplyr::case_when(
          is_special ~ treatment,
          is_combo   ~ paste0(pick_base_by_priority_vec(treatment), "-based combinations"),
          TRUE       ~ treatment
        ),
        combo_sig  = ifelse(is_combo, norm_combo_vec(treatment), NA_character_)
      ) %>%
      distinct(recordid, node_name, combo_sig) %>%
      mutate(rx = 1L) %>%
      group_by(node_name) %>%
      summarise(rx = sum(rx), .groups = "drop")
  }
  
  emp_rx_tbl <- mk_rx_counts_side(emp_pairs$recordid, emp_pairs$emp_treatment)
  def_rx_tbl <- mk_rx_counts_side(def_pairs$recordid, def_pairs$def_treatment)
  
  
  # Collect into a global list for exporting to Excel later
  cat("Empirical special node counts:\n")
  print(emp_rx_tbl %>% dplyr::filter(node_name %in% special))
  
  cat("Definitive special node counts:\n")
  print(def_rx_tbl %>% dplyr::filter(node_name %in% special))
  
  special_counts_global[[paste0(subset_key, "_Empirical")]]  <<- emp_rx_tbl %>%
    dplyr::filter(node_name %in% special) %>%
    dplyr::mutate(window = "Empirical", subset = subset_key)
  
  special_counts_global[[paste0(subset_key, "_Definitive")]] <<- def_rx_tbl %>%
    dplyr::filter(node_name %in% special) %>%
    dplyr::mutate(window = "Definitive", subset = subset_key)
  
  # Maps for labels (non-special nodes use "(n)"; special nodes show plain text)
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
  
  # SPECIAL NODES: show without counts
  label_emp <- function(t) {
    if (t %in% special) return(t)
    if (grepl(",", t))
      return(safe_pick(emp_based_map, pick_base_by_priority_vec(t),
                       paste0(pick_base_by_priority_vec(t), "-based combinations (0)")))
    safe_pick(emp_rx_map, t, paste0(t, " (0)"))
  }
  label_def <- function(t) {
    if (t %in% special) return(t)
    if (grepl(",", t))
      return(safe_pick(def_based_map, pick_base_by_priority_vec(t),
                       paste0(pick_base_by_priority_vec(t), "-based combinations (0)")))
    safe_pick(def_rx_map, t, paste0(t, " (0)"))
  }
  
  links_named <- joined %>%
    transmute(
      emp_node = vapply(emp_treatment, label_emp, character(1)),
      def_node = vapply(def_treatment, label_def, character(1))
    ) %>%
    count(emp_node, def_node, name = "n")
  
  links_rep <- tidyr::uncount(links_named, weights = n) %>%
    mutate(emp_node = as.character(emp_node), def_node = as.character(def_node))
  
  df_long <- ggsankey::make_long(links_rep, emp_node, def_node)
  
  strip_paren <- function(s) sub(" \\(.*\\)$", "", s)
  visible_emp_nodes <- unique(strip_paren(as.character(df_long$node[df_long$x == "emp_node"])))
  visible_def_nodes <- unique(strip_paren(as.character(df_long$node[df_long$x == "def_node"])))
  
  emp_node_rx_map <- setNames(emp_rx_tbl$rx, emp_rx_tbl$node_name)
  def_node_rx_map <- setNames(def_rx_tbl$rx, def_rx_tbl$node_name)
  
  # Totals EXCLUDING special nodes
  N_emp_rx <- sum(emp_node_rx_map[setdiff(visible_emp_nodes, special)], na.rm = TRUE)
  N_def_rx <- sum(def_node_rx_map[setdiff(visible_def_nodes, special)], na.rm = TRUE)
  
  # Title with n patients
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
    sprintf("<b>%s (%s episodes)</b>", title_core, format(n, big.mark = "", scientific = FALSE))
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
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor"
  ))
  tgcre_order <- rev(c(
    "Carbapenem", "Carbapenem-based combinations",
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor-based combinations",
    "Third-generation cephalosporin", "Third-generation cephalosporin-based combinations",
    "Aminoglycoside-based combinations"
  ))
  cre_order <- rev(c(
    "Polymyxin-based combinations",
    "Carbapenem", "Carbapenem-based combinations",
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor-based combinations",
    "Aminoglycoside-based combinations"
  ))
  crp_order <- rev(c(
    "Polymyxin-based combinations",
    "Carbapenem", "Carbapenem-based combinations",
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
    "Anti-pseudomonal penicillin/beta-lactamase inhibitor-based combinations",
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
  
  #
  node_base <- function(s) { s <- sub(" \\(.*\\)$", "", s); sub("-based combinations$", "", s) }
  df_long$node_base <- node_base(as.character(df_long$node))
  
  # Colors
  col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual', ]
  mycol <- unlist(mapply(RColorBrewer::brewer.pal, col_pals$maxcolors, rownames(col_pals)))
  vibrant_first <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
                     "#FFD92F", "#A65628", "#F781BF", "#66C2A5", "#E7298A", 
                     "#1B9E77", "#D95F02", "#7570B3", "#E6AB02", "#A6761D" )
  mycol <- c(vibrant_first, setdiff(mycol, vibrant_first))
  uniq_bases <- unique(df_long$node_base)
  color_bases <- setdiff(uniq_bases, c("None","Died","Discharged"))
  if (length(color_bases) > length(mycol)) mycol <- rep(mycol, length.out = length(color_bases))
  base_colors <- setNames(mycol[seq_along(color_bases)], color_bases)
  base_colors[c("None","Died","Discharged")] <- "gray80"
  
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
             label = paste0("Empirical prescriptions\n(total = ", N_emp_rx, ")"),
             hjust = 1, vjust = 1, size = 7, fontface = "bold", family = "Times New Roman") +
    annotate("text", x = 2.05, y = Inf,
             label = paste0("Definitive prescriptions\n(total = ", N_def_rx, ")"),
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

labels      <- letters[seq_along(plots)]
x_positions <- rep(c(0.086, 0.084, 0.086, 0.086, 0.087, 0.086), length.out = length(plots))

label_and_save <- function(p, label, x_pos, file, width = 47, height = 20) {
  lab_g <- grid::textGrob(
    paste0("(", label, ")"),
    x = grid::unit(x_pos, "npc"),
    y = grid::unit(0.9, "npc"),
    just = c("left", "top"),
    gp = grid::gpar(fontsize = 20, fontface = "bold", family = "Times New Roman")
  )
  g <- gridExtra::arrangeGrob(p, top = lab_g)
  cairo_pdf(filename = file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.draw(g)
}

for (i in seq_along(plots)) {
  fn <- sprintf("output/figure/sankey_based_p%s.pdf", i)
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
  pdfs <- list.files("output/figure", pattern = "^sankey_based_p\\d+\\.pdf$", full.names = TRUE)
  for (f in pdfs) system2(pdfcrop_bin, c(shQuote(f), shQuote(f)))
} else {
  message("pdfcrop not found; skipping cropping.")
}


# Combine and export special node counts to Excel
if (length(special_counts_global)) {
  dir.create("output/table", recursive = TRUE, showWarnings = FALSE)
  special_counts_all <- dplyr::bind_rows(special_counts_global)
  wb3 <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb3, "special_nodes")
  openxlsx::writeData(wb3, "special_nodes", special_counts_all)
  openxlsx::saveWorkbook(wb3, "output/table/special_node_counts.xlsx", overwrite = TRUE)
}

# =======================
# Distinct combinations inside each "-based combinations" node
# =======================
# 1) Detailed combos per window
detailed_combos_for_ids <- function(id_vec, window = c("emp","def"), priority_vec){
  window <- match.arg(window)
  
  norm_combo_vec <- norm_combo_vec_factory(priority_vec)
  pick_base_vec  <- pick_base_by_priority_factory(priority_vec)
  
  ws <- build_window_sets_for_ids(id_vec, priority_vec)
  emp_pairs_raw <- ws$emp_final %>% dplyr::distinct(recordid, emp_treatment)
  def_pairs_raw <- ws$def_final %>% dplyr::distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  side_pairs <- if (window == "emp") {
    filt$emp %>% dplyr::rename(treatment = emp_treatment)
  } else {
    filt$def %>% dplyr::rename(treatment = def_treatment)
  }
  
  special <- c("None","Died","Discharged")
  
  combos_df <- side_pairs %>%
    dplyr::filter(!is.na(treatment),
                  !treatment %in% special,
                  grepl(",", treatment)) %>%     # true combinations only
    dplyr::mutate(
      combination_list = norm_combo_vec(treatment), # normalize order
      base             = pick_base_vec(combination_list) # assign base
    ) %>%
    dplyr::distinct(recordid, base, combination_list)    # one patient once per combo
  
  out <- combos_df %>%
    dplyr::group_by(base, combination_list) %>%
    dplyr::summarise(detailed_rx = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(detailed_rx), base, combination_list)
  
  return(out)
}


# 2) Base-level summary
mk_rx_counts_with_distinct_local <- function(recordids, trts, priority_vec){
  pick_base_by_priority_vec <- pick_base_by_priority_factory(priority_vec)
  norm_combo_vec <- norm_combo_vec_factory(priority_vec)
  
  tibble(recordid = recordids, treatment = trts) %>%
    dplyr::mutate(
      is_special = treatment %in% c("None","Died","Discharged"),
      is_combo   = !is_special & grepl(",", treatment),
      node_name  = dplyr::case_when(
        is_special ~ treatment,
        is_combo   ~ paste0(pick_base_by_priority_vec(treatment), "-based combinations"),
        TRUE       ~ treatment
      ),
      combo_sig  = dplyr::if_else(is_combo, norm_combo_vec(treatment), NA_character_)
    ) %>%
    dplyr::distinct(recordid, node_name, combo_sig) -> tmp
  
  rx_tbl <- tmp %>% dplyr::group_by(node_name) %>% dplyr::summarise(rx = dplyr::n(), .groups = "drop")
  combo_tbl <- tmp %>%
    dplyr::filter(!is.na(combo_sig)) %>%
    dplyr::group_by(node_name) %>%
    dplyr::summarise(distinct_combos = dplyr::n_distinct(combo_sig), .groups = "drop")
  
  dplyr::full_join(rx_tbl, combo_tbl, by = "node_name") %>%
    dplyr::mutate(distinct_combos = dplyr::coalesce(distinct_combos, 0L))
}

combo_base_summary_for_ids <- function(id_vec, window = c("emp","def"), priority_vec){
  window <- match.arg(window)
  ws <- build_window_sets_for_ids(id_vec, priority_vec)
  emp_pairs_raw <- ws$emp_final %>% dplyr::distinct(recordid, emp_treatment)
  def_pairs_raw <- ws$def_final %>% dplyr::distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  side_pairs <- if (window == "emp") {
    filt$emp %>% dplyr::rename(treatment = emp_treatment)
  } else {
    filt$def %>% dplyr::rename(treatment = def_treatment)
  }
  
  counts <- mk_rx_counts_with_distinct_local(side_pairs$recordid, side_pairs$treatment, priority_vec)
  special   <- c("None","Died","Discharged")
  denom_all <- sum(counts$rx[!counts$node_name %in% special], na.rm = TRUE)
  
  counts %>%
    dplyr::filter(grepl("-based combinations$", node_name)) %>%
    dplyr::mutate(
      base = sub("-based combinations$", "", node_name),
      denom_all  = denom_all,
      pct_of_all = dplyr::if_else(denom_all > 0, round(100 * rx / denom_all, 1), NA_real_)
    ) %>%
    dplyr::arrange(dplyr::desc(rx), base) %>%
    dplyr::select(base, node_name, rx, distinct_combos, denom_all, pct_of_all)
}


# 3) Build nested summary + details 

one_subset_all_with_detail <- function(s) {
  priority_vec <- priority_all_for(s$key, df_data$anti_used)
  
  # --- Empirical window ---
  emp_summary <- combo_base_summary_for_ids(s$ids, "emp", priority_vec = priority_vec) %>%
    dplyr::mutate(window = "Empirical") %>%
    dplyr::select(-node_name)
  
  emp_detail_raw <- detailed_combos_for_ids(s$ids, "emp", priority_vec = priority_vec)
  
  emp_detail_nested <- emp_detail_raw %>%
    dplyr::select(base, combination_list, detailed_rx) %>%
    dplyr::group_by(base) %>%
    tidyr::nest(details = c(combination_list, detailed_rx)) %>%
    dplyr::ungroup()
  
  emp_tab <- emp_summary %>%
    dplyr::left_join(emp_detail_nested, by = "base")
  
  # --- Definitive window ---
  def_summary <- combo_base_summary_for_ids(s$ids, "def", priority_vec = priority_vec) %>%
    dplyr::mutate(window = "Definitive") %>%
    dplyr::select(-node_name)
  
  def_detail_raw <- detailed_combos_for_ids(s$ids, "def", priority_vec = priority_vec)
  
  def_detail_nested <- def_detail_raw %>%
    dplyr::select(base, combination_list, detailed_rx) %>%
    dplyr::group_by(base) %>%
    tidyr::nest(details = c(combination_list, detailed_rx)) %>%
    dplyr::ungroup()
  
  def_tab <- def_summary %>%
    dplyr::left_join(def_detail_nested, by = "base")
  
  # --- Combine and sanitize list-column ---
  combined_summary <- dplyr::bind_rows(emp_tab, def_tab) %>%
    dplyr::mutate(
      subset_key   = s$key,
      subset_title = s$title,
      node_name    = paste0(base, "-based combinations")
    ) %>%
    dplyr::relocate(subset_key, subset_title, window, node_name, base)
  
  # Ensure a well-formed list-column even when empty/NULL
  combined_summary <- combined_summary %>%
    dplyr::mutate(
      details = purrr::map(
        details,
        ~{
          if (is.null(.x)) {
            tibble::tibble(combination_list = character(), detailed_rx = integer())
          } else {
            x <- tibble::as_tibble(.x)
            req <- c("combination_list","detailed_rx")
            if (!all(req %in% names(x))) tibble::tibble(combination_list = character(), detailed_rx = integer()) else x
          }
        }
      )
    )
  
  return(combined_summary)
}


# 4) Build all and export to Excel
combos_all_with_detail <- purrr::map(subset_list, one_subset_all_with_detail) %>%
  dplyr::bind_rows()

cat("\nSummary Table with Nested Detail Column (combos_all_with_detail):\n")
print(utils::head(combos_all_with_detail, 20))

wb_combined <- openxlsx::createWorkbook()
dir.create("output/table", recursive = TRUE, showWarnings = FALSE)

for (s in subset_list) {
  sh_name <- s$key
  openxlsx::addWorksheet(wb_combined, sh_name)
  
  summary_data <- combos_all_with_detail %>%
    dplyr::filter(subset_key == s$key)
  
  # PART 1: Summary (exclude list-column)
  summary_output <- summary_data %>%
    dplyr::select(window, base, node_name, rx, distinct_combos, pct_of_all, denom_all)
  
  openxlsx::writeData(wb_combined, sh_name, x = s$title, startRow = 1, startCol = 1)
  openxlsx::writeData(wb_combined, sh_name, x = "Treatment Combination Summary", startRow = 2, startCol = 1)
  openxlsx::writeData(wb_combined, sh_name, x = summary_output, startRow = 3, startCol = 1, withFilter = TRUE)
  
  # PART 2: Detailed Breakdown
  start_row_detail <- 5 + nrow(summary_output)
  
  detailed_output <- summary_data %>%
    dplyr::filter(rx > 0) %>%
    tidyr::unnest(details, keep_empty = TRUE) %>%
    dplyr::mutate(
      node_name = paste0(base, "-based combinations (", rx, ")")
    ) %>%
    dplyr::rename(
      Prescriptions =  window,
      `Combination names` = node_name,
      `Distinct combinations` = combination_list,
      Counts = detailed_rx
    ) %>%
    dplyr::select(Prescriptions, `Combination names`, `Distinct combinations`, Counts) %>%
    # >>> KEY CHANGE: within each Prescription + Combination group,
    # show the group name only for the first row; blank for subsequent rows
    dplyr::group_by(Prescriptions, `Combination names`) %>%
    dplyr::mutate(
      `Combination names` = dplyr::if_else(dplyr::row_number() == 1, `Combination names`, "")
    ) %>%
    dplyr::ungroup()
  
  if (nrow(detailed_output) > 0) {
    openxlsx::writeData(
      wb_combined, sh_name,
      x = "--- Detailed Combination Breakdown ---",
      startRow = start_row_detail, startCol = 1
    )
    openxlsx::writeData(
      wb_combined, sh_name,
      x = detailed_output,
      startRow = start_row_detail + 2, startCol = 1, withFilter = TRUE
    )
    
  } else {
    openxlsx::writeData(
      wb_combined, sh_name,
      x = "--- Detailed Combination Breakdown: No combinations found. ---",
      startRow = start_row_detail, startCol = 1
    )
  }
}

openxlsx::saveWorkbook(
  wb_combined,
  file = "output/table/summary_and_detailed_combinations_combined.xlsx",
  overwrite = TRUE
)


mk_rx_counts_with_distinct_local <- function(recordids, trts, priority_vec){
  pick_base_by_priority_vec <- pick_base_by_priority_factory(priority_vec)
  norm_combo_vec <- norm_combo_vec_factory(priority_vec)
  
  tibble(recordid = recordids, treatment = trts) %>%
    mutate(
      is_special = treatment %in% c("None","Died","Discharged"),
      is_combo   = !is_special & grepl(",", treatment),
      node_name  = dplyr::case_when(
        is_special ~ treatment,
        is_combo   ~ paste0(pick_base_by_priority_vec(treatment), "-based combinations"),
        TRUE       ~ treatment
      ),
      combo_sig  = ifelse(is_combo, norm_combo_vec(treatment), NA_character_)
    ) %>%
    distinct(recordid, node_name, combo_sig) -> tmp
  
  rx_tbl <- tmp %>% group_by(node_name) %>% summarise(rx = n(), .groups = "drop")
  combo_tbl <- tmp %>%
    filter(!is.na(combo_sig)) %>%
    group_by(node_name) %>%
    summarise(distinct_combos = n_distinct(combo_sig), .groups = "drop")
  
  full_join(rx_tbl, combo_tbl, by = "node_name") %>%
    mutate(distinct_combos = dplyr::coalesce(distinct_combos, 0L))
}

combo_base_summary_for_ids <- function(id_vec, window = c("emp","def"), priority_vec){
  window <- match.arg(window)
  ws <- build_window_sets_for_ids(id_vec, priority_vec)
  
  emp_pairs_raw <- ws$emp_final %>% dplyr::distinct(recordid, emp_treatment)
  def_pairs_raw <- ws$def_final %>% dplyr::distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  side_pairs <- if (window == "emp") {
    filt$emp %>% dplyr::rename(treatment = emp_treatment)
  } else {
    filt$def %>% dplyr::rename(treatment = def_treatment)
  }
  
  counts <- mk_rx_counts_with_distinct_local(side_pairs$recordid, side_pairs$treatment, priority_vec)
  
  special   <- c("None","Died","Discharged")
  denom_all <- sum(counts$rx[!counts$node_name %in% special], na.rm = TRUE)
  
  counts %>%
    dplyr::filter(grepl("-based combinations$", node_name)) %>%
    dplyr::mutate(
      base = sub("-based combinations$", "", node_name),
      denom_all = denom_all,
      pct_of_all = ifelse(denom_all > 0, round(100 * rx / denom_all, 1), NA_real_)
    ) %>%
    dplyr::arrange(dplyr::desc(rx), base) %>%
    dplyr::select(base, node_name, rx, distinct_combos, denom_all, pct_of_all)
}

one_subset_all <- function(s) {
  priority_vec <- priority_all_for(s$key, df_data$anti_used)
  emp_tab <- combo_base_summary_for_ids(s$ids, "emp", priority_vec = priority_vec) %>% mutate(window = "Empirical")
  def_tab <- combo_base_summary_for_ids(s$ids, "def", priority_vec = priority_vec) %>% mutate(window = "Definitive")
  bind_rows(emp_tab, def_tab) %>%
    mutate(subset_key = s$key, subset_title = s$title) %>%
    relocate(subset_key, subset_title, window)
}
combos_all <- purrr::map_dfr(subset_list, one_subset_all)
print(combos_all, n = 100)

# Export distinct-combo summaries
wb <- createWorkbook()
for (s in subset_list) {
  sh <- s$key
  addWorksheet(wb, sh)
  tab <- combos_all %>%
    dplyr::filter(subset_key == s$key) %>%
    dplyr::select(window, base, node_name, rx, distinct_combos, pct_of_all, denom_all)
  writeData(wb, sh, x = s$title, startRow = 1, startCol = 1)
  writeData(wb, sh, x = tab, startRow = 3, startCol = 1, withFilter = TRUE)
}
saveWorkbook(wb, file = "output/table/distinct_combinations_by_base.xlsx", overwrite = TRUE)


# =========================
# Combination vs Single — NODE-LEVEL (matches Sankey denominator)
#   Denominator = number of (recordid, node_name) pairs excluding special nodes
#   Node construction is identical to the Sankey code.
# =========================

# Build node instances exactly like the plot (per side)
.build_side_node_instances <- function(side_pairs, priority_vec){
  special <- c("None","Died","Discharged")
  pick_base_by_priority_vec <- pick_base_by_priority_factory(priority_vec)
  norm_combo_vec            <- norm_combo_vec_factory(priority_vec)
  
  tibble::tibble(recordid = side_pairs$recordid, treatment = side_pairs$treatment) %>%
    dplyr::mutate(
      is_special = treatment %in% special,
      is_combo   = !is_special & grepl(",", treatment),
      node_name  = dplyr::case_when(
        is_special ~ treatment,
        is_combo   ~ paste0(pick_base_by_priority_vec(treatment), "-based combinations"),
        TRUE       ~ treatment
      ),
      combo_sig  = dplyr::if_else(is_combo, norm_combo_vec(treatment), NA_character_)
    ) %>%
    dplyr::distinct(recordid, node_name, combo_sig)   # <- node-level instances
}

combo_single_for_ids <- function(id_vec, window = c("emp","def"), priority_vec){
  window <- match.arg(window)
  ws <- build_window_sets_for_ids(id_vec, priority_vec)
  
  emp_pairs_raw <- ws$emp_final %>% dplyr::distinct(recordid, emp_treatment)
  def_pairs_raw <- ws$def_final %>% dplyr::distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  side_pairs <- if (window == "emp") {
    filt$emp %>% dplyr::rename(treatment = emp_treatment)
  } else {
    filt$def %>% dplyr::rename(treatment = def_treatment)
  }
  
  # Build node instances like the Sankey
  inst <- .build_side_node_instances(side_pairs, priority_vec)
  
  # Exclude special nodes (same as the plot when summing N)
  special <- c("None","Died","Discharged")
  inst_ns <- inst %>% dplyr::filter(!node_name %in% special)
  
  # Empty guard
  if (nrow(inst_ns) == 0) {
    return(tibble::tibble(
      window       = ifelse(window == "emp", "Empirical", "Definitive"),
      type         = character(),
      counts       = integer(),
      total_counts = integer(),
      percentage   = numeric()
    ))
  }
  
  # Node-level classification: "-based combinations" -> Combination, otherwise Single
  out <- inst_ns %>%
    dplyr::mutate(type = ifelse(grepl("-based combinations$", node_name), "Combination", "Single")) %>%
    dplyr::count(type, name = "counts") %>%
    dplyr::mutate(
      total_counts = sum(counts),
      percentage   = round(100 * counts / total_counts, 1),
      window       = ifelse(window == "emp", "Empirical", "Definitive")
    ) %>%
    dplyr::select(window, type, counts, total_counts, percentage)
  
  out
}

# ---- Run for all subsets
combo_single_all <- purrr::map_dfr(subset_list, function(s){
  priority_vec <- priority_all_for(s$key, df_data$anti_used)
  dplyr::bind_rows(
    combo_single_for_ids(s$ids, "emp", priority_vec) %>% dplyr::mutate(subset = s$key, .before = 1),
    combo_single_for_ids(s$ids, "def", priority_vec) %>% dplyr::mutate(subset = s$key, .before = 1)
  )
})

cat("\n=== Combination vs Single (node-level; Sankey-matched denominators) ===\n")
print(combo_single_all, n = 200)

# ---- Cross-check: reproduce the exact Sankey N and verify equality
sankey_N_for_ids <- function(id_vec, window = c("emp","def"), priority_vec){
  window <- match.arg(window)
  ws <- build_window_sets_for_ids(id_vec, priority_vec)
  emp_pairs_raw <- ws$emp_final %>% dplyr::distinct(recordid, emp_treatment)
  def_pairs_raw <- ws$def_final %>% dplyr::distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  side_pairs <- if (window == "emp") {
    filt$emp %>% dplyr::rename(treatment = emp_treatment)
  } else {
    filt$def %>% dplyr::rename(treatment = def_treatment)
  }
  
  inst <- .build_side_node_instances(side_pairs, priority_vec)
  special <- c("None","Died","Discharged")
  # Sankey's n = sum over node rx excluding specials = number of instances excluding specials
  nrow(inst %>% dplyr::filter(!node_name %in% special))
}

denom_compare <- purrr::map_dfr(subset_list, function(s){
  priority_vec <- priority_all_for(s$key, df_data$anti_used)
  
  emp_tbl <- combo_single_for_ids(s$ids, "emp", priority_vec)  %>% dplyr::mutate(subset = s$key)
  def_tbl <- combo_single_for_ids(s$ids, "def", priority_vec)  %>% dplyr::mutate(subset = s$key)
  
  dplyr::bind_rows(emp_tbl, def_tbl) %>%
    dplyr::distinct(subset, window, total_counts) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      sankey_n = ifelse(window == "Empirical",
                        sankey_N_for_ids(s$ids, "emp", priority_vec),
                        sankey_N_for_ids(s$ids, "def", priority_vec))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::rename(total_counts_used = total_counts)
})

cat("\n[CHECK] Denominator alignment (total_counts_used vs sankey_n):\n")
print(denom_compare, n = 200)



# ---------------------------------
# Excel: daily details per organism (using per-subset window sets)
summarize_phase <- function(df_sub, col) {
  df_sub %>%
    dplyr::transmute(treatment = .data[[col]]) %>%
    dplyr::mutate(treatment = ifelse(is.na(treatment) | treatment == "", "None", as.character(treatment))) %>%
    dplyr::count(treatment, name = "n") %>%
    dplyr::arrange(dplyr::desc(n)) %>%
    dplyr::mutate(
      percentage = round(100 * n / sum(n), 1),
      total_n = sum(n)
    )
}

if (!dir.exists("output/table")) dir.create("output/table", recursive = TRUE)
wb <- createWorkbook()
safe_name <- function(s) {
  s <- gsub("[\\/:?*\\[\\]]", "_", s)
  s <- gsub("\\s+", "_", s)
  substr(s, 1, 31)
}

for (x in subset_list) {
  key <- x$key; ids <- x$ids
  priority_vec <- priority_all_for(key, df_data$anti_used)
  ws <- build_window_sets_for_ids(ids, priority_vec)
  emp_sub <- ws$emp_final
  def_sub <- ws$def_final
  
  emp_pairs_raw <- emp_sub %>% dplyr::distinct(recordid, emp_treatment)
  def_pairs_raw <- def_sub %>% dplyr::distinct(recordid, def_treatment)
  filt <- apply_pair_exclusion(emp_pairs_raw, def_pairs_raw)
  
  emp_sub_f <- filt$emp %>% dplyr::rename(treatment = emp_treatment)
  def_sub_f <- filt$def %>% dplyr::rename(treatment = def_treatment)
  
  emp_tab <- emp_sub_f %>%
    dplyr::rename(emp_treatment = treatment) %>%
    summarize_phase(col = "emp_treatment")
  
  def_tab <- def_sub_f %>%
    dplyr::rename(def_treatment = treatment) %>%
    summarize_phase(col = "def_treatment")
  
  addWorksheet(wb, sheetName = safe_name(paste0(key, "_Empirical_daily")))
  writeData(wb, sheet = safe_name(paste0(key, "_Empirical_daily")), emp_tab)
  
  addWorksheet(wb, sheetName = safe_name(paste0(key, "_Definitive_daily")))
  writeData(wb, sheet = safe_name(paste0(key, "_Definitive_daily")), def_tab)
}

saveWorkbook(wb, file = "output/table/daily_details_by_organism.xlsx", overwrite = TRUE)


# ------------------------------------------------------------
# Newer broad-spectrum combo counts (uses anti_names; Sankey-consistent)
# ------------------------------------------------------------
cre_df <- df_data %>%
  dplyr::filter(ent_car == 1) %>%
  dplyr::select(recordid, anti_names)


newer_pat_strict <- "(ceftazidime\\s*[-/]?\\s*avibactam|meropenem\\s*[-/]?\\s*vaborbactam|imipenem.*relebactam|cefiderocol)"
newer_pat_broad  <- "(ceftazidime\\s*[-/]?\\s*avibactam|meropenem|imipenem|cefiderocol)"


by_patient <- cre_df %>%
  dplyr::mutate(
    hit_strict = grepl(newer_pat_strict, anti_names %||% "", ignore.case = TRUE),
    hit_broad  = grepl(newer_pat_broad,  anti_names %||% "", ignore.case = TRUE)
  ) %>%
  dplyr::group_by(recordid) %>%
  dplyr::summarise(
    used_strict = any(hit_strict, na.rm = TRUE),
    used_broad  = any(hit_broad,  na.rm = TRUE),
    .groups = "drop"
  )

N_cre <- dplyr::n_distinct(by_patient$recordid)


summary_tbl <- tibble::tibble(
  metric = c("CRE patients total",
             "Used STRICT (new combos only)",
             "Used BROAD (any mero/imipenem/cefiderocol or CZA)",
             "Overlap (both strict & broad)",
             "Strict only",
             "Broad only",
             "Neither"),
  n = c(
    N_cre,
    sum(by_patient$used_strict),
    sum(by_patient$used_broad),
    sum(by_patient$used_strict & by_patient$used_broad),
    sum(by_patient$used_strict & !by_patient$used_broad),
    sum(!by_patient$used_strict & by_patient$used_broad),
    sum(!by_patient$used_strict & !by_patient$used_broad)
  )
) %>%
  dplyr::mutate(pct = round(100 * n / N_cre, 1))

print(summary_tbl)

########
# =========================================================
# Overall episode-level ever use
# Denominator: df_baseline episodes
# Numerator: each episode counted once if ever used the drug
# =========================================================

# master denominator
episode_master <- df_baseline %>%
  dplyr::distinct(recordid)

N_total_episode <- nrow(episode_master)

# ever use by episode
overall_episode_use <- episode_master %>%
  dplyr::left_join(
    df_data %>%
      dplyr::mutate(
        anti_group2 = dplyr::case_when(
          anti_used == "Other anti-pseudomonal penicillin/beta-lactamase inhibitor" ~
            "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
          TRUE ~ anti_used
        )
      ) %>%
      dplyr::distinct(recordid, anti_group2),
    by = "recordid"
  ) %>%
  dplyr::filter(!is.na(anti_group2)) %>%
  dplyr::count(anti_group2, name = "n_episode") %>%
  dplyr::mutate(
    total_episode = N_total_episode,
    pct = round(100 * n_episode / total_episode, 1)
  ) %>%
  dplyr::arrange(dplyr::desc(n_episode), anti_group2)

overall_episode_use

# =========================================================
# CRA + CRE + CRP episode-level empirical / definitive use
# Denominator: treated episodes in each window
# Numerator: each episode counted once if ever used the drug in that window
# =========================================================

# target episodes from df_baseline / episode_num
ccp_ids <- union(union(aci_car_id, ent_car_id), pse_car_id)

# use your existing window builder
priority_vec_ccp <- c(
  "Polymyxin", "Carbapenem",
  "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
  "Other anti-pseudomonal penicillin/beta-lactamase inhibitor",
  "Other beta-lactam/beta-lactamase inhibitor",
  "Aminoglycoside", "Third-generation cephalosporin",
  "Glycylcycline", "Sulfonamide-trimethoprim-combination"
)

ws_ccp <- build_window_sets_for_ids(ccp_ids, priority_vec_ccp)

emp_cols <- c("emp_m1","emp_0","emp_p1","emp_p2")
def_cols <- c("def_d3","def_d4","def_d5")

# master denominator from baseline
ccp_master <- df_baseline %>%
  dplyr::filter(recordid %in% ccp_ids) %>%
  dplyr::distinct(recordid) %>%
  dplyr::left_join(ws_ccp$emp_def_wide, by = "recordid")

# helper: collapse "Other anti-pseudomonal..." into main class when searching
contains_class <- function(x, drug) {
  ifelse(
    is.na(x),
    FALSE,
    stringr::str_detect(x, stringr::fixed(drug))
  )
}

ccp_master2 <- ccp_master %>%
  dplyr::mutate(
    # treated denominator by window
    emp_any_treated = if_any(
      all_of(emp_cols),
      ~ !is.na(.) & . != "None"
    ),
    def_any_treated = if_any(
      all_of(def_cols),
      ~ !is.na(.) & . != "None"
    ),
    any_treated_either = emp_any_treated | def_any_treated,
    
    # ever used carbapenem / polymyxin in empirical window
    emp_carbapenem = if_any(
      all_of(emp_cols),
      ~ contains_class(., "Carbapenem")
    ),
    emp_polymyxin = if_any(
      all_of(emp_cols),
      ~ contains_class(., "Polymyxin")
    ),
    
    # ever used carbapenem / polymyxin in definitive window
    def_carbapenem = if_any(
      all_of(def_cols),
      ~ contains_class(., "Carbapenem")
    ),
    def_polymyxin = if_any(
      all_of(def_cols),
      ~ contains_class(., "Polymyxin")
    )
  )

# denominators
emp_den <- ccp_master2 %>%
  dplyr::filter(emp_any_treated) %>%
  nrow()

def_den <- ccp_master2 %>%
  dplyr::filter(def_any_treated) %>%
  nrow()

# numerators
emp_carb_num <- ccp_master2 %>%
  dplyr::filter(emp_any_treated, emp_carbapenem) %>%
  nrow()

emp_poly_num <- ccp_master2 %>%
  dplyr::filter(emp_any_treated, emp_polymyxin) %>%
  nrow()

def_carb_num <- ccp_master2 %>%
  dplyr::filter(def_any_treated, def_carbapenem) %>%
  nrow()

def_poly_num <- ccp_master2 %>%
  dplyr::filter(def_any_treated, def_polymyxin) %>%
  nrow()

result_ccp <- tibble::tibble(
  metric = c(
    "Empirical denominator",
    "Empirical carbapenem",
    "Empirical polymyxin",
    "Definitive denominator",
    "Definitive carbapenem",
    "Definitive polymyxin"
  ),
  n = c(
    emp_den,
    emp_carb_num,
    emp_poly_num,
    def_den,
    def_carb_num,
    def_poly_num
  ),
  pct = c(
    NA,
    round(100 * emp_carb_num / emp_den, 1),
    round(100 * emp_poly_num / emp_den, 1),
    NA,
    round(100 * def_carb_num / def_den, 1),
    round(100 * def_poly_num / def_den, 1)
  )
)

print(result_ccp)
