# Clear
rm(list = ls())

# Load packages
suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr, dplyr, tidyr, lubridate, purrr,
    ggplot2, ggtext, ggsankey, RColorBrewer,
    grid, gridExtra, Cairo, openxlsx, stringr,
    patchwork
  )
})

# Load data
wd <- "./"; setwd(wd)
df_data <- readRDS("data/clean_data/anti_treat_index.RData")
df_baseline <- readRDS("data/clean_data/baseline_outcomes_index.RData")
episode_num <- readRDS("data/clean_data/episode_num_for_treatment.RData")

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
        anti_group == "Anti-pseudomonal penicillin/beta-lactamase inhibitor" ~ 
        "Other anti-pseudomonal penicillin/beta-lactamase inhibitor",
      
      # Default
      TRUE ~ anti_group
    )
  )

# Per-pathogen priority for anti_used (full names, no abbreviations)
priority_map <- list(
  CRA = c("Sulbactam","Polymyxin","Carbapenem",
          "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
          "Other beta-lactam/beta-lactamase inhibitor",
          "Aminoglycoside","Third-generation cephalosporin",
          "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  `3GCRE` = c("Carbapenem","Polymyxin",
              "Anti-pseudomonal penicillin/beta-lactamase inhibitor",
              "Other beta-lactam/beta-lactamase inhibitor",
              "Aminoglycoside","Third-generation cephalosporin",
              "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  CRE = c("Ceftazidime/avibactam","Polymyxin","Carbapenem",
          "Other anti-pseudomonal penicillin/beta-lactamase inhibitor",
          "Other beta-lactam/beta-lactamase inhibitor",
          "Aminoglycoside","Third-generation cephalosporin",
          "Glycylcycline","Sulfonamide-trimethoprim-combination"),
  CRP = c("Anti-pseudomonal penicillin/beta-lactamase inhibitor",
          "Polymyxin","Carbapenem",
          "Other beta-lactam/beta-lactamase inhibitor","Aminoglycoside"),
  VRE = c("Lipopeptide","Oxazolidinone","Penicillin",
          "Third-generation cephalosporin","Phosphonic"),
  MRSA = c("Glycopeptide","Lipopeptide","Oxazolidinone",
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
# daily_antibiotics <- function(df,
#                               id_col="recordid",
#                               group_col="anti_used",
#                               start_col="anti_start",
#                               end_col  ="anti_end",
#                               windows=list(
#                                 emp=list(ref_col="inf_onset", offsets=c(-1,0,1,2)),
#                                 def=list(ref_col="spec_date", offsets=c(3,4,5))
#                               ),
#                               same_day_priority=NULL){
#   stopifnot(all(c(id_col,group_col,start_col,end_col) %in% names(df)))
#   to_date <- function(x) as.Date(x)
#   df <- df %>% mutate(
#     !!start_col := to_date(.data[[start_col]]),
#     !!end_col   := to_date(.data[[end_col]])
#   )
#   rx_base <- df %>% distinct(.data[[id_col]], .data[[group_col]], .data[[start_col]], .data[[end_col]])
#   anti_by_day <- rx_base %>%
#     filter(!is.na(.data[[start_col]]), !is.na(.data[[end_col]]), .data[[end_col]] >= .data[[start_col]]) %>%
#     rowwise() %>%
#     mutate(day_date = list(seq(.data[[start_col]], .data[[end_col]], by="day"))) %>%
#     unnest(day_date) %>% ungroup() %>%
#     group_by(.data[[id_col]], day_date) %>%
#     summarise(anti = combine_same_day(.data[[group_col]], priority = same_day_priority), .groups="drop")
#   
#   mk_grid <- function(ref_col, offsets, prefix){
#     ids <- df %>%
#       select(all_of(c(id_col, ref_col))) %>% distinct() %>%
#       mutate(!!ref_col := to_date(.data[[ref_col]])) %>% filter(!is.na(.data[[ref_col]]))
#     ids %>%
#       rowwise() %>%
#       mutate(.dates = list(.data[[ref_col]] + lubridate::days(offsets)), .offs = list(offsets)) %>%
#       unnest(c(.dates, .offs)) %>%
#       ungroup() %>%
#       mutate(win = prefix, rel = .offs,
#              label = paste0(prefix, ifelse(.offs >= 0, paste0("+", .offs), .offs))) %>%
#       rename(day_date = .dates, ref_date = !!ref_col)
#   }
#   grids <- purrr::imap(windows, ~ mk_grid(.x$ref_col, .x$offsets, .y)) %>% bind_rows()
#   daily_long <- grids %>%
#     left_join(anti_by_day, by = setNames(c(id_col,"day_date"), c(id_col,"day_date"))) %>%
#     mutate(anti = tidyr::replace_na(anti, "None")) %>%
#     select(all_of(c(id_col)), day_date, win, label, rel, anti)
#   
#   col_name <- function(win, rel){
#     if (win=="emp"){
#       if (rel<0) paste0("emp_m", abs(rel))
#       else if (rel==0) "emp_0" else paste0("emp_p", rel)
#     } else paste0("def_d", rel)
#   }
#   daily_wide <- daily_long %>%
#     mutate(col = pmap_chr(list(win,rel), col_name)) %>%
#     select(all_of(c(id_col)), col, anti) %>%
#     distinct() %>%
#     pivot_wider(names_from = col, values_from = anti)
#   
#   list(long=daily_long, wide=daily_wide)
# }
daily_antibiotics <- function(df,
                              master_ids_df,
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
  
  df <- df %>%
    mutate(
      !!start_col := to_date(.data[[start_col]]),
      !!end_col   := to_date(.data[[end_col]])
    )
  
  master_ids_df <- master_ids_df %>%
    mutate(across(any_of(c("inf_onset", "spec_date")), to_date))
  
  rx_base <- df %>%
    distinct(.data[[id_col]], .data[[group_col]], .data[[start_col]], .data[[end_col]])
  
  anti_by_day <- rx_base %>%
    filter(!is.na(.data[[start_col]]), !is.na(.data[[end_col]]), .data[[end_col]] >= .data[[start_col]]) %>%
    rowwise() %>%
    mutate(day_date = list(seq(.data[[start_col]], .data[[end_col]], by = "day"))) %>%
    unnest(day_date) %>%
    ungroup() %>%
    group_by(.data[[id_col]], day_date) %>%
    summarise(anti = combine_same_day(.data[[group_col]], priority = same_day_priority), .groups = "drop")
  
  mk_grid <- function(ref_col, offsets, prefix){
    ids <- master_ids_df %>%
      select(all_of(c(id_col, ref_col))) %>%
      distinct() %>%
      filter(!is.na(.data[[ref_col]]))
    
    ids %>%
      rowwise() %>%
      mutate(.dates = list(.data[[ref_col]] + lubridate::days(offsets)),
             .offs  = list(offsets)) %>%
      unnest(c(.dates, .offs)) %>%
      ungroup() %>%
      mutate(
        win   = prefix,
        rel   = .offs,
        label = paste0(prefix, ifelse(.offs >= 0, paste0("+", .offs), .offs))
      ) %>%
      rename(day_date = .dates, ref_date = !!ref_col)
  }
  
  grids <- purrr::imap(windows, ~ mk_grid(.x$ref_col, .x$offsets, .y)) %>%
    bind_rows()
  
  daily_long <- grids %>%
    left_join(anti_by_day, by = setNames(c(id_col, "day_date"), c(id_col, "day_date"))) %>%
    mutate(anti = tidyr::replace_na(anti, "None")) %>%
    select(all_of(c(id_col)), day_date, win, label, rel, anti)
  
  col_name <- function(win, rel){
    if (win == "emp") {
      if (rel < 0) paste0("emp_m", abs(rel))
      else if (rel == 0) "emp_0"
      else paste0("emp_p", rel)
    } else {
      paste0("def_d", rel)
    }
  }
  
  daily_wide <- daily_long %>%
    mutate(col = pmap_chr(list(win, rel), col_name)) %>%
    select(all_of(c(id_col)), col, anti) %>%
    distinct() %>%
    pivot_wider(names_from = col, values_from = anti)
  
  list(long = daily_long, wide = daily_wide)
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
# aci_car_id  <- df_data$recordid[which(df_data$aci_car == 1)]
# ent_thir_id <- df_data$recordid[which(df_data$ent_thir == 1)]
# ent_car_id  <- df_data$recordid[which(df_data$ent_car == 1)]
# pse_car_id  <- df_data$recordid[which(df_data$pse_car == 1)]
# entc_van_id <- df_data$recordid[which(df_data$entc_van == 1)]
# sa_meth_id  <- df_data$recordid[which(df_data$sa_meth == 1)]

# Organism recordid subsets: use episode_num as the master denominator
aci_car_id  <- episode_num %>% filter(CRA   == TRUE) %>% pull(recordid) %>% unique()
ent_thir_id <- episode_num %>% filter(`3GCRE` == TRUE) %>% pull(recordid) %>% unique()
ent_car_id  <- episode_num %>% filter(CRE   == TRUE) %>% pull(recordid) %>% unique()
pse_car_id  <- episode_num %>% filter(CRP   == TRUE) %>% pull(recordid) %>% unique()
entc_van_id <- episode_num %>% filter(VRE   == TRUE) %>% pull(recordid) %>% unique()
sa_meth_id  <- episode_num %>% filter(MRSA  == TRUE) %>% pull(recordid) %>% unique()

subset_list <- list(
  list(key="CRA", ids = aci_car_id, title = "Carbapenem-resistant Acinetobacter spp.", space=2, height=8),
  list(key="3GCRE", ids = ent_thir_id, title = "Third-generation cephalosporin-resistant Enterobacterales", space=16, height=13),
  list(key="CRE",   ids = ent_car_id,  title = "Carbapenem-resistant Enterobacterales", space=2, height=8),
  list(key="CRP",   ids = pse_car_id,  title = "Carbapenem-resistant Pseudomonas spp.",  space=1, height=8),
  list(key="VRE",   ids = entc_van_id, title = "Vancomycin-resistant Enterococcus spp.",  space=1,  height=8),
  list(key="MRSA",  ids = sa_meth_id,  title = "Methicillin-resistant Staphylococcus aureus", space=1, height=8)
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
# build_window_sets_for_ids <- function(id_vec, priority_vec){
#   res_local <- daily_antibiotics(
#     df = df_data %>% dplyr::filter(recordid %in% id_vec),
#     id_col="recordid",
#     group_col="anti_used",
#     start_col="anti_start",
#     end_col  ="anti_end",
#     windows = list(
#       emp=list(ref_col="inf_onset", offsets=c(-1,0,1,2)),
#       def=list(ref_col="spec_date", offsets=c(3,4,5))
#     ),
#     same_day_priority = priority_vec
#   )
#   emp_def_wide_local <- res_local$wide
#   emp_cols <- intersect(c("emp_m1","emp_0","emp_p1","emp_p2"), names(emp_def_wide_local))
#   def_cols <- intersect(c("def_d3","def_d4","def_d5"), names(emp_def_wide_local))
#   emp_def_wide_local <- emp_def_wide_local %>%
#     mutate(across(all_of(c(emp_cols, def_cols)), as.character))
#   normalize_daily_local <- normalize_daily_local_factory(priority_vec)
#   
#   # Empirical window set
#   emp_final_local <- {
#     long <- emp_def_wide_local %>%
#       select(recordid, all_of(emp_cols)) %>%
#       pivot_longer(cols = all_of(emp_cols), names_to="day", values_to="treatment") %>%
#       mutate(treatment = vapply(treatment, normalize_daily_local, character(1)))
#     all_none <- long %>% group_by(recordid) %>% summarise(all_none = all(treatment=="None"), .groups="drop")
#     keep_non_none <- long %>% inner_join(all_none, by="recordid") %>%
#       filter(!all_none, treatment!="None") %>%
#       distinct(recordid, treatment)
#     keep_none <- all_none %>% filter(all_none) %>%
#       transmute(recordid, treatment="None")
#     bind_rows(keep_non_none, keep_none) %>% rename(emp_treatment = treatment)
#   }
#   
#   # Definitive window set (None -> outcome)
#   def_final_local <- {
#     long <- emp_def_wide_local %>%
#       select(recordid, all_of(def_cols)) %>%
#       pivot_longer(cols = all_of(def_cols), names_to="day", values_to="treatment") %>%
#       mutate(treatment = vapply(treatment, normalize_daily_local, character(1)))
#     all_none <- long %>% group_by(recordid) %>% summarise(all_none = all(treatment=="None"), .groups="drop")
#     keep_non_none <- long %>% inner_join(all_none, by="recordid") %>%
#       filter(!all_none, treatment!="None") %>%
#       distinct(recordid, treatment) %>% rename(def_treatment = treatment)
#     keep_outcome <- all_none %>% filter(all_none) %>%
#       left_join(df_outcomes, by="recordid") %>%
#       mutate(def_treatment = ifelse(is.na(outcomes), "None", outcomes)) %>%
#       select(recordid, def_treatment) %>% distinct()
#     bind_rows(keep_non_none, keep_outcome)
#   }
#   
#   list(emp_final = emp_final_local,
#        def_final = def_final_local,
#        emp_def_wide = emp_def_wide_local)
# }
build_window_sets_for_ids <- function(id_vec, priority_vec){
  
  master_ids_df <- df_baseline %>%
    select(recordid, inf_onset, spec_date) %>%
    filter(recordid %in% id_vec) %>%
    distinct()
  
  res_local <- daily_antibiotics(
    df = df_data %>% dplyr::filter(recordid %in% id_vec),
    master_ids_df = master_ids_df,
    id_col = "recordid",
    group_col = "anti_used",
    start_col = "anti_start",
    end_col = "anti_end",
    windows = list(
      emp = list(ref_col = "inf_onset", offsets = c(-1, 0, 1, 2)),
      def = list(ref_col = "spec_date", offsets = c(3, 4, 5))
    ),
    same_day_priority = priority_vec
  )
  
  emp_def_wide_local <- res_local$wide
  
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
      pivot_longer(cols = all_of(emp_cols), names_to = "day", values_to = "treatment") %>%
      mutate(
        treatment = ifelse(is.na(treatment) | treatment == "", "None", treatment),
        treatment = vapply(treatment, normalize_daily_local, character(1))
      )
    
    all_none <- long %>%
      group_by(recordid) %>%
      summarise(all_none = all(treatment == "None"), .groups = "drop")
    
    keep_non_none <- long %>%
      inner_join(all_none, by = "recordid") %>%
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
      pivot_longer(cols = all_of(def_cols), names_to = "day", values_to = "treatment") %>%
      mutate(
        treatment = ifelse(is.na(treatment) | treatment == "", "None", treatment),
        treatment = vapply(treatment, normalize_daily_local, character(1))
      )
    
    all_none <- long %>%
      group_by(recordid) %>%
      summarise(all_none = all(treatment == "None"), .groups = "drop")
    
    keep_non_none <- long %>%
      inner_join(all_none, by = "recordid") %>%
      filter(!all_none, treatment != "None") %>%
      distinct(recordid, treatment) %>%
      rename(def_treatment = treatment)
    
    keep_outcome <- all_none %>%
      filter(all_none) %>%
      left_join(df_outcomes, by = "recordid") %>%
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

# Pair-level exclusion
# apply_pair_exclusion <- function(emp_df, def_df) {
#   special <- c("None","Died","Discharged")
#   joined <- emp_df %>%
#     full_join(def_df, by = "recordid", relationship = "many-to-many") %>%
#     mutate(
#       emp_treatment = if_else(is.na(emp_treatment), "None", emp_treatment),
#       def_treatment = if_else(is.na(def_treatment), "None", def_treatment)
#     ) %>%
#     # Drop Emp=None with Def in {None, Died, Discharged}
#     filter(!(emp_treatment == "None" & def_treatment %in% special))
#   
#   list(
#     joined = joined,
#     emp = joined %>% distinct(recordid, emp_treatment),
#     def = joined %>% distinct(recordid, def_treatment)
#   )
# }
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


# Global container for special nodes export
special_counts_global <- list()

# Sankey plotting with <5% -> Others + symmetric exception, and fixed totals
plot_sankey_based_for_ids <- function(id_vec, main_title, subset_key, space_value=40){
  special <- c("None","Died","Discharged")
  cazavi_node <- "Ceftazidime/avibactam"
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
          is_special ~ treatment,
          subset_key == "CRE" & grepl("Ceftazidime/avibactam", treatment, ignore.case = TRUE) ~ cazavi_node,
          is_combo ~ pick_base_by_priority_vec(treatment),
          TRUE ~ treatment
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
  
  # Keep only overall top 3 treatment nodes; all others -> Others
  top_n_keep <- 3
  
  overall_keep <- bind_rows(
    emp_rx_tbl0 %>% select(node_name, rx) %>% mutate(side = "emp"),
    def_rx_tbl0 %>% select(node_name, rx) %>% mutate(side = "def")
  ) %>%
    filter(!node_name %in% special) %>%
    group_by(node_name) %>%
    summarise(rx = sum(rx), .groups = "drop") %>%
    arrange(desc(rx), node_name) %>%
    slice_head(n = top_n_keep) %>%
    pull(node_name)
  
  emp_keep <- overall_keep
  def_keep <- overall_keep
  
  #
  remap_emp <- function(nm) {
    if (nm %in% special) return(nm)
    
    # For CRE: always keep CAZ-AVI as its own node
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", nm, ignore.case = TRUE)) {
      return(cazavi_node)
    }
    
    if (nm %in% emp_keep) nm else "Others"
  }
  
  remap_def <- function(nm) {
    if (nm %in% special) return(nm)
    
    # For CRE: always keep CAZ-AVI as its own node
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", nm, ignore.case = TRUE)) {
      return(cazavi_node)
    }
    
    if (nm %in% def_keep) nm else "Others"
  }
  
  # Build raw node names from treatments (before adding counts)
  build_node_name <- function(t, side = c("emp","def")){
    side <- match.arg(side)
    
    if (t %in% special) return(t)
    
    # CRE: any treatment containing CAZ-AVI -> same base node
    if (subset_key == "CRE" && grepl("Ceftazidime/avibactam", t, ignore.case = TRUE)) {
      return(cazavi_node)
    }
    
    # other combinations use base drug by priority
    if (grepl(",", t)) return(pick_base_by_priority_vec(t))
    
    t
  }
  
  joined_nodes <- joined %>%
    transmute(
      emp_node_raw = vapply(emp_treatment, build_node_name, character(1), side="emp"),
      def_node_raw = vapply(def_treatment, build_node_name, character(1), side="def")
    )
  
  # --- Build links ---
  joined_collapsed <- joined_nodes %>%
    mutate(
      emp_node = vapply(emp_node_raw, remap_emp, character(1)),
      def_node = vapply(def_node_raw, remap_def, character(1))
    ) %>%
    count(emp_node, def_node, name = "n")
  
  # --- SIDE COUNTS ---
  make_node_label <- function(node_name, n_mono, n_combo, rx, special){
    if (node_name %in% special) return(node_name)
    if (node_name == "Others")  return(paste0("Others (", rx, ")"))
    paste0(node_name, "-based (", rx, ")")
  }
  
  make_side_counts <- function(pairs_df, side = c("emp","def")) {
    side <- match.arg(side)
    
    raw_treat <- if (side == "emp") pairs_df$emp_treatment else pairs_df$def_treatment
    
    # node raw: exactly same logic as links
    node_raw <- vapply(raw_treat, build_node_name, character(1), side = side)
    
    # combo is determined by original treatment string
    is_combo_raw <- grepl(",", raw_treat)
    
    # combo signature (order-stabilised)
    combo_sig <- ifelse(is_combo_raw, norm_combo_vec(raw_treat), NA_character_)
    
    # apply remap (<5% -> Others) consistently
    node_final <- if (side == "emp") vapply(node_raw, remap_emp, character(1)) else vapply(node_raw, remap_def, character(1))
    
    dfu <- tibble(
      recordid   = pairs_df$recordid,
      node_name  = node_final,
      node_raw   = node_raw,
      is_combo   = is_combo_raw,
      combo_sig  = combo_sig
    ) %>%
      distinct(recordid, node_name, node_raw, is_combo, combo_sig)
    
    # ---- MONO counts ----
    # non-Others mono: per patient once
    mono_non_others <- dfu %>%
      filter(!is_combo, node_name != "Others") %>%
      group_by(node_name) %>%
      summarise(n_mono = n_distinct(recordid), .groups = "drop")
    
    # Others mono: keep "S:node_raw" so one patient can contribute multiple mono classes into Others (old behaviour)
    mono_others <- dfu %>%
      filter(!is_combo, node_name == "Others") %>%
      mutate(mono_sig = paste0("S:", node_raw)) %>%
      summarise(node_name = "Others",
                n_mono = n_distinct(paste(recordid, mono_sig, sep = "||")),
                .groups = "drop")
    
    mono_tbl <- bind_rows(mono_non_others, mono_others)
    
    # ---- COMBO counts ----
    # combo: per patient × combo signature (this matches your old “based combinations” magnitude)
    combo_tbl <- dfu %>%
      filter(is_combo) %>%
      mutate(combo_sig2 = paste0("C:", combo_sig)) %>%
      group_by(node_name) %>%
      summarise(n_combo = n_distinct(paste(recordid, combo_sig2, sep = "||")), .groups = "drop")
    
    full_join(mono_tbl, combo_tbl, by = "node_name") %>%
      mutate(
        n_mono  = tidyr::replace_na(n_mono, 0L),
        n_combo = tidyr::replace_na(n_combo, 0L),
        rx      = n_mono + n_combo
      )
  }
  
  # ---- build side tables (must have n_mono/n_combo/rx) ----
  emp_rx_tbl <- make_side_counts(emp_pairs, "emp")
  def_rx_tbl <- make_side_counts(def_pairs, "def")
  
  # ---- maps for label lookup ----
  emp_rx_map <- setNames(
    mapply(
      make_node_label,
      emp_rx_tbl$node_name,
      emp_rx_tbl$n_mono,
      emp_rx_tbl$n_combo,
      emp_rx_tbl$rx,
      MoreArgs = list(special = special)
    ),
    emp_rx_tbl$node_name
  )
  
  def_rx_map <- setNames(
    mapply(
      make_node_label,
      def_rx_tbl$node_name,
      def_rx_tbl$n_mono,
      def_rx_tbl$n_combo,
      def_rx_tbl$rx,
      MoreArgs = list(special = special)
    ),
    def_rx_tbl$node_name
  )
  
  safe_pick <- function(map, key, fallback) {
    out <- unname(map[key])
    if (is.null(out) || is.na(out)) fallback else out
  }
  
  label_emp <- function(t) {
    nm_raw <- build_node_name(t, side = "emp")
    nm_fin <- remap_emp(nm_raw)
    safe_pick(emp_rx_map, nm_fin, paste0(nm_fin, " (0)"))
  }
  
  label_def <- function(t) {
    nm_raw <- build_node_name(t, side = "def")
    nm_fin <- remap_def(nm_raw)
    safe_pick(def_rx_map, nm_fin, paste0(nm_fin, " (0)"))
  }
  
  #
  links_named <- joined_collapsed %>%
    transmute(
      emp_node = vapply(emp_node, label_emp, character(1)),
      def_node = vapply(def_node, label_def, character(1)),
      n
    )
  links_rep <- tidyr::uncount(links_named, weights = pmax(1, round(n/30))) %>%
    mutate(emp_node = as.character(emp_node), def_node = as.character(def_node))
  
  df_long <- ggsankey::make_long(links_rep, emp_node, def_node)
  
  # manual line breaks for long labels
  manual_wrap_label <- function(x) {
    x <- gsub(
      "Other anti-pseudomonal penicillin/beta-lactamase inhibitor-based",
      "Other anti-pseudomonal penicillin/\nbeta-lactamase inhibitor-based",
      x,
      fixed = TRUE
    )
    x <- gsub(
      "Anti-pseudomonal penicillin/beta-lactamase inhibitor-based",
      "Anti-pseudomonal penicillin/\nbeta-lactamase inhibitor-based",
      x,
      fixed = TRUE
    )
    x <- gsub(
      "Other beta-lactam/beta-lactamase inhibitor-based",
      "Other beta-lactam/\nbeta-lactamase inhibitor-based",
      x,
      fixed = TRUE
    )
    x <- gsub(
      "Sulfonamide-trimethoprim-combination-based",
      "Sulfonamide-trimethoprim/\ncombination-based",
      x,
      fixed = TRUE
    )
    x <- gsub(
      "Third-generation cephalosporin-based",
      "Third-generation\ncephalosporin-based",
      x,
      fixed = TRUE
    )
    x
  }
  
  df_long$label_wrap <- manual_wrap_label(as.character(df_long$node))
  
  # Visible node bases
  strip_paren <- function(s) sub(" \\(.*\\)$", "", s)
  node_base <- function(s) {
    s <- sub(" \\(.*\\)$", "", s)
    s <- sub("-based$", "", s)
    s
  }
  
  df_long$node_base <- node_base(as.character(df_long$node))
  
  # Totals for side headers (exclude special nodes)
  emp_node_rx_map <- setNames(emp_rx_tbl$rx, emp_rx_tbl$node_name)
  def_node_rx_map <- setNames(def_rx_tbl$rx, def_rx_tbl$node_name)
  
  visible_emp_nodes <- unique(node_base(as.character(df_long$node[df_long$x == "emp_node"])))
  visible_def_nodes <- unique(node_base(as.character(df_long$node[df_long$x == "def_node"])))
  
  keep_emp <- setdiff(visible_emp_nodes, special)
  keep_def <- setdiff(visible_def_nodes, special)
  
  keep_emp <- intersect(keep_emp, names(emp_node_rx_map))
  keep_def <- intersect(keep_def, names(def_node_rx_map))
  
  N_emp_rx <- sum(emp_node_rx_map[keep_emp], na.rm = TRUE)
  N_def_rx <- sum(def_node_rx_map[keep_def], na.rm = TRUE)
  
  # Unique treatments = distinct treatment names in daily table, excluding special nodes
  N_emp_unique <- emp_pairs %>%
    dplyr::rename(treatment = emp_treatment) %>%
    dplyr::filter(!treatment %in% special) %>%
    dplyr::distinct(treatment) %>%
    nrow()
  
  N_def_unique <- def_pairs %>%
    dplyr::rename(treatment = def_treatment) %>%
    dplyr::filter(!treatment %in% special) %>%
    dplyr::distinct(treatment) %>%
    nrow()
  
  # Title (n patients)
  N_patient <- length(unique(id_vec))
  
  make_title_html <- function(main_title, n, subset_key){
    title_core <- main_title
    
    if (grepl("Acinetobacter", main_title, ignore.case = TRUE)) {
      title_core <- "Carbapenem-resistant <i>Acinetobacter</i> spp."
    } else if (grepl("Pseudomonas", main_title, ignore.case = TRUE)) {
      title_core <- "Carbapenem-resistant <i>Pseudomonas</i> spp."
    } else if (grepl("Enterococcus", main_title, ignore.case = TRUE)) {
      title_core <- "Vancomycin-resistant <i>Enterococcus</i> spp."
    } else if (grepl("Staphylococcus aureus", main_title, ignore.case = TRUE)) {
      title_core <- "Methicillin-resistant <i>Staphylococcus aureus</i>"
    }
    
    sprintf(
      "<b>%s (%s; %s episodes)</b>",
      title_core,
      subset_key,
      format(n, big.mark = "", scientific = FALSE)
    )
  }
  
  title_html <- make_title_html(main_title, N_patient, subset_key)
  
  # Node ordering
  custom_orders <- list(
    CRA   = rev(c("Sulbactam","Polymyxin","Carbapenem","Anti-pseudomonal penicillin/beta-lactamase inhibitor")),
    `3GCRE`= rev(c("Carbapenem","Anti-pseudomonal penicillin/beta-lactamase inhibitor","Third-generation cephalosporin","Aminoglycoside")),
    CRE   = rev(c("Ceftazidime/avibactam","Polymyxin","Carbapenem","Anti-pseudomonal penicillin/beta-lactamase inhibitor","Aminoglycoside")),
    CRP   = rev(c("Anti-pseudomonal penicillin/beta-lactamase inhibitor","Polymyxin","Carbapenem","Aminoglycoside")),
    VRE   = rev(c("Lipopeptide","Oxazolidinone","Penicillin","Third-generation cephalosporin","Phosphonic")),
    MRSA  = rev(c("Glycopeptide","Lipopeptide","Oxazolidinone","Sulfonamide-trimethoprim-combination","Lincosamide"))
  )
  
  node_chr <- unique(as.character(df_long$node))
  
  special_nodes <- c(
    grep("^None( \\(|$)", node_chr, value = TRUE),
    grep("^Died( \\(|$)", node_chr, value = TRUE),
    grep("^Discharged( \\(|$)", node_chr, value = TRUE),
    grep("^Others( \\(|$)", node_chr, value = TRUE)
  )
  
  # base of each displayed node (drop counts + mono/combo suffix)
  node_base_chr <- node_base(node_chr)
  
  # custom base order for this organism
  custom_base_order <- custom_orders[[subset_key]]
  custom_base_order <- custom_base_order[!is.na(custom_base_order)]
  
  # nodes in custom order (keep the exact label variants present)
  ordered_nodes <- c()
  if (!is.null(custom_base_order)) {
    for (b in custom_base_order) {
      ordered_nodes <- c(ordered_nodes, node_chr[node_base_chr == b])
    }
  }
  
  # remaining nodes not in special/custom: alphabetical by base
  remaining_nodes <- setdiff(node_chr, c(special_nodes, ordered_nodes))
  remaining_bases <- sort(unique(node_base(remaining_nodes)))
  remaining_ordered <- c()
  for (b in remaining_bases) {
    remaining_ordered <- c(remaining_ordered, remaining_nodes[node_base(remaining_nodes) == b])
  }
  
  final_levels <- c(special_nodes, remaining_ordered, ordered_nodes)
  
  df_long$node <- factor(df_long$node, levels = final_levels)
  df_long$node_base <- node_base(as.character(df_long$node))
  
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
    
    nm0 <- sub(" \\(.*\\)$", "", unique_nodes[i]) 
    
    unname(col)
  }, character(1))
  names(node_colors) <- unique_nodes
  
  # Plot
  ggplot(
    df_long,
    aes(x = x, next_x = next_x, node = node, next_node = next_node,
        fill = node, label = label_wrap)
  ) +
    geom_sankey(flow.alpha = 0.75, smooth = 9, width = 0.05, alpha = 1, space = space_value) +
    geom_sankey_label(
      aes(
        x = stage(x, after_stat = x + 0.1 * dplyr::case_when(x == 1 ~ -0.27, x == 2 ~ 0.27, .default = 0)),
        hjust = dplyr::case_when(x == "emp_node" ~ 1, x == "def_node" ~ 0, .default = 0.5)
      ),
      size = 6, color = "black", fill = NA, fontface = "bold",
      label.r = unit(0.3, "mm"), space = space_value,
      family = "Times New Roman"
    ) +
    theme_void() +
    annotate("text", x = 0.97, y = Inf,
             label = paste0(N_emp_rx, " empirical prescriptions\n(",
                            N_emp_unique, " unique regimens)"),
             hjust = 1, vjust = 0.5, size = 6.5, fontface = "bold", family = "Times New Roman") +
    annotate("text", x = 2.03, y = Inf,
             label = paste0(N_def_rx, " definitive prescriptions\n(",
                            N_def_unique, " unique regimens)"),
             hjust = 0, vjust = 0.5, size = 6.5, fontface = "bold", family = "Times New Roman") +
    labs(title = title_html) +
    theme(
      plot.title = ggtext::element_markdown(hjust = 0.5, size = 20, 
                                            family = "Times New Roman", margin = margin(b = 40)),
      legend.position = "none",
      plot.margin = margin(t = 30, r = 10, b = 30, l = 10)
    ) +
    scale_fill_manual(values = node_colors) +
    coord_cartesian(clip = "off")
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

# heights <- vapply(subset_list, function(s) s$height, numeric(1))
# dir.create("output/figure", recursive = TRUE, showWarnings = FALSE)
# 
# labels <- rep(letters[1:3], length.out = length(plots)) 
# x_positions <- rep(c(0.095, 0.06, 0.06, 0.062, 0.062, 0.062), length.out = length(plots))
# 
# label_and_save <- function(p, label, x_pos, file, width = 16, height = 8) {
#   lab_g <- grid::textGrob(
#     paste0("(", label, ")"),
#     x = grid::unit(x_pos, "npc"),
#     y = grid::unit(0.9, "npc"),
#     just = c("left", "top"),
#     gp = grid::gpar(fontsize = 18, fontface = "bold", family = "Times New Roman")
#   )
#   g <- gridExtra::arrangeGrob(p, top = lab_g)
#   
#   # Use grDevices::cairo_pdf
#   grDevices::cairo_pdf(filename = file, width = width, height = height)
#   on.exit(grDevices::dev.off(), add = TRUE)
#   grid::grid.draw(g)
# }
# 
# for (i in seq_along(plots)) {
#   fn <- sprintf("output/figure/sankey_others_p%s.pdf", i)
#   label_and_save(
#     p = plots[[i]],
#     label = labels[i],
#     x_pos = x_positions[i],
#     file = fn,
#     width = 16,
#     height = heights[i]
#   )
# }
# 
# pdfcrop_bin <- Sys.which("pdfcrop")
# if (nzchar(pdfcrop_bin)) {
#   pdfs <- list.files("output/figure", pattern = "^sankey_others_p\\d+\\.pdf$", full.names = TRUE)
#   for (f in pdfs) system2(pdfcrop_bin, c(shQuote(f), shQuote(f)))
# } else {
#   message("pdfcrop not found; skipping cropping.")
# }

# Save
combined_all <- (plots[[1]] | plots[[2]]) /
  (plots[[3]] | plots[[4]]) /
  (plots[[5]] | plots[[6]])

cairo_pdf(
  "output/figure/sankey_others_all_combined.pdf",
  width = 32,
  height = 24
)

print(combined_all)
dev.off()



# # 3*2
# combined_all <- (plots[[1]] | plots[[2]] | plots[[3]]) /
#   (plots[[4]] | plots[[5]] | plots[[6]])
# 
# cairo_pdf(
#   "output/figure/sankey_others_all_combined.pdf",
#   width = 48,
#   height = 20
# )
# 
# print(combined_all)
# dev.off()
# 
# 
# 1-3
combined_1_3 <- plots[[1]] / plots[[2]] / plots[[3]]

cairo_pdf(
  "output/figure/sankey_others_p1_3_combined.pdf",
  width = 16,
  height = 24
)

print(combined_1_3)
dev.off()

