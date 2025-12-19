# Clear & packages
rm(list = ls())

suppressPackageStartupMessages({
  require(pacman)
  pacman::p_load(
    magrittr, dplyr, tidyr, lubridate, purrr,
    ggplot2, ggtext, RColorBrewer, stringr, scales,
    Cairo, openxlsx
  )
})

# Working directory & data
wd <- "./"; setwd(wd)
df_data     <- readRDS("data/clean data/anti_treat_index.RData")
df_baseline <- readRDS("data/clean data/baseline_outcomes_index.RData")

# Keep rows with known infection type
df_baseline <- df_baseline[!is.na(df_baseline$infection_types), ]

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

priority_all_for <- function(key, pool_values) {
  base <- priority_map[[key]]
  others <- setdiff(sort(unique(pool_values)), base)
  c(base, others)
}

# Helpers (merge/signature/windows)
normalize_group <- function(x) gsub("\\s*-\\s*", "-", x)

combine_same_day <- function(v, priority=NULL){
  v <- unique(na.omit(v))
  if (!length(v)) return(NA_character_)
  v <- normalize_group(v)
  if (is.null(priority) || !length(priority)) return(paste(sort(v), collapse=", "))
  ord <- match(v, priority)
  max_ord <- suppressWarnings(max(ord, na.rm = TRUE)); if(!is.finite(max_ord)) max_ord <- 0
  ord[is.na(ord)] <- max_ord + seq_len(sum(is.na(ord)))
  paste(v[order(ord)], collapse=", ")
}

union_names <- function(x){
  x <- x[!is.na(x) & x != ""]
  if (!length(x)) return("")
  toks <- unique(unlist(strsplit(x, "\\s*\\|\\s*")))
  toks <- sort(unique(trimws(toks[toks!=""])))
  if (!length(toks)) "" else paste(toks, collapse=" | ")
}

canonical_class_sig <- function(s, priority){
  s <- ifelse(is.na(s) | s=="", "None", s)
  if (s == "None") return("None")
  toks <- strsplit(s, "\\s*,\\s*")[[1]]
  toks <- normalize_group(toks)
  ord  <- match(toks, priority)
  maxo <- suppressWarnings(max(ord, na.rm = TRUE)); if(!is.finite(maxo)) maxo <- 0
  ord[is.na(ord)] <- maxo + seq_len(sum(is.na(ord)))
  paste(toks[order(ord)], collapse = ", ")
}

daily_class_and_names <- function(df,
                                  id_col="recordid",
                                  class_col="anti_used",
                                  name_col="anti_names",
                                  start_col="anti_start",
                                  end_col  ="anti_end",
                                  windows=list(
                                    emp=list(ref_col="inf_onset", offsets=c(-1,0,1,2)),
                                    def=list(ref_col="spec_date", offsets=c(3,4,5))
                                  ),
                                  class_priority=NULL){
  stopifnot(all(c(id_col,class_col,name_col,start_col,end_col) %in% names(df)))
  to_date <- function(x) as.Date(x)
  
  df <- df %>% mutate(
    !!start_col := to_date(.data[[start_col]]),
    !!end_col   := to_date(.data[[end_col]])
  )
  
  rx_base <- df %>%
    distinct(.data[[id_col]], .data[[class_col]], .data[[name_col]], .data[[start_col]], .data[[end_col]])
  
  per_day <- rx_base %>%
    filter(!is.na(.data[[start_col]]), !is.na(.data[[end_col]]), .data[[end_col]] >= .data[[start_col]]) %>%
    rowwise() %>%
    mutate(day_date = list(seq(.data[[start_col]], .data[[end_col]], by="day"))) %>%
    unnest(day_date) %>% ungroup() %>%
    group_by(.data[[id_col]], day_date) %>%
    summarise(
      class_sig = combine_same_day(.data[[class_col]], priority = class_priority),
      name_set  = {
        nm <- unique(na.omit(.data[[name_col]]))
        nm <- nm[nm != ""]
        paste(sort(nm), collapse=" | ")
      },
      .groups="drop"
    )
  
  mk_grid <- function(ref_col, offsets, prefix){
    ids <- df %>% select(all_of(c(id_col, ref_col))) %>% distinct() %>%
      mutate(!!ref_col := to_date(.data[[ref_col]])) %>% filter(!is.na(.data[[ref_col]]))
    ids %>% rowwise() %>%
      mutate(.dates = list(.data[[ref_col]] + lubridate::days(offsets)), .offs = list(offsets)) %>%
      unnest(c(.dates, .offs)) %>% ungroup() %>%
      mutate(win = prefix, rel = .offs) %>%
      rename(day_date = .dates, ref_date = !!ref_col)
  }
  grids <- purrr::imap(windows, ~ mk_grid(.x$ref_col, .x$offsets, .y)) %>% bind_rows()
  
  daily_long <- grids %>%
    left_join(per_day, by = setNames(c(id_col,"day_date"), c(id_col,"day_date"))) %>%
    mutate(
      class_sig = tidyr::replace_na(class_sig, "None"),
      name_set  = tidyr::replace_na(name_set, "")
    ) %>%
    select(all_of(c(id_col)), day_date, win, rel, class_sig, name_set)
  
  daily_long
}

make_phase_pairs <- function(long_df, phase = c("emp","def"), class_priority){
  phase <- match.arg(phase)
  dfp <- long_df %>% filter(win == phase)
  dfp <- dfp %>% mutate(class_sig_c = vapply(class_sig, canonical_class_sig, character(1), priority = class_priority))
  
  out <- dfp %>%
    group_by(recordid) %>%
    group_modify(~{
      dd <- .x
      if (all(dd$class_sig_c == "None")) {
        tibble(treatment_class = "None", names_union = "")
      } else {
        dd2 <- dd %>% filter(class_sig_c != "None")
        dd2 %>%
          group_by(class_sig_c) %>%
          summarise(names_union = union_names(name_set), .groups="drop") %>%
          rename(treatment_class = class_sig_c)
      }
    }) %>%
    ungroup()
  out
}

# Pairwise exclusion: drop records where Emp=None and Def in {None, Died, Discharged}
apply_pair_exclusion_pairs <- function(emp_df, def_df) {
  special <- c("None","Died","Discharged")
  joined <- emp_df %>%
    dplyr::rename(emp_treatment = treatment_class,
                  emp_names     = names_union) %>%
    dplyr::full_join(
      def_df %>% dplyr::rename(def_treatment = treatment_class,
                               def_names     = names_union),
      by = "recordid",
      relationship = "many-to-many"
    ) %>%
    dplyr::mutate(
      emp_treatment = dplyr::if_else(is.na(emp_treatment), "None", emp_treatment),
      def_treatment = dplyr::if_else(is.na(def_treatment), "None", def_treatment),
      emp_names     = tidyr::replace_na(emp_names, ""),
      def_names     = tidyr::replace_na(def_names, "")
    ) %>%
    dplyr::filter(!(emp_treatment == "None" & def_treatment %in% special))
  
  list(
    joined   = joined,
    emp_keep = joined %>% dplyr::select(recordid, treatment_class = emp_treatment, names_union = emp_names) %>% dplyr::distinct(),
    def_keep = joined %>% dplyr::select(recordid, treatment_class = def_treatment, names_union = def_names) %>% dplyr::distinct()
  )
}

# anti_name -> AWaRe (pick the "strictest")
aware_priority <- c("Not recommended", "Reserve","Watch","Access","Unknown")

name_aware_tbl <- df_data %>%
  transmute(
    anti_name  = anti_names,
    WHO_AWaRe  = ifelse(is.na(WHO_AWaRe), "Unknown", WHO_AWaRe)
  ) %>%
  filter(!is.na(anti_name) & anti_name != "") %>%
  distinct() %>%
  group_by(anti_name) %>%
  summarise(
    WHO_AWaRe = {
      v <- unique(WHO_AWaRe)
      idx <- match(v, aware_priority)
      v[which.min(ifelse(is.na(idx), length(aware_priority)+1, idx))]
    },
    .groups="drop"
  )
aware_by_name <- setNames(name_aware_tbl$WHO_AWaRe, name_aware_tbl$anti_name)

names_to_aware <- function(names_union){
  if (is.na(names_union) || names_union=="") return("Unknown")
  toks <- unique(unlist(strsplit(names_union, "\\s*\\|\\s*")))
  toks <- toks[toks != ""]
  if (!length(toks)) return("Unknown")
  v <- unname(aware_by_name[toks])
  v[is.na(v)] <- "Unknown"
  u <- unique(v)
  idx <- match(u, aware_priority)
  u[which.min(ifelse(is.na(idx), length(aware_priority)+1, idx))]
}


# Define pathogen subsets
aci_car_id  <- df_data$recordid[df_data$aci_car == 1]
ent_thir_id <- df_data$recordid[df_data$ent_thir == 1]
ent_car_id  <- df_data$recordid[df_data$ent_car == 1]
pse_car_id  <- df_data$recordid[df_data$pse_car == 1]
entc_van_id <- df_data$recordid[df_data$entc_van == 1]
sa_meth_id  <- df_data$recordid[df_data$sa_meth == 1]

subset_list <- list(
  list(key="CRA",   ids = aci_car_id,  title = "Carbapenem-resistant Acinetobacter spp."),
  list(key="3GCRE", ids = ent_thir_id, title = "Third-generation cephalosporin-resistant Enterobacterales"),
  list(key="CRE",   ids = ent_car_id,  title = "Carbapenem-resistant Enterobacterales"),
  list(key="CRP",   ids = pse_car_id,  title = "Carbapenem-resistant Pseudomonas spp."),
  list(key="VRE",   ids = entc_van_id, title = "Vancomycin-resistant Enterococcus spp."),
  list(key="MRSA",  ids = sa_meth_id,  title = "Methicillin-resistant Staphylococcus aureus")
)
plot_order <- c("MRSA","VRE","CRP","CRE","3GCRE","CRA")


# For each pathogen: run EMP/DEF with its priority → pairwise exclusion → map AWaRe → count
make_awr_counts <- function(pairs_df, aware_col, class_col, ids, side_label, pathogen_key){
  tab <- pairs_df %>%
    mutate(across(all_of(c(aware_col,class_col)), as.character)) %>%
    filter(recordid %in% ids) %>%
    filter(.data[[class_col]] != "None") %>%
    transmute(AWaRe = .data[[aware_col]]) %>%
    count(AWaRe, name = "counts") %>%
    arrange(desc(counts))
  denom <- sum(tab$counts[tab$AWaRe != "Unknown"])
  tab %>%
    mutate(
      total_counts = denom,
      percentage   = ifelse(AWaRe=="Unknown", NA_real_, 100*counts/denom),
      df_name      = side_label,
      id_list_name = pathogen_key
    )
}

all_chunks <- list()

for (s in subset_list) {
  ids <- s$ids
  key <- s$key
  prio_vec <- priority_all_for(key, pool_values = df_data$anti_used)
  
  daily_long_sub <- daily_class_and_names(
    df = df_data %>% dplyr::filter(recordid %in% ids),
    id_col="recordid",
    class_col="anti_used",
    name_col="anti_names",
    start_col="anti_start",
    end_col  ="anti_end",
    windows = list(
      emp=list(ref_col="inf_onset", offsets=c(-1,0,1,2)),
      def=list(ref_col="spec_date", offsets=c(3,4,5))
    ),
    class_priority = prio_vec
  )
  
  emp_pairs <- make_phase_pairs(daily_long_sub, phase = "emp", class_priority = prio_vec)
  def_pairs <- make_phase_pairs(daily_long_sub, phase = "def", class_priority = prio_vec)
  
  filt <- apply_pair_exclusion_pairs(emp_pairs, def_pairs)
  emp_pairs_filt <- filt$emp_keep %>% mutate(emp_AWaRe = vapply(names_union, names_to_aware, character(1)))
  def_pairs_filt <- filt$def_keep %>% mutate(def_AWaRe = vapply(names_union, names_to_aware, character(1)))
  
  all_chunks[[length(all_chunks)+1]] <- make_awr_counts(emp_pairs_filt, "emp_AWaRe", "treatment_class", ids, "Empirical prescriptions",  key)
  all_chunks[[length(all_chunks)+1]] <- make_awr_counts(def_pairs_filt, "def_AWaRe", "treatment_class", ids, "Definitive prescriptions", key)
}

# Combine data for plotting
# (drop "Unknown" to keep it out of the legend)
all_combine <- bind_rows(all_chunks) %>%
  filter(AWaRe != "Unknown") %>%                      
  mutate(
    WHO_AWaRe    = factor(AWaRe, levels = c("Not recommended","Access","Watch","Reserve")),
    id_list_name = factor(id_list_name, levels = plot_order),
    df_name      = factor(df_name, levels = c("Empirical prescriptions","Definitive prescriptions"))
  ) %>%
  select(WHO_AWaRe, counts, percentage, total_counts, df_name, id_list_name)


# Plot AWaRe stacked bars (Unknown removed; legend excludes it)
p <- ggplot(all_combine,
            aes(y = id_list_name, x = percentage, fill = WHO_AWaRe)) +
  geom_bar(stat = "identity", width = 0.7, na.rm = TRUE) +
  facet_grid( ~ df_name, scales = "free", space = "free") +
  scale_fill_manual(
    values = c("Access" = "#319539",
               "Watch"  = "#f4b122",
               "Reserve"= "#f10c2d",
               "Not recommended" = "#956131"),
    breaks = c("Reserve","Watch","Access","Not recommended")    # <-- no "Unknown" here
  ) +
  labs(y = "", x = "", fill = "WHO AWaRe category") +
  geom_text(aes(label = ifelse(is.na(percentage), "", counts)),
            position = position_stack(vjust = 0.5),
            colour = "black", size = 3.5, family = "Times New Roman", na.rm = TRUE) +
  scale_x_continuous(expand = c(0.005,0), labels = scales::percent_format(scale = 1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.background  = element_blank(),
        panel.border      = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black", 
                                   family = "Times New Roman", margin = margin(t = 4, unit = "pt")),
        axis.text.y = element_text(size = 12, colour = "black", family = "Times New Roman"),
        strip.text.x = element_text(size = 12, colour = "black", margin = margin(0.35,0,0.35,0, "cm")),
        legend.position = "bottom",
        legend.box.spacing = unit(0.1, "cm"),
        legend.margin = margin(t = -5, b = 2),
        legend.title = element_text(size = 12, family = "Times New Roman",
                                    colour = "black", hjust = 0.5, vjust = 1,
                                    margin = margin(r = 8, t = 4, unit = "pt")),
        legend.text = element_text(size = 12, family = "Times New Roman",
                                   colour = "black",
                                   margin = margin(l = 3, r = 6, unit = "pt")),
        legend.spacing = unit(0.5, "cm"),
        plot.margin = margin(0, 30, 0, 0),
        panel.spacing = unit(1, "cm"),
        strip.background = element_rect(color = "black", fill = "lightgrey", linewidth = NA))

# Save
CairoPDF("output/figure/AWR_index.pdf", width = 17, height = 5)
print(p)
dev.off()
