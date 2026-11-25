# =============================================================================
# Microbial Co-occurrence Network Reconstruction (SPIEC-EASI protocol)
# MiMB Chapter — Bacterial Pangenomics
#
# Demonstration dataset: two functionally distinct tanks from a WWTP
#   - Liquid waste
#   - Denitrification
# Restricted to ONE season (Winter) to keep sample sizes adequate for
# SPIEC-EASI / StARS stability selection (see MIN_SAMPLES below).
#
# Pipeline: phyloseq subsetting 
#           -> tax_glom (Genus) -> prevalence filtering
#           -> SPIEC-EASI network fitting per tank (method = "mb")
#           -> degree extraction -> Spearman correlation
#           -> Plot A (node centrality), Plot B (degree vs. abundance),
#              Plot C (ggraph network topology)
#
# All user-adjustable settings are grouped in Section 2 ("Parameters"),
# so the protocol can be re-run on a different dataset, season, tank pair,
# or genus subset without editing code further down.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load required packages
# -----------------------------------------------------------------------------
setwd("path/to/your/project/")

library(SpiecEasi)

# Data handling / phyloseq
library(phyloseq)
library(tidyverse)

# Network objects and plotting
library(igraph)
library(ggraph)
library(tidygraph)
library(patchwork)

# Plot formatting helpers
library(scales)
library(ggrepel)

# Utility
library(tictoc)

# -----------------------------------------------------------------------------
# 2. Parameters (edit these to adapt the protocol to a new run)
# -----------------------------------------------------------------------------

# Load the pyloseq object for the tutorial
load("network_tutorial_phyloseq_object.RData")

# --- Subsetting ---
SEASON <- "Winter"                              # season to analyze
tanks  <- c("Liquid waste", "Denitrification")  # tanks to compare

# --- Filtering thresholds ---
PREVALENCE_THRESHOLD <- 0.10   # keep genera present in >= 10% of samples
MIN_SAMPLES <- 10              # minimum samples per tank required to fit a
# network; below this, StARS stability
# selection becomes unreliable and the tank
# is skipped rather than fit anyway

# --- SPIEC-EASI fitting ---
SE_METHOD <- "mb"              # Meinshausen-Buhlmann neighborhood selection
SE_LAMBDA_MIN_RATIO <- 1e-2#3
SE_NLAMBDA <- 20#50
SE_REP_NUM <- 50                # StARS subsampling replicates


# --- Plot aesthetics ---
tank_colors <- c("Liquid waste" = "#9A32CD", "Denitrification" = "#79CDCD")


# --- Genus selection for Plots A and B ---
# Choose whether these two plots are restricted to the smallest set of genera that jointly account for AB_CUMULATIVE_THRESHOLD of total
# relative abundance, or show every genus retained after filtering.

FOCAL_ONLY <- TRUE
AB_CUMULATIVE_THRESHOLD <- 0.6


# -----------------------------------------------------------------------------
# 2. Subset to the tanks/season of interest
# -----------------------------------------------------------------------------
# Reducing the phyloseq object BEFORE tax_glom() keeps the agglomeration step fast and avoids collapsing taxa that will be discarded anyway.

keep_samples <- sample_data(ps_tutorial)$Season == SEASON

ps_tut <- prune_samples(keep_samples, ps_tutorial)

# Drop taxa that are all-zero after removing the other samples — this is what actually shrinks the object before tax_glom.
ps_tut <- prune_taxa(taxa_sums(ps_tut) > 0, ps_tut)

cat("Samples kept:", nsamples(ps_tut), "\n")
cat("Taxa kept (pre-glom):", ntaxa(ps_tut), "\n")
print(table(data.frame(sample_data(ps_tut))$Tank))


# -----------------------------------------------------------------------------
# 4. Aggregate to genus level
# -----------------------------------------------------------------------------

ps_genus <- tax_glom(ps_tut, taxrank = "Genus")


# -----------------------------------------------------------------------------
# 5. Filter low-prevalence taxa
# -----------------------------------------------------------------------------
# Genera present in fewer than PREVALENCE_THRESHOLD of samples are dropped.
# This reduces the compositional sparsity that otherwise destabilizes SPIEC-EASI's neighborhood selection.

ps_filtered <- filter_taxa(ps_genus,
                           function(x) sum(x > 0) > (PREVALENCE_THRESHOLD * length(x)),
                           prune = TRUE)

cat("Taxa after filtering:", ntaxa(ps_filtered), "\n")
cat("Samples:", nsamples(ps_filtered), "\n")

# Confirm per-tank sample sizes meet MIN_SAMPLES before proceeding.
sampling_info_numbers <- data.frame(sample_data(ps_filtered)) %>%
  filter(Season == SEASON) %>%
  group_by(Tank, Season) %>%
  summarise(n = n(), .groups = "drop")

print(sampling_info_numbers)


# -----------------------------------------------------------------------------
# 6. Compute mean relative abundance per genus/tank
# -----------------------------------------------------------------------------
# Relative abundance is calculated on the full filtered object (before splitting by tank) so normalization is consistent across tanks.

ps_rel <- transform_sample_counts(ps_filtered, function(x) x / sum(x))

abund_df <- psmelt(ps_rel) %>%
  group_by(Genus, Phylum, Tank) %>%
  summarise(MeanRelAbund = mean(Abundance), .groups = "drop")

toc()

# -----------------------------------------------------------------------------
# 7. Network reconstruction per tank
# -----------------------------------------------------------------------------
# Each tank is fit as an independent SPIEC-EASI network (rather than one joint network) so that per-tank topology can be compared directly. 
# Tanks with fewer than MIN_SAMPLES samples are skipped rather than fit anyway, since StARS stability selection is unreliable at very small n.

degree_results <- list()
graph_results  <- list()

tic("Total network reconstruction")

for (t in tanks) {
  
  tic(t)
  
  # --- Fit the network ---
  keep <- sample_data(ps_filtered)$Tank == t & sample_data(ps_filtered)$Season == SEASON
  ps_sub <- prune_samples(keep, ps_filtered)
  ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
  
  n_samp <- nsamples(ps_sub)
  n_taxa <- ntaxa(ps_sub)
  cat("  |", n_samp, "samples,", n_taxa, "taxa")
  
  if (n_samp < MIN_SAMPLES) {
    cat("  --> SKIPPED (too few samples)\n")
    next
  }
  
  se_net <- tryCatch(
    spiec.easi(ps_sub,
               method = SE_METHOD,
               lambda.min.ratio = SE_LAMBDA_MIN_RATIO,
               nlambda = SE_NLAMBDA,
               pulsar.params = list(rep.num = SE_REP_NUM)),
    error = function(e) {
      cat("  --> SpiecEasi ERROR:", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (is.null(se_net)) next
  
  ig <- adj2igraph(getRefit(se_net))
  
  # --- Extract and label results ---
  # Attach genus/phylum labels to nodes. make.unique() guards against duplicate genus names within a tank (e.g. paraphyletic genera split
  # across multiple ASV clusters after tax_glom).
  genus_names  <- as.character(tax_table(ps_sub)[, "Genus"])
  phylum_names <- as.character(tax_table(ps_sub)[, "Phylum"])
  genus_names_unique <- make.unique(genus_names, sep = "_dup")
  V(ig)$name <- genus_names_unique
  
  graph_results[[t]] <- ig
  
  node_degrees <- degree(ig)
  network_density <- edge_density(ig)
  
  degree_df <- data.frame(Genus = genus_names,
                          GenusUnique = genus_names_unique,
                          Phylum = phylum_names,
                          Degree = node_degrees,
                          RelDegree = node_degrees / (n_taxa - 1),  # normalized 0-1, comparable across tanks of different size
                          Tank = t,
                          Season = SEASON,
                          NetDensity = network_density,
                          N_samples = n_samp,
                          stringsAsFactors = FALSE)
  
  degree_results[[t]] <- degree_df
  cat("  --> OK | network density:", round(network_density, 3), "\n")
  toc()
  
}

toc()


# -----------------------------------------------------------------------------
# 8. Combine network results with abundance data
# -----------------------------------------------------------------------------
tic()
final_degree_df <- bind_rows(degree_results)
cat("\nTotal genus-network entries:", nrow(final_degree_df), "\n")

final_df <- left_join(final_degree_df, abund_df, by = c("Genus", "Phylum", "Tank")) %>%
  mutate(LogAbund = log10(MeanRelAbund + 1e-6),
         Tank = factor(Tank, levels = tanks))

# write.csv(final_df, "input_files/degree_vs_abundance_twotank_winter.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# 9. Statistical test: degree vs. abundance correlation
# -----------------------------------------------------------------------------
# Spearman correlation per tank, testing whether highly connected genera also tend to be highly abundant. BH-adjusted across tanks.

cor_results <- final_df %>%
  group_by(Tank) %>%
  filter(!is.na(MeanRelAbund), !is.na(RelDegree), n() >= 5) %>%
  summarise(rho = cor(MeanRelAbund, RelDegree, method = "spearman", use = "complete.obs"),
            pval = cor.test(MeanRelAbund, RelDegree, method = "spearman", exact = FALSE)$p.value,
            n_taxa = n(),
            .groups = "drop") %>%
  mutate(padj = p.adjust(pval, method = "BH"),
         sig = case_when(padj < 0.001 ~ "***",
                         padj < 0.01  ~ "**",
                         padj < 0.05  ~ "*",
                         TRUE ~ "ns"),
         label = paste0("\u03c1=", round(rho, 2), ", n=", n_taxa, "\n", sig))

print(cor_results)
# write.csv(cor_results, "input_files/spearman_twotank_winter.csv", row.names = FALSE)
toc()

# -----------------------------------------------------------------------------
# 10. Plot 1 — Node centrality (relative degree) per genus, faceted by Tank
# -----------------------------------------------------------------------------
# Rank genera by relative abundance WITHIN EACH TANK, and keep the minimal top set that reaches the threshold in that tank. Ranking per tank (rather
# than on a combined/pooled total) matters here because Liquid waste and Denitrification have distinct dominant taxa — a combined ranking would let
# whichever tank has the larger raw abundances crowd out genera that are locally dominant but globally modest in the other tank. The two per-tank
# sets are then unified so a genus important in either tank is retained in both facets for comparison
tic()
# --- Liquid waste: genera needed to reach 95% ---
genera_liquidwaste_rank <- abund_df %>% filter(Tank == "Liquid waste") %>%
  arrange(desc(MeanRelAbund)) %>% mutate(CumFrac = cumsum(MeanRelAbund) / sum(MeanRelAbund))

focal_genera_liquidwaste <- genera_liquidwaste_rank %>% 
  filter(row_number() <= (which(CumFrac >= AB_CUMULATIVE_THRESHOLD)[1])) %>%
  pull(Genus)

cat("Liquid waste:",length(unique(focal_genera_liquidwaste)), "genera account for",
    AB_CUMULATIVE_THRESHOLD * 100, "% of tank abundance\n")

# --- Denitrification: genera needed to reach 95% ---
genera_denitrification_rank <- abund_df %>% filter(Tank == "Denitrification") %>%
  arrange(desc(MeanRelAbund)) %>% mutate(CumFrac = cumsum(MeanRelAbund) / sum(MeanRelAbund))

focal_genera_denitrification <- genera_denitrification_rank %>% 
  filter(row_number() <= (which(CumFrac >= AB_CUMULATIVE_THRESHOLD)[1])) %>%
  pull(Genus)

cat("Denitrification:",length(unique(focal_genera_denitrification)), "genera account for",
    AB_CUMULATIVE_THRESHOLD * 100, "% of tank abundance\n")

# --- Labels ---
n_labels <- data.frame(Tank = tanks,
                       n_taxa = c(length(focal_genera_liquidwaste), length(focal_genera_denitrification)),
                       label = paste0("n = ",c(length(focal_genera_liquidwaste), length(focal_genera_denitrification))))


# Build a lookup table instead
focal_lookup <- bind_rows(tibble(Tank = "Liquid waste",     Genus = focal_genera_liquidwaste),
                          tibble(Tank = "Denitrification",  Genus = focal_genera_denitrification))
# Restricted to focal_genera_abundant when FOCAL_ONLY = TRUE.
plot_a_data <- if (FOCAL_ONLY) {
  final_df %>% semi_join(focal_lookup, by = c("Tank", "Genus"))
} else {
  final_df
}


ggplot(plot_a_data, aes(x = RelDegree, y = reorder(Genus, RelDegree,decreasing = T), colour = Tank)) +
  geom_point(size = 3, alpha = 0.85, na.rm = TRUE) +
  scale_colour_manual(values = tank_colors) +
  facet_wrap(~Tank, ncol = 2, scales = "free_y") +
  geom_text(data = n_labels, aes(label = label),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.7,
            inherit.aes = FALSE, size = 3.2, color = "grey20") +
  theme_minimal(base_size = 11) +
  labs(title = "Node centrality",
       subtitle = paste0("Genera accounting for ",
                         if (FOCAL_ONLY) paste0(AB_CUMULATIVE_THRESHOLD * 100, "%") else "all"," of tank abundance in ", SEASON),
       x = "Relative degree", y = NULL, colour = "Tank") +
  theme(axis.text.y = element_text(face = "italic", size = 8),
        strip.text = element_text(size = 9, face = "bold"),
        legend.position = "bottom")



# -----------------------------------------------------------------------------
# 11. Plot 2 — Network topology (ggraph), Liquid waste vs. Denitrification
# -----------------------------------------------------------------------------
# Converts each tank's igraph object to a tidygraph, attaches node-level attributes (degree, relative degree, phylum) from final_df so
# color/size mappings stay consistent with Plots 1 and 3, and optionally restricts nodes to focal_genera_wwtp via the same FOCAL_ONLY option.

ig_liquidwaste     <- graph_results[["Liquid waste"]]
ig_denitrification <- graph_results[["Denitrification"]]

make_tidygraph <- function(ig, tank_name, focal_only = FALSE, focal_lookup = NULL) {
  node_attrs <- final_df %>%
    filter(Tank == tank_name) %>%
    select(GenusUnique, Genus, Phylum, Degree, RelDegree)
  
  tg <- as_tbl_graph(ig) %>% activate(nodes) %>%
    left_join(node_attrs, by = c("name" = "GenusUnique"))
  
  if (focal_only) { 
    tank_focal_genera <- focal_lookup %>% filter(Tank == tank_name) %>% pull(Genus)
    tg <- tg %>% filter(Genus %in% tank_focal_genera) }
  
  tg
}

tg_liquidwaste <- make_tidygraph(ig_liquidwaste, "Liquid waste",
                                 focal_only = FOCAL_ONLY, focal_lookup = focal_lookup)
tg_denitrification <- make_tidygraph(ig_denitrification, "Denitrification",
                                     focal_only = FOCAL_ONLY, focal_lookup = focal_lookup)

# Shared layout seed so both panels are visually comparable.
set.seed(42)

# Node labels: when showing every genus (FOCAL_ONLY = FALSE), labeling all of them gets cluttered, so only the top 10% by RelDegree are labeled.
# When FOCAL_ONLY = TRUE the node set is already small and curated, so every node is labeled.
plot_network <- function(tg, title, fill_color, label_all) {
  ggraph(tg, layout = "fr") +
    geom_edge_link(colour = "grey75", edge_width = 0.3, alpha = 0.6) +
    geom_node_point(aes(size = RelDegree), colour = fill_color, alpha = 0.85) +
    {if (label_all) { geom_node_text(aes(label = Genus), repel = TRUE, size = 2.8, fontface = "italic",
                     segment.size = 0.2, max.overlaps = 20)
    }else{geom_node_text(aes(label = ifelse(RelDegree > quantile(RelDegree, 0.90, na.rm = TRUE),Genus, NA)),
                     repel = TRUE, size = 2.6, fontface = "italic",segment.size = 0.2, max.overlaps = 15)}} +
    scale_size_continuous(range = c(1, 8), name = "Relative\ndegree") +
    theme_graph(base_family = "sans") +
    labs(title = title) +
    theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9))
}

p_liquidwaste <- plot_network(
  tg_liquidwaste,
  paste0("Liquid waste | n = ", gorder(tg_liquidwaste), " nodes, ", gsize(tg_liquidwaste), " edges"),
  tank_colors["Liquid waste"],
  label_all = FOCAL_ONLY
)

p_denitrification <- plot_network(
  tg_denitrification,
  paste0("Denitrification | n = ", gorder(tg_denitrification), " nodes, ", gsize(tg_denitrification), " edges"),
  tank_colors["Denitrification"],
  label_all = FOCAL_ONLY
)


plot_2 <- p_liquidwaste + p_denitrification +
  plot_annotation(title = "Co-occurrence network topology by tank",
                  subtitle = paste0("Genera accounting for ",
                                    if (FOCAL_ONLY) paste0(AB_CUMULATIVE_THRESHOLD * 100, "%") else "all"," of tank abundance in ", SEASON),
                  theme = theme(plot.title = element_text(size = 13, face = "bold")))

plot_2


# -----------------------------------------------------------------------------
# 12. Plot 3 — Degree vs. abundance scatter, faceted by Tank
# -----------------------------------------------------------------------------
# Plot 3 always uses the full filtered genus set (final_df), since it's a global abundance-degree relationship rather than a per-genus display.

ggplot(final_df, aes(x = MeanRelAbund, y = RelDegree)) +
  geom_point(alpha = 0.6, size = 1.8, color = "gray70") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              linewidth = 0.5, linetype = "dashed") +
  geom_text(data = cor_results, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
            inherit.aes = FALSE, size = 3.2, color = "grey10") +
  facet_wrap(~Tank, ncol = 2, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(title = paste0("Network degree vs. Mean relative abundance in ", SEASON),
       x = "Mean relative abundance", y = "Relative network degree")

toc()
