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
