# Analysis registry ---------------------------------------------------------------
# One entry per analysis offered in the app. To add an analysis, add an entry:
#   label     : name shown in the menu and in plot titles
#   group     : menu section
#   transform : key of TRANSFORM_LABELS (fct_transform.R)
#   distance  : key of compute_distance() (fct_distance.R); also used by PERMANOVA
#   dist_label: how the distance is named in the summary
#   engine    : pcoa | pca | ca | nmds | umap  (fct_ordination.R)
#   prefix    : axis name prefix
#   drivers_transform (optional): transformation used for the Drivers tab when
#               the distance does its own transformation (e.g. Aitchison)

METHODS <- list(
  pcoa_bray = list(
    label = "PCoA – Bray-Curtis", group = "PCoA",
    transform = "total", distance = "bray", dist_label = "Bray-Curtis",
    engine = "pcoa", prefix = "PCoA"
  ),
  pcoa_jaccard = list(
    label = "PCoA – Jaccard", group = "PCoA",
    transform = "pa", distance = "jaccard", dist_label = "Jaccard",
    engine = "pcoa", prefix = "PCoA"
  ),
  # Aitchison distances are computed by vegdist() directly on raw counts;
  # `drivers_transform` is the transformation used to rank the driving taxa.
  pcoa_aitchison = list(
    label = "PCoA – Aitchison", group = "PCoA",
    transform = "none", distance = "aitchison", dist_label = "Aitchison",
    engine = "pcoa", prefix = "PCoA", drivers_transform = "clr"
  ),
  pcoa_robust_aitchison = list(
    label = "PCoA – robust Aitchison", group = "PCoA",
    transform = "none", distance = "robust.aitchison", dist_label = "Robust Aitchison",
    engine = "pcoa", prefix = "PCoA", drivers_transform = "rclr"
  ),
  pca_clr = list(
    label = "PCA – CLR", group = "PCA",
    transform = "clr", distance = "euclidean", dist_label = "Aitchison",
    engine = "pca", prefix = "PC"
  ),
  pca_rclr = list(
    label = "PCA – robust CLR", group = "PCA",
    transform = "rclr", distance = "euclidean", dist_label = "Robust Aitchison",
    engine = "pca", prefix = "PC"
  ),
  pca_hellinger = list(
    label = "PCA – Hellinger", group = "PCA",
    transform = "hellinger", distance = "euclidean", dist_label = "Hellinger",
    engine = "pca", prefix = "PC"
  ),
  ca = list(
    label = "Correspondence analysis", group = "Other",
    transform = "none", distance = "chisq", dist_label = "Chi-square",
    engine = "ca", prefix = "CA"
  ),
  nmds_bray = list(
    label = "NMDS – Bray-Curtis", group = "Other",
    transform = "total", distance = "bray", dist_label = "Bray-Curtis",
    engine = "nmds", prefix = "NMDS"
  ),
  umap_bray = list(
    label = "UMAP – Bray-Curtis", group = "UMAP",
    transform = "total", distance = "bray", dist_label = "Bray-Curtis",
    engine = "umap", prefix = "UMAP"
  ),
  umap_clr = list(
    label = "UMAP – CLR", group = "UMAP",
    transform = "clr", distance = "euclidean", dist_label = "Aitchison",
    engine = "umap", prefix = "UMAP"
  ),
  umap_rclr = list(
    label = "UMAP – robust CLR", group = "UMAP",
    transform = "rclr", distance = "euclidean", dist_label = "Robust Aitchison",
    engine = "umap", prefix = "UMAP"
  )
)

# Engines that need the distance matrix to compute the ordination itself
ENGINES_NEED_DIST <- c("pcoa", "nmds")

# Choices for selectInput(), grouped by menu section
method_choices <- function() {
  groups <- unique(vapply(METHODS, `[[`, "", "group"))
  stats::setNames(lapply(groups, function(g) {
    ids <- names(METHODS)[vapply(METHODS, `[[`, "", "group") == g]
    stats::setNames(ids, vapply(METHODS[ids], `[[`, "", "label"))
  }), groups)
}
