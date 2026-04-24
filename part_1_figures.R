# GAI Declaration:
# I, [201945294], declare that i have not employed ChatGPT in the creation of this script 
# to help structure and refine R and shell scripts used in the analysis, 
# debug code and resolve errors during implementation. 
# All outputs generated evaluated, tested, and modified where necessary to ensure correctness 
# and suitability for this project.

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
theme_set(theme_minimal(base_size = 12))

# fugure 1- PacBio HiFi assembly benchmarking
assembly_df <- data.frame(
  Sample = c("GN3", "GN6", "GN9"),
  Contigs = c(4, 4, 4),
  Total_Length = c(5243268, 5481946, 5320437),
  Largest_Contig = c(4947120, 5248224, 5028452),
  N50 = c(4947120, 5248224, 5028452),
  L50 = c(1, 1, 1),
  GC = c(50.66, 50.66, 50.82)
)

assembly_long <- assembly_df %>%
  select(Sample, Contigs, Total_Length, N50) %>%
  pivot_longer(cols = -Sample, names_to = "Metric", values_to = "Value")

p1 <- ggplot(assembly_long, aes(x = Sample, y = Value, fill = Metric)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(
    ~Metric,
    scales = "free_y",
    nrow = 1,
    labeller = as_labeller(c(
      Contigs = "Contigs (count)",
      N50 = "N50 (bp)",
      Total_Length = "Genome size (bp)"
    ))
  ) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Assembly benchmarking across PacBio HiFi clinical isolates",
    x = "Sample",
    y = "Metric value"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p1
ggsave("Figure1_assembly_metrics.png", p1, width = 10, height = 4, dpi = 300)


# sequencing technology comparison
tech_df <- data.frame(
  Dataset = c("GN3 PacBio HiFi 30x",
              "S HiFi long 30x",
              "S ONT long 10x",
              "S ONT long 30x",
              "S ONT long 100x",
              "S ONT raw long 30x",
              "S ONT short 30x"),
  Contigs = c(4, 3, 2, 4, 3, 2, 1471),
  N50 = c(4947120, 4879401, 4874661, 4875538, 4874929, 4874956, 2239)
)


# Contigs on log scale
p2a <- ggplot(tech_df, aes(x = reorder(Dataset, Contigs), y = Contigs)) +
  geom_col() +
  coord_flip() +
  scale_y_log10(labels = comma) +
  labs(
    title = "Contig count across sequencing datasets",
    x = "Dataset",
    y = "Contigs (log10 scale)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p2a
ggsave("Figure2a_contig_comparison.png", p2a, width = 9, height = 5, dpi = 300)

# N50 as separate figure
p2b <- ggplot(tech_df, aes(x = reorder(Dataset, N50), y = N50)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(labels = comma) +
  labs(
    title = "N50 across sequencing datasets",
    x = "Dataset",
    y = "N50 (bp)"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p2b
ggsave("Figure2b_n50_comparison.png", p2b, width = 9, height = 5, dpi = 300)



# Annotation comparison- Prokka vs Bakta

annot_df <- data.frame(
  Tool = c("Prokka", "Bakta"),
  CDS = c(4848, 4841),
  rRNA = c(22, 22),
  tRNA = c(90, 87),
  tmRNA = c(1, 1),
  Hypothetical_Proteins = c(1227, 228),
  Functional_CDS = c(3621, 4613)
)

annot_df$Tool <- factor(annot_df$Tool, levels = c("Prokka", "Bakta"))

annot_long <- annot_df %>%
  pivot_longer(cols = -Tool, names_to = "Metric", values_to = "Count")

p3 <- ggplot(annot_long, aes(x = Tool, y = Count, fill = Tool)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(
    ~Metric,
    scales = "free_y",
    nrow = 2,
    labeller = as_labeller(c(
      CDS = "CDS",
      Functional_CDS = "Functionally annotated CDS",
      Hypothetical_Proteins = "Hypothetical proteins",
      rRNA = "rRNA",
      tRNA = "tRNA",
      tmRNA = "tmRNA"
    ))
  ) +
  labs(
    title = "Comparison of GN3 annotation outputs from Prokka and Bakta",
    x = "Tool",
    y = "Count"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

p3
ggsave("Figure3_annotation_comparison.png", p3, width = 10, height = 6, dpi = 300)


# functional annotation percentage
func_df <- data.frame(
  Tool = c("Prokka", "Bakta"),
  Functional_Annotation_Percent = c(74.7, 95.3)
)

ggplot(func_df, aes(x = Tool, y = Functional_Annotation_Percent, fill = Tool)) +
  geom_col(width = 0.6) +
  ylim(0, 100) +
  labs(
    title = "Functional annotation coverage of GN3 CDS predictions",
    x = "Tool",
    y = "Functionally annotated CDS (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

    

    
    
    
    
    
    
    
    