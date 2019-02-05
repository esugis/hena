# This script draws the comparison plots 
#between genemania and integrated ds co-expression and ppis

library("UpSetR")

print("Reading in individual counts")
coexp_int <- read.table(file="results/comparisons/integrated_int_coexp_noduplicated_count.txt")
coexp_int <- coexp_int[1,1]

coexp_gm <- read.table(file="results/comparisons/gm_coexp_count.txt")
coexp_gm <- coexp_gm[1,1]

coexp_string <- read.table(file="results/comparisons/string_coexp_count.txt")
coexp_string <- coexp_string[1,1]

ppi_int <- read.table(file="results/comparisons/integrated_int_ppi_noduplicated_count.txt")
ppi_int <- ppi_int[1,1]

ppi_gm <- read.table(file="results/comparisons/gm_ppi_count.txt")
ppi_gm <- ppi_gm[1,1]

ppi_string <- read.table(file="results/comparisons/string_ppi_count.txt")
ppi_string <- ppi_string[1,1]


genint_gm <- read.table(file="results/comparisons/gm_genint_count.txt")
genint_gm <- genint_gm[1,1]

genint_int <- read.table(file="results/comparisons/integrated_int_epi_noduplicated_count.txt")
genint_int <- genint_int[1,1]


print("Reading in the numbers of overlapping pairs")

# Overlapping co-expression
coexp_int_gm <- read.table(file="results/comparisons/gm_intds_coexp_overlap_sort_count.txt")
coexp_int_gm <- coexp_int_gm[1,1]

coexp_int_string <- read.table(file="results/comparisons/string_intds_coexp_overlap_sort_count.txt")
coexp_int_string <- coexp_int_string[1,1]

coexp_gm_string <- read.table(file="results/comparisons/gm_string_coexp_overlap_sort_count.txt")
coexp_gm_string <- coexp_gm_string[1,1]

# Overlapping PPIs
ppi_int_gm <- read.table(file="results/comparisons/gm_intds_ppi_overlap_sort_count.txt")
ppi_int_gm <- ppi_int_gm[1,1] 

ppi_int_string <- read.table(file="results/comparisons/string_intds_ppi_overlap_sort_count.txt")
ppi_int_string <- ppi_int_string[1,1]

ppi_gm_string <- read.table(file="results/comparisons/gm_string_ppi_overlap_sort_count.txt")
ppi_gm_string <- ppi_gm_string[1,1]


# Overlapping genetic interactions
genint_int_gm <- read.table(file="results/comparisons/gm_intds_genint_overlap_sort_count.txt")
genint_int_gm <- genint_int_gm[1,1]

# Coexpression
expressionInput <- c(
                     "HENA" = coexp_int,
                     "GeneMania" = coexp_gm,
                     "STRING" = coexp_string,
                     "GeneMania&STRING" = coexp_gm_string,
                     "GeneMania&HENA" = coexp_int_gm,
                     "STRING&HENA"  = coexp_int_string,
                     "HENA&STRING&GeneMania" = 0)

# Order the sets 
pdf("results/comparisons/coexpression.pdf")
upset(fromExpression(expressionInput), sets= c("GeneMania","STRING","HENA"),
      keep.order = T,
      main.bar.color = "gray23",
      sets.bar.color = "#56B4E9",
      #sets.x.label = "Total number of interactions in the sets",
      text.scale=1.5,
    set_size.angles=0
)
dev.off()



################ PPI

######## PPI no scale
expressionInput <- c(
  "HENA" = ppi_int,
  "GeneMania" = ppi_gm,
  "STRING" = ppi_string,
  "GeneMania&STRING" = ppi_gm_string,
  "GeneMania&HENA" = ppi_int_gm,
  "STRING&HENA"  = ppi_int_string,
  "HENA&STRING&GeneMania" = 0)


# Order the sets 
pdf("results/comparisons/ppi.pdf")
upset(fromExpression(expressionInput),sets= c("GeneMania","STRING","HENA"),
      keep.order = T,
      main.bar.color = "gray23",
      sets.bar.color = "#56B4E9",
      #sets.x.label = "Total number of interactions in the sets",
# mainbar.y.label = "Number of intersections",
#     mb.ratio = c(0.7, 0.3),
      text.scale=1.5,
      set_size.angles=0
)

dev.off()

################# Genetic interactions
expressionInput <- c(
  "HENA" = genint_int,
  "GeneMania" = genint_gm,
  "GeneMania&HENA" = genint_int_gm)

# Order the sets 
pdf("results/comparisons/genetic_interactions.pdf")
upset(fromExpression(expressionInput), sets= c("GeneMania","HENA"),
      keep.order = T,
      main.bar.color = "gray23",
      sets.bar.color = "#56B4E9",
      #sets.x.label = "Total number of interactions in the sets",
#mainbar.y.label = "Number of intersections",
      text.scale=1.5,
      set_size.angles=0
)
dev.off()





