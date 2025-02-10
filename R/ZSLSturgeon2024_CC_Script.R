# ---
# Title: "ZSL Sturgeon eDNA Metabarcoding Results -  2024 Samples"
# Author: "Clare E. Collins using some previous scripts from Graham Sellers and Nathan P. Griffiths as a starter"
# Date: "20 December 2024"
# ---
# 
# ## Prepare working environment
# 
# Clear R memory <<rm(list=ls())>>, set working directory (check <<getwd()>>; Set <<setwd()>>) and load required packages.
# 

# To ensure reproducibility, print details about the version of R being used for analysis.
sessionInfo()
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19045)

## Install packages (including versions used for reproducibility)
library(plyr) #v1.8.9
library(dplyr) #V1.1.2
library(reshape) #v0.8.9
library(ggplot2) #v3.4.2
library(tidyverse) #v2.0.0


## Import raw Tapirs outputs for ZSL samples
ZSL2024 = read.csv2("data/ZSLsevern24_mifish_blast98_denoise.tsv", 
                     sep = "\t", header = TRUE)

## make '.' a '-'
colnames(ZSL2024) = gsub('\\.', '-', colnames(ZSL2024))

## Rename first column
colnames(ZSL2024)[1] = "Assignment"

## Remove underscore from species names
ZSL2024$Assignment = gsub("_", " ", ZSL2024$Assignment)

### Controls ####

## Create separate dataframe of controls
ZSLcontrols2024 = ZSL2024[grep('Assignment|NEG|POS|EB|FB', colnames(ZSL2024))]

## Remove species that do not feature in any of the controls
ZSLcontrols2024_species = ZSLcontrols2024[rowSums(ZSLcontrols2024[(2:14)]) > 0,]

## Turn table into long format
ZSLcontrols2024_species = cbind(ZSLcontrols2024_species[1], stack(ZSLcontrols2024_species[-1]))

## Change column names and add a column naming the run
colnames(ZSLcontrols2024_species)[2:3] = c("Reads","Controls") 
ZSLcontrols2024_species$Run = "Severn ZSL 2024"

## Order data by taxonomic assignment
ZSLcontrols2024_species = ZSLcontrols2024_species[order(ZSLcontrols2024_species$Assignment),]

## Create column specifying type of negative control
ZSLcontrols2024_species$Type = ifelse(grepl("FB", ZSLcontrols2024_species$Controls), "Field",
                                   ifelse(grepl("EB", ZSLcontrols2024_species$Controls), "Extraction",
                                          ifelse(grepl("NEG", ZSLcontrols2024_species$Controls), "Neg",
                                                               "Pos")))

## Create factor to order heatmap by
ZSLcontrols2024_species$fType = factor(ZSLcontrols2024_species$Type, 
                                    levels=c("Field","Extraction",
                                             "Neg", "Pos"))

## Plot contamination found in negative controls
options(scipen=999)
plot_control = ggplot(ZSLcontrols2024_species, aes(x=Controls, 
                                                y=fct_rev(as_factor(Assignment)), 
                                                fill=Reads))
plot_control = plot_control + geom_tile(colour="black")
plot_control = plot_control + scale_fill_gradientn(name="Total read counts", 
                                                   limits=c(0,100000),
                                                   breaks=c(0,50000, 100000),
                                                   colours=c("white","grey","black"), 
                                                   values=c(0,0.1,1))
plot_control = plot_control + labs(x="Process controls", y="Taxonomic assignment")
plot_control = plot_control + theme_bw()
plot_control = plot_control + theme(panel.grid.major = element_line(colour="white"),
                                    panel.grid.minor = element_line(colour="white"), 
                                    axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                                    axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                                    axis.ticks.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.text.y = element_text(colour="black", face = "italic"),
                                    text = element_text(size=10),
                                    legend.key.size = unit(1, 'lines'))
plot_control = plot_control + facet_grid(Run ~ fType, scale = "free", space = "free")
plot_control

png("output/ZSLSturgen2024_Unfiltered_Contamination.png", width = 5000, height = 2500, units = "px", res = 345)
plot_control
dev.off()


### Samples ####

## Make Assignment column the row names for these transposed data
rownames(ZSL2024) = ZSL2024$Assignment
ZSL2024 = ZSL2024[, -grep('Assignment', colnames(ZSL2024))]

## Sort by taxonomy
ZSL2024 = ZSL2024[order(ZSL2024$taxonomy),]

## isolate fish and others to continue working on fish in ZSL2024
ZSL2024_nonfish = ZSL2024[-grep('Actinopteri|Petromyzontiformes', ZSL2024$taxonomy),]
ZSL2024 = ZSL2024[grep('Actinopteri|Petromyzontiformes', ZSL2024$taxonomy),]

## Grab taxonomy, before removing column below
ZSL2024_tax = ZSL2024$taxonomy

## Remove negatives, positives, filtration and extraction blanks, and taxonomy
ZSL2024 = ZSL2024[, -grep('taxonomy|NEG|POS|EB|FB', colnames(ZSL2024))]

## Check positive doesn't appear in samples and if 0 remove positive control (Maylandia zebra)
ZSL2024_POSITIVEWARNING = NULL
if (all(ZSL2024['Maylandia zebra',] <= 0)) { 
  ## Remove the positive control row
  ZSL2024 = ZSL2024[!(rownames(ZSL2024) == 'Maylandia zebra'),]
} else {
  ZSL2024_POSITIVEWARNING = subset(ZSL2024, rownames(ZSL2024) == 'Maylandia zebra')
}


## Now remove positive row if still there
ZSL2024 = ZSL2024[!(rownames(ZSL2024) == 'Maylandia zebra'),]

## Noise filter data at 0.1% (but record data removed in a df ZSL2024_noisefiltered):
## This function transposes (t) df (swapping columns and rownames) and divides each element transposed in ZSL2024 by the corresponding column sum, normalising the columns of the matrix and then checks whether each normalised element is less (or more) than 0.001 returning logical (TRUE/FALSE), the true values are returned as 0 and the false values are left as before in the corresponding df
ZSL2024_noisefiltered = ZSL2024
ZSL2024_noisefiltered[t(t(ZSL2024)/colSums(ZSL2024)) > 0.001] = 0
ZSL2024[t(t(ZSL2024)/colSums(ZSL2024)) < 0.001] = 0
write.csv(ZSL2024_noisefiltered, "output/ZSL2024_noisefiltered.csv")

## Drop (but record - species in these dropped samples only) samples with less than 500 total (i.e. sum across all species) reads:
ZSL2024_samplesdropped = ZSL2024[!colSums(ZSL2024) > 500]
ZSL2024_samplesdropped = ZSL2024_samplesdropped[rowSums(ZSL2024_samplesdropped[])>0,]
write.csv(ZSL2024_samplesdropped,'output/ZSL2024_droppedless500reads.csv')

ZSL2024 = ZSL2024[!colSums(ZSL2024) < 500]

## Remove assignments to Cyprinidae, Percidae and unassigned
ZSL2024 = ZSL2024[-grep('Cyprinidae|Percidae|Leuciscidae|unassigned', rownames(ZSL2024)),]

## Change Alosa_alosa to Alosa spp.
rownames(ZSL2024) = sub('Alosa alosa', 'Alosa spp.', rownames(ZSL2024))

## Change Chelon to Chelon spp.
rownames(ZSL2024) = sub('Chelon', 'Chelon spp.', rownames(ZSL2024))

## Merge Leuciscus to Leuciscus leuciscus and change to Leuciscus spp.
Leuciscus_leuciscus = ZSL2024['Leuciscus leuciscus',]
Leuciscus = ZSL2024['Leuciscus',]
ZSL2024['Leuciscus leuciscus',] = Leuciscus_leuciscus + Leuciscus
ZSL2024 = ZSL2024[!(rownames(ZSL2024) == 'Leuciscus'),]
rownames(ZSL2024) = sub('Leuciscus leuciscus', 'Leuciscus spp.', rownames(ZSL2024))

## change Lampetra fluviatilis to Lampetra_spp. as could be L. fluviatilis or L. planeri
rownames(ZSL2024) = sub('Lampetra fluviatilis', 'Lampetra spp.', rownames(ZSL2024))

## remove (but record) species with no assignments
ZSL2024_speciesnoassignments = ZSL2024[rowSums(ZSL2024) == 0,]
write.csv(ZSL2024_speciesnoassignments, 'output/ZSL2024_speciesnoassignments.csv')
ZSL2024 = ZSL2024[rowSums(ZSL2024) > 0,]

### Bubble Plot ####

## Transpose and sort (create new df, copy rownames to a column, move that column to the start and remove rownames, then transpose df)
bubble_df = ZSL2024
bubble_df$assignment = rownames(bubble_df)
r = dim(bubble_df) # get dimensions
bubble_df = bubble_df[,c((r[2]),1:(r[2]-1))] # move last column to start of df using the dimensions r
rownames(bubble_df) = NULL
bubble_df = data.frame(t(bubble_df))

## Assign assignment row as column names and delete assignment row
colnames(bubble_df) = bubble_df['assignment',]
bubble_df = bubble_df[!(rownames(bubble_df) == 'assignment'),]

## Assign rownames column as site column, move it to the start and delete the rownames
bubble_df$site = rownames(bubble_df)
r = dim(bubble_df) # get dimensions
bubble_df = bubble_df[,c((r[2]),1:(r[2]-1))] # move last column to start of df using the dimensions r
rownames(bubble_df) = NULL

## Remove sample number from site ID
bubble_df$site = gsub("^(.{4}).", "\\1", bubble_df$site)

## Convert characters to numeric for all columns except the first
bubble_df[,-1] = lapply(bubble_df[,-1], function(x) as.numeric(as.character(x)))

## Merge (sum) reads by site (e.g. D1S1A and D1S1B will be merged as now both labelled as D1S1) PACKAGE plyr
bubble_df = ddply(bubble_df, .(site), numcolwise(sum))

## Add assigned reads column (get dimensions (rows and columns) of the df; print the number of columns and use this to set the larger number needed to add reads after rowSums)
r=dim(bubble_df)
bubble_df$assigned = rowSums(bubble_df[2:(r[2])])

## Change site names so samples run downstream in river Severn and easily identify Wye:
bubble_df$site = ifelse(grepl("D3", bubble_df$site), paste0("S1", bubble_df$site),
                         ifelse(grepl("D2", bubble_df$site), paste0("S2", bubble_df$site),
                                ifelse(grepl("D1", bubble_df$site), paste0("S3", bubble_df$site),
                                       ifelse(grepl("D4", bubble_df$site), paste0("W1", bubble_df$site),
                                              bubble_df$site))))

## Re-sort site column and Add site as rownames (and then remove site column)
bubble_df = bubble_df[order(bubble_df$site),]
rownames(bubble_df) = bubble_df$site
bubble_df = bubble_df[, -grep('site', colnames(bubble_df))]


## Make proportional reads eDNA dataset 
ZSL_prop_reads_df = bubble_df # create new df
ZSL_prop_reads_df = ZSL_prop_reads_df/ZSL_prop_reads_df$assigned # make proportions
r = dim(ZSL_prop_reads_df) # get df dimensions
ZSL_prop_reads_df = ZSL_prop_reads_df[, 1:(r[2]-1)] # use r above to remove end column
ZSL_prop_reads_df$ZSL_site = rownames(ZSL_prop_reads_df) # move rownames to column
r = dim(ZSL_prop_reads_df) # get df dimensions
ZSL_prop_reads_df = ZSL_prop_reads_df[,c((r[2]),1:(r[2]-1))] # use dimensions r to move end column to front
rownames(ZSL_prop_reads_df) = NULL
dim(ZSL_prop_reads_df)
ZSL_prop_reads_df = melt(ZSL_prop_reads_df, id=(1)) ## melt from reshape PACKAGE to transform df


## Change variable to taxonomic
colnames(ZSL_prop_reads_df) = sub('variable', 'taxonomic', colnames(ZSL_prop_reads_df))


## Make bubble plot
### Make Bubble Plot Figure
##remove values of 0 and then order correctly (should already be, but just in case)
ZSL_set1 = subset(ZSL_prop_reads_df, ZSL_prop_reads_df$value>0)
# ZSL_set1$ZSL_site = as.factor(ZSL_set1$ZSL_site)

## Add common names, importing them from dictionary file, matching to taxonomic names and then inserting into this df
common_names_df = read.csv('data/severn2024_CommonNames.txt', sep = '\t', header = F, col.names = c('taxonomic', 'common'))
taxonomic_common = match(ZSL_set1$taxonomic, common_names_df$taxonomic)
ZSL_set1$common = common_names_df$common[taxonomic_common]

## Add sampling events column to ZSL set and column order in both
ZSL_set1$event = ZSL_set1$ZSL_site
ZSL_set1 = ZSL_set1[,c(1, 5, 2, 4, 3)]

## Change site from sample ID/event to Severn or Wye
ZSL_set1$ZSL_site = ifelse(grepl("S[1-3]D", ZSL_set1$ZSL_site), "Severn",
                         ifelse(grepl("W1D", ZSL_set1$ZSL_site), "Wye",
                                ZSL_set1$ZSL_site))

## Create column with day sampled
ZSL_set1$day = ifelse(grepl("S3D", ZSL_set1$event), "25 June 2024",
                      ifelse(grepl("S2D", ZSL_set1$event), "26 June 2024",
                             ifelse(grepl("S1D", ZSL_set1$event), "27 June 2024",
                                    ifelse(grepl("W1D", ZSL_set1$event), "28 June 2024",
                                           ZSL_set1$event
                             ))))

## Create simple sample names that match the map
ZSL_set1$sample <- ifelse(grepl("S1D3S[1-8]", ZSL_set1$event), 
                          paste0("S3-", sub(".*S([1-8])", "\\1", ZSL_set1$event)),
                          ifelse(grepl("S2D2S[1-8]", ZSL_set1$event), 
                                 paste0("S2-", sub(".*S([1-8])", "\\1", ZSL_set1$event)),
                                 ifelse(grepl("S3D1S[1-8]", ZSL_set1$event), 
                                        paste0("S1-", sub(".*S([1-8])", "\\1", ZSL_set1$event)),
                                        ifelse(grepl("W1D4S[1-5]", ZSL_set1$event), 
                                               paste0("S4-", sub(".*S([1-5])", "\\1", ZSL_set1$event)),
                                               ZSL_set1$event))))


## Check species richness per event (alpha diversity) and save to file for report
ZSL_sr = table(ZSL_set1$event)
ZSL_sr = as.data.frame(ZSL_sr)
write.csv(ZSL_sr, file = "output/ZSLlocationSR.csv", row.names = FALSE)

## Check species richness per river

# Aggregate to count the number of unique taxonomic values per ZSL_site
ZSL_river_sr = aggregate(taxonomic ~ ZSL_site, data = ZSL_set1, FUN = function(x) length(unique(x)))

# Save to file for report
write.csv(ZSL_river_sr, file = "output/ZSLriverSR.csv", row.names = FALSE)

## Common name spatial bubbleplot

## Using PACKAGE library(dplyr) try to get the coordinates of the river Severn and Wye data to identify them in the plot
# Ensure that sample is treated as an ordered factor across all data
ZSL_set1$sample = factor(ZSL_set1$sample, levels = unique(ZSL_set1$sample))

# Calculate the xmin and xmax for each river based on the full sample order
river_shading = ZSL_set1 %>%
  group_by(ZSL_site) %>%
  reframe(  ## Keeps ZSL_site automatically!
    xmin = min(as.numeric(sample)) - 0.5,
    xmax = max(as.numeric(sample)) + 0.5,
    xmid = (min(as.numeric(sample)) + max(as.numeric(sample))) / 2
  )

## Add fill colours for rivers
river_shading = river_shading %>%
  mutate(fill_colour = case_when(
    ZSL_site == "Severn" ~ "lightblue",
    ZSL_site == "Wye" ~ "lightgrey",
    TRUE ~ "white" ## fallback colour
  ))

## Convert ZSL_set1 to a factor with levels based on unique values in the common field
ZSL_set1$common = factor(ZSL_set1$common, levels = unique(ZSL_set1$common))
## Reorder 'sample' based on 'event'
ZSL_set1$sample = factor(ZSL_set1$sample, levels = unique(ZSL_set1$sample[order(ZSL_set1$event)]))


## Ensure the factor levels are sorted alphabetically by common name
ZSL_set1$common = factor(ZSL_set1$common, levels = sort(levels(ZSL_set1$common)))
ZSL_set1 = ZSL_set1[order(ZSL_set1$common), ]

## Create the plot
ZSL_com_bubble_plot = ggplot (ZSL_set1, aes(x=sample, y=common, fill=day, size=value*100)) +

## Add background shading for River Severn and River Wye based on calculated xmin and xmax
  geom_rect(data = river_shading, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = river_shading$fill_colour, alpha = 0.2, inherit.aes = FALSE) +

## Add the bubbles for fish data    
  geom_point(alpha=0.9, shape=21, color="black") +
  
  ## Set bubble size to be relative to % of sample
  scale_size(range = c(1, 10), name="Relative fish reads (%)") +
  
  ## Set bubble fill to match sampling days and also river colour
  scale_fill_manual(values = c(     "25 June 2024" = "#440154",
                                    "26 June 2024" = "#2a768e",
                                    "27 June 2024" = "#70cf57",
                                    "28 June 2024" = "#fde725"
                                    ), 
                                    guide = guide_legend(override.aes = list(size = 5)))+  ## Increase dot size in the legend
  
  ## Set theme and customize appearance
  theme_bw() +
  labs(fill = "Sampling Date") + ## Change legend title
  theme(legend.position = "right", plot.margin = margin(t = 17, r = 5, b = 5, l = 5))+  ## Align legend right and add extra space at the top for River names
  xlab("Sampling Site") +
  ylab("") +
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(size = 9),
        axis.title.x = element_text(size = 10))+
  scale_y_discrete(limits = rev(levels(ZSL_set1$common))) +
  scale_x_discrete(limits = levels(ZSL_set1$sample)) +
  
  ## Add titles for each river
  annotate("text", 
           x = river_shading$xmid[river_shading$ZSL_site == "Severn"], 
           y = 29.1,
             # max(as.numeric(factor(ZSL_set1$common))) + 5, 
           label = "River Severn", 
           size = 4, 
           hjust = 0.5 
           # vjust = 1, 
           # fontface = "bold",
           # colour = "midnightblue"
           ) +
  annotate("text", 
           x = river_shading$xmid[river_shading$ZSL_site == "Wye"], 
           y = 29.1,
             # max(as.numeric(factor(ZSL_set1$common))) + 5,
           label = "River Wye", 
           size = 4, 
           hjust = 0.5 
           # vjust = 1, 
           # fontface = "bold",
           # colour = "darkgrey"
           ) +  # River Wye title
  ## Ensure annotations are not clipped
  coord_cartesian(ylim = c(1,28), clip = "off")


png("output/ZSL_Fish2024_Summary_common.png", width = 4500, height = 2500, units = "px", res = 345)
ZSL_com_bubble_plot
dev.off()

