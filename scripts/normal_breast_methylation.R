
library(minfi)
library(ComplexHeatmap)

# Loading data (GSE69914)
my_files <- list.files("/Volumes/Data/Project_3/normal_breast_methylation/GSE69914_RAW/", pattern = "Raw.txt$")

# Loadk Illumina 450K annotation
illumina450K_annotation <- read.csv("/Volumes/Data/Project_3/normal_breast_methylation/humanmethylation450_15017482_v1-2.csv",
                                    skip = 6)

# Filtering files. Get only normal tissue
files_to_keep <- paste0("XX0XX0XXXXXXXXXXXX0XXXXXX0XXXXX0XXXXXXXXXXXXX0XXXX",
        "XX0XXX0XXXXXXXXXXX0XXX00XXXXXXXX0XXXXXXXXXXXXXX0XX",
        "XXXXX0XXXXXXXXXXXXXXXX00XXXXXXXXXXXXX0XXXXXXXXXXXX",
        "XXXXXXXXXXXXXXXXXX0XXXXX0XXXXXXXXX0XXXXXXXXXXXXX00",
        "XXXXXXX0XXX0XXXX0XXXXXXXXXXXXXX0XXXXXXXXXX0XXXXX0X",
        "XXX0XXXX0XXXXXXX0XXXXXXXXXXXXXX00XXXXXXXXXX00XX0XX",
        "XXXXXXXXXXX0XXX00XXXXXXXXX0XXXXXX0XXXXXXXXXXXX0XXX",
        "XXXX0X0XXXXXXXXX0XXXXXXXXXXX0XXXXXXXXX0XXX0XXXXXX0",
        "XXX0XXX")
files_to_keep <- strsplit(files_to_keep, split="")[[1]]

# samples selection within a series
indices_files <- which(files_to_keep != "X")  # eliminate samples marked as "X
my_files <- my_files[indices_files]

# Readong files and generating matrix

# Define colnames based on a matrix
example <- read.table(paste0("/Volumes/Data/Project_3/normal_breast_methylation/GSE69914_RAW/", "GSM1712369_BCFD3_Raw.txt"))

normal_tissue_methylation <- data.frame(row.names = unname(example$V1))

for(file in my_files) {
  
  my_file <- read.table(paste0("/Volumes/Data/Project_3/normal_breast_methylation/GSE69914_RAW/", file))
  column_name <- strsplit(file, split = "_")[[1]][1]
  
  normal_tissue_methylation[unname(my_file$V1), column_name] <- unname(my_file$V2)                 
}

normal_tissue_methylation <- normal_tissue_methylation[-1,]


gbp4_cpgs <- names(genes)[genes == "GBP4"]

gbp4_cpgs <- gbp4_cpgs[gbp4_cpgs %in% rownames(normal_tissue_methylation)]

data_subset <- data.frame(lapply(normal_tissue_methylation[gbp4_cpgs,], function(x) as.numeric(as.character(x))))
rownames(data_subset) <- gbp4_cpgs


Heatmap(data_subset,
        show_row_names = TRUE,       # Display row names
        show_column_names = FALSE,   # Hide column names
        col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))  # Define color gradient
)

# Adjust margins 
par(mar = c(6, 4, 2, 2)) 

# Create the boxplot

# Compute median values for each category
medians <- apply(data_subset, 1, median, na.rm = TRUE)

# Define a color gradient: blue (0), white (0.5), red (1)
color_palette <- colorRampPalette(c("blue", "white", "red"))(100)

# Map medians to the color palette
median_scaled <- medians * 100 # Rescale to range 1-100
box_colors <- color_palette[round(median_scaled)] # Assign colors based on median

# Plot with category-specific colors
boxplot(t(data_subset), 
        las = 2, 
        xlab = "", 
        ylab = "Beta Value", 
        col = box_colors,  # Assign colors based on median
        ylim = c(0,1), 
        cex.axis = 0.8, 
        cex.lab = 1.2,  
        cex.main = 1.4, 
        border = "darkblue", 
        notch = FALSE, 
        frame = FALSE)



