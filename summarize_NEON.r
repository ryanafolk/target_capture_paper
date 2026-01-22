# Load CSVs
extractions <- read.csv("./metadata/intbio_metadata_draft.csv")
sitemap <- read.csv("./metadata/NEON-SiteMap-Table.csv")
ontarget <- read.csv("./metadata/capture_qc.csv")

# Change working directory to base level directory of GitHub repository local copy

##########
# Extraction map
##########
names(extractions)[12] <- "siteCode"

# Suppress non-sample designations
#extractions$siteCode <- gsub(" .*","",gsub("-.*", "", extractions$Sample.ID, perl = TRUE))


# Full list of site codes (look for errors)
unique(extractions$siteCode)

library(dplyr)

# Add site information from NEON website
extractions <- left_join(extractions, sitemap, by="siteCode")


library(ggplot2)
library(maps)

world_map <- map_data("world")

# Summarize counts
summary_df <- extractions %>% count(longitude, latitude, siteCode)
summary_df <- na.omit(summary_df)

# Prepare map
ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray60") +
  coord_cartesian(xlim = c(-160, -50), ylim = c(10, 90)) +
  geom_point(data = summary_df, aes(x = longitude, y = latitude, size = n), color = "steelblue", alpha = 0.7) +
  scale_size_continuous(range = c(3, 10)) +  # Adjust size scale
  theme_minimal() +
  labs(title = "Extractions by NEON site",
       x = "Longitude", y = "Latitude", size = "Occurrences")
  
         
##########
# Qubit results
##########

# Plot
extractions_clean <- extractions
extractions_clean$Total.DNA[extractions_clean$Total.DNA == "#DIV/0!"] <- 0
extractions_clean$Total.DNA <- as.numeric(extractions_clean$Total.DNA)
                 
ggplot(extractions_clean, aes(x = Domain, y = as.numeric(Total.DNA))) +
  #geom_violin(fill = "lightblue", color = "gray40", na.rm = TRUE) + scale_y_continuous(trans=scales::pseudo_log_trans(base = 10)) + 
  geom_boxplot(fill = "white", outlier.shape = NA) +  # inset boxplot
  theme_minimal() +
  labs(title = "DNA extraction result vs. domain", x = "Domain", y = "Total DNA yield (ng)") + 
  coord_cartesian(ylim = c(0, 1000))

##########
# On-target results
##########

# Add sample metadata.
names(ontarget)[1] <- "Sample.ID"
ontarget <- left_join(ontarget, extractions, by="Sample.ID")

ontarget <- ontarget[ontarget$Type %in% c("nodule", "root", "rhizosphere"),]
ontarget$X16S.on.target.. <- as.numeric(gsub("%", "", ontarget$X16S.on.target..))
ontarget <- ontarget[ontarget$X16S.on.target.. < 100, ]
ontarget$ITS.on.target.. <- as.numeric(gsub("%", "", ontarget$ITS.on.target..))
ontarget <- ontarget[ontarget$ITS.on.target.. < 100, ]
ontarget$Symbiosis.on.target.. <- as.numeric(gsub("%", "", ontarget$Symbiosis.on.target..))
ontarget <- ontarget[ontarget$Symbiosis.on.target.. < 100, ]

# Box plots 16S
ggplot(ontarget, aes(x = Type, y = X16S.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(aes(alpha = Family == "Fabaceae"), position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20") +
  	scale_alpha_manual(values = c(`TRUE` = 0.03, `FALSE` = 0.15), guide = "none") +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of 16S\nvs sample type",
		x = "Sample type", y = "% on target"
	)

ggplot(ontarget, aes(x = Family, y = X16S.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(aes(alpha = Family == "Fabaceae"), position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20") +
  	scale_alpha_manual(values = c(`TRUE` = 0.03, `FALSE` = 0.15), guide = "none") +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of 16S vs.\nhost family",
		x = "Family", y = "% on target"
	)
	
ggplot(ontarget, aes(x = Domain, y = X16S.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20", alpha = 0.05) +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of 16S vs. NEON domain",
		x = "Domain", y = "% on target"
	)  

 
  
# Box plots ITS
ggplot(ontarget, aes(x = Type, y = ITS.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(aes(alpha = Family == "Fabaceae"), position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20") +
  	scale_alpha_manual(values = c(`TRUE` = 0.03, `FALSE` = 0.15), guide = "none") +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of ITS\nvs sample type",
		x = "Sample type", y = "% on target"
	)

ggplot(ontarget, aes(x = Family, y = ITS.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(aes(alpha = Family == "Fabaceae"), position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20") +
  	scale_alpha_manual(values = c(`TRUE` = 0.03, `FALSE` = 0.15), guide = "none") +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of ITS vs.\nhost family",
		x = "Family", y = "% on target"
	)

ggplot(ontarget, aes(x = Domain, y = ITS.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20", alpha = 0.05) +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of ITS vs. NEON domain",
		x = "Domain", y = "% on target"
	)  

# Box plots symbiosis genes
ggplot(ontarget, aes(x = Type, y = Symbiosis.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(aes(alpha = Family == "Fabaceae"), position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20") +
  	scale_alpha_manual(values = c(`TRUE` = 0.03, `FALSE` = 0.15), guide = "none") +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of\nsymbiosis genes vs sample type",
		x = "Sample type", y = "% on target"
	)
	
ggplot(ontarget[ontarget$Type == "nodule", ], aes(x = Family, y = Symbiosis.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(aes(alpha = Family == "Fabaceae"), position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20") +
  	scale_alpha_manual(values = c(`TRUE` = 0.1, `FALSE` = 0.25), guide = "none") +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of symbiosis genes vs.\nhost family (nodules only)",
		x = "Family", y = "% on target"
	)

  
ggplot(ontarget[ontarget$Type == "nodule", ], aes(x = Domain, y = Symbiosis.on.target..)) +
	geom_boxplot(width = 0.5, fill = "white", outlier.shape = NA, color = "gray20") +
	geom_point(position = position_jitter(width = 0.12, height = 0), size = 1.0, color = "gray20", alpha = 0.1) +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of symbiosis genes vs. NEON domain (nodules only)",
		x = "Domain", y = "% on target"
	)



# Scatter plots of different on-target results
# Also investigated nodule only. Looks the same so the relationship does not depend on sample type

m <- lm(X16S.on.target.. ~ ITS.on.target.., data = ontarget)
pval <- summary(m)$coefficients[2, 4]
r2   <- summary(m)$r.squared

ggplot(ontarget, aes(x = X16S.on.target.., y = ITS.on.target..)) +
	geom_point(size = 1, alpha = 0.15) +
	geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "gray30") +
	annotate(
		"text", x = 2, y = 98, hjust = 0, vjust = 1,
		label = sprintf("p = %.2g, R² = %.3f", pval, r2),
		size = 4) +
	coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of 16S vs. ITS",
		x = "% on target (16S)",
		y = "% on target (ITS)")

m <- lm(Symbiosis.on.target.. ~ X16S.on.target.., data = ontarget)
pval <- summary(m)$coefficients[2, 4]
r2   <- summary(m)$r.squared

ggplot(ontarget, aes(x = Symbiosis.on.target.., y = X16S.on.target..)) +
	geom_point(size = 1, alpha = 0.15) +
	geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "gray30") +
	annotate(
		"text", x = 2, y = 98, hjust = 0, vjust = 1,
		label = sprintf("p = %.2g, R² = %.3f", pval, r2),
		size = 4) +
	coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of 16S vs. symbiosis genes",
		x = "% on target (symbiosis genes)",
		y = "% on target (16S)")
	
m <- lm(ITS.on.target.. ~ Symbiosis.on.target.., data = ontarget)
pval <- summary(m)$coefficients[2, 4]
r2   <- summary(m)$r.squared

ggplot(ontarget, aes(x = Symbiosis.on.target.., y = ITS.on.target..)) +
	geom_point(size = 1, alpha = 0.15) +
	geom_smooth(method = "lm", se = TRUE, linewidth = 0.8, color = "gray30") +
	annotate(
		"text", x = 2, y = 98, hjust = 0, vjust = 1,
		label = sprintf("p = %.2g, R² = %.3f", pval, r2),
		size = 4) +
	coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
	theme_minimal() +
	labs(
		title = "Hybridization efficiency of ITS vs.\nsymbiosis genes",
		x = "% on target (symbiosis genes)",
		y = "% on target (ITS)")

