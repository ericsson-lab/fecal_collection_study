{
  library(tidyverse)
  library(janitor)
  library(readxl)
  
  select <- dplyr::select
  std_error <- function(x) sd(x)/sqrt(length(x))
  
}

# Load metadata
metadata <- read_excel('metadata.xlsx') %>% 
  clean_names() %>% 
  filter(!cage %in% c(143, 172143)) # same cage.

# Load rarefied table
orig_table <- read_tsv('data/core-metrics-results_10734/rarefied-feature-table_taxa.tsv', 
                       skip = 1) %>% 
  rename(featureid = `#OTU ID`)

# Pull all samples, remove taxa column
table <- orig_table %>% 
  select(featureid, any_of(metadata$sampleid))

# Load and pull taxa. 
# fill down taxa to prevent NA
taxonomy <- orig_table %>% 
  select(featureid, Taxon) %>% 
  separate(Taxon, 
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), 
           sep = '; ') %>% 
  pivot_longer(-featureid, names_to = 'level', values_to = 'taxon') %>% 
  group_by(featureid) %>% 
  fill(taxon) %>% 
  mutate(taxon = str_sub(taxon, 4, nchar(taxon)))

