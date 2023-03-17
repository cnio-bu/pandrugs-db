rm(list = ls()) # R version 4.1.3
library(XML) # XML_3.99-0.9
library(dbparser) # dbparser_1.2.0
library(tidyverse) # tidyverse_1.3.1
setwd(".")
outdir <- "results/"

# --- Data ---
read_drugbank_xml_db("data/drugbank.xml")

# --- Code ---
# Outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Run parsers
drugBank <- run_all_parsers()

# Drug names
dnames <- drugBank$general_information %>% 
  rename(preferred.DB = name, DB.IDs = primary_key) %>% 
  group_by(DB.IDs) %>% 
  mutate(n.preferred = length(unique(preferred.DB))) %>% 
  select(preferred.DB, DB.IDs, n.preferred) %>% unique
table(dnames$n.preferred) # One term per drug

# Targets
targets <- drugBank$targets_polypeptides %>% select(gene_name, parent_id) %>%
  rename(id = parent_id) %>% unique %>% 
  merge(drugBank$targets, by = "id") %>% filter(organism == "Humans") %>%
  rename(targets = gene_name, DB.IDs = parent_key) %>% 
  select(targets, DB.IDs) %>% unique

# Drug-gene interactions
drug_gene <- merge(dnames, targets, by = "DB.IDs", all = TRUE) %>%
  select(-n.preferred) %>% na.omit() %>% filter(targets != "") %>% unique()

# Save
write.table(drug_gene, file = paste0(outdir, "DrugBank.tsv"), col.names = TRUE, 
            row.names = FALSE, sep = "\t")
