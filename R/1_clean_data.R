#
# Clean ArHA dataset
#

# Latest ArHa dataset
dat <- readRDS("data/Project_ArHa_database_2026-08-17.rds")

# SDM estimates from Gonzalo
sdm_dat <- read.csv("data/sdm_estimates.csv")

  # Inclusion criteria:
# Pathogen species, assay type, host species, >1 tested, non-NA values, 

# Cleaning function: pathogen data
clean_pathogen_data <- function(data){
  
  cleaned_pathogen_data <- data %>% 
    filter(
      taxonomic_level == "species",  # pathogen species, not family-level
      assay != "Missing", assay != "Other",  # assay type defined
      number_tested > 0,
      !is.na(number_tested), !is.na(number_positive), !is.na(pathogen_species_cleaned), 
      !is.na(host_record_id), !is.na(pathogen_record_id), !is.na(study_id))
  
  return(cleaned_pathogen_data)
}

# Cleaning function: host data
clean_host_data <- function(data){
  
  cleaned_host_data <- data %>% 
    filter(
      taxon_rank == "species",  # host species, not higher level
      coord_status == "valid",  # only valid spatial coordinates
      coordinate_resolution_processed %in% c("site", "village", "town", "city", "adm3"), # coordinate resolution at some meaningful spatial scale
      temporal_resolution %in% c("full_date", "day_range_resolution", "month_year"),  # temporal resolution at some meaningful scale to capture seasonal effects
      !is.na(host_species), 
      !is.na(host_record_id), !is.na(study_id))
  
  return(cleaned_host_data)
}

# Clean data by applying cleaning functions
path_data <- clean_pathogen_data(data = dat$pathogen)
host_data <- clean_host_data(data = dat$host)

# Combine
names(host_data)
names(path_data)
intersect(names(host_data), names(path_data))
# full_data <- full_join(path_data, host_data)
full_data <- inner_join(path_data, host_data)

View(full_data)

full_data$host_species %>% n_distinct()


# Which host_record_id/study_id pairs in path_data have no match in host_data?
missing_hosts <- anti_join(path_data, host_data, by = c("host_record_id", "study_id"))
missing_hosts %>% distinct(host_record_id, study_id) %>% nrow()

# Check them against the RAW host table (before clean_host_data filtering)
dat$host %>%
  filter(host_record_id %in% missing_hosts$host_record_id) %>%
  distinct(host_record_id, taxon_rank, coord_status, coordinate_resolution_processed, temporal_resolution) %>% 
  View()
