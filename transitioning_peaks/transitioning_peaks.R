# Load packages
library(stringr)
library(tidyverse)
library(gplots)

# Locate files - patient csvs with annotated peaks
infiles <- list.files('Downloads/patient_csvs', pattern = '*.csv', full.names = T)

# Create vector of sample names
samples <- gsub('Downloads/patient_csvs/', replacement = '', infiles)
samples <- substr(samples, 1, 3)

# Create vector of states
states <- str_sub(infiles, -6, -5)

# Create vector of sample and state 
sample_states <- paste0(samples, '_' , states)

# Create vector of unique sample names for later
samples_unique <- unique(samples)

# Load data
dat <- lapply(infiles, function(x){
  read.csv(x)
})

# Set names of list
names(dat) <- sample_states

# Save dat as RDS 
saveRDS(dat, file = 'dat.RDS')
#dat <- readRDS('dat.RDS')

# Create a peak name in each dataframe and subset down to select columns
for (x in names(dat)) {
  dat[[x]] <- dat[[x]] %>%
    mutate(peak = paste(seqnames, start, end, sep = '-'),
           celltype = substr(x, 1, 1),
           patientID = substr(x, 2, 3),
           peak_patient =  paste(peak, patientID, sep = '-')) %>%
    select(peak, celltype, patientID, V4, peak_patient)
}

# ARE WE INCLUDING NON_TRANSITIONING PEAKS?? - YES 
# Remove NON-transitioning peaks
# Merge together GIC and NSC peaks into 2 different dataframes
g_dat <- dat[grep('G', names(dat))] %>%
  reduce(full_join) 
colnames(g_dat)
colnames(g_dat) <- c('peak', 'celltype', 'patientID', 'G_State', 'peak_patient')
g_dat <- g_dat %>%
  select(peak_patient, G_State)

n_dat <- dat[grep('N', names(dat))] %>%
  reduce(full_join)
colnames(n_dat)
colnames(n_dat) <- c('peak', 'celltype', 'patientID', 'N_State', 'peak_patient')
n_dat <- n_dat %>%
  select(peak_patient, N_State) 

# Then merge g_dat and n_dat for each peak have a column for GIC and NSC state
g_n_dat <- full_join(g_dat, n_dat, by = 'peak_patient')

# Save g_dat, n_dat and g_n_dat as RDS 
saveRDS(g_dat, file = 'g_dat.RDS')
saveRDS(n_dat, file = 'n_dat.RDS')
saveRDS(g_n_dat, file = 'g_n_dat.RDS')
#readRDS('g_dat.RDS')
#readRDS('n_dat.RDS')
#readRDS('g_n_dat.RDS')

# Find all rows where GIC state and NSC state are the same i.e. non-transitioning peaks
# This will be used later to remove non-transitioning peaks
non_trans <- g_n_dat %>%
  filter(N_State == G_State)

# Merge together all dataframes - 1 big dataframe containing all peaks
dat_join <- dat %>%
  reduce(full_join)
colnames(dat_join)
colnames(dat_join) <- c('peak', 'celltype', 'patientID', 'State', 'peak_patient')

# Save dat_join as RDS 
saveRDS(dat_join, file = 'dat_join.RDS')
#readRDS('dat_join.RDS')

# Remove non-transitioning peaks from dat_join
dat_trans <- dat_join %>%
  filter(!peak_patient %in% non_trans$peak_patient)

# Next we want peaks common between GIC and NSC in each patientID i.e. not just unique to GIC or NSC
# Find number of times each peak is seen in each patientID group 
# (should be 2 for peaks in GIC and NSC)
# Create a new column to add the patient ID to the peak name - we will use this to filter
patientID_occurences <- dat_trans %>%
  count(peak, patientID, sort = TRUE) %>%
  filter(n == 2)

patientID_occurences <- patientID_occurences %>%
  mutate(peak_patient = paste(peak, patientID, sep = '-'))

# Subset down to the above peaks
dat_join_filt <- dat_trans %>%
  filter(peak_patient %in% patientID_occurences$peak_patient) 

# Then count occurences of peak
# We want peaks that occur in at least 2 patients
# They should occur in GIC and NSC in each patient
# Therefore they should occur at least 4 times 
# n occurences should all be even (use to check)
peaks_2more_patients <- dat_join_filt %>%
  count(peak) %>%
  filter(n >= 4)

# Now filter dat_join_filt further to the above peaks
dat_join_filt <- dat_join_filt %>%
  filter(peak %in% peaks_2more_patients$peak)

# See how many peaks change from 'state' in NSC to 'state' in GIC
nsc_state <- 'E2'
gic_state <- 'E3'
test <- dat_join_filt %>%
  filter((celltype == 'N' & State == nsc_state)|
           (celltype == 'G' & State == gic_state)) %>%
  group_by(peak_patient) %>%
  count(patientID) %>%
  filter(n == 2)

# Print result
print(
  paste0(
    'Number of peaks transitioning between NSC-',
    nsc_state,
    ' and GIC-',
    gic_state,
    ' = ',
    length(unique(test$peak_patient))
  )
)

# Now run a nested for loop to get number of transitioning states for all comparisons 

# Make a matrix for results to go in
n_trans_peaks <- matrix(nrow = length(unique(states)),
                        ncol = length(unique(states)))
rownames(n_trans_peaks) <- paste0('GIC-', unique(states))
colnames(n_trans_peaks) <- paste0('NSC-', unique(states))

# Loop through all potential changes from NSC > GIC
for (nsc_state in unique(states)) {
  for (gic_state in unique(states)) {
    # See how many peaks change from 'state' in NSC to 'state' in GIC
    test <- dat_join_filt %>%
      filter((celltype == 'N' & State == nsc_state)|
               (celltype == 'G' & State == gic_state)) %>%
      group_by(peak_patient) %>%
      count(patientID) %>%
      filter(n == 2)
    # Pull out integer of peaks
    npeaks <- length(unique(test$peak_patient))
    # Add the result matrix
    n_trans_peaks[paste0('GIC-', gic_state), paste0('NSC-', nsc_state)] <- npeaks
  }
}

# Plot heatmap
heatmap.2(
  n_trans_peaks,
  Rowv = F,
  Colv = F,
  trace = 'none',
  col = bluered,
  density.info = 'none',
  cellnote = n_trans_peaks,
  notecol = 'black'
)

