################### Simulation in RNA-seq ######################

####################### 1 ##############################################

####### Install polyester package ##########
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('polyester')

library(polyester)
library(Biostrings)

####### Upload CDS fasta ############
CDS_fasta = readDNAStringSet(filepath = 'FASTAPATH')
writeXStringSet(CDS_fasta, 'CDS.fasta')

####### Read table #######
tabla_epi = read.csv('TABLEPATH.CSV', sep = ',')

####### Generate sub-tables with the genes and associated est_counts. ########
conteos_epi = cbind.data.frame(tabla_epi$target_id, tabla_epi$est_counts)
colnames(conteos_epi) = c('gene_ID', 'est_counts')

#######  Add a column with the reading values I want for the simulation (I simply multiply the RPM values by a factor).######
conteos_epi$cant_lecturas = round(conteos_epi$est_counts, 0) 

####### Generate a matrix that will contain the desired reading values.
countmat_epis = matrix(conteos_epi$cant_lecturas, nrow=length(CDS_fasta), ncol=1)
write.csv(x = conteos_epi, 'NEWTABLEPATH')

####### Run simulation
simulate_experiment_countmat('CDS.fasta', readmat=countmat_epis, 
                             outdir='OUTPUTPATH.FASTA')


