library(devtools)
install_github("PRIDE-R/prideR")
library(prideR)

archive.df <- as.data.frame(search.list.ProjectSummary("",0,20000))
archive.df <- subset(archive.df, project.tags!="Not available")

species.counts <- table(archive.df$species)
low.freq.species <- as.factor(names(which(species.counts<3)))
archive.df[archive.df$species %in% low.freq.species,]$species <- "Other"
archive.df$species <- as.factor(archive.df$species)

tissues.counts <- table(archive.df$tissues)
low.freq.tissues <- as.factor(names(which(tissues.counts<3)))
archive.df[archive.df$tissues %in% low.freq.tissues,]$tissues <- "Other"
archive.df$tissues <- as.factor(archive.df$tissues)

ptm.counts <- table(archive.df$ptm.names)
low.freq.ptm <- as.factor(names(which(ptm.counts<3)))
archive.df[archive.df$ptm.names %in% low.freq.ptm,]$ptm.names <- "Other"
archive.df$ptm.names <- as.factor(archive.df$ptm.names)

instrument.counts <- table(archive.df$instrument.names)
low.freq.instrument <- as.factor(names(which(instrument.counts<3)))
archive.df[archive.df$instrument.names %in% low.freq.instrument,]$instrument.names <- "Other"
archive.df$instrument.names <- as.factor(archive.df$instrument.names)

archive.df$project.tags <- as.factor(archive.df$project.tags)
archive.df$submissionType <- as.factor(archive.df$submissionType)
