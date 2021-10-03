library(tidyverse)
library(shiny)
library(RColorBrewer)

fun.read.count <- function(file.bion = ""){
	if(nchar(file.bion) == 0){bion.orig <- read.table("BION.taxmap.species.tsv", sep = "\t", header=T, comment.char="")} else {bion.orig <- read.table(paste0(file.bion), sep = "\t", header=T, comment.char="")}
	bion.orig.del <- bion.orig[-c(1, 2, 3), ]		
	bion.orig.del[,1] <-sub("\\#", "",bion.orig.del[,1]) #eemaldab reast „#“
	bion.orig.del[,1] = as.integer(bion.orig.del[,1])		
	bion.orig.del = bion.orig.del[,!grepl("Row|Sim|Fav|Mark",colnames(bion.orig.del))] #eemaldab üleliigsed veerud MERGE TABELID
	names(bion.orig.del)[1] <-sub("X..", "", names(bion.orig.del)[1]) #eemaldab reast „X..“
	names(bion.orig.del)[dim(bion.orig.del)[2]] <- paste("taxa") #paneb viimasele veerule nimeks „taxa“
	input.reads = colSums(bion.orig.del[,-ncol(bion.orig.del)])
	return(input.reads)
}

get.data <- function(level, longnames = F, cutoff = "", file.bion = ""){
	options(warn=-1)
	library(tidyverse)
	library(vegan)
	#library(reshape2)
	#library(data.table)
	if("mean" %in% ls()){rm(mean)}


	if(nchar(file.bion) == 0){bion.orig <- read.table("BION.taxmap.species.tsv", sep = "\t", header=T, comment.char="")} else {bion.orig <- read.table(file.bion, sep = "\t", header=T, comment.char="")}
	bion.orig.del <- bion.orig[-c(1, 2, 3), ]
	bion.orig.del[,1] <-sub("\\#", "",bion.orig.del[,1]) #eemaldab reast „#“
	bion.orig.del[,1] = as.integer(bion.orig.del[,1])
	bion.orig.del = bion.orig.del[,!grepl("Row|Sim|Fav|Mark",colnames(bion.orig.del))] #eemaldab üleliigsed veerud MERGE TABELID
	
	names(bion.orig.del)[1] <-sub("X..", "", names(bion.orig.del)[1]) #eemaldab reast „X..“
	colnames(bion.orig.del) <-gsub("^X", "", colnames(bion.orig.del)) #eemaldab reast „X“, lõpus on uuesti vaja see rida, sest X-d tulevad mingis andmeraamiks tegemise etapis tagasi, kui veerunimed algavad numbritega

	names(bion.orig.del)[dim(bion.orig.del)[2]] <-paste("taxa") #paneb viimasele veerule nimeks „taxa“

	bion.orig.del[,-ncol(bion.orig.del)] = t(decostand(t(bion.orig.del[,-ncol(bion.orig.del)]), "total")) 
	if(nchar(cutoff) > 0){bion.orig.del[,-ncol(bion.orig.del)][bion.orig.del[,-ncol(bion.orig.del)] < cutoff] <- 0}  #kui soovin alla mingi piiri väärtused 0-ga asendada
	bion.orig.del = bion.orig.del[rowSums(bion.orig.del[-ncol(bion.orig.del)]) > 0,] # eemaldame nullsummaga veerud

	data.st = bion.orig.del[-ncol(bion.orig.del)] #eemaldan viimase veeru
	
	taxa = as.data.frame(bion.orig.del$taxa) #tõstame taksonoomia eraldi
	colnames(taxa) = "taxa"

	taxa.s = separate(taxa, taxa,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ";[ ]")
	#taxa.s = separate(taxa, taxa,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), "\\;[:space:]")
	taxa.s[is.na(taxa.s)] <- "unclassified"
	taxa.s =data.frame(taxa.s, taxa)

	genus = data.frame(taxa.s[,6],data.st)
	taxa.s[grep("/|^unclassified$",taxa.s[,2]),2] = "p__unclassified" #kaotan ära “/”ga perekonnanimed
	taxa.s[grep("/|^unclassified$",taxa.s[,3]),3] = "c__unclassified" #kaotan ära “/”ga liiginimed
	taxa.s[grep("/|^unclassified$",taxa.s[,4]),4] = "o__unclassified" #kaotan ära “/”ga perekonnanimed
	taxa.s[grep("/|^unclassified$",taxa.s[,5]),5] = "f__unclassified" #kaotan ära “/”ga liiginimed
	taxa.s[grep("/|^unclassified$",taxa.s[,6]),6] = "g__unclassified" #kaotan ära “/”ga perekonnanimed
	taxa.s[grep("/|^unclassified$",taxa.s[,7]),7] = "s__unclassified" #kaotan ära “/”ga liiginimed
	group.genus_species = data.frame(taxa.s[,6], taxa.s[,7], koos=paste0(taxa.s[,6]," ",taxa.s[,7]))

	if(level == "kingdom"){data = data.frame(taxa.s[,1],data.st)}; if(level == "kingdom" & longnames == T){data = data.frame(tax = gsub("; p__.*", "", taxa.s[,"taxa"]),data.st); data = data[grepl("d__", data$tax),]}		
	if(level == "phylum"){data = data.frame(taxa.s[,2],data.st)}; if(level == "phylum" & longnames == T){data = data.frame(tax = gsub("; c__.*", "", taxa.s[,"taxa"]),data.st); data = data[grepl("p__", data$tax),]}	
	if(level == "class"){data = data.frame(taxa.s[,3],data.st)}; if(level == "class" & longnames == T){data = data.frame(tax = gsub("; o__.*", "", taxa.s[,"taxa"]),data.st); data = data[grepl("c__", data$tax),]}	
	if(level == "order"){data = data.frame(taxa.s[,4],data.st)}; if(level == "order" & longnames == T){data = data.frame(tax = gsub("; f__.*", "", taxa.s[,"taxa"]),data.st); data = data[grepl("o__", data$tax),]}	
	if(level == "family"){data = data.frame(taxa.s[,5],data.st)}; if(level == "family" & longnames == T){data = data.frame(tax = gsub("; g__.*", "", taxa.s[,"taxa"]),data.st); data = data[grepl("f__", data$tax),]}	
	if(level == "genus"){data = data.frame(taxa.s[,6],data.st)}; if(level == "genus" & longnames == T){data = data.frame(tax = gsub("; s__.*", "", taxa.s[,"taxa"]),data.st); data = data[grepl("g__", data$tax),]}	
	if(level == "species"){data = data.frame(group.genus_species[,3], data.st)}; if(level == "species" & longnames == T){data = data.frame(tax = taxa.s[,"taxa"],data.st); data = data[grepl("s__", data$tax),]}	

	colnames(data)[1]= c("tax")
	data = group_by(data, tax) %>%  summarise_all(funs(sum)) %>% as.data.frame()
	colnames(data) <-gsub("^X", "", colnames(data)) #eemaldab veerunimede algusest „X“ 
	options(warn=0)
	return(data)
}

# species = get.data(level = "species", file.bion = "raw_data/BION.taxmap.species.tsv")
# species.long <- species %>% gather(-tax, key = "sample", value = "proportion") #alternatiiv: species.long = melt(species) #pakett - reshape
# species.long$sample = factor(species.long$sample, levels = unique(species.long$sample))
# un.genera = unique(gsub("g__| s__.*$|_..*","",species.long$tax))
# un.species = unique(gsub("g__|s__","",species.long$tax))


