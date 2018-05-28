args <- commandArgs(trailingOnly = TRUE)
options(scipen=999)

inFile = args[1]
outDir = args[2]
logfile = args[3]
min_freq = as.numeric(args[4])
min_cells = as.numeric(args[5])
mergeOn = args[6]

cat("<html><table><tr><td>Starting analysis</td></tr>", file=logfile, append=F)

library(ggplot2)
library(reshape2)
library(data.table)
library(grid)
library(parallel)
#require(xtable)
cat("<tr><td>Reading input</td></tr>", file=logfile, append=T)
dat = read.table(inFile, header=T, sep="\t", dec=".", fill=T, stringsAsFactors=F)

needed_cols = c("Patient",  "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")
if(!all(needed_cols %in% names(dat))){
	cat("Missing column(s):<br />", file=logfile, append=F)
	missing_cols = needed_cols[!(needed_cols %in% names(dat))]
	for(col in missing_cols){
		cat(paste(col, "<br />"), file=logfile, append=T)
	}
	stop("Not all columns are present")
}

if(!("Total_Read_Count" %in% names(dat))){
	dat$Total_Read_Count = 0
}

dat = dat[,c("Patient",  "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")] 
dat = dat[!(nchar(dat$Clone_Sequence) < 2),]

dat$dsPerM = 0
dat = dat[!is.na(dat$Patient),]
dat$Related_to_leukemia_clone = F

setwd(outDir)
cat("<tr><td>Selecting first V/J Genes</td></tr>", file=logfile, append=T)
dat$V_Segment_Major_Gene = as.factor(as.character(lapply(strsplit(as.character(dat$V_Segment_Major_Gene), "; "), "[[", 1)))
dat$J_Segment_Major_Gene = as.factor(as.character(lapply(strsplit(as.character(dat$J_Segment_Major_Gene), "; "), "[[", 1)))

cat("<tr><td>Calculating Frequency</td></tr>", file=logfile, append=T)

dat$Frequency = ((10^dat$Log10_Frequency)*100)
dat = dat[dat$Frequency <= 100,] #just in case?

dat = dat[dat$Frequency >= min_freq,]

patient.sample.counts = data.frame(data.table(dat)[, list(count=.N), by=c("Patient", "Sample")])
patient.sample.counts = data.frame(data.table(patient.sample.counts)[, list(count=.N), by=c("Patient")])

print("Found the following patients with number of samples:")
print(patient.sample.counts)

patient.sample.counts.pairs = patient.sample.counts[patient.sample.counts$count %in% 1:2,]
patient.sample.counts.triplets = patient.sample.counts[patient.sample.counts$count == 3,]



triplets = dat[dat$Patient %in% patient.sample.counts.triplets$Patient,]
dat = dat[dat$Patient %in% patient.sample.counts.pairs$Patient,]

cat("<tr><td>Normalizing to lowest cell count within locus</td></tr>", file=logfile, append=T)

dat$locus_V = substring(dat$V_Segment_Major_Gene, 0, 4)
dat$locus_J = substring(dat$J_Segment_Major_Gene, 0, 4)
min_cell_count = data.frame(data.table(dat)[, list(min_cell_count=min(.SD$Cell_Count)), by=c("Patient", "locus_V", "locus_J")])

dat$min_cell_paste = paste(dat$Patient, dat$locus_V, dat$locus_J)
min_cell_count$min_cell_paste = paste(min_cell_count$Patient, min_cell_count$locus_V, min_cell_count$locus_J)

min_cell_count = min_cell_count[,c("min_cell_paste", "min_cell_count")]
print(paste("rows:", nrow(dat)))
dat = merge(dat, min_cell_count, by="min_cell_paste")
print(paste("rows:", nrow(dat)))
dat$normalized_read_count = round(dat$Clone_Molecule_Count_From_Spikes / dat$Cell_Count * dat$min_cell_count / 2, digits=2) #??????????????????????????????????? wel of geen / 2

dat = dat[dat$normalized_read_count >= min_cells,]

dat$paste = paste(dat$Sample, dat$Clone_Sequence)

patients = split(dat, dat$Patient, drop=T)
intervalReads = rev(c(0,10,25,50,100,250,500,750,1000,10000))
intervalFreq = rev(c(0,0.01,0.05,0.1,0.5,1,5))
V_Segments = c(".*", "IGHV", "IGHD", "IGKV", "IGKV", "IgKINTR", "TRGV", "TRDV", "TRDD" , "TRBV")
J_Segments = c(".*", ".*", ".*", "IGKJ", "KDE", ".*", ".*", ".*", ".*", ".*")
Titles = c("Total", "IGH-Vh-Jh", "IGH-Dh-Jh", "Vk-Jk", "Vk-Kde" , "Intron-Kde", "TCRG", "TCRD-Vd-Dd", "TCRD-Dd-Dd", "TCRB-Vb-Jb")
Titles = factor(Titles, levels=Titles)
TitlesOrder = data.frame("Title"=Titles, "TitlesOrder"=1:length(Titles))

single_patients = dat[NULL,]

patient.merge.list = list() #cache the 'both' table, 2x speedup for more memory...
patient.merge.list.second = list()
  scatter_locus_data_list = list()
cat(paste("<table border='0' style='font-family:courier;'>", sep=""), file="multiple_matches.html", append=T)
cat(paste("<table border='0' style='font-family:courier;'>", sep=""), file="single_matches.html", append=T)
patientCountOnColumn <- function(x, product, interval, on, appendtxt=F){
  if (!is.data.frame(x) & is.list(x)){
    x = x[[1]]
  }
  #x$Sample = factor(x$Sample, levels=unique(x$Sample))
  x = data.frame(x,stringsAsFactors=F)
  onShort = "reads"
  if(on == "Frequency"){
    onShort = "freq"
  }
  onx = paste(on, ".x", sep="")
  ony = paste(on, ".y", sep="")
  splt = split(x, x$Sample, drop=T)
  type="pair"
  if(length(splt) == 1){
    print(paste(paste(x[1,which(colnames(x) == "Patient")]), "has one sample"))
    splt[[2]] = splt[[1]][NULL,]
    type="single"
  }
  patient1 = splt[[1]]
  patient2 = splt[[2]]
  
  oneSample = patient1[1,"Sample"]
  twoSample = patient2[1,"Sample"]
  patient = x[1,"Patient"]
  
  switched = F
  if(length(grep(".*_Right$", twoSample)) == 1 || length(grep(".*_Dx_BM$", twoSample)) == 1 || length(grep(".*_Dx$", twoSample)) == 1 ){
    tmp = twoSample
    twoSample = oneSample
    oneSample = tmp
    tmp = patient1
    patient1 = patient2
    patient2 = tmp
    switched = T
  }
  if(appendtxt){
    cat(paste(patient, oneSample, twoSample, type, sep="\t"), file="patients.txt", append=T, sep="", fill=3)
  }
  cat(paste("<tr><td>", patient, "</td>", sep=""), file=logfile, append=T)
  
  if(mergeOn == "Clone_Sequence"){
    patient1$merge = paste(patient1$Clone_Sequence)
    patient2$merge = paste(patient2$Clone_Sequence)
  } else {
    patient1$merge = paste(patient1$V_Segment_Major_Gene, patient1$J_Segment_Major_Gene, patient1$CDR3_Sense_Sequence)
    patient2$merge = paste(patient2$V_Segment_Major_Gene, patient2$J_Segment_Major_Gene, patient2$CDR3_Sense_Sequence)
  }
  
  scatterplot_data_columns = c("Patient", "Sample", "Frequency", "normalized_read_count", "V_Segment_Major_Gene", "J_Segment_Major_Gene", "merge")
  scatterplot_data = patient1[NULL,scatterplot_data_columns]
  scatterplot.data.type.factor = c(oneSample, twoSample, paste(c(oneSample, twoSample), "In Both"))
  scatterplot_data$type = character(0)
  scatterplot_data$link = numeric(0)
  scatterplot_data$on = character(0)
  
  patientMerge = merge(patient1, patient2, by.x="merge", by.y="merge")[NULL,] #blegh
  
  cs.exact.matches = patient1[patient1$Clone_Sequence %in% patient2$Clone_Sequence,]$Clone_Sequence
  
  start.time = proc.time()
  merge.list = c()
  
  if(patient %in% names(patient.merge.list)){
    patientMerge = patient.merge.list[[patient]]
    merge.list[["second"]] = patient.merge.list.second[[patient]]
    scatterplot_data = scatter_locus_data_list[[patient]]
    cat(paste("<td>", nrow(patient1), " in ", oneSample, " and ", nrow(patient2), " in ", twoSample, ", ", nrow(patientMerge), " in both (fetched from cache)</td></tr>", sep=""), file=logfile, append=T)
    
    #print(names(patient.merge.list))
  } else {
    #fuzzy matching here...
    
    patient1.fuzzy = patient1
    patient2.fuzzy = patient2
    
    patient1.fuzzy$merge = paste(patient1.fuzzy$locus_V, patient1.fuzzy$locus_J)
    patient2.fuzzy$merge = paste(patient2.fuzzy$locus_V, patient2.fuzzy$locus_J)
    
    patient.fuzzy = rbind(patient1.fuzzy, patient2.fuzzy)
    patient.fuzzy = patient.fuzzy[order(nchar(patient.fuzzy$Clone_Sequence)),]
    
    merge.list = list()
    
    merge.list[["second"]] = vector()
    
    link.count = 1
    
    while(nrow(patient.fuzzy) > 1){
      first.merge = patient.fuzzy[1,"merge"]
      first.clone.sequence = patient.fuzzy[1,"Clone_Sequence"]
      first.sample = patient.fuzzy[1,"Sample"]
      merge.filter = first.merge == patient.fuzzy$merge
      
      #length.filter = nchar(patient.fuzzy$Clone_Sequence) - nchar(first.clone.sequence) <= 9
      
      first.sample.filter = first.sample == patient.fuzzy$Sample
      second.sample.filter = first.sample != patient.fuzzy$Sample
      
      #first match same sample, sum to a single row, same for other sample
      #then merge rows like 'normal'
      
      sequence.filter = grepl(paste("^", first.clone.sequence, sep=""), patient.fuzzy$Clone_Sequence)
      
      
      
      #match.filter = merge.filter & grepl(first.clone.sequence, patient.fuzzy$Clone_Sequence) & length.filter & sample.filter
      first.match.filter = merge.filter & sequence.filter & first.sample.filter
      second.match.filter = merge.filter & sequence.filter & second.sample.filter
      
      first.rows = patient.fuzzy[first.match.filter,]
      second.rows = patient.fuzzy[second.match.filter,]
      
      first.rows.v = table(first.rows$V_Segment_Major_Gene)
      first.rows.v = names(first.rows.v[which.max(first.rows.v)])
      first.rows.j = table(first.rows$J_Segment_Major_Gene)
      first.rows.j = names(first.rows.j[which.max(first.rows.j)])
      
      first.sum = data.frame(merge = first.clone.sequence,
                             Patient = patient,
                             Receptor = first.rows[1,"Receptor"],
                             Sample = first.rows[1,"Sample"],
                             Cell_Count = first.rows[1,"Cell_Count"],
                             Clone_Molecule_Count_From_Spikes = sum(first.rows$Clone_Molecule_Count_From_Spikes),
                             Log10_Frequency = log10(sum(first.rows$Frequency)),
                             Total_Read_Count = sum(first.rows$Total_Read_Count),
                             dsPerM = sum(first.rows$dsPerM),
                             J_Segment_Major_Gene = first.rows.j,
                             V_Segment_Major_Gene = first.rows.v,
                             Clone_Sequence = first.clone.sequence,
                             CDR3_Sense_Sequence = first.rows[1,"CDR3_Sense_Sequence"],
                             Related_to_leukemia_clone = F,
                             Frequency = sum(first.rows$Frequency),
                             locus_V = first.rows[1,"locus_V"],
                             locus_J = first.rows[1,"locus_J"],
                             min_cell_count = first.rows[1,"min_cell_count"],
                             normalized_read_count = sum(first.rows$normalized_read_count),
                             paste = first.rows[1,"paste"],
                             min_cell_paste = first.rows[1,"min_cell_paste"])
      
      if(nrow(second.rows) > 0){
        second.rows.v = table(second.rows$V_Segment_Major_Gene)
        second.rows.v = names(second.rows.v[which.max(second.rows.v)])
        second.rows.j = table(second.rows$J_Segment_Major_Gene)
        second.rows.j = names(second.rows.j[which.max(second.rows.j)])
        
        second.sum = data.frame(merge = first.clone.sequence,
                                Patient = patient,
                                Receptor = second.rows[1,"Receptor"],
                                Sample = second.rows[1,"Sample"],
                                Cell_Count = second.rows[1,"Cell_Count"],
                                Clone_Molecule_Count_From_Spikes = sum(second.rows$Clone_Molecule_Count_From_Spikes),
                                Log10_Frequency = log10(sum(second.rows$Frequency)),
                                Total_Read_Count = sum(second.rows$Total_Read_Count),
                                dsPerM = sum(second.rows$dsPerM),
                                J_Segment_Major_Gene = second.rows.j,
                                V_Segment_Major_Gene = second.rows.v,
                                Clone_Sequence = first.clone.sequence,
                                CDR3_Sense_Sequence = second.rows[1,"CDR3_Sense_Sequence"],
                                Related_to_leukemia_clone = F,
                                Frequency = sum(second.rows$Frequency),
                                locus_V = second.rows[1,"locus_V"],
                                locus_J = second.rows[1,"locus_J"],
                                min_cell_count = second.rows[1,"min_cell_count"],
                                normalized_read_count = sum(second.rows$normalized_read_count),
                                paste = second.rows[1,"paste"],
                                min_cell_paste = second.rows[1,"min_cell_paste"])
        
        #print(names(patientMerge))
        #print(merge(first.sum, second.sum, by="merge"))		   
        patientMerge = rbind(patientMerge, merge(first.sum, second.sum, by="merge"))
        #print("test2")
        patient.fuzzy = patient.fuzzy[!(first.match.filter | second.match.filter),]
        
        hidden.clone.sequences = c(first.rows[-1,"Clone_Sequence"], second.rows[second.rows$Clone_Sequence != first.clone.sequence,"Clone_Sequence"])
        merge.list[["second"]] = append(merge.list[["second"]], hidden.clone.sequences)
        
        tmp.rows = rbind(first.rows, second.rows)
        #print("test3")
        tmp.rows = tmp.rows[order(nchar(tmp.rows$Clone_Sequence)),]
        
        
        #add to the scatterplot data
        scatterplot.row = first.sum[,scatterplot_data_columns]
        scatterplot.row$type = paste(first.sum[,"Sample"], "In Both")
        scatterplot.row$link = link.count
        scatterplot.row$on = onShort
        
        scatterplot_data = rbind(scatterplot_data, scatterplot.row)
        
        scatterplot.row = second.sum[,scatterplot_data_columns]
        scatterplot.row$type = paste(second.sum[,"Sample"], "In Both")
        scatterplot.row$link = link.count
        scatterplot.row$on = onShort
        
        scatterplot_data = rbind(scatterplot_data, scatterplot.row)    
        
        #write some information about the match to a log file
        if (nrow(first.rows) > 1 | nrow(second.rows) > 1) {
          cat(paste("<tr><td>", patient, " row ", 1:nrow(tmp.rows), "</td><td>", tmp.rows$Sample, ":</td><td>", tmp.rows$Clone_Sequence, "</td><td>", tmp.rows$normalized_read_count, "</td></tr>", sep=""), file="multiple_matches.html", append=T)
        } else {
          second.clone.sequence = second.rows[1,"Clone_Sequence"]
          if(nchar(first.clone.sequence) != nchar(second.clone.sequence)){
            cat(paste("<tr bgcolor='#DDD'><td>", patient, " row ", 1:nrow(tmp.rows), "</td><td>", tmp.rows$Sample, ":</td><td>", tmp.rows$Clone_Sequence, "</td><td>", tmp.rows$normalized_read_count, "</td></tr>", sep=""), file="single_matches.html", append=T)
          } else {
            #cat(paste("<tr><td>", patient, " row ", 1:nrow(tmp.rows), "</td><td>", tmp.rows$Sample, ":</td><td>", tmp.rows$Clone_Sequence, "</td><td>", tmp.rows$normalized_read_count, "</td></tr>", sep=""), file="single_matches.html", append=T)
          }
        }
        
      } else if(nrow(first.rows) > 1) {
        if(patient1[1,"Sample"] == first.sample){
          patient1 = patient1[!(patient1$Clone_Sequence %in% first.rows$Clone_Sequence),]
          patient1 = rbind(patient1, first.sum)
        } else {
          patient2 = patient2[!(patient2$Clone_Sequence %in% first.rows$Clone_Sequence),]
          patient2 = rbind(patient2, first.sum)
        }
        
        hidden.clone.sequences = c(first.rows[-1,"Clone_Sequence"])
        merge.list[["second"]] = append(merge.list[["second"]], hidden.clone.sequences)
        
        patient.fuzzy = patient.fuzzy[-first.match.filter,]
        
        #add to the scatterplot data
        scatterplot.row = first.sum[,scatterplot_data_columns]
        scatterplot.row$type = first.sum[,"Sample"]
        scatterplot.row$link = link.count
        scatterplot.row$on = onShort
        
        scatterplot_data = rbind(scatterplot_data, scatterplot.row)
        
        cat(paste("<tr bgcolor='#DDF'><td>", patient, " row ", 1:nrow(first.rows), "</td><td>", first.rows$Sample, ":</td><td>", first.rows$Clone_Sequence, "</td><td>", first.rows$normalized_read_count, "</td></tr>", sep=""), file="single_matches.html", append=T)
      } else {
        patient.fuzzy = patient.fuzzy[-1,]
        
        #add to the scatterplot data
        scatterplot.row = first.sum[,scatterplot_data_columns]
        scatterplot.row$type = first.sum[,"Sample"]
        scatterplot.row$link = link.count
        scatterplot.row$on = onShort
        
        scatterplot_data = rbind(scatterplot_data, scatterplot.row)
      }
      link.count = link.count + 1    
    }
    patient.merge.list[[patient]] <<- patientMerge
    patient.merge.list.second[[patient]] <<- merge.list[["second"]]
    
    sample.order = data.frame(type = c(oneSample, twoSample, paste(c(oneSample, twoSample), "In Both")),type.order = 1:4)
    scatterplot_data = merge(scatterplot_data, sample.order, by="type")
    
    scatter_locus_data_list[[patient]] <<- scatterplot_data
    cat(paste("<td>", nrow(patient1), " in ", oneSample, " and ", nrow(patient2), " in ", twoSample, ", ", nrow(patientMerge), " in both (finding both took ", (proc.time() - start.time)[[3]], "s)</td></tr>", sep=""), file=logfile, append=T)
  }
  
  patient1 = patient1[!(patient1$Clone_Sequence %in% patient.merge.list.second[[patient]]),]
  patient2 = patient2[!(patient2$Clone_Sequence %in% patient.merge.list.second[[patient]]),]
  
  
  patientMerge$thresholdValue = pmax(patientMerge[,onx], patientMerge[,ony])
  #patientMerge$thresholdValue = pmin(patientMerge[,onx], patientMerge[,ony])
  res1 = vector()
  res2 = vector()
  resBoth = vector()
  read1Count = vector()
  read2Count = vector()
  locussum1 = vector()
  locussum2 = vector()
  
  #for(iter in 1){
  for(iter in 1:length(product[,1])){
    threshhold = product[iter,"interval"]
    V_Segment = paste(".*", as.character(product[iter,"V_Segments"]), ".*", sep="")
    J_Segment = paste(".*", as.character(product[iter,"J_Segments"]), ".*", sep="")
    #both = (grepl(V_Segment, patientMerge$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge$J_Segment_Major_Gene.x) & patientMerge[,onx] > threshhold & patientMerge[,ony] > threshhold) #both higher than threshold
    both = (grepl(V_Segment, patientMerge$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge$J_Segment_Major_Gene.x) & patientMerge$thresholdValue > threshhold) #highest of both is higher than threshold
    one = (grepl(V_Segment, patient1$V_Segment_Major_Gene) & grepl(J_Segment, patient1$J_Segment_Major_Gene) & patient1[,on] > threshhold & !(patient1$merge %in% patientMerge[both,]$merge))
    two = (grepl(V_Segment, patient2$V_Segment_Major_Gene) & grepl(J_Segment, patient2$J_Segment_Major_Gene) & patient2[,on] > threshhold & !(patient2$merge %in% patientMerge[both,]$merge))
    read1Count = append(read1Count, sum(patient1[one,]$normalized_read_count))
    read2Count = append(read2Count, sum(patient2[two,]$normalized_read_count))
    res1 = append(res1, sum(one))
    res2 = append(res2, sum(two))
    resBoth = append(resBoth, sum(both))
    locussum1 = append(locussum1, sum(patient1[(grepl(V_Segment, patient1$V_Segment_Major_Gene) & grepl(J_Segment, patient1$J_Segment_Major_Gene)),]$normalized_read_count))
    locussum2 = append(locussum2, sum(patient2[(grepl(V_Segment, patient2$V_Segment_Major_Gene) & grepl(J_Segment, patient2$J_Segment_Major_Gene)),]$normalized_read_count))
    #threshhold = 0
    if(threshhold != 0 | T){
      if(sum(one) > 0){
        dfOne = patient1[one,c("V_Segment_Major_Gene", "J_Segment_Major_Gene", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")]
        colnames(dfOne) = c("Proximal segment", "Distal segment", "normalized_read_count", "Frequency", "Clone Sequence", "Related_to_leukemia_clone")
        filenameOne = paste(oneSample, "_", product[iter, "Titles"], "_", threshhold, sep="")
        write.table(dfOne, file=paste(filenameOne, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
      }
      if(sum(two) > 0){
        dfTwo = patient2[two,c("V_Segment_Major_Gene", "J_Segment_Major_Gene", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")]
        colnames(dfTwo) = c("Proximal segment", "Distal segment", "normalized_read_count", "Frequency", "Clone Sequence", "Related_to_leukemia_clone")
        filenameTwo = paste(twoSample, "_", product[iter, "Titles"], "_", threshhold, sep="")
        write.table(dfTwo, file=paste(filenameTwo, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
      }
    }
	scatterplot_locus_data = scatterplot_data[grepl(V_Segment, scatterplot_data$V_Segment_Major_Gene) & grepl(J_Segment, scatterplot_data$J_Segment_Major_Gene),]
	if(nrow(scatterplot_locus_data) > 0){
		scatterplot_locus_data$Rearrangement = product[iter, "Titles"]
	}
	p = NULL
    #print(paste("nrow scatterplot_locus_data", nrow(scatterplot_locus_data)))
    if(nrow(scatterplot_locus_data) != 0){
      if(on == "normalized_read_count"){
        write.table(scatterplot_locus_data, file=paste(oneSample, twoSample, product[iter, "Titles"], "scatterplot_locus_data.txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
        scales = 10^(0:6) #(0:ceiling(log10(max(scatterplot_locus_data$normalized_read_count))))
        p = ggplot(scatterplot_locus_data, aes(factor(reorder(type, type.order)), normalized_read_count, group=link)) + geom_line() + scale_y_log10(breaks=scales,labels=scales, limits=c(1,1e6)) + scale_x_discrete(breaks=levels(scatterplot_data$type), labels=levels(scatterplot_data$type), drop=FALSE)
      } else {
        p = ggplot(scatterplot_locus_data, aes(factor(reorder(type, type.order)), Frequency, group=link)) + geom_line() + scale_y_log10(limits=c(0.0001,100), breaks=c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), labels=c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100")) + scale_x_discrete(breaks=levels(scatterplot_data$type), labels=levels(scatterplot_data$type), drop=FALSE)
      }
      p = p + geom_point(aes(colour=type), position="dodge")
      p = p + xlab("In one or both samples") + ylab(onShort) + ggtitle(paste(patient1[1,"Patient"], patient1[1,"Sample"], patient2[1,"Sample"], onShort, product[iter, "Titles"]))
    } else {
      p = ggplot(NULL, aes(x=c("In one", "In Both"),y=0)) + geom_blank(NULL) + xlab("In one or both of the samples") + ylab(onShort) + ggtitle(paste(patient1[1,"Patient"], patient1[1,"Sample"], patient2[1,"Sample"], onShort, product[iter, "Titles"]))
    }
	file_name = paste(patient1[1,"Patient"], "_", patient1[1,"Sample"], "_", patient2[1,"Sample"], "_", onShort, "_", product[iter, "Titles"],"_scatter.png", sep="")
	print(paste("Writing figure:", file_name))
    png(file_name)
    print(p)
    dev.off()
    if(sum(both) > 0){
      dfBoth = patientMerge[both,c("V_Segment_Major_Gene.x", "J_Segment_Major_Gene.x", "normalized_read_count.x", "Frequency.x", "Related_to_leukemia_clone.x", "Clone_Sequence.x", "V_Segment_Major_Gene.y", "J_Segment_Major_Gene.y", "normalized_read_count.y", "Frequency.y", "Related_to_leukemia_clone.y")]
      colnames(dfBoth) = c(paste("Proximal segment", oneSample), paste("Distal segment", oneSample), paste("Normalized_Read_Count", oneSample), paste("Frequency", oneSample), paste("Related_to_leukemia_clone", oneSample),"Clone Sequence", paste("Proximal segment", twoSample), paste("Distal segment", twoSample), paste("Normalized_Read_Count", twoSample), paste("Frequency", twoSample), paste("Related_to_leukemia_clone", twoSample))
      filenameBoth = paste(oneSample, "_", twoSample, "_", product[iter, "Titles"], "_", threshhold, sep="")
      write.table(dfBoth, file=paste(filenameBoth, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
    } 
  }
  patientResult = data.frame("Locus"=product$Titles, "J_Segment"=product$J_Segments, "V_Segment"=product$V_Segments, "cut_off_value"=paste(">", product$interval, sep=""), "Both"=resBoth, "tmp1"=res1, "read_count1" = round(read1Count), "tmp2"=res2, "read_count2"= round(read2Count), "Sum"=res1 + res2 + resBoth, "percentage" = round((resBoth/(res1 + res2 + resBoth)) * 100, digits=2), "Locus_sum1"=locussum1, "Locus_sum2"=locussum2)
  if(sum(is.na(patientResult$percentage)) > 0){
    patientResult[is.na(patientResult$percentage),]$percentage = 0
  }
  colnames(patientResult)[6] = oneSample
  colnames(patientResult)[8] = twoSample
  colnamesBak = colnames(patientResult)
  colnames(patientResult) = c("Ig/TCR gene rearrangement type", "Distal Gene segment", "Proximal gene segment", "cut_off_value", paste("Number of sequences ", patient, "_Both", sep=""), paste("Number of sequences", oneSample, sep=""), paste("Normalized Read Count", oneSample), paste("Number of sequences", twoSample, sep=""), paste("Normalized Read Count", twoSample), paste("Sum number of sequences", patient), paste("Percentage of sequences ", patient, "_Both", sep=""), paste("Locus Sum", oneSample), paste("Locus Sum", twoSample))
  write.table(patientResult, file=paste(patient, "_", onShort, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
  colnames(patientResult) = colnamesBak
  
  patientResult$Locus = factor(patientResult$Locus, Titles)
  patientResult$cut_off_value = factor(patientResult$cut_off_value, paste(">", interval, sep=""))
  
  plt = ggplot(patientResult[,c("Locus", "cut_off_value", "Both")])
  plt = plt + geom_bar( aes( x=factor(cut_off_value), y=Both), stat='identity', position="dodge", fill="#79c36a")
  plt = plt + facet_grid(.~Locus) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plt = plt + geom_text(aes(ymax=max(Both), x=cut_off_value,y=Both,label=Both), angle=90, hjust=0)
  plt = plt + xlab("Reads per locus") + ylab("Count") + ggtitle("Number of clones in both")
  plt = plt + theme(plot.margin = unit(c(1,8.8,0.5,1.5), "lines"))
  file_name = paste(patient, "_", onShort, ".png", sep="")
  print(paste("Writing figure:", file_name))
  png(file_name, width=1920, height=1080)
  print(plt)
  dev.off()
  #(t,r,b,l)
  plt = ggplot(patientResult[,c("Locus", "cut_off_value", "percentage")])
  plt = plt + geom_bar( aes( x=factor(cut_off_value), y=percentage), stat='identity', position="dodge", fill="#79c36a")
  plt = plt + facet_grid(.~Locus) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plt = plt + geom_text(aes(ymax=max(percentage), x=cut_off_value,y=percentage,label=percentage), angle=90, hjust=0)
  plt = plt + xlab("Reads per locus") + ylab("Count") + ggtitle("% clones in both left and right")
  plt = plt + theme(plot.margin = unit(c(1,8.8,0.5,1.5), "lines"))
  
  file_name = paste(patient, "_percent_", onShort, ".png", sep="")
  print(paste("Writing figure:", file_name))
  png(file_name, width=1920, height=1080)
  print(plt)
  dev.off()
  
  patientResult = melt(patientResult[,c('Locus','cut_off_value', oneSample, twoSample)] ,id.vars=1:2)
  patientResult$relativeValue = patientResult$value * 10
  patientResult[patientResult$relativeValue == 0,]$relativeValue = 1
  plt = ggplot(patientResult)
  plt = plt + geom_bar( aes( x=factor(cut_off_value), y=relativeValue, fill=variable), stat='identity', position="dodge")
  plt = plt + facet_grid(.~Locus) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  plt = plt + scale_y_continuous(trans="log", breaks=10^c(0:10), labels=c(0, 10^c(0:9)))
  plt = plt + geom_text(data=patientResult[patientResult$variable == oneSample,], aes(ymax=max(value), x=cut_off_value,y=relativeValue,label=value), angle=90, position=position_dodge(width=0.9), hjust=0, vjust=-0.2)
  plt = plt + geom_text(data=patientResult[patientResult$variable == twoSample,], aes(ymax=max(value), x=cut_off_value,y=relativeValue,label=value), angle=90, position=position_dodge(width=0.9), hjust=0, vjust=0.8)
  plt = plt + xlab("Reads per locus") + ylab("Count") + ggtitle(paste("Number of clones in only ", oneSample, " and only ", twoSample, sep=""))
  file_name = paste(patient, "_", onShort, "_both.png", sep="")
  print(paste("Writing figure:", file_name))
  png(file_name, width=1920, height=1080)
  print(plt)
  dev.off()
}

if(length(patients) > 0){
	cat("<tr><td>Starting Frequency analysis</td></tr>", file=logfile, append=T)

	interval = intervalFreq
	intervalOrder = data.frame("interval"=paste(">", interval, sep=""), "intervalOrder"=1:length(interval))
	product = data.frame("Titles"=rep(Titles, each=length(interval)), "interval"=rep(interval, times=10), "V_Segments"=rep(V_Segments, each=length(interval)), "J_Segments"=rep(J_Segments, each=length(interval)))
	for (current_patient in patients){
		print(paste("Started working", unique(current_patient$Patient), "Frequency analysis"))
		patientCountOnColumn(current_patient, product=product, interval=interval, on="Frequency", appendtxt=T)
	}

	cat("<tr><td>Starting Cell Count analysis</td></tr>", file=logfile, append=T)

	interval = intervalReads
	intervalOrder = data.frame("interval"=paste(">", interval, sep=""), "intervalOrder"=1:length(interval))
	product = data.frame("Titles"=rep(Titles, each=length(interval)), "interval"=rep(interval, times=10), "V_Segments"=rep(V_Segments, each=length(interval)), "J_Segments"=rep(J_Segments, each=length(interval)))
	for (current_patient in patients){
		print(paste("Started working on ", unique(current_patient$Patient), "Read Count analysis"))
		patientCountOnColumn(current_patient, product=product, interval=interval, on="normalized_read_count")
	}
}
if(nrow(single_patients) > 0){
	scales = 10^(0:6) #(0:ceiling(log10(max(scatterplot_locus_data$normalized_read_count))))
	p = ggplot(single_patients, aes(Rearrangement, normalized_read_count)) + scale_y_log10(breaks=scales,labels=as.character(scales)) + expand_limits(y=c(0,1000000))
	p = p + geom_point(aes(colour=type), position="jitter")
	p = p + xlab("In one or both samples") + ylab("Reads")
	p = p + facet_grid(.~Patient) + ggtitle("Scatterplot of the reads of the patients with a single sample")
	file_name = "singles_reads_scatterplot.png"
  	print(paste("Writing figure:", file_name))
	png(file_name, width=640 * length(unique(single_patients$Patient)) + 100, height=1080)
	print(p)
	dev.off()

	#p = ggplot(single_patients, aes(Rearrangement, Frequency)) + scale_y_continuous(limits = c(0, 100)) + expand_limits(y=c(0,100))
	p = ggplot(single_patients, aes(Rearrangement, Frequency)) + scale_y_log10(limits=c(0.0001,100), breaks=c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), labels=c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100")) + expand_limits(y=c(0,100))
	p = p + geom_point(aes(colour=type), position="jitter")
	p = p + xlab("In one or both samples") + ylab("Frequency")
	p = p + facet_grid(.~Patient) + ggtitle("Scatterplot of the frequency of the patients with a single sample")
	file_name = "singles_freq_scatterplot.png"
  	print(paste("Writing figure:", file_name))
	png(file_name, width=640 * length(unique(single_patients$Patient)) + 100, height=1080)
	print(p)
	dev.off()
} else {
	empty <- data.frame()
	p = ggplot(empty) + geom_point() + xlim(0, 10) + ylim(0, 100) + xlab("In one or both samples") + ylab("Frequency") + ggtitle("Scatterplot of the frequency of the patients with a single sample")
	
	png("singles_reads_scatterplot.png", width=400, height=300)
	print(p)
	dev.off()	
	
	png("singles_freq_scatterplot.png", width=400, height=300)
	print(p)
	dev.off()
}

patient.merge.list = list() #cache the 'both' table, 2x speedup for more memory...
patient.merge.list.second = list()

tripletAnalysis <- function(patient1, label1, patient2, label2, patient3, label3, product, interval, on, appendTriplets= FALSE){
	onShort = "reads"
	if(on == "Frequency"){
	onShort = "freq"
	}
	onx = paste(on, ".x", sep="")
	ony = paste(on, ".y", sep="")
	onz = paste(on, ".z", sep="")
	type="triplet"

	threshholdIndex = which(colnames(product) == "interval")
	V_SegmentIndex = which(colnames(product) == "V_Segments")
	J_SegmentIndex = which(colnames(product) == "J_Segments")
	titleIndex = which(colnames(product) == "Titles")
	sampleIndex = which(colnames(patient1) == "Sample")
	patientIndex = which(colnames(patient1) == "Patient")
	oneSample = paste(patient1[1,sampleIndex], sep="")
	twoSample = paste(patient2[1,sampleIndex], sep="")
	threeSample = paste(patient3[1,sampleIndex], sep="")

	if(mergeOn == "Clone_Sequence"){
	patient1$merge = paste(patient1$Clone_Sequence)
		patient2$merge = paste(patient2$Clone_Sequence)
		patient3$merge = paste(patient3$Clone_Sequence)

	} else {
		patient1$merge = paste(patient1$V_Segment_Major_Gene, patient1$J_Segment_Major_Gene, patient1$CDR3_Sense_Sequence)
		patient2$merge = paste(patient2$V_Segment_Major_Gene, patient2$J_Segment_Major_Gene, patient2$CDR3_Sense_Sequence)
		patient3$merge = paste(patient3$V_Segment_Major_Gene, patient3$J_Segment_Major_Gene, patient3$CDR3_Sense_Sequence)
	}

	#patientMerge = merge(patient1, patient2, by="merge")[NULL,]
	patient1.fuzzy = patient1
	patient2.fuzzy = patient2
	patient3.fuzzy = patient3

	cat(paste("<tr><td>", label1, "</td>", sep=""), file=logfile, append=T)

	patient1.fuzzy$merge = paste(patient1.fuzzy$locus_V, patient1.fuzzy$locus_J)
	patient2.fuzzy$merge = paste(patient2.fuzzy$locus_V, patient2.fuzzy$locus_J)
	patient3.fuzzy$merge = paste(patient3.fuzzy$locus_V, patient3.fuzzy$locus_J)

	patient.fuzzy = rbind(patient1.fuzzy ,patient2.fuzzy, patient3.fuzzy)
	patient.fuzzy = patient.fuzzy[order(nchar(patient.fuzzy$Clone_Sequence)),]

	other.sample.list = list()
	other.sample.list[[oneSample]] = c(twoSample, threeSample)
	other.sample.list[[twoSample]] = c(oneSample, threeSample)
	other.sample.list[[threeSample]] = c(oneSample, twoSample)

	patientMerge = merge(patient1, patient2, by="merge")
	patientMerge = merge(patientMerge, patient3, by="merge")
	colnames(patientMerge)[which(!grepl("(\\.x$)|(\\.y$)|(merge)", names(patientMerge)))] = paste(colnames(patientMerge)[which(!grepl("(\\.x$)|(\\.y$)|(merge)", names(patientMerge), perl=T))], ".z", sep="")
	#patientMerge$thresholdValue = pmax(patientMerge[,onx], patientMerge[,ony], patientMerge[,onz])
	patientMerge = patientMerge[NULL,]

	duo.merge.list = list()

	patientMerge12 = merge(patient1, patient2, by="merge")
	#patientMerge12$thresholdValue = pmax(patientMerge12[,onx], patientMerge12[,ony])
	patientMerge12 = patientMerge12[NULL,]
	duo.merge.list[[paste(oneSample, twoSample)]] = patientMerge12
	duo.merge.list[[paste(twoSample, oneSample)]] = patientMerge12

	patientMerge13 = merge(patient1, patient3, by="merge")
	#patientMerge13$thresholdValue = pmax(patientMerge13[,onx], patientMerge13[,ony])
	patientMerge13 = patientMerge13[NULL,]
	duo.merge.list[[paste(oneSample, threeSample)]] = patientMerge13
	duo.merge.list[[paste(threeSample, oneSample)]] = patientMerge13

	patientMerge23 = merge(patient2, patient3, by="merge")
	#patientMerge23$thresholdValue = pmax(patientMerge23[,onx], patientMerge23[,ony])
	patientMerge23 = patientMerge23[NULL,]
	duo.merge.list[[paste(twoSample, threeSample)]] = patientMerge23
	duo.merge.list[[paste(threeSample, twoSample)]] = patientMerge23

	merge.list = list()
	merge.list[["second"]] = vector()
	
	#print(paste(nrow(patient1), nrow(patient2), nrow(patient3), label1, label2, label3))
	
	start.time = proc.time()
	if(paste(label1, "123") %in% names(patient.merge.list)){
		patientMerge = patient.merge.list[[paste(label1, "123")]]
		patientMerge12 = patient.merge.list[[paste(label1, "12")]]
		patientMerge13 = patient.merge.list[[paste(label1, "13")]]
		patientMerge23 = patient.merge.list[[paste(label1, "23")]]

		#merge.list[["second"]] = patient.merge.list.second[[label1]]

		cat(paste("<td>", nrow(patient1), " in ", label1, " and ", nrow(patient2), " in ", label2, nrow(patient3), " in ", label3, ", ", nrow(patientMerge), " in both (fetched from cache)</td></tr>", sep=""), file=logfile, append=T)
	} else {
		while(nrow(patient.fuzzy) > 0){
			first.merge = patient.fuzzy[1,"merge"]
			first.clone.sequence = patient.fuzzy[1,"Clone_Sequence"]
			first.sample = paste(patient.fuzzy[1,"Sample"], sep="")
			
			merge.filter = first.merge == patient.fuzzy$merge
			
			second.sample = other.sample.list[[first.sample]][1]
			third.sample = other.sample.list[[first.sample]][2]

			sample.filter.1 = first.sample == patient.fuzzy$Sample
			sample.filter.2 = second.sample == patient.fuzzy$Sample
			sample.filter.3 = third.sample == patient.fuzzy$Sample

			sequence.filter = grepl(paste("^", first.clone.sequence, sep=""), patient.fuzzy$Clone_Sequence)

			match.filter.1 = sample.filter.1 & sequence.filter & merge.filter
			match.filter.2 = sample.filter.2 & sequence.filter & merge.filter
			match.filter.3 = sample.filter.3 & sequence.filter & merge.filter

			matches.in.1 = any(match.filter.1)
			matches.in.2 = any(match.filter.2)
			matches.in.3 = any(match.filter.3)

			rows.1 = patient.fuzzy[match.filter.1,]

			sum.1 = data.frame(merge = first.clone.sequence,
								 Patient = label1,
								 Receptor = rows.1[1,"Receptor"],
								 Sample = rows.1[1,"Sample"],
								 Cell_Count = rows.1[1,"Cell_Count"],
								 Clone_Molecule_Count_From_Spikes = sum(rows.1$Clone_Molecule_Count_From_Spikes),
								 Log10_Frequency = log10(sum(rows.1$Frequency)),
								 Total_Read_Count = sum(rows.1$Total_Read_Count),
								 dsPerM = sum(rows.1$dsPerM),
								 J_Segment_Major_Gene = rows.1[1,"J_Segment_Major_Gene"],
								 V_Segment_Major_Gene = rows.1[1,"V_Segment_Major_Gene"],
								 Clone_Sequence = first.clone.sequence,
								 CDR3_Sense_Sequence = rows.1[1,"CDR3_Sense_Sequence"],
								 Related_to_leukemia_clone = F,
								 Frequency = sum(rows.1$Frequency),
								 locus_V = rows.1[1,"locus_V"],
								 locus_J = rows.1[1,"locus_J"],
								 uniqueID = rows.1[1,"uniqueID"],
								 normalized_read_count = sum(rows.1$normalized_read_count))
			sum.2 = sum.1[NULL,]
			rows.2 = patient.fuzzy[match.filter.2,]
			if(matches.in.2){
				sum.2 = data.frame(merge = first.clone.sequence,
								   Patient = label1,
								   Receptor = rows.2[1,"Receptor"],
								   Sample = rows.2[1,"Sample"],
								   Cell_Count = rows.2[1,"Cell_Count"],
								   Clone_Molecule_Count_From_Spikes = sum(rows.2$Clone_Molecule_Count_From_Spikes),
								   Log10_Frequency = log10(sum(rows.2$Frequency)),
								   Total_Read_Count = sum(rows.2$Total_Read_Count),
								   dsPerM = sum(rows.2$dsPerM),
								   J_Segment_Major_Gene = rows.2[1,"J_Segment_Major_Gene"],
								   V_Segment_Major_Gene = rows.2[1,"V_Segment_Major_Gene"],
								   Clone_Sequence = first.clone.sequence,
								   CDR3_Sense_Sequence = rows.2[1,"CDR3_Sense_Sequence"],
								   Related_to_leukemia_clone = F,
								   Frequency = sum(rows.2$Frequency),
								   locus_V = rows.2[1,"locus_V"],
								   locus_J = rows.2[1,"locus_J"],
								   uniqueID = rows.2[1,"uniqueID"],
								   normalized_read_count = sum(rows.2$normalized_read_count))
			}

		sum.3 = sum.1[NULL,]
		rows.3 = patient.fuzzy[match.filter.3,]
		if(matches.in.3){
			sum.3 = data.frame(merge = first.clone.sequence,
							   Patient = label1,
							   Receptor = rows.3[1,"Receptor"],
							   Sample = rows.3[1,"Sample"],
							   Cell_Count = rows.3[1,"Cell_Count"],
							   Clone_Molecule_Count_From_Spikes = sum(rows.3$Clone_Molecule_Count_From_Spikes),
							   Log10_Frequency = log10(sum(rows.3$Frequency)),
							   Total_Read_Count = sum(rows.3$Total_Read_Count),
							   dsPerM = sum(rows.3$dsPerM),
							   J_Segment_Major_Gene = rows.3[1,"J_Segment_Major_Gene"],
							   V_Segment_Major_Gene = rows.3[1,"V_Segment_Major_Gene"],
							   Clone_Sequence = first.clone.sequence,
							   CDR3_Sense_Sequence = rows.3[1,"CDR3_Sense_Sequence"],
							   Related_to_leukemia_clone = F,
							   Frequency = sum(rows.3$Frequency),
							   locus_V = rows.3[1,"locus_V"],
							   locus_J = rows.3[1,"locus_J"],
							   uniqueID = rows.3[1,"uniqueID"],
							   normalized_read_count = sum(rows.3$normalized_read_count))
		}

	  if(matches.in.2 & matches.in.3){
			merge.123 = merge(sum.1, sum.2, by="merge")
			merge.123 = merge(merge.123, sum.3, by="merge")
			colnames(merge.123)[which(!grepl("(\\.x$)|(\\.y$)|(merge)", names(merge.123)))] = paste(colnames(merge.123)[which(!grepl("(\\.x$)|(\\.y$)|(merge)", names(merge.123), perl=T))], ".z", sep="")
			#merge.123$thresholdValue = pmax(merge.123[,onx], merge.123[,ony], merge.123[,onz])

			patientMerge = rbind(patientMerge, merge.123)
			patient.fuzzy = patient.fuzzy[!(match.filter.1 | match.filter.2 | match.filter.3),]

			hidden.clone.sequences = c(rows.1[-1,"Clone_Sequence"], rows.2[rows.2$Clone_Sequence != first.clone.sequence,"Clone_Sequence"], rows.3[rows.3$Clone_Sequence != first.clone.sequence,"Clone_Sequence"])
			merge.list[["second"]] = append(merge.list[["second"]], hidden.clone.sequences)

		} else if (matches.in.2) {
			#other.sample1 = other.sample.list[[first.sample]][1]
			#other.sample2 = other.sample.list[[first.sample]][2]

			second.sample = sum.2[,"Sample"]

			current.merge.list = duo.merge.list[[paste(first.sample, second.sample)]]

			merge.12 = merge(sum.1, sum.2, by="merge")

			current.merge.list = rbind(current.merge.list, merge.12)
			duo.merge.list[[paste(first.sample, second.sample)]] = current.merge.list

			patient.fuzzy = patient.fuzzy[!(match.filter.1 | match.filter.2),]

			hidden.clone.sequences = c(rows.1[-1,"Clone_Sequence"], rows.2[rows.2$Clone_Sequence != first.clone.sequence,"Clone_Sequence"])
			merge.list[["second"]] = append(merge.list[["second"]], hidden.clone.sequences)

		} else if (matches.in.3) {

			#other.sample1 = other.sample.list[[first.sample]][1]
			#other.sample2 = other.sample.list[[first.sample]][2]

			second.sample = sum.3[,"Sample"]

			current.merge.list = duo.merge.list[[paste(first.sample, second.sample)]]

			merge.13 = merge(sum.1, sum.3, by="merge")

			current.merge.list = rbind(current.merge.list, merge.13)
			duo.merge.list[[paste(first.sample, second.sample)]] = current.merge.list

			patient.fuzzy = patient.fuzzy[!(match.filter.1 | match.filter.3),]

			hidden.clone.sequences = c(rows.1[-1,"Clone_Sequence"], rows.3[rows.3$Clone_Sequence != first.clone.sequence,"Clone_Sequence"])
			merge.list[["second"]] = append(merge.list[["second"]], hidden.clone.sequences)

		} else if(nrow(rows.1) > 1){
			patient1 = patient1[!(patient1$Clone_Sequence %in% rows.1$Clone_Sequence),]
			#print(names(patient1)[names(patient1) %in% sum.1])
			#print(names(patient1)[!(names(patient1) %in% sum.1)])
			#print(names(patient1))
			#print(names(sum.1))
			#print(summary(sum.1))
			#print(summary(patient1))
			#print(dim(sum.1))
			#print(dim(patient1))
			#print(head(sum.1[,names(patient1)]))
			patient1 = rbind(patient1, sum.1[,names(patient1)])
			patient.fuzzy = patient.fuzzy[-match.filter.1,]
		} else {
			patient.fuzzy = patient.fuzzy[-1,]
		}

		tmp.rows = rbind(rows.1, rows.2, rows.3)
		tmp.rows = tmp.rows[order(nchar(tmp.rows$Clone_Sequence)),]

		if (sum(match.filter.1) > 1 | sum(match.filter.2) > 1 | sum(match.filter.1) > 1) {
		cat(paste("<tr><td>", label1, " row ", 1:nrow(tmp.rows), "</td><td>", tmp.rows$Sample, ":</td><td>", tmp.rows$Clone_Sequence, "</td><td>", tmp.rows$normalized_read_count, "</td></tr>", sep=""), file="multiple_matches.html", append=T)
		} else {
		}

	}
		patient.merge.list[[paste(label1, "123")]] = patientMerge

		patientMerge12 = duo.merge.list[[paste(oneSample, twoSample)]]
		patientMerge13 = duo.merge.list[[paste(oneSample, threeSample)]]
		patientMerge23 = duo.merge.list[[paste(twoSample, threeSample)]]

		patient.merge.list[[paste(label1, "12")]] = patientMerge12
		patient.merge.list[[paste(label1, "13")]] = patientMerge13
		patient.merge.list[[paste(label1, "23")]] = patientMerge23

		#patient.merge.list.second[[label1]] = merge.list[["second"]]
	}
	cat(paste("<td>", nrow(patient1), " in ", label1, " and ", nrow(patient2), " in ", label2, nrow(patient3), " in ", label3, ", ", nrow(patientMerge), " in both (finding both took ", (proc.time() - start.time)[[3]], "s)</td></tr>", sep=""), file=logfile, append=T)
	patientMerge$thresholdValue = pmax(patientMerge[,onx], patientMerge[,ony], patientMerge[,onz])
	patientMerge12$thresholdValue = pmax(patientMerge12[,onx], patientMerge12[,ony])
	patientMerge13$thresholdValue = pmax(patientMerge13[,onx], patientMerge13[,ony])
	patientMerge23$thresholdValue = pmax(patientMerge23[,onx], patientMerge23[,ony])

	#patientMerge$thresholdValue = pmin(patientMerge[,onx], patientMerge[,ony], patientMerge[,onz])
	#patientMerge12$thresholdValue = pmin(patientMerge12[,onx], patientMerge12[,ony])
	#patientMerge13$thresholdValue = pmin(patientMerge13[,onx], patientMerge13[,ony])
	#patientMerge23$thresholdValue = pmin(patientMerge23[,onx], patientMerge23[,ony])

	patient1 = patient1[!(patient1$Clone_Sequence %in% merge.list[["second"]]),]
	patient2 = patient2[!(patient2$Clone_Sequence %in% merge.list[["second"]]),]
	patient3 = patient3[!(patient3$Clone_Sequence %in% merge.list[["second"]]),]

	if(F){
		patientMerge = merge(patient1, patient2, by="merge")
		patientMerge = merge(patientMerge, patient3, by="merge")
		colnames(patientMerge)[which(!grepl("(\\.x$)|(\\.y$)|(merge)", names(patientMerge)))] = paste(colnames(patientMerge)[which(!grepl("(\\.x$)|(\\.y$)|(merge)", names(patientMerge), perl=T))], ".z", sep="")
		patientMerge$thresholdValue = pmax(patientMerge[,onx], patientMerge[,ony], patientMerge[,onz])
		patientMerge12 = merge(patient1, patient2, by="merge")
		patientMerge12$thresholdValue = pmax(patientMerge12[,onx], patientMerge12[,ony])
		patientMerge13 = merge(patient1, patient3, by="merge")
		patientMerge13$thresholdValue = pmax(patientMerge13[,onx], patientMerge13[,ony])
		patientMerge23 = merge(patient2, patient3, by="merge")
		patientMerge23$thresholdValue = pmax(patientMerge23[,onx], patientMerge23[,ony])
	}

	scatterplot_data_columns = c("Clone_Sequence", "Frequency", "normalized_read_count", "V_Segment_Major_Gene", "J_Segment_Major_Gene", "merge")
	scatterplot_data = rbind(patient1[,scatterplot_data_columns], patient2[,scatterplot_data_columns], patient3[,scatterplot_data_columns])
	scatterplot_data = scatterplot_data[!duplicated(scatterplot_data$merge),]

	scatterplot_data$type = factor(x="In one", levels=c("In one", "In two", "In three", "In multiple"))

	res1 = vector()
	res2 = vector()
	res3 = vector()
	res12 = vector()
	res13 = vector()
	res23 = vector()
	resAll = vector()
	read1Count = vector()
	read2Count = vector()
	read3Count = vector()

	if(appendTriplets){
		cat(paste(label1, label2, label3, sep="\t"), file="triplets.txt", append=T, sep="", fill=3)
	}
	for(iter in 1:length(product[,1])){
	threshhold = product[iter,threshholdIndex]
	V_Segment = paste(".*", as.character(product[iter,V_SegmentIndex]), ".*", sep="")
	J_Segment = paste(".*", as.character(product[iter,J_SegmentIndex]), ".*", sep="")
	#all = (grepl(V_Segment, patientMerge$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge$J_Segment_Major_Gene.x) & patientMerge[,onx] > threshhold & patientMerge[,ony] > threshhold & patientMerge[,onz] > threshhold) 
	all = (grepl(V_Segment, patientMerge$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge$J_Segment_Major_Gene.x) & patientMerge$thresholdValue > threshhold)

	one_two = (grepl(V_Segment, patientMerge12$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge12$J_Segment_Major_Gene.x) & patientMerge12$thresholdValue > threshhold & !(patientMerge12$merge %in% patientMerge[all,]$merge))
	one_three = (grepl(V_Segment, patientMerge13$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge13$J_Segment_Major_Gene.x) & patientMerge13$thresholdValue > threshhold & !(patientMerge13$merge %in% patientMerge[all,]$merge))
	two_three = (grepl(V_Segment, patientMerge23$V_Segment_Major_Gene.x) & grepl(J_Segment, patientMerge23$J_Segment_Major_Gene.x) & patientMerge23$thresholdValue > threshhold & !(patientMerge23$merge %in% patientMerge[all,]$merge))

	one = (grepl(V_Segment, patient1$V_Segment_Major_Gene) & grepl(J_Segment, patient1$J_Segment_Major_Gene) & patient1[,on] > threshhold & !(patient1$merge %in% patientMerge[all,]$merge) & !(patient1$merge %in% patientMerge12[one_two,]$merge) & !(patient1$merge %in% patientMerge13[one_three,]$merge))
	two = (grepl(V_Segment, patient2$V_Segment_Major_Gene) & grepl(J_Segment, patient2$J_Segment_Major_Gene) & patient2[,on] > threshhold & !(patient2$merge %in% patientMerge[all,]$merge) & !(patient2$merge %in% patientMerge12[one_two,]$merge) & !(patient2$merge %in% patientMerge23[two_three,]$merge))
	three = (grepl(V_Segment, patient3$V_Segment_Major_Gene) & grepl(J_Segment, patient3$J_Segment_Major_Gene) & patient3[,on] > threshhold & !(patient3$merge %in% patientMerge[all,]$merge) & !(patient3$merge %in% patientMerge13[one_three,]$merge) & !(patient3$merge %in% patientMerge23[two_three,]$merge))

	read1Count = append(read1Count, sum(patient1[one,]$normalized_read_count) + sum(patientMerge[all,]$normalized_read_count.x))
	read2Count = append(read2Count, sum(patient2[two,]$normalized_read_count) + sum(patientMerge[all,]$normalized_read_count.y))
	read3Count = append(read3Count, sum(patient3[three,]$normalized_read_count) + sum(patientMerge[all,]$normalized_read_count.z))
	res1 = append(res1, sum(one))
	res2 = append(res2, sum(two))
	res3 = append(res3, sum(three))
	resAll = append(resAll, sum(all))
	res12 = append(res12, sum(one_two))
	res13 = append(res13, sum(one_three))
	res23 = append(res23, sum(two_three))
	#threshhold = 0
	if(threshhold != 0){
		if(sum(one) > 0){
			dfOne = patient1[one,c("V_Segment_Major_Gene", "J_Segment_Major_Gene", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")]
			colnames(dfOne) = c("Proximal segment", "Distal segment", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")
			filenameOne = paste(label1, "_", product[iter, titleIndex], "_", threshhold, sep="")
			write.table(dfOne, file=paste(filenameOne, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
		}
		if(sum(two) > 0){
			dfTwo = patient2[two,c("V_Segment_Major_Gene", "J_Segment_Major_Gene", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")]
			colnames(dfTwo) = c("Proximal segment", "Distal segment", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")
			filenameTwo = paste(label2, "_", product[iter, titleIndex], "_", threshhold, sep="")
			write.table(dfTwo, file=paste(filenameTwo, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
		}
		if(sum(three) > 0){
			dfThree = patient3[three,c("V_Segment_Major_Gene", "J_Segment_Major_Gene", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")]
			colnames(dfThree) = c("Proximal segment", "Distal segment", "normalized_read_count", "Frequency", "Clone_Sequence", "Related_to_leukemia_clone")
			filenameThree = paste(label3, "_", product[iter, titleIndex], "_", threshhold, sep="")
			write.table(dfThree, file=paste(filenameThree, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
		}
		if(sum(one_two) > 0){
			dfOne_two = patientMerge12[one_two,c("V_Segment_Major_Gene.x", "J_Segment_Major_Gene.x", "normalized_read_count.x", "Frequency.x", "Related_to_leukemia_clone.x", "Clone_Sequence.x", "V_Segment_Major_Gene.y", "J_Segment_Major_Gene.y", "normalized_read_count.y", "Frequency.y", "Related_to_leukemia_clone.y")]
			colnames(dfOne_two) = c(paste("Proximal segment", oneSample), paste("Distal segment", oneSample), paste("Normalized_Read_Count", oneSample), paste("Frequency", oneSample), paste("Related_to_leukemia_clone", oneSample),"Clone_Sequence", paste("Proximal segment", twoSample), paste("Distal segment", twoSample), paste("Normalized_Read_Count", twoSample), paste("Frequency", twoSample), paste("Related_to_leukemia_clone", twoSample))
			filenameOne_two = paste(label1, "_", label2, "_", product[iter, titleIndex], "_", threshhold, onShort, sep="")
			write.table(dfOne_two, file=paste(filenameOne_two, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
		}
		if(sum(one_three) > 0){
			dfOne_three = patientMerge13[one_three,c("V_Segment_Major_Gene.x", "J_Segment_Major_Gene.x", "normalized_read_count.x", "Frequency.x", "Related_to_leukemia_clone.x", "Clone_Sequence.x", "V_Segment_Major_Gene.y", "J_Segment_Major_Gene.y", "normalized_read_count.y", "Frequency.y", "Related_to_leukemia_clone.y")]
			colnames(dfOne_three) = c(paste("Proximal segment", oneSample), paste("Distal segment", oneSample), paste("Normalized_Read_Count", oneSample), paste("Frequency", oneSample), paste("Related_to_leukemia_clone", oneSample),"Clone_Sequence", paste("Proximal segment", threeSample), paste("Distal segment", threeSample), paste("Normalized_Read_Count", threeSample), paste("Frequency", threeSample), paste("Related_to_leukemia_clone", threeSample))
			filenameOne_three = paste(label1, "_", label3, "_", product[iter, titleIndex], "_", threshhold, onShort, sep="")
			write.table(dfOne_three, file=paste(filenameOne_three, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
		}
		if(sum(two_three) > 0){
			dfTwo_three = patientMerge23[two_three,c("V_Segment_Major_Gene.x", "J_Segment_Major_Gene.x", "normalized_read_count.x", "Frequency.x", "Related_to_leukemia_clone.x", "Clone_Sequence.x", "V_Segment_Major_Gene.y", "J_Segment_Major_Gene.y", "normalized_read_count.y", "Frequency.y", "Related_to_leukemia_clone.y")]
			colnames(dfTwo_three) = c(paste("Proximal segment", twoSample), paste("Distal segment", twoSample), paste("Normalized_Read_Count", twoSample), paste("Frequency", twoSample), paste("Related_to_leukemia_clone", twoSample),"Clone_Sequence", paste("Proximal segment", threeSample), paste("Distal segment", threeSample), paste("Normalized_Read_Count", threeSample), paste("Frequency", threeSample), paste("Related_to_leukemia_clone", threeSample))
			filenameTwo_three = paste(label2, "_", label3, "_", product[iter, titleIndex], "_", threshhold, onShort, sep="")
			write.table(dfTwo_three, file=paste(filenameTwo_three, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
		}
	} else { #scatterplot data
		scatterplot_locus_data = scatterplot_data[grepl(V_Segment, scatterplot_data$V_Segment_Major_Gene) & grepl(J_Segment, scatterplot_data$J_Segment_Major_Gene),]
		scatterplot_locus_data = scatterplot_locus_data[!(scatterplot_locus_data$merge %in% merge.list[["second"]]),]
		in_two = (scatterplot_locus_data$merge %in% patientMerge12[one_two,]$merge) | (scatterplot_locus_data$merge %in% patientMerge13[one_three,]$merge) | (scatterplot_locus_data$merge %in% patientMerge23[two_three,]$merge)
		if(sum(in_two) > 0){
			scatterplot_locus_data[in_two,]$type = "In two"
		}
		in_three = (scatterplot_locus_data$merge %in% patientMerge[all,]$merge)
		if(sum(in_three)> 0){
			scatterplot_locus_data[in_three,]$type = "In three"
		}
		not_in_one = scatterplot_locus_data$type != "In one"
		if(sum(not_in_one) > 0){
			#scatterplot_locus_data[not_in_one,]$type = "In multiple"
		}
		p = NULL
		if(nrow(scatterplot_locus_data) != 0){
			filename.scatter = paste(label1, "_", label2, "_", label3, "_", product[iter, titleIndex], "_scatter_", threshhold, sep="")
			write.table(scatterplot_locus_data, file=paste(filename.scatter, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
			if(on == "normalized_read_count"){
			scales = 10^(0:6) #(0:ceiling(log10(max(scatterplot_locus_data$normalized_read_count))))
			p = ggplot(scatterplot_locus_data, aes(type, normalized_read_count)) + scale_y_log10(breaks=scales,labels=scales, limits=c(1, 1e6))
			} else {
			p = ggplot(scatterplot_locus_data, aes(type, Frequency)) + scale_y_log10(limits=c(0.0001,100), breaks=c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100), labels=c("0.0001", "0.001", "0.01", "0.1", "1", "10", "100")) + expand_limits(y=c(0,100))
			#p = ggplot(scatterplot_locus_data, aes(type, Frequency)) + scale_y_continuous(limits = c(0, 100)) + expand_limits(y=c(0,100))
			}
			p = p + geom_point(aes(colour=type), position="jitter")
			p = p + xlab("In one or in multiple samples") + ylab(onShort) + ggtitle(paste(label1, label2, label3, onShort, product[iter, titleIndex]))
		} else {
			p = ggplot(NULL, aes(x=c("In one", "In multiple"),y=0)) + geom_blank(NULL) + xlab("In two or in three of the samples") + ylab(onShort) + ggtitle(paste(label1, label2, label3, onShort, product[iter, titleIndex]))
		}
		file_name = paste(label1, "_", label2, "_", label3, "_", onShort, "_", product[iter, titleIndex],"_scatter.png", sep="")
  		print(paste("Writing figure:", file_name))
		png(file_name)
		print(p)
	  dev.off()
	} 
	if(sum(all) > 0){
		dfAll = patientMerge[all,c("V_Segment_Major_Gene.x", "J_Segment_Major_Gene.x", "normalized_read_count.x", "Frequency.x", "Related_to_leukemia_clone.x", "Clone_Sequence.x", "V_Segment_Major_Gene.y", "J_Segment_Major_Gene.y", "normalized_read_count.y", "Frequency.y", "Related_to_leukemia_clone.y", "V_Segment_Major_Gene.z", "J_Segment_Major_Gene.z", "normalized_read_count.z", "Frequency.z", "Related_to_leukemia_clone.z")]
		colnames(dfAll) = c(paste("Proximal segment", oneSample), paste("Distal segment", oneSample), paste("Normalized_Read_Count", oneSample), paste("Frequency", oneSample), paste("Related_to_leukemia_clone", oneSample),"Clone_Sequence", paste("Proximal segment", twoSample), paste("Distal segment", twoSample), paste("Normalized_Read_Count", twoSample), paste("Frequency", twoSample), paste("Related_to_leukemia_clone", twoSample), paste("Proximal segment", threeSample), paste("Distal segment", threeSample), paste("Normalized_Read_Count", threeSample), paste("Frequency", threeSample), paste("Related_to_leukemia_clone", threeSample))
		filenameAll = paste(label1, "_", label2, "_", label3, "_", product[iter, titleIndex], "_", threshhold, sep="")
		write.table(dfAll, file=paste(filenameAll, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
	}
	}
	#patientResult = data.frame("Locus"=product$Titles, "J_Segment"=product$J_Segments, "V_Segment"=product$V_Segments, "cut_off_value"=paste(">", product$interval, sep=""), "All"=resAll, "tmp1"=res1, "read_count1" = round(read1Count), "tmp2"=res2, "read_count2"= round(read2Count), "tmp3"=res3, "read_count3"=round(read3Count))
	patientResult = data.frame("Locus"=product$Titles, "J_Segment"=product$J_Segments, "V_Segment"=product$V_Segments, "cut_off_value"=paste(">", product$interval, sep=""), "All"=resAll, "tmp1"=res1, "tmp2"=res2, "tmp3"=res3, "tmp12"=res12, "tmp13"=res13, "tmp23"=res23)
	colnames(patientResult)[6] = oneSample
	colnames(patientResult)[7] = twoSample
	colnames(patientResult)[8] = threeSample
	colnames(patientResult)[9] = paste(oneSample, twoSample, sep="_")
	colnames(patientResult)[10] = paste(oneSample, twoSample, sep="_")
	colnames(patientResult)[11] = paste(oneSample, twoSample, sep="_")

	colnamesBak = colnames(patientResult)
	colnames(patientResult) = c("Ig/TCR gene rearrangement type", "Distal Gene segment", "Proximal gene segment", "cut_off_value", "Number of sequences All", paste("Number of sequences", oneSample), paste("Number of sequences", twoSample), paste("Number of sequences", threeSample), paste("Number of sequences", oneSample, twoSample), paste("Number of sequences", oneSample, threeSample), paste("Number of sequences", twoSample, threeSample))
	write.table(patientResult, file=paste(label1, "_", label2, "_", label3, "_", onShort, ".txt", sep=""), quote=F, sep="\t", dec=",", row.names=F, col.names=T)
	colnames(patientResult) = colnamesBak

	patientResult$Locus = factor(patientResult$Locus, Titles)
	patientResult$cut_off_value = factor(patientResult$cut_off_value, paste(">", interval, sep=""))

	plt = ggplot(patientResult[,c("Locus", "cut_off_value", "All")])
	plt = plt + geom_bar( aes( x=factor(cut_off_value), y=All), stat='identity', position="dodge", fill="#79c36a")
	plt = plt + facet_grid(.~Locus) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	plt = plt + geom_text(aes(ymax=max(All), x=cut_off_value,y=All,label=All), angle=90, hjust=0)
	plt = plt + xlab("Reads per locus") + ylab("Count") + ggtitle("Number of clones in All")
	plt = plt + theme(plot.margin = unit(c(1,8.8,0.5,1.5), "lines"))
	file_name = paste(label1, "_", label2, "_", label3, "_", onShort, "_total_all.png", sep="")
	print(paste("Writing figure:", file_name))
	png(file_name, width=1920, height=1080)
	print(plt)
	dev.off()

	fontSize = 4

	bak = patientResult
	patientResult = melt(patientResult[,c('Locus','cut_off_value', oneSample, twoSample, threeSample)] ,id.vars=1:2)
	patientResult$relativeValue = patientResult$value * 10
	patientResult[patientResult$relativeValue == 0,]$relativeValue = 1
	plt = ggplot(patientResult)
	plt = plt + geom_bar( aes( x=factor(cut_off_value), y=relativeValue, fill=variable), stat='identity', position="dodge")
	plt = plt + facet_grid(.~Locus) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
	plt = plt + scale_y_continuous(trans="log", breaks=10^c(0:10), labels=c(0, 10^c(0:9)))
	plt = plt + geom_text(data=patientResult[patientResult$variable == oneSample,], aes(ymax=max(value), x=cut_off_value,y=relativeValue,label=value), angle=90, position=position_dodge(width=0.9), hjust=0, vjust=-0.7, size=fontSize)
	plt = plt + geom_text(data=patientResult[patientResult$variable == twoSample,], aes(ymax=max(value), x=cut_off_value,y=relativeValue,label=value), angle=90, position=position_dodge(width=0.9), hjust=0, vjust=0.4, size=fontSize)
	plt = plt + geom_text(data=patientResult[patientResult$variable == threeSample,], aes(ymax=max(value), x=cut_off_value,y=relativeValue,label=value), angle=90, position=position_dodge(width=0.9), hjust=0, vjust=1.5, size=fontSize)
	plt = plt + xlab("Reads per locus") + ylab("Count") + ggtitle("Number of clones in only one sample")
	file_name = paste(label1, "_", label2, "_", label3, "_", onShort, "_indiv_all.png", sep="")
	print(paste("Writing figure:", file_name))
	png(file_name, width=1920, height=1080)
	print(plt)
	dev.off()
}

if(nrow(triplets) != 0){

	cat("<tr><td>Starting triplet analysis</td></tr>", file=logfile, append=T)

	triplets$uniqueID = paste(triplets$Patient, triplets$Sample, sep="_")

	cat("<tr><td>Normalizing to lowest cell count within locus</td></tr>", file=logfile, append=T)

	triplets$locus_V = substring(triplets$V_Segment_Major_Gene, 0, 4)
	triplets$locus_J = substring(triplets$J_Segment_Major_Gene, 0, 4)
	min_cell_count = data.frame(data.table(triplets)[, list(min_cell_count=min(.SD$Cell_Count)), by=c("uniqueID", "locus_V", "locus_J")])

	triplets$min_cell_paste = paste(triplets$uniqueID, triplets$locus_V, triplets$locus_J)
	min_cell_count$min_cell_paste = paste(min_cell_count$uniqueID, min_cell_count$locus_V, min_cell_count$locus_J)

	min_cell_count = min_cell_count[,c("min_cell_paste", "min_cell_count")]

	triplets = merge(triplets, min_cell_count, by="min_cell_paste")

	triplets$normalized_read_count = round(triplets$Clone_Molecule_Count_From_Spikes / triplets$Cell_Count * triplets$min_cell_count / 2, digits=2)

	triplets = triplets[triplets$normalized_read_count >= min_cells,]

	column_drops = c("min_cell_count", "min_cell_paste")

	triplets = triplets[,!(colnames(triplets) %in% column_drops)]

	cat("<tr><td>Starting Cell Count analysis</td></tr>", file=logfile, append=T)

	interval = intervalReads
	intervalOrder = data.frame("interval"=paste(">", interval, sep=""), "intervalOrder"=1:length(interval))
	product = data.frame("Titles"=rep(Titles, each=length(interval)), "interval"=rep(interval, times=10), "V_Segments"=rep(V_Segments, each=length(interval)), "J_Segments"=rep(J_Segments, each=length(interval)))

	triplets = split(triplets, triplets$Patient, drop=T)
	#print(nrow(triplets))
	for(triplet in triplets){
		samples = unique(triplet$Sample)
		one = triplet[triplet$Sample == samples[1],]
		two = triplet[triplet$Sample == samples[2],]
		three = triplet[triplet$Sample == samples[3],]
		
		print(paste(nrow(triplet), nrow(one), nrow(two), nrow(three)))
		tripletAnalysis(one, one[1,"uniqueID"], two, two[1,"uniqueID"], three, three[1,"uniqueID"], product=product, interval=interval, on="normalized_read_count", T)
	}

	cat("<tr><td>Starting Frequency analysis</td></tr>", file=logfile, append=T)

	interval = intervalFreq
	intervalOrder = data.frame("interval"=paste(">", interval, sep=""), "intervalOrder"=1:length(interval))
	product = data.frame("Titles"=rep(Titles, each=length(interval)), "interval"=rep(interval, times=10), "V_Segments"=rep(V_Segments, each=length(interval)), "J_Segments"=rep(J_Segments, each=length(interval)))

	for(triplet in triplets){
		samples = unique(triplet$Sample)
		one = triplet[triplet$Sample == samples[1],]
		two = triplet[triplet$Sample == samples[2],]
		three = triplet[triplet$Sample == samples[3],]
		tripletAnalysis(one, one[1,"uniqueID"], two, two[1,"uniqueID"], three, three[1,"uniqueID"], product=product, interval=interval, on="Frequency", F)
	}
} else {
  cat("", file="triplets.txt")
}
cat("</table></html>", file=logfile, append=T)
