library(ggplot2)
library(data.table)
windowsFonts(Arial=windowsFont("TT Arial"))

needed_patients = c("5784", "6449", "5257", "6991", "8503", "8478", "10666", "10891", "11255", "11288", "4714", "5117", "5756", "5766", "5864", "6190", "6386", "6417", "6558", "6603", "6874", "6892", "6974", "7170", "8099", "8549", "9974", "12401", "14209", "15582", "9563", "9850", "10151", "10265", "10275", "10944", "11441", "11938", "12210", "14559", "5930", "5456", "5657", "261195", "130890", "30390")
needed_patients = c("5784", "6449", "5257", "6991", "8503", "8478", "10666", "10891", "11255", "11288", "4714", "5117", "5756", "5766", "5864", "6190", "6386", "6417", "6558", "6603", "6874", "6892", "6974", "7170", "8099", "8549", "9974", "12401", "14209", "15582", "9563", "9850", "10151", "10265", "10275", "11441", "11938", "12210", "14559", "5930", "5456", "5657", "261195", "130890", "30390")

#input_file = "D:/wd/prisca/heatmaps_vincent/2015-11-20-sorted.txt"

#input_file = "D:/wd/prisca/heatmaps_vincent/vs_2015_fltr.txt"
#input_data1 = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

#names(input_data1)[names(input_data1) == "dsPerm"] = "dsPerM"


#input_file = "D:/wd/prisca/heatmaps_vincent/VanDongen_clones_evol_may2015_complete_2015_8_20.txt"
#input_data2 = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

#input_data2 = input_data2[,names(input_data1)]

#input_data = rbind(input_data1, input_data2)

#input_data = input_data[grepl(paste("VanDongen_cALL_",needed_patients, sep="", collapse="|"), input_data$Patient),]


#write.table(x = input_data, file = "Dx_R_patients_2017.txt", quote = F, append = F, sep = "\t", row.names = F, col.names = T)
#input_file = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017.txt"
#input_data = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

input_file = "D:/wd/_oud/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original.txt"
input_data = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

#input_file = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original_filter.txt"
#input_data = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

#output_dir = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original_output"
#output_dir = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original_output_only_read_filter"
#output_dir = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original_output_freq_read_filter"
#output_dir = "D:/wd/prisca/heatmaps_vincent/2018-01-22"

read_thresholds = c(0,10,25,50,100,250,500,750,1000,10000)
read_thresholds = c(0,25,100,1000,10000)
read_threshold_labels = c("0-25", "25-100", "100-1000", "1000-10000", ">10000")
freq_thresholds_2 = c(0,0.01,0.05,0.1,0.5,1,5)
freq_thresholds = c(0,0.01,0.1,1,5)
freq_threshold_labels = c("0-0.01", "0.01-0.1", "0.1-1", "1-5", ">5")

rearrangement_types = c("Total", "IGH-Vh-Jh", "IGH-Dh-Jh", "Vk-Jk", "Vk-Kde", "Intron-Kde", "TCRG", "TCRD-Vd-Dd", "TCRD-Dd-Dd", "TCRB-Vb-Jb", "IGK", "TCRD")
#rearrangement_types = c("Total", "IGH-Vh-Jh", "IGH-Dh-Jh", "IGK", "TCRG", "TCRD", "TCRB-Vb-Jb")

#patient_samples = unique(input_data[,c("Patient", "Sample")])

#patient_samples = patient_samples[order(patient_samples$Patient),]

#for vs_2015_fltr data
#patients = unique(patient_samples$Patient)
#input_dir = "ziv_rosen_vanDongen_report_2015/"

#for may2015 data
#patients = c("VanDongen_cALL_10666", "VanDongen_cALL_10891", "VanDongen_cALL_11255", "VanDongen_cALL_11288", "VanDongen_cALL_12401", "VanDongen_cALL_14209", "VanDongen_cALL_15582", "VanDongen_cALL_4714", "VanDongen_cALL_5117", "VanDongen_cALL_5257", "VanDongen_cALL_5756", "VanDongen_cALL_5766", "VanDongen_cALL_5784", "VanDongen_cALL_5864", "VanDongen_cALL_6190", "VanDongen_cALL_6386", "VanDongen_cALL_6417", "VanDongen_cALL_6449", "VanDongen_cALL_6558", "VanDongen_cALL_6603", "VanDongen_cALL_6874", "VanDongen_cALL_6892", "VanDongen_cALL_6974", "VanDongen_cALL_6991", "VanDongen_cALL_7170", "VanDongen_cALL_8099", "VanDongen_cALL_8478", "VanDongen_cALL_8503", "VanDongen_cALL_8549", "VanDongen_cALL_9974")
#input_dir = "VanDongen_clones_evol_may2015_complete_2015_8_20/"

#for 40 patients from both
#input_dir = "may2015_ziv_rosen_combined/"

input_data$Frequency = ((10^input_data$Log10_Frequency)*100)

result_list = list()

#merge the files of multiple locus into new files for the a combined locus
merge.locus = function(patient.list, locus.list, new.locus.name, work_dir){ 
  for(patient in patient.list){
    new_patient_data = F
    new_sample1_data = F
    new_sample2_data = F
    this_patient_samples = patient_samples[patient_samples$Patient == patient,]
    sample1 = as.character(this_patient_samples[1, "Sample"])
    sample2 = as.character(this_patient_samples[2, "Sample"])
    
    if(grepl("_Dx", sample2)){
      print(paste("Switching", sample1, "with", sample2))
      tmp = sample1
      sample1 = sample2
      sample2 = tmp
    }
    
    #combined.patient.data = read.table(file = file.path(work_dir, paste(sample1, sample2, "Total", "0.txt", sep="_"), sep=""), header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")[NULL,]
    #combined.patient.data = read.table(file = file.path(work_dir, paste(sample2, sample1, "Total", "0.txt", sep="_"), sep=""), header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")[NULL,]
    combined.patient.data = combined.patient.data[NULL,]
    
    #combined.sample1.data = read.table(file = file.path(work_dir, paste(sample1, "Total", "0.txt", sep="_"), sep=""), header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
    combined.sample1.data = combined.sample1.data[NULL,] # file.path(work_dir, paste(sample1, "Total", "0.txt", sep="_"), sep="")[NULL,]

    #combined.sample2.data = read.table(file = file.path(work_dir, paste(sample2, "Total", "0.txt", sep="_"), sep=""), header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
    combined.sample2.data = combined.sample2.data[NULL,]
    
    for(locus in locus.list){
      print(paste(patient, locus))
      
      patient_data_file = file.path(work_dir, paste(sample1, sample2, locus, "0.txt", sep="_"), sep="")
      if(!file.exists(patient_data_file)){
        patient_data_file = file.path(work_dir, paste(sample2, sample1, locus, "0.txt", sep="_"), sep="")
      }
      
      if(file.exists(patient_data_file)){
        patient_data = read.table(file = patient_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
        new_patient_data = T
        print(paste(patient, "patient", nrow(patient_data)))
        if(is.null(combined.patient.data)){
          combined.patient.data = patient_data
        } else {
          combined.patient.data = rbind(combined.patient.data, patient_data)
        }
      } else {
        print(paste(patient, "patient", "doesn't exist"))
      }

      sample1_data_file = file.path(work_dir, paste(sample1, locus, "0.txt", sep="_"), sep="")
      if(file.exists(sample1_data_file)){
        sample1_data = read.table(file = sample1_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
        new_sample1_data = T
        print(paste(patient, "sample1", nrow(sample1_data)))
        if(is.null(combined.sample1.data)){
          combined.sample1.data = sample1_data
        } else {
          combined.sample1.data = rbind(combined.sample1.data, sample1_data)
        }
      } else {
        print(paste(patient, "sample1", "doesn't exist:", sample1_data_file))
      }
      
      sample2_data_file = file.path(work_dir, paste(sample2, locus, "0.txt", sep="_"), sep="")
      if(file.exists(sample2_data_file)){
        sample2_data = read.table(file = sample2_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
        new_sample2_data = T
        print(paste(patient, "sample2", nrow(sample2_data)))
        if(is.null(combined.sample2.data)){
          combined.sample2.data = sample2_data
        } else {
          combined.sample2.data = rbind(combined.sample2.data, sample2_data)
        }
      } else {
        print(paste(patient, "sample2", "doesn't exist"))
      }
    }
    
    if(new_patient_data){
      print(paste(patient, "combined patient", nrow(combined.patient.data)))
      combined.patient.data.file = file.path(work_dir, paste(sample1, sample2, new.locus.name, "0.txt", sep="_"))
      write.table(x = combined.patient.data, file = combined.patient.data.file, quote = F, row.names = F, col.names = T, sep = "\t", dec = ",")
    } else {
      print(paste("No new data for", sample1, "and", sample2, "would have written:", names(combined.patient.data)))
    }
    
    if(new_sample1_data){
      print(paste(patient, "combined sample1", nrow(combined.sample1.data)))
      combined.sample1.data.file = file.path(work_dir, paste(sample1, new.locus.name, "0.txt", sep="_"))
      write.table(x = combined.sample1.data, file = combined.sample1.data.file, quote = F, row.names = F, col.names = T, sep = "\t", dec = ",")
    } else {
      print(paste("No new data for", sample1, "would have written:", names(combined.sample1.data)))
    }
    
    if(new_sample2_data){
      print(paste(patient, "combined sample2", nrow(combined.sample2.data)))
      combined.sample2.data.file = file.path(work_dir, paste(sample2, new.locus.name, "0.txt", sep="_"))
      write.table(x = combined.sample2.data, file = combined.sample2.data.file, quote = F, row.names = F, col.names = T, sep = "\t", dec = ",")
    } else {
      print(paste("No new data for", sample2, "would have written:", names(combined.sample2.data)))
    }
  }
}

merge.locus(unique(patient_samples$Patient), c("TCRD-Vd-Dd", "TCRD-Dd-Dd"), "TCRD", "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original")
merge.locus(unique(patient_samples$Patient), c("Vk-Jk", "Vk-Kde", "Intron-Kde"), "IGK", "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original")

write_fasta = function(id, seq, fasta_file){
  cat(c(paste(">", id, sep=""), seq), sep="\n", append=TRUE, file=fasta_file)
}

write_patient_to_fasta = function(dat, patient, fasta_file){
  if(nrow(dat) > 0){
    for(i in 1:nrow(dat)){
      write_fasta(
        paste(
          patient, 
          "Dx/R", 
          dat[i,"Proximal.segment.Dx"],
          dat[i,"Distal.segment.Dx"],
          dat[i,"Normalized_Read_Count.Dx"],
          dat[i,"Frequency.Dx"],
          dat[i,"Proximal.segment.R"],
          dat[i,"Distal.segment.R"],
          dat[i,"Normalized_Read_Count.R"],
          dat[i,"Frequency.R"]
        ),
        dat[i,"Clone.Sequence"],
        fasta_file
      )
    }
  }
}

write_sample_to_fasta = function(dat, patient, sample, fasta_file){
  if(nrow(dat) > 0){
    for(i in 1:nrow(dat)){
      write_fasta(
        paste(
          patient, 
          sample, 
          dat[i,"Proximal.segment"],
          dat[i,"Distal.segment"],
          dat[i,"Normalized_Read_Count"],
          dat[i,"Frequency"]
        ),
        dat[i,"Clone.Sequence"],
        fasta_file
      )
    }
  }
}

filter_df_for_paper_threshold = function(df, on="Normalized_Read_Count", on.v="V_Segment_Major_Gene", on.j="J_Segment_Major_Gene"){
  minimum_read_threshold = data.frame(
    V_segments = as.character(c("IGHV", "IGHD", "IGKV", "IGKV", "IgKINTR", "TRGV", "TRDV", "TRDD" , "TRBV")),
    J_segments = as.character(c(".*", ".*", "IGKJ", "KDE", ".*", ".*", ".*", ".*", ".*")),
    read_thresholds = c(30, 30, 70, 40, 170, 50, 20, 30, 30)
  )
  for(i in 1:nrow(minimum_read_threshold)){
    V_segment = minimum_read_threshold[i, "V_segments"]
    J_segment = minimum_read_threshold[i, "J_segments"]
    read_threshold = minimum_read_threshold[i, "read_thresholds"]
    clones_in_locus = grepl(V_segment, df[,on.v]) & grepl(J_segment, df[,on.j])
    fltr = df[,on] >= read_threshold
    df = df[fltr | !clones_in_locus,]
  }
  return(df)
}

#patient_ids = unique(patient_samples$Patient)
#input_dir = "D:/wd/_oud/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original"
#rearrangement_type = "Total"
#output_dir = "D:/wd/_oud/prisca/heatmaps_vincent/tmp"
#only_both = F
#diag_freq=0
#recid_freq=0

outside_tmp = NULL

sum_heatmaps_data = list()

post_prisca_patient_analysis = function(patient_ids, input_dir, rearrangement_type, output_dir, diag_freq=0, recid_freq=0, only_both=F) {
  if(!dir.exists(input_dir)){
    print(paste("Input dir doesn't exist:", input_dir))
    return()
  }
  if(!dir.exists(output_dir)){
    print(paste("Output dir doesn't exist, creating:", output_dir))
    dir.create(output_dir)
  }
  
  large_clones_fasta = file.path(output_dir, "large_clones.fasta")
  
  if(file.exists(large_clones_fasta)){
    file.remove(large_clones_fasta)
  }
  
  sample1.bins = NULL
  sample2.bins = NULL
  all.bins.samples = NULL
  all.bins.both = NULL
  all_patient_data = NULL
  all_sample_data = NULL
  
  rearrangment_list = list()
  
  for(patient in patient_ids){
    this_patient_samples = patient_samples[patient_samples$Patient == patient,]
    sample1 = as.character(this_patient_samples[1, "Sample"])
    sample2 = as.character(this_patient_samples[2, "Sample"])
    
    if(grepl("_Dx", sample2)){ #switch samples so sample 1 is diagnosis
      print(paste("Switching", sample1, "with", sample2))
      tmp = sample1
      sample1 = sample2
      sample2 = tmp
    }
    
    #find the file with all the sequences for this patient and 
    patient_data_file = file.path(input_dir, paste(sample1, sample2, rearrangement_type, "0.txt", sep="_"), sep="")
    if(!file.exists(patient_data_file)){
      patient_data_file = file.path(input_dir, paste(sample2, sample1, rearrangement_type, "0.txt", sep="_"), sep="")
    }
    
    #create an empty patient file as a fallback
    patient_data = data.frame(Clone.Sequence="")
    patient_data[,paste(rep(c("Proximal.segment", "Distal.segment", "Related_to_leukemia_clone"), each=2), c(sample1, sample2), sep=".")] = ""
    patient_data[,paste(rep(c("Normalized_Read_Count", "Frequency"), each=2), c(sample1, sample2), sep=".")] = 0
    patient_data$Normalized_Read_Count = 0
    patient_data = patient_data[NULL,]
    
    #if the patient file exists, read it
    if(file.exists(patient_data_file)){
      patient_data = read.table(file = patient_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
      patient_data$V_Segment_Major_Gene = patient_data[,paste("Proximal.segment", sample1, sep=".")]
      patient_data$J_Segment_Major_Gene = patient_data[,paste("Distal.segment", sample1, sep=".")]
      patient_data$Normalized_Read_Count = pmax(
        patient_data[,paste("Normalized_Read_Count", sample1, sep=".")],
        patient_data[,paste("Normalized_Read_Count", sample2, sep=".")]
      )
      patient_data = filter_df_for_paper_threshold(patient_data, on="Normalized_Read_Count", on.v = "V_Segment_Major_Gene", on.j = "J_Segment_Major_Gene")
    }
    patient_data = patient_data[,names(patient_data)[!(names(patient_data) %in% c("V_Segment_Major_Gene", "J_Segment_Major_Gene", "Normalized_Read_Count"))]]
    patient_data = patient_data[
      patient_data[,paste("Frequency", sample1, sep=".")] >= diag_freq &
      patient_data[,paste("Frequency", sample2, sep=".")] >= recid_freq,
    ]
    
    if(any(patient_data[,paste("Frequency", sample1, sep=".")] > 100 | patient_data[,paste("Frequency", sample2, sep=".")] > 100)){
      #print(paste("--------------------------------", patient, "-----------------------------------"))
      patient_data = patient_data[patient_data[,paste("Frequency", sample1, sep=".")] <= 100 & patient_data[,paste("Frequency", sample2, sep=".")] <= 100,]
    }
    
    #make a new data.frame with Dx/R colnames instead of sample1/sample2
    tmp = data.frame(patient_data)
    names(tmp)[grepl(sample1, names(tmp))] = gsub(sample1, "Dx", names(tmp)[grepl(sample1, names(tmp))])
    names(tmp)[grepl(sample2, names(tmp))] = gsub(sample2, "R", names(tmp)[grepl(sample2, names(tmp))])
    
    #write biggest clones to fasta
    if(rearrangement_type == "Total"){
      write_patient_to_fasta(tmp[tmp$Frequency.Dx > 5,], patient, large_clones_fasta)
    }
    
    #rbind it to all previous patients
    if(nrow(tmp) > 0){
      tmp$Patient = patient
      if(is.null(all_patient_data)){
        all_patient_data = tmp
      } else {
        all_patient_data = rbind(all_patient_data, tmp)
      }
    }
    
    #read the data for sample1 from the file
    sample1_data_file = file.path(input_dir, paste(sample1, rearrangement_type, "0.txt", sep="_"), sep="")
    sample1_data = data.frame("Proximal segment"=character(0), "Distal segment"=character(0), "normalized_read_count"=numeric(0), "Frequency"=numeric(0), "Clone Sequence"=character(0), "Related_to_leukemia_clone"=character(0))
    if(file.exists(sample1_data_file) & !only_both){
      sample1_data = read.table(file = sample1_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
      
      if(rearrangement_type == "Total"){
        write_sample_to_fasta(sample1_data[sample1_data$Frequency > 5,], patient, sample1, large_clones_fasta)
      }
      
      #just the read count and frequency
      sample1_data = sample1_data[,c("normalized_read_count", "Frequency", "Proximal.segment", "Distal.segment")]
      patient_sample1_data = patient_data[,paste(c("Normalized_Read_Count", "Frequency"), sample1, sep=".")]
      names(patient_sample1_data) = c("normalized_read_count", "Frequency")
      #sample1_data = rbind(sample1_data, patient_sample1_data)
      
      sample1_data$normalized_read_count = as.numeric(sample1_data$normalized_read_count)
      sample1_data[is.na(sample1_data$normalized_read_count), "normalized_read_count"] = 0
      sample1_data$Frequency = as.numeric(sample1_data$Frequency)
      sample1_data[is.na(sample1_data$Frequency), "Frequency"] = 0
      sample1_data = filter_df_for_paper_threshold(sample1_data, on="normalized_read_count", on.v = "Proximal.segment", on.j = "Distal.segment")
      sample1_data = sample1_data[sample1_data$Frequency >= diag_freq,]
    }
    
    #collect the rearrangements + recid freq/RC with freq>5 and Rc>10k
    if(nrow(sample1_data) > 0){
      sample1_rearrangements = data.frame(
        Proximal.segment.Dx = sample1_data$Proximal.segment,
        Distal.segment.Dx = sample1_data$Distal.segment,
        Normalized_Read_Count.Dx = sample1_data$normalized_read_count, 
        Frequency.Dx = sample1_data$Frequency,
        Normalized_Read_Count.R = NA, 
        Frequency.R = NA
      )
      rearrangment_list[[patient]] = rbind(
        tmp[, c("Normalized_Read_Count.Dx", "Frequency.Dx", "Proximal.segment.Dx", "Distal.segment.Dx", "Normalized_Read_Count.R", "Frequency.R")],
        sample1_rearrangements
      )
    } else {
      rearrangment_list[[patient]] = tmp[, 
        c("Normalized_Read_Count.Dx", "Frequency.Dx", "Proximal.segment.Dx", "Distal.segment.Dx", "Normalized_Read_Count.R", "Frequency.R")
      ]
    }
    
    if(nrow(sample1_data) > 0){
      tmp = data.frame(sample1_data)
      tmp$type = "Diagnose"
      
      if(is.null(all_sample_data)){
        all_sample_data = tmp
      } else {
        all_sample_data = rbind(all_sample_data, tmp)
      }
    }
    
    #read the data for sample2 from the file
    sample2_data_file = file.path(input_dir, paste(sample2, rearrangement_type, "0.txt", sep="_"), sep="")
    sample2_data = data.frame("Proximal segment"=character(0), "Distal segment"=character(0), "normalized_read_count"=numeric(0), "Frequency"=numeric(0), "Clone Sequence"=character(0), "Related_to_leukemia_clone"=character(0))
    if(file.exists(sample2_data_file) & !only_both){
      sample2_data = read.table(file = sample2_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
      
      if(rearrangement_type == "Total"){
        write_sample_to_fasta(sample2_data[sample2_data$Frequency > 5,], patient, sample2, large_clones_fasta)
      }
      
      #just the read count and frequency
      sample2_data = sample2_data[,c("normalized_read_count", "Frequency", "Proximal.segment", "Distal.segment")]
      patient_sample2_data = patient_data[,paste(c("Normalized_Read_Count", "Frequency"), sample2, sep=".")]
      names(patient_sample2_data) = c("normalized_read_count", "Frequency")
      #sample2_data = rbind(sample2_data, patient_sample2_data)
      
      sample2_data$normalized_read_count = as.numeric(sample2_data$normalized_read_count)
      sample2_data[is.na(sample2_data$normalized_read_count), "normalized_read_count"] = 0
      sample2_data$Frequency = as.numeric(sample2_data$Frequency)
      sample2_data[is.na(sample2_data$Frequency), "Frequency"] = 0
      sample2_data = filter_df_for_paper_threshold(sample2_data, on="normalized_read_count", on.v = "Proximal.segment", on.j = "Distal.segment")
      sample2_data = sample2_data[sample2_data$Frequency >= recid_freq,]
    }
    if(nrow(sample2_data) > 0){
      tmp = data.frame(sample2_data)
      tmp$type = "Recidief"
      
      if(is.null(all_sample_data)){
        all_sample_data = tmp
      } else {
        all_sample_data = rbind(all_sample_data, tmp)
      }
    }
  
  
    print(paste("Both: ", nrow(patient_data), ", ", sample1, ": ", nrow(sample1_data), ", ", sample2, ": ", nrow(sample2_data)))
    
    #bin the data
    sample1_data$read_bin = cut(
      sample1_data$normalized_read_count, 
      breaks=c(read_thresholds, Inf), 
      labels=read_threshold_labels #c("0-10", "10-25", "25-50", "50-100", "100-250", "250-500", "500-750", "750-1000", "1000-10000", ">10000")
    )
    sample1_data$freq_bin = cut(
      sample1_data$Frequency, 
      breaks=c(freq_thresholds, Inf), 
      labels=freq_threshold_labels #c("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.5", "0.5-1", "1-5", ">5")
    )
    
    sample2_data$read_bin = cut(
      sample2_data$normalized_read_count,
      breaks=c(read_thresholds, Inf),
      labels=read_threshold_labels #c("0-10", "10-25", "25-50", "50-100", "100-250", "250-500", "500-750", "750-1000", "1000-10000", ">10000")
    )
    sample2_data$freq_bin = cut(
      sample2_data$Frequency, 
      breaks=c(freq_thresholds, Inf), 
      labels=freq_threshold_labels #c("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.5", "0.5-1", "1-5", ">5")
    )
    
    sample1_bins = data.frame(table(sample1_data$read_bin, sample1_data$freq_bin))
    names(sample1_bins) = c("Reads", "Freq", "dx.count")
    
    sample1_bins$log.dx.count = log(sample1_bins$dx.count)
    sample1_bins[sample1_bins$log.dx.count == -Inf, "log.dx.count"] = 0
    
    sample2_bins = data.frame(table(sample2_data$read_bin, sample2_data$freq_bin))
    names(sample2_bins) = c("Reads", "Freq", "r.count")
    
    sample2_bins$log.r.count = log(sample2_bins$r.count)
    sample2_bins[sample2_bins$log.r.count == -Inf, "log.r.count"] = 0
    
    sample1_bins$r.count = sample2_bins$r.count
    sample1_bins$r.percentage = sample1_bins$r.count / sample1_bins$dx.count * 100
    sample1_bins[sample1_bins$r.percentage %in% c(NaN, Inf),"r.percentage"] = 0
    
    patient_bins = data.frame( # both reads vs both freq
      table(
        cut(
          pmax(patient_data[,paste("Normalized_Read_Count", sample1, sep=".")], patient_data[,paste("Normalized_Read_Count", sample2, sep=".")]),
          breaks=c(read_thresholds, Inf), 
          labels=c("0-25", "25-100", "100-1000", "1000-10000", ">10000") #c("0-10", "10-25", "25-50", "50-100", "100-250", "250-500", "500-750", "750-1000", "1000-10000", ">10000")
        ),
        cut(
          pmax(patient_data[,paste("Frequency", sample1, sep=".")], patient_data[,paste("Frequency", sample2, sep=".")]),
          breaks=c(freq_thresholds, Inf), 
          labels=c("0-0.01", "0.01-0.1", "0.1-1", "1-5", ">5") #c("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.5", "0.5-1", "1-5", ">5")
        )
      )
    )
    names(patient_bins) = c("Reads", "Freq", "both.count")

    bins = merge(patient_bins, sample1_bins[c("Reads", "Freq", "dx.count", "r.count")], by=c("Reads", "Freq"))
    bins$dx.count.log = log(bins$dx.count)
    bins[bins$dx.count.log == -Inf,"dx.count.log"] = 0
    
    bins$total.count = bins$dx.count + bins$both.count
    
    if(is.null(all.bins.both)){
      tmp.bins = bins
      names(tmp.bins)[3:length(names(tmp.bins))] = gsub("$", patient, names(tmp.bins)[3:length(names(tmp.bins))])
      all.bins.both = tmp.bins
    } else {
      tmp.bins = bins
      names(tmp.bins)[3:length(names(tmp.bins))] = gsub("$", patient, names(tmp.bins)[3:length(names(tmp.bins))])
      all.bins.both = merge(all.bins.both, tmp.bins, by=c("Reads", "Freq"))
    }
    
    if(is.null(all.bins.samples)){
      tmp.bins = sample1_bins[,c("Reads", "Freq", "dx.count", "r.count")]
      names(tmp.bins)[3:length(names(tmp.bins))] = gsub("$", patient, names(tmp.bins)[3:length(names(tmp.bins))])
      all.bins.samples = tmp.bins
    } else {
      tmp.bins = sample1_bins[,c("Reads", "Freq", "dx.count", "r.count")]
      names(tmp.bins)[3:length(names(tmp.bins))] = gsub("$", patient, names(tmp.bins)[3:length(names(tmp.bins))])
      all.bins.samples = merge(all.bins.samples, tmp.bins, by=c("Reads", "Freq"))
    }
    
    # stability_per_locus_patient_frequency
    
    stability_per_locus_patient_frequency = data.frame(patient="", diag_freq=0, recid_freq=0, type="", locus="", count=0, stringsAsFactors = F)
    stability_per_locus_patient_frequency_file = file.path(output_dir, "stability_per_locus_patient_frequency_long.txt")
    if(!file.exists(stability_per_locus_patient_frequency_file)){ #make sure column names are there
      write.table(x = stability_per_locus_patient_frequency[-1,], file = stability_per_locus_patient_frequency_file, quote = F, sep = "\t", row.names = F, col.names = T)
    }
    stability_per_locus_patient_frequency = rbind(stability_per_locus_patient_frequency, c(patient, diag_freq, recid_freq, "diag", rearrangement_type, nrow(sample1_data)))
    stability_per_locus_patient_frequency = rbind(stability_per_locus_patient_frequency, c(patient, diag_freq, recid_freq, "recid", rearrangement_type, nrow(sample2_data)))
    stability_per_locus_patient_frequency = rbind(stability_per_locus_patient_frequency, c(patient, diag_freq, recid_freq, "both", rearrangement_type, nrow(patient_data)))
    
    write.table(x = stability_per_locus_patient_frequency[-1,], file = stability_per_locus_patient_frequency_file, quote = F, sep = "\t", row.names = F, col.names = F, append = T)
  }
  print("sum/mean plots...")
  
  if(diag_freq == 0 & recid_freq == 0){ #collect the rearrangements + recid freq/RC with freq>5 and Rc>10k
    for(rearrangement_freq_filter in c(5)){
      for(rearrangement_rc_filter in c(1000, 10000)){
        max.nrow = max(sapply( #find the maximum rows needed after filtering for freq/rc
          rearrangment_list, 
          function(tmp_df) sum(tmp_df$Normalized_Read_Count.Dx >= rearrangement_rc_filter & tmp_df$Frequency.Dx >= rearrangement_freq_filter)
        ))
        rearrangment_df = data.frame(empty=1:max.nrow)
        for(patient_id in names(rearrangment_list)){
          patient_rearrangement = rearrangment_list[[patient_id]]
          patient_rearrangement = patient_rearrangement[
            patient_rearrangement$Normalized_Read_Count.Dx >= rearrangement_rc_filter & patient_rearrangement$Frequency.Dx >= rearrangement_freq_filter,
          ]
          rearrangment_df[,paste(patient_id, "Proximal.segment.Dx", sep=".")] = patient_rearrangement[1:max.nrow,"Proximal.segment.Dx"]
          rearrangment_df[,paste(patient_id, "Distal.segment.Dx", sep=".")] = patient_rearrangement[1:max.nrow,"Distal.segment.Dx"]
          rearrangment_df[,paste(patient_id, "Normalized_Read_Count.R", sep=".")] = patient_rearrangement[1:max.nrow,"Normalized_Read_Count.R"]
          rearrangment_df[,paste(patient_id, "Frequency.R", sep=".")] = patient_rearrangement[1:max.nrow,"Frequency.R"]
        }
        
        rearrangement_data_path = file.path(output_dir, paste("rearrangement_per_patient_freq_", rearrangement_freq_filter, "_rc_", rearrangement_rc_filter, "_", rearrangement_type, ".txt", sep=""))
        
        write.table(x = rearrangment_df, file = rearrangement_data_path, quote = F, row.names = F, col.names = T, sep = "\t", dec = ",")
      }
    }
  }
  
  all.bins.both[is.na(all.bins.both)] = 0
  all.bins.both$Reads = as.character(all.bins.both$Reads)
  all.bins.both$Freq = as.character(all.bins.both$Freq)
  
  all.bins.both$sum.count = rowSums(all.bins.both[,names(all.bins.both)[grepl("^dx.countV", names(all.bins.both))]])
  all.bins.both$sum.both.count = rowSums(all.bins.both[,names(all.bins.both)[grepl("^both.count", names(all.bins.both))]])
  all.bins.both$sum.total.count = rowSums(all.bins.both[,names(all.bins.both)[grepl("^total.count", names(all.bins.both))]])
  
  all.bins.both$mean.count = rowMeans(all.bins.both[,names(all.bins.both)[grepl("^dx.countV", names(all.bins.both))]])
  all.bins.both$mean.both.count = rowMeans(all.bins.both[,names(all.bins.both)[grepl("^both.count", names(all.bins.both))]])
  all.bins.both$mean.total.count = rowMeans(all.bins.both[,names(all.bins.both)[grepl("^total.count", names(all.bins.both))]])
  
  all.bins.both$sum.count.log = log(all.bins.both$sum.count)
  
  all.bins.both$mean.count.log = log(all.bins.both$mean.count)
  
  cols = c("Reads", "Freq", "sum.count", "mean.count", "sum.both.count", "mean.both.count", "sum.total.count", "mean.total.count", "sum.count.log", "mean.count.log")
  cols.sub = c("sum.count", "mean.count", "sum.both.count", "mean.both.count", "sum.total.count", "mean.total.count", "sum.count.log", "mean.count.log")
  all.bins.both.tmp = all.bins.both[,cols]
  for(freq_threshold_label in freq_threshold_labels){ #create the sums for the heatmap
    tmp = data.frame(
      Freq=freq_threshold_label,
      Reads = "Sum"
    )
    tmp[,cols.sub] = colSums(all.bins.both[all.bins.both$Freq == freq_threshold_label, cols.sub])
    all.bins.both.tmp = rbind(all.bins.both.tmp, tmp)
  }
  for(read_threshold_label in read_threshold_labels){
    tmp = data.frame(
      Freq="Sum",
      Reads = read_threshold_label
    )
    tmp[,cols.sub] = colSums(all.bins.both[all.bins.both$Reads == read_threshold_label, cols.sub])
    all.bins.both.tmp = rbind(all.bins.both.tmp, tmp)
  }
  
  sum.sum = data.frame(t(colSums(all.bins.both.tmp[all.bins.both.tmp$Reads == "Sum", cols.sub])))
  sum.sum$Reads = "Sum"
  sum.sum$Freq = "Sum"
  all.bins.both.tmp = rbind(all.bins.both.tmp, sum.sum)
  
  all.bins.both = all.bins.both.tmp
  
  all.bins.both$perc.in.both.sum = all.bins.both$sum.both.count /  all.bins.both$sum.total.count * 100
  all.bins.both$perc.in.both.mean = all.bins.both$mean.both.count /  all.bins.both$mean.total.count * 100
  
  all.bins.both$heatmap.text.sum = paste(round(all.bins.both$perc.in.both.sum, 1), "\n", all.bins.both$sum.both.count, "/", all.bins.both$sum.total.count, sep="")
  all.bins.both[is.nan(all.bins.both$perc.in.both.sum), "heatmap.text.sum"] = ""
  
  all.bins.both$heatmap.text.mean = paste(round(all.bins.both$perc.in.both.mean, 1), "\n", round(all.bins.both$mean.both.count, 1), " / ", round(all.bins.both$mean.total.count, 1), sep="")
  all.bins.both[is.nan(all.bins.both$perc.in.both.mean), "heatmap.text.mean"] = ""
  
  #ordering of axis
  all.bins.both$Reads = factor(x = all.bins.both$Reads, levels = c("0-25", "25-100", "100-1000", "1000-10000", ">10000", "Sum"), ordered = T)
  all.bins.both$Freq = factor(x = all.bins.both$Freq, levels = c("0-0.01", "0.01-0.1", "0.1-1", "1-5", ">5", "Sum"), ordered = T)
  
  print(ggplot(all.bins.both, aes(Reads, Freq)) + 
          geom_tile(aes(fill=perc.in.both.sum), colour = NA) + 
          geom_text(aes(label=heatmap.text.sum), size=6) +
          scale_fill_gradient(low = "white",high = "steelblue", na.value="transparent", limits=c(0,100),name="Percentage\nIn Both") +
          theme(
            text=element_text(size=30, family="Arial"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.text=element_text(size=17),
            legend.title=element_text(size=20)
          ) + xlab("Absolute read count at diagnosis") + ylab("Clone frequency at diagnosis (%)") +
          ggtitle(paste("Sum of", rearrangement_type, "patients")))
  
  #ggsave(file.path(output_dir, paste("sum_of_patients_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".pdf", sep="")), width = 14, height = 12)
  ggsave(file.path(output_dir, paste("sum_of_patients_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".tiff", sep="")), device = "tiff", width = 14, height = 12)
  
  write.table(x = all.bins.both[,c("Reads", "Freq","sum.count", "sum.both.count", "sum.count.log")], file = file.path(output_dir, paste("sum_of_patients_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".txt", sep="")), quote = F, row.names = F, col.names = T, sep = "\t")
  
  print(ggplot(all.bins.both, aes(Reads, Freq)) + 
          geom_tile(aes(fill=perc.in.both.mean), colour = NA) + 
          geom_text(aes(label=heatmap.text.mean), size=6) +
          scale_fill_gradient(low = "white",high = "steelblue", na.value="transparent", limits=c(0,100),name="Percentage\nIn Both") +
          theme(
            text=element_text(size=30, family="Arial"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.text=element_text(size=17),
            legend.title=element_text(size=20)
          ) + xlab("Absolute read count at diagnosis") + ylab("Clone frequency at diagnosis (%)") +
          ggtitle(paste("Mean of", rearrangement_type, "patients")))
  
  #ggsave(file.path(output_dir, paste("mean_of_patients_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".pdf", sep="")), width = 14, height = 12)
  ggsave(file.path(output_dir, paste("mean_of_patients_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".tiff", sep="")), device = "tiff", width = 14, height = 12)
  
  write.table(x = all.bins.both[,c("Reads", "Freq","mean.count", "mean.both.count", "mean.count.log")], file = file.path(output_dir, paste("mean_of_patients_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".txt", sep="")), quote = F, row.names = F, col.names = T, sep = "\t")
  
  sum_heatmaps_data[[rearrangement_type]] <<- all.bins.both
  if(!is.null(all_sample_data) && nrow(all_sample_data) > 0){
    # dx freq vs dx reads heatmap
    all_sample_data$read_bin = cut(
      all_sample_data$normalized_read_count, 
      breaks=c(read_thresholds, Inf), 
      labels=read_threshold_labels
    )
    all_sample_data$freq_bin = cut(
      all_sample_data$Frequency, 
      breaks=c(freq_thresholds, Inf), 
      labels=freq_threshold_labels
    )
    
    read.vs.freq.count = data.frame(table(all_sample_data$read_bin, all_sample_data$freq_bin, all_sample_data$type))
    names(read.vs.freq.count) = c("Reads", "Frequency", "type", "count")
    
    read.vs.freq.count = dcast(read.vs.freq.count, Reads + Frequency ~ type, value.var="count")
    
    read.vs.freq.count$perc.diag = round(read.vs.freq.count$Diagnose / max(read.vs.freq.count$Diagnose) * 100, 2)
    
    for(freq_threshold_label in freq_threshold_labels){ #create the sums for the heatmap
      tmp = data.frame(
        Reads = "Sum",
        Frequency=freq_threshold_label
      )
      tmp$Diagnose = sum(read.vs.freq.count[read.vs.freq.count$Frequency == freq_threshold_label,"Diagnose"])
      tmp$Recidief = sum(read.vs.freq.count[read.vs.freq.count$Frequency == freq_threshold_label,"Recidief"])
      tmp$perc.diag = NA
      read.vs.freq.count = rbind(read.vs.freq.count, tmp)
    }
    for(read_threshold_label in read_threshold_labels){
      tmp = data.frame(
        Reads = read_threshold_label,
        Frequency="Sum"
      )
      tmp$Diagnose = sum(read.vs.freq.count[read.vs.freq.count$Reads == read_threshold_label,"Diagnose"])
      tmp$Recidief = sum(read.vs.freq.count[read.vs.freq.count$Reads == read_threshold_label,"Recidief"])
      tmp$perc.diag = NA
      read.vs.freq.count = rbind(read.vs.freq.count, tmp)
    }
    sum.sum = data.frame(
      Reads="Sum", 
      Frequency="Sum",
      Diagnose=sum(read.vs.freq.count[read.vs.freq.count$Reads == "Sum", "Diagnose"]),
      Recidief=sum(read.vs.freq.count[read.vs.freq.count$Frequency == "Sum", "Recidief"]),
      perc.diag=NA
    )
    read.vs.freq.count = rbind(read.vs.freq.count, sum.sum)
    
    read.vs.freq.count$Reads = factor(x = read.vs.freq.count$Reads, levels = c("0-25", "25-100", "100-1000", "1000-10000", ">10000", "Sum"), ordered = T)
    read.vs.freq.count$Frequency = factor(x = read.vs.freq.count$Frequency, levels = c("0-0.01", "0.01-0.1", "0.1-1", "1-5", ">5", "Sum"), ordered = T)
    
    print(ggplot(read.vs.freq.count, aes(Reads, Frequency)) + 
            geom_tile(aes(fill=perc.diag), colour = NA) + 
            geom_text(aes(label=Diagnose), size=6) +
            scale_fill_gradient(low = "white",high = "steelblue", na.value="transparent", limits=c(0,100),name="Percentage") +
            theme(
              text=element_text(size=30, family="Arial"),
              axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            ) +
            ggtitle(paste("Sum of", rearrangement_type, "Dx")))
    
    #ggsave(file.path(output_dir, paste("sum_of_dx_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".pdf", sep="")), width = 14, height = 12)
    ggsave(file.path(output_dir, paste("sum_of_dx_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".tiff", sep="")), device = "tiff", width = 14, height = 12)
    
    write.table(x = read.vs.freq.count, file = file.path(output_dir, paste("sum_of_dx_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".txt", sep="")), quote = F, row.names = F, col.names = T, sep = "\t")
    
  }
  #all_patient_data frequency scatter plot for Dx/R with linked dots
  if(!is.null(all_patient_data)){
    all_patient_data_fltr = all_patient_data
    #all_patient_data_fltr = all_patient_data[
    #  pmax(all_patient_data$Frequency.Dx, all_patient_data$Frequency.R) > 5 & pmax(all_patient_data$Normalized_Read_Count.Dx, all_patient_data$Normalized_Read_Count.R) > 10000,] # 
    if(nrow(all_patient_data_fltr) > 0){  
      all_patient_data_fltr$link = 1:nrow(all_patient_data_fltr)
      
      all_patient_data_fltr_D = all_patient_data_fltr[,c("Normalized_Read_Count.Dx", "Frequency.Dx", "link")]
      all_patient_data_fltr_D$type = "Diagnose Stable"
      all_patient_data_fltr_R = all_patient_data_fltr[,c("Normalized_Read_Count.R", "Frequency.R", "link")]
      all_patient_data_fltr_R$type = "Recidief Stable"
      
      names(all_patient_data_fltr_D) = c("Normalized_Read_Count", "Frequency", "link", "type")
      names(all_patient_data_fltr_R) = c("Normalized_Read_Count", "Frequency", "link", "type")
      
      all_patient_data_rbind = rbind(all_patient_data_fltr_D, all_patient_data_fltr_R)
      
      #all_sample_data = all_sample_data[all_sample_data$Frequency > 5 & all_sample_data$Normalized_Read_count > 10000,] # 
      if(!is.null(all_sample_data)){
        if(nrow(all_sample_data) > 0){
          all_sample_data$link = (max(all_patient_data_fltr$link) + 1):(max(all_patient_data_fltr$link) + nrow(all_sample_data))
          all_sample_data = all_sample_data[,c("Normalized_Read_Count", "Frequency", "link", "type")]
          names(all_sample_data) = c("Normalized_Read_Count", "Frequency", "link", "type")
          all_patient_data_rbind = rbind(all_patient_data_rbind, all_sample_data)
        }
      }
      
      type.count = data.frame(table(all_patient_data_rbind$type))
      names(type.count) = c("type", "Freq")
      type.count$type.label = paste(type.count$type, " (", type.count$Freq, ")", sep="")
      
      all_patient_data_rbind$type = factor(x = all_patient_data_rbind$type, levels = type.count$type, labels=type.count$type.label)
      
      write.table(x = all_patient_data_rbind, file = file.path(output_dir, paste("linked_scatter_dx_r_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".txt", sep="")), quote = F, row.names = F, col.names = T, sep = "\t")
      
      all_patient_data_rbind$type = factor(
        x = all_patient_data_rbind$type,
        levels = (c(
          paste("Diagnose (", type.count[type.count$type == "Diagnose", "Freq"], ")", sep=""),
          paste("Diagnose Stable (", type.count[type.count$type == "Diagnose Stable", "Freq"], ")", sep=""),
          paste("Recidief Stable (", type.count[type.count$type == "Recidief Stable", "Freq"], ")", sep=""),
          paste("Recidief (", type.count[type.count$type == "Recidief", "Freq"], ")", sep="")
        )),
        labels = (c(
          paste("Diagnosis (", type.count[type.count$type == "Diagnose", "Freq"], ")", sep=""),
          paste("Diagnosis Stable (", type.count[type.count$type == "Diagnose Stable", "Freq"], ")", sep=""),
          paste("Relapse Stable (", type.count[type.count$type == "Recidief Stable", "Freq"], ")", sep=""),
          paste("Relapse (", type.count[type.count$type == "Recidief", "Freq"], ")", sep="")
        )),
        ordered = T
      )
      
      p = ggplot(all_patient_data_rbind, aes(type, Frequency, group=link)) + 
        geom_line() + theme(
          text=element_text(size=30,  family="Arial"),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size=20)
        ) +
        scale_x_discrete(all_patient_data_rbind$type, drop=FALSE)
      p = p + geom_point(aes(colour=type), position=position_dodge(width = 0)) + ggtitle(rearrangement_type) + xlab("Type") + ylab("Frequency")
      
      print(p)
      #ggsave(file.path(output_dir, paste("linked_frequency_scatter_dx_r_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".pdf", sep="")), width = 14, height = 12)
      ggsave(file.path(output_dir, paste("linked_frequency_scatter_dx_r_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".tiff", sep="")), device = "tiff", width = 14, height = 12)
      write.table(x = all_patient_data_rbind, file = file.path(output_dir, paste("linked_scatter_dx_r_diag_", diag_freq, "_recid_", recid_freq, "_", rearrangement_type, ".txt", sep="")), quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }
}

outside_tmp = NULL

sum_heatmaps_data = list()
windowsFonts(Arial=windowsFont("TT Arial"))

post_prisca_patient_analysis(
  unique(patient_samples$Patient), 
  "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original", 
  "Total", 
  "D:/wd/prisca/heatmaps_vincent/tmp",
  diag_freq = 0,
  recid_freq = 0
)

patient_samples = unique(input_data[,c("Patient", "Sample")])
patient_samples = patient_samples[order(patient_samples$Patient),]
patient_ids = unique(patient_samples$Patient)
patient_ids = patient_ids[grepl(paste(needed_patients, collapse="|"), patient_ids)]




for(rearrangement_type in rearrangement_types){ #everything
  for(diag_freq in c(0, 0.01, 0.1, 1, 5)){
    for(recid_freq in c(0, 1, 5)){
      print(paste("Running", rearrangement_type, diag_freq, recid_freq))
      post_prisca_patient_analysis(
        patient_ids, 
        "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original", 
        rearrangement_type, 
        "D:/wd/prisca/heatmaps_vincent/2018-4-10_result",
        diag_freq = diag_freq,
        recid_freq = recid_freq
      )
    }
  }
}

for(rearrangement_type in rearrangement_types){ #just for the relapse
  for(diag_freq in c(0.01, 0.1, 1, 5)){
    for(recid_freq in c(1, 5)){
      print(paste("Running", rearrangement_type, diag_freq, recid_freq))
      post_prisca_patient_analysis(
        unique(patient_samples$Patient), 
        "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original", 
        rearrangement_type, 
        "D:/wd/prisca/heatmaps_vincent/2018-1-9_relapse_result",
        diag_freq = diag_freq,
        recid_freq = recid_freq,
        only_both = T
      )
    }
  }
}

create_stability_wide_from_short = function(input_file, output_file){
  input = read.table(file=input_file, header = T, sep = "\t", quote = "")
  input$patient.type = paste(input$patient, input$type, sep="_")
  
  output = dcast(input, locus + diag_freq + recid_freq ~ patient.type, value.var = "count")
  
  output$only_diag_sum = rowSums(output[,names(output)[grepl("diag$", names(output))]])
  output$only_recid_sum = rowSums(output[,names(output)[grepl("recid$", names(output))]])
  output$only_both_sum = rowSums(output[,names(output)[grepl("both$", names(output))]])
  
  output$diag_sum = rowSums(output[,names(output)[grepl("diag$|both$", names(output))]])
  output$recid_sum = rowSums(output[,names(output)[grepl("recid$|both$", names(output))]])
  
  output$only_diag_mean = rowMeans(output[,names(output)[grepl("diag$", names(output))]])
  output$only_recid_mean = rowMeans(output[,names(output)[grepl("recid$", names(output))]])
  output$only_both_mean = rowMeans(output[,names(output)[grepl("both$", names(output))]])
  
  output$diag_positive = rowSums(output[,names(output)[grepl("diag$", names(output))]] > 0)
  output$recid_positive = rowSums(output[,names(output)[grepl("recid$", names(output))]] > 0)
  output$both_positve = rowSums(output[,names(output)[grepl("both$", names(output))]] > 0)
  
  output$diag_positve_perc = output$diag_positive / length(unique(input$patient)) * 100
  output$recid_positve_perc = output$diag_positive / length(unique(input$patient)) * 100
  output$both_positve_perc = output$both_positve / length(unique(input$patient)) * 100
  
  output = output[,c(1:3, (ncol(output) - 13):ncol(output),which(grepl("^VanDongen", names(output))))]
  
  write.table(x = output, file = output_file, quote = F, row.names = F, col.names = T, sep = "\t")
}

stability_long = file.path("D:/wd/prisca/heatmaps_vincent/2018-2-27_result/", "stability_per_locus_patient_frequency_long.txt")
stability_wide = file.path("D:/wd/prisca/heatmaps_vincent/2018-2-27_result/", "stability_per_locus_patient_frequency.txt")

#input_file = stability_long
#output_file = stability_wide

create_stability_wide_from_short(stability_long, stability_wide)


#generate the stable/total sequence counts per loci/patient
all.bins[,c("Reads", "Freq","mean.count", "mean.both.count", "mean.count.log")]

dat = sum_heatmaps_data[["Total"]][1,]
dat$locus = "Total"

for(rearrangement_type in rearrangement_types[2:length(rearrangement_types)]){
  tmp = sum_heatmaps_data[[rearrangement_type]][1,]
  tmp$locus = rearrangement_type
  dat = rbind(dat, tmp)
}

dat2 = data.frame(locus=dat$locus, stable=dat$sum.both.count, total=dat$sum.total.count)

for(patient in unique(patient_samples$Patient)){
  patient.cols = dat[,grepl(patient, names(dat))]
  dat2[,patient] = paste(patient.cols[,1], "/", patient.cols[,4], sep="")
}

write.table(x = dat2, file = file.path(output_dir, "stable_by_patient_locus.txt"), quote = F, row.names = F, col.names = T, sep = "\t")


#pdf of single plots
library(ggplot2)
library(data.table)

input_file = "D:/wd/prisca/heatmaps_vincent/single_samples_paper.txt"

input_data = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T)


if(!("Total_Read_Count" %in% names(input_data))){
  input_data$Total_Read_Count = 0
}
input_data = input_data[,c("Patient",  "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")] 

input_data$dsPerM = 0
input_data = input_data[!is.na(input_data$Patient),]
input_data$Related_to_leukemia_clone = F

input_data$V_Segment_Major_Gene = as.factor(as.character(lapply(strsplit(as.character(input_data$V_Segment_Major_Gene), "; "), "[[", 1)))
input_data$J_Segment_Major_Gene = as.factor(as.character(lapply(strsplit(as.character(input_data$J_Segment_Major_Gene), "; "), "[[", 1)))

input_data$Frequency = ((10^input_data$Log10_Frequency)*100)

input_data$locus_V = substring(input_data$V_Segment_Major_Gene, 0, 4)
input_data$locus_J = substring(input_data$J_Segment_Major_Gene, 0, 4)

min_cell_count = data.frame(data.table(input_data)[, list(min_cell_count=min(.SD$Cell_Count)), by=c("Patient", "locus_V", "locus_J")])

input_data$min_cell_paste = paste(input_data$Patient, input_data$locus_V, input_data$locus_J)
min_cell_count$min_cell_paste = paste(min_cell_count$Patient, min_cell_count$locus_V, min_cell_count$locus_J)

min_cell_count = min_cell_count[,c("min_cell_paste", "min_cell_count")]
input_data = merge(input_data, min_cell_count, by="min_cell_paste")

input_data$normalized_read_count = round(input_data$Clone_Molecule_Count_From_Spikes / input_data$Cell_Count * input_data$min_cell_count / 2, digits=2)

bak = input_data
#input_data = bak

input_data = input_data[grepl("IGHV", input_data$locus_V),]

input_data[input_data$normalized_read_count < 1,"normalized_read_count"] = 1
input_data$x = ""

scales = c(1,10,30,100)

p = ggplot(input_data, aes(x, normalized_read_count)) + 
  scale_y_log10(breaks=scales,labels=scales, limits=c(1,100)) + 
  geom_jitter() + facet_grid(Patient ~ .)


#making new dataset with just HD_050515_unsorted, 11444 and 12553
dongen_file = "D:/wd/prisca/heatmaps_vincent/VanDongen_clones_evol_may2015_complete_2015_8_20.txt"
dongen = read.table(file=dongen_file, header = T, sep = "\t", quote = "",fill = T)

singles2015_file = "D:/wd/prisca/heatmaps_vincent/2015-11-20-singles.txt" 
#singles2015_file = "D:/wd/prisca/heatmaps_vincent/2015-11-20-sorted.txt" 
singles2015 = read.table(file=singles2015_file, header = T, sep = "\t", quote = "",fill = T)
table(singles2015$Patient)

dongen = dongen[grepl("11440|12553", dongen$Patient),]
singles2015 = singles2015[singles2015$Patient == "Patient_HD_050515_unsorted",]

names(dongen)[!(names(dongen) %in% names(singles2015))]
names(singles2015)[!(names(singles2015) %in% names(dongen))]
unique(dongen$Patient)
unique(singles2015$Patient)

singles2015$Related_to_leukemia_clone = F

combined = rbind(dongen, singles2015)
unique(combined$Patient)

combined_file = "D:/wd/prisca/heatmaps_vincent/single_samples_paper.txt"
write.table(x = combined, file = combined_file, quote = F, row.names = F, col.names = T, sep = "\t")




#making new dataset with just the 16278_Right_Trio_26402_Right_Trio_26759_Right_Trio triple
dongen_file = "D:/wd/prisca/heatmaps_vincent/VanDongen_clones_evol_may2015_complete_2015_8_20.txt"
dongen = read.table(file=dongen_file, header = T, sep = "\t", quote = "",fill = T)

#bak = dongen
#dongen = bak

dongen = dongen[grepl("16278_Right|26402_Right|26759_Right", dongen$Patient),]
dongen$Patient = as.character(dongen$Patient)
dongen$Sample = as.character(dongen$Sample)
dongen$uniqueID = ""

dongen[grepl("16278_Right", dongen$Sample),"uniqueID"] = "16278_26402_26759_Right"
dongen[grepl("26402_Right", dongen$Sample),"uniqueID"] = "16278_26402_26759_Right"
dongen[grepl("26759_Right", dongen$Sample),"uniqueID"] = "16278_26402_26759_Right" 

dongen = dongen[dongen$uniqueID == "16278_26402_26759_Right",]
dongen$Patient = "16278_26402_26759_Right"

complete_triple = "D:/wd/prisca/heatmaps_vincent/16278_26402_26759_Right/complete_16278_26402_26759_Right.txt"
write.table(x = dongen, file = complete_triple, quote = F, row.names = F, col.names = T, sep = "\t")



#make scatter from the loci_scatter.txt
triple_file = "D:/wd/prisca/heatmaps_vincent/16278_26402_26759_Right/16278_26402_26759_Right_16278_IGH-Vh-Jh_scatter_0.txt"
triple = read.table(file=triple_file, header = T, sep = "\t", quote = "",fill = T, dec = ",")
triple[triple$normalized_read_count < 1,"normalized_read_count"] = 1
triple$type = factor(x = as.numeric(triple$type), levels = c("In one", "In two", "In three"))
scales=c()

p = ggplot(triple, aes(type, normalized_read_count)) + 
  scale_y_log10(breaks=scales,labels=scales, limits=c(1, 1e6)) +
  geom_jitter()




#collect all the total >0 files into a single data.frame
library(data.table)

complete_data = NULL #all the sequence frequency/read counts with an extra columns for patient, sample and diag/recid/both

for(rearrangement_type in rearrangement_types){
  for(patient in patients){
    this_patient_samples = patient_samples[patient_samples$Patient == patient,]
    sample1 = as.character(this_patient_samples[1, "Sample"])
    sample2 = as.character(this_patient_samples[2, "Sample"])
    
    if(grepl("_Dx", sample2)){ #switch samples so sample 1 is diagnosis
      print(paste("Switching", sample1, "with", sample2))
      tmp = sample1
      sample1 = sample2
      sample2 = tmp
    }
    
    print(paste(rearrangement_type, patient, sample1, sample2))
    
    #find the file with all the sequences for this patient and 
    patient_data_file = file.path(input_dir, paste(sample1, sample2, rearrangement_type, "0.txt", sep="_"))
    if(!file.exists(patient_data_file)){
      patient_data_file = file.path(input_dir, paste(sample2, sample1, rearrangement_type, "0.txt", sep="_"))
    }
    
    #create an empty patient file as a fallback
    patient_data = data.frame(Clone.Sequence="")
    patient_data[,paste(rep(c("Proximal.segment", "Distal.segment", "Related_to_leukemia_clone"), each=2), c(sample1, sample2), sep=".")] = ""
    patient_data[,paste(rep(c("Normalized_Read_Count", "Frequency"), each=2), c(sample1, sample2), sep=".")] = 0
    patient_data = patient_data[NULL,]
    
    #if the patient file exists, read it
    if(file.exists(patient_data_file)){
      patient_data = read.table(file = patient_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
      
      if(nrow(patient_data) > 0){
        patient_data = data.frame(
          patient=patient,
          Proximal.segment=patient_data[,paste("Proximal.segment", sample1, sep=".")],
          Distal.segment=patient_data[,paste("Distal.segment", sample1, sep=".")],
          Normalized_Read_count_diag=patient_data[,paste("Normalized_Read_Count", sample1, sep=".")],
          Normalized_Read_count_rec=patient_data[,paste("Normalized_Read_Count", sample2, sep=".")],
          Frequency_diag=patient_data[,paste("Frequency", sample1, sep=".")],
          Frequency_rec=patient_data[,paste("Frequency", sample2, sep=".")],
          type="both",
          locus=rearrangement_type
        )
        
        if(is.null(complete_data)){
          complete_data = patient_data
        } else {
          complete_data = rbind(complete_data, patient_data)
        }
      }
    }
    
    sample1_data_file = file.path(input_dir, paste(sample1, rearrangement_type, "0.txt", sep="_"))
    sample1_data = sample1_data[NULL,]
    if(file.exists(sample1_data_file)){
      sample1_data = read.table(file = sample1_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
      
      if(nrow(patient_data) > 0){
        sample1_data = data.frame(
          patient=patient,
          Proximal.segment=sample1_data$Proximal.segment,
          Distal.segment=sample1_data$Distal.segment,
          Normalized_Read_count_diag=sample1_data$normalized_read_count,
          Normalized_Read_count_rec=NA,
          Frequency_diag=sample1_data$Frequency,
          Frequency_rec=NA,
          type="diag",
          locus=rearrangement_type
        )
        
        if(is.null(complete_data)){
          complete_data = sample1_data
        } else {
          complete_data = rbind(complete_data, sample1_data)
        }
      }
    }
    
    sample2_data_file = file.path(input_dir, paste(sample2, rearrangement_type, "0.txt", sep="_"))
    sample2_data = sample2_data[NULL,]
    if(file.exists(sample2_data_file)){
      sample2_data = read.table(file = sample2_data_file, header = T, sep = "\t", quote = "", fill = T, stringsAsFactors = F, dec=",")
      
      if(nrow(patient_data) > 0){
        sample2_data = data.frame(
          patient=patient,
          Proximal.segment=sample2_data$Proximal.segment,
          Distal.segment=sample2_data$Distal.segment,
          Normalized_Read_count_diag=NA,
          Normalized_Read_count_rec=sample2_data$normalized_read_count,
          Frequency_diag=NA,
          Frequency_rec=sample2_data$Frequency,
          type="recid",
          locus=rearrangement_type
        )
        
        if(is.null(complete_data)){
          complete_data = sample2_data
        } else {
          complete_data = rbind(complete_data, sample2_data)
        }
      }
    }
  }
}
write.table(x = complete_data, file = "D:/wd/prisca/heatmaps_vincent/dx_rec_from_orig_complete_data.txt", quote = F, sep = "\t", row.names = F, col.names = T)
complete_data_fltr = complete_data[
  pmax(complete_data$Normalized_Read_count_diag, complete_data$Normalized_Read_count_rec) > 100 &&
  pmax(complete_data$Frequency_diag, complete_data$Frequency_rec, na.rm = T) > 0.01,
]
complete_data = read.table(file = "D:/wd/prisca/heatmaps_vincent/dx_rec_from_orig_complete_data.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F)



#mean number of clones per patient and per frequency (table1 in Tables.xlsx)
result = data.frame(locus=rearrangement_types)

complete_data = data.table(complete_data)

count_by_patient_locus = complete_data[type %in% c("both", "diag"), list(count=.N), by=list(patient, locus)]
mean_by_locus = count_by_patient_locus[,list(Mean=mean(count), Min=min(count), Max=max(count)), by=list(locus)]

for(freq_threshold in freq_thresholds_2){
  count_by_patient = complete_data[type %in% c("both", "diag") & Frequency > freq_threshold, list(count=.N), by=list(patient, locus)]
  count_by_locus = count_by_patient[, list(count=mean(count)), by=list(locus)]
  names(count_by_locus) = c("locus", paste(">", freq_threshold, sep=""))
  mean_by_locus = merge(mean_by_locus, count_by_locus, by="locus", all=T)
}

write.table(mean_by_locus, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)




#stable/stable+diag (table2)

overall_count_by_locus = merge(
  complete_data[type %in% c("both"),list(stable=.N),by=list(locus)],
  complete_data[type %in% c("both", "diag"),list(total=.N),by=list(locus)]
)

write.table(overall_count_by_locus, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)









#create a table with per patient/locus a total number of clones (total, >5%, >1%, %0.1%, >0.01%)
library(reshape2)

complete_data = read.table(file = "D:/wd/prisca/heatmaps_vincent/dx_rec_from_orig_complete_data.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F)

#filter complete_data for the thresholds from the paper
complete_data$Normalized_Read_Count = pmax(
  complete_data$Normalized_Read_count_diag,
  complete_data$Normalized_Read_count_rec,
  na.rm = T
)

complete_data$V_Segment_Major_Gene = complete_data$Proximal.segment
complete_data$J_Segment_Major_Gene = complete_data$Distal.segment

complete_data = filter_df_for_paper_threshold(complete_data)

complete_data = complete_data[,names(complete_data)[!(names(complete_data) %in% c("Normalized_Read_Count", "V_Segment_Major_Gene", "J_Segment_Major_Gene"))]]

complete_data$Frequency = pmax(
  complete_data$Frequency_diag,
  complete_data$Frequency_rec,
  na.rm = T
)

all_clones_per_patient_locus = clones_per_patient_locus[NULL,]

for(frq in c(0, 0.01, 0.1, 1, 5)){
  fltr = complete_data$Frequency >= frq
  clones_per_patient_locus = data.frame(table(complete_data[fltr,]$patient, complete_data[fltr,]$type, complete_data[fltr,]$locus))
  names(clones_per_patient_locus) = c("Patient", "type", "locus", "count")
  clones_per_patient_locus$patient_type = paste(clones_per_patient_locus$Patient, clones_per_patient_locus$type, sep="_")
  clones_per_patient_locus$Frequency = frq
  
  clones_per_patient_locus = dcast(clones_per_patient_locus, locus + Frequency ~ patient_type, value.var="count")
  all_clones_per_patient_locus = rbind(all_clones_per_patient_locus, clones_per_patient_locus)
}

all_clones_per_patient_locus$diag_mean = rowMeans(all_clones_per_patient_locus[,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_diag", names(all_clones_per_patient_locus))]])
all_clones_per_patient_locus$recid_mean = rowMeans(all_clones_per_patient_locus[,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_recid", names(all_clones_per_patient_locus))]])
all_clones_per_patient_locus$both_mean = rowMeans(all_clones_per_patient_locus[,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_both", names(all_clones_per_patient_locus))]])

all_clones_per_patient_locus$diag_positive = rowSums(all_clones_per_patient_locus[,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_diag", names(all_clones_per_patient_locus))]] > 0)
all_clones_per_patient_locus$recid_positive = rowSums(all_clones_per_patient_locus[,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_recid", names(all_clones_per_patient_locus))]] > 0)
all_clones_per_patient_locus$both_positive = rowSums(all_clones_per_patient_locus[,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_both", names(all_clones_per_patient_locus))]] > 0)

length(all_clones_per_patient_locus[1,names(all_clones_per_patient_locus)[grepl("^VanDongen_cALL_[0-9]+_diag", names(all_clones_per_patient_locus))]] > 0)

all_clones_per_patient_locus$diag_positive_perc = all_clones_per_patient_locus$diag_positive / 43 * 100
all_clones_per_patient_locus$recid_positive_perc = all_clones_per_patient_locus$recid_positive / 43 * 100
all_clones_per_patient_locus$both_positive_perc = all_clones_per_patient_locus$both_positive / 43 * 100

all_clones_per_patient_locus = all_clones_per_patient_locus[,c(
  "locus", "Frequency", 
  "diag_mean", "recid_mean", "both_mean", 
  "diag_positive", "recid_positive", "both_positive",
  "diag_positive_perc", "recid_positive_perc", "both_positive_perc",
  names(all_clones_per_patient_locus)[grepl("VanDongen", names(all_clones_per_patient_locus))]
)]

write.table(all_clones_per_patient_locus, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T) # stability_per_locus_patient_frequency


#compare the dx/r file with the 2 seperate files
V_Segments = c(".*", "IGHV", "IGHD", "IGKV", "IGKV", "IgKINTR", "TRGV", "TRDV", "TRDD" , "TRBV")
J_Segments = c(".*", ".*", ".*", "IGKJ", "KDE", ".*", ".*", ".*", ".*", ".*")
Titles = c("Total", "IGH-Vh-Jh", "IGH-Dh-Jh", "Vk-Jk", "Vk-Kde" , "Intron-Kde", "TCRG", "TCRD-Vd-Dd", "TCRD-Dd-Dd", "TCRB-Vb-Jb")

targets = data.frame(V=V_Segments, J=J_Segments, label=Titles)

targets$count_may2015 = 0
targets$count_nov2015 = 0
targets$count_dx.r = 0

for(i in 1:nrow(targets)){
  v = targets[i,"V"]
  j = targets[i,"J"]
  targets[i,"count_may2015"] = sum(grepl(v, input_data1$V_Segment_Major_Gene) & grepl(j, input_data1$J_Segment_Major_Gene))
  targets[i,"count_nov2015"] = sum(grepl(v, input_data2$V_Segment_Major_Gene) & grepl(j, input_data2$J_Segment_Major_Gene))
  targets[i,"count_dx.r"] = sum(grepl(v, input_data$V_Segment_Major_Gene) & grepl(j, input_data$J_Segment_Major_Gene))
}

write.table(targets, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

may2015.subset = input_data2[input_data2$Patient %in% input_data$Patient,]
all(may2015.subset$Clone_Sequence %in% input_data$Clone_Sequence)


#make new dx/r file from the original data
input_file1 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/VanDongen_clones_evol.txt"
input_file2 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/VanDongen_clones_evol_may2015.txt"
input_file3 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/VanDongen_clones_evol_may2015_2.txt"
input_file4 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/VanDongen_clones_evol_oct2014.txt"

input_file5 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/Ziv Rosen - vanDongen_report_6_23_2015_May_June_Only.csv"
input_file6 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/Ziv Rosen - vanDongen_report_7_16_2015_July_Only.csv"
input_file7 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/Ziv Rosen - vanDongen_report_7_31_2015.csv"
input_file8 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/Ziv Rosen - vanDongen_report_8_28_2015.csv"
input_file9 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/Ziv Rosen - vanDongen_report_9_9_2015.csv"
input_file10 = "D:/wd/prisca/heatmaps_vincent/dx_r_from_original/VanDongen_clones_evol_may2015_complete_2015_8_20.txt"

input_data1 = read.table(file=input_file1, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)
input_data2 = read.table(file=input_file2, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)
input_data3 = read.table(file=input_file3, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)
input_data4 = read.table(file=input_file4, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

input_data5 = read.table(file=input_file5, header = T, sep = ",", quote = "",fill = T, stringsAsFactors = F)
input_data6 = read.table(file=input_file6, header = T, sep = ",", quote = "",fill = T, stringsAsFactors = F)
input_data7 = read.table(file=input_file7, header = T, sep = ",", quote = "",fill = T, stringsAsFactors = F)
input_data8 = read.table(file=input_file8, header = T, sep = ",", quote = "",fill = T, stringsAsFactors = F)
input_data9 = read.table(file=input_file9, header = T, sep = ",", quote = "",fill = T, stringsAsFactors = F)
input_data10 = read.table(file=input_file10, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

new_col_names = c("Patient",  "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")

input_data1 = input_data1[,c("Patient", "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")]
names(input_data1) = new_col_names

input_data2 = input_data2[,c("patient", "Receptor", "sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data2) = new_col_names

input_data3$Receptor = ""
input_data3 = input_data3[,c("patient", "Receptor", "sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data3) = new_col_names

input_data4$Receptor = ""
input_data4 = input_data4[,c("Patient", "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")]
names(input_data4) = new_col_names

input_data5 = input_data5[,c("patient", "Receptor", "Sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data5) = new_col_names

input_data6 = input_data6[,c("patient", "Receptor", "Sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data6) = new_col_names

input_data7 = input_data7[,c("patient", "Receptor", "Sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data7) = new_col_names

input_data8 = input_data8[,c("patient", "Receptor", "Sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data8) = new_col_names

input_data9 = input_data9[,c("patient", "Receptor", "Sample", "nCells", "nMol", "Log10_Frequency", "nRead", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "cSeq")]
names(input_data9) = new_col_names

input_data10 = input_data10[,c("Patient", "Receptor", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene", "CDR3_Sense_Sequence", "Clone_Sequence")]
names(input_data10) = new_col_names

#test = rbind(head(input_data1), head(input_data2), head(input_data3), head(input_data4))
#test = rbind(head(input_data5), head(input_data6), head(input_data7), head(input_data8), head(input_data9))

#input_data = rbind(input_data1, input_data2, input_data3, input_data4, input_data5, input_data6, input_data7, input_data8, input_data9)

#input_data = input_data[grepl(paste(needed_patients, collapse="|"), input_data$Patient),]

#write.table(input_data, file = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#length(unique(input_data$Patient))

#compare the old file with the new file made from the original files

old_table = data.frame(table(input_data$Patient))
names(old_table) = c("Patient", "old_count")
new_table = data.frame(table(input_data_orig$Patient))
names(new_table) = c("Patient", "new_count")
mrg = merge(old_table, new_table, by="Patient")

input_files = c(input_file1, input_file2, input_file3, input_file4, input_file5, input_file6, input_file7, input_file8, input_file9, input_file10)
input_datas = list(input_data1, input_data2, input_data3, input_data4, input_data5, input_data6, input_data7, input_data8, input_data9, input_data10)

mrg = NULL

for(i in 1:length(input_files)){
  current_input_file = input_files[i]
  current_input_data = input_datas[[i]]
  
  tmp = data.frame(table(current_input_data$Patient))
  names(tmp) = c("Patient", basename(current_input_file))
  
  if(is.null(mrg)){
    mrg = tmp
  } else {
    mrg = merge(mrg, tmp, by="Patient", all=T)
  }
}

mrg = mrg[grepl(paste(needed_patients, collapse="|"), mrg$Patient),]
mrg[is.na(mrg)] = 0

write.table(mrg, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)


#compare a patient from both 'VanDongen_clones_evol_may2015.txt' (input_data2) and 'VanDongen_clones_evol_oct2014.txt' (input_data4) 

input_data2.patient = input_data2[input_data2$Patient == "VanDongen_cALL_15582",]
input_data4.patient = input_data4[input_data4$Patient == "VanDongen_cALL_15582",]

input_data2.patient = input_data2.patient[ordered(-input_data2.patient$Total_Read_Count),]
input_data4.patient = input_data4.patient[ordered(-input_data4.patient$Total_Read_Count),]

tmp = rbind(head(input_data2.patient), head(input_data4.patient))

no_clone_sequence_columns = c("Patient", "Sample", "Cell_Count", "Clone_Molecule_Count_From_Spikes", "Log10_Frequency", "Total_Read_Count", "J_Segment_Major_Gene", "V_Segment_Major_Gene")

all(input_data2.patient[,no_clone_sequence_columns] == input_data4.patient[,no_clone_sequence_columns])

#check if the clones in 'VanDongen_clones_evol.txt' (input_data1) are just the sequences from 'VanDongen_clones_evol_may2015.txt' (input_data2) twice

input_data1.patient = input_data1[input_data1$Patient == "VanDongen_cALL_15582",]
input_data2.patient = input_data2[input_data2$Patient == "VanDongen_cALL_15582",]

input_data1.patient = input_data1.patient[order(-input_data1.patient$Total_Read_Count),]
input_data2.patient = input_data2.patient[order(-input_data2.patient$Total_Read_Count),]

#make the new Dx/R file from:
#'VanDongen_clones_evol_may2015.txt' (input_data2)
#'VanDongen_clones_evol_may2015_2.txt' (input_data3)
#'Ziv Rosen - vanDongen_report_6_23_2015_May_June_Only.csv' (input_data5)
#'Ziv Rosen - vanDongen_report_7_16_2015_July_Only.csv' (input_data6)
#'Ziv Rosen - vanDongen_report_7_31_2015.csv' (input_data7)
#'Ziv Rosen - vanDongen_report_8_28_2015.csv' (input_data8)
#'Ziv Rosen - vanDongen_report_9_9_2015.csv' (input_data9)

input_data = rbind(input_data2, input_data3, input_data5, input_data6, input_data7, input_data8, input_data9)

input_data = input_data[grepl(paste(needed_patients, collapse="|"), input_data$Patient),]

#swap some Dx/R samples

test = input_data[grepl("5930|5456|5657|261195|130890|30390", input_data$Sample),]
df.t = data.frame(table(test$Patient, test$Sample))
df.t = df.t[df.t$Freq > 0,]
df.t = df.t[order(df.t$Var1),]

input_data[input_data$Patient == "VanDongen_cALL_5930" & input_data$Sample == "5930_Dx", "Sample"] = "5930_R"
input_data[input_data$Patient == "VanDongen_cALL_261195" & input_data$Sample == "261195_dx", "Patient"] = "VanDongen_cALL_5930"

input_data[input_data$Patient == "VanDongen_cALL_5456" & input_data$Sample == "5456_Dx", "Sample"] = "5456_R"
input_data[input_data$Patient == "VanDongen_cALL_130890" & input_data$Sample == "130890_dx", "Patient"] = "VanDongen_cALL_5456"

input_data[input_data$Patient == "VanDongen_cALL_5657" & input_data$Sample == "5657_Dx", "Sample"] = "5657_R"
input_data[input_data$Patient == "VanDongen_cALL_030390" & input_data$Sample == "30390_dx", "Patient"] = "VanDongen_cALL_5657"

test = input_data[grepl("5930|5456|5657|261195|130890|30390", input_data$Sample),]
df.t.new = data.frame(table(test$Patient, test$Sample))
df.t.new = df.t.new[df.t.new$Freq > 0,]
df.t.new = df.t.new[order(df.t.new$Var1),]

#check that everything transfered ok
mrg = merge(df.t, df.t.new, by="Var2", all=TRUE)


#fix other samples from the master file
old_new_samples = read.table(file="dx_r_from_original/new_sample_names.txt", header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)
old_new_samples$modified = 0

for(i in 1:nrow(old_new_samples)){
  old_sample = old_new_samples[i,"old"]
  new_sample = old_new_samples[i,"new"]
  
  fltr = input_data$Sample == old_sample
  print(paste(old_sample, "has", sum(fltr), "clones"))
  input_data[fltr,"Sample"] = new_sample
  old_new_samples[i,"modified"] = sum(fltr)
}


write.table(input_data, file = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

length(unique(input_data$Patient))


#make a frequency count of the patient/samples and compare it to the new 
input_data = read.table(file="D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original.txt", header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

patient_sample_tbl = data.frame(table(input_data$Patient, input_data$Sample))
names(patient_sample_tbl) = c("Patient", "Sample", "n_before")
patient_sample_tbl = patient_sample_tbl[patient_sample_tbl$n_before > 0,]

patient_sample_tbl_after = data.frame(table(input_data$Patient, input_data$Sample))
names(patient_sample_tbl_after) = c("Patient", "Sample", "n_after")
patient_sample_tbl_after = patient_sample_tbl_after[patient_sample_tbl_after$n_after > 0,]

patient_sample_tbl_mrg = merge(patient_sample_tbl, patient_sample_tbl_after, by=c("Patient", "Sample"), all=T)

write.table(patient_sample_tbl_mrg, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)





# compare 
#'VanDongen_clones_evol_may2015' (input_data2) 
#'VanDongen_clones_evol_may2015_2' (input_data3)
#'VanDongen_clones_evol_may2015_complete_2015_8_20' (input_file10)


nrow(input_data2) + nrow(input_data3)
nrow(input_data10)

input_data10 = input_data10[order(-input_data10$Total_Read_Count),]



#create some smaller patient files for testing
patient_clone_count = data.frame(table(input_data$Patient))
names(patient_clone_count) = c("Patient", "count")
patient_clone_count = patient_clone_count[order(patient_clone_count$count),]

for(i in 1:4){
  patient = patient_clone_count[i,"Patient"]
  input_data_sub = input_data[input_data$Patient == patient_clone_count[i,"Patient"],]
  write.table(input_data_sub, file = file.path("D:/wd/prisca/heatmaps_vincent/", paste(patient,".txt", sep="")), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
}






#create new dataset with all clones filtered on the minimum read threshold from the paper

minimum_read_threshold = data.frame(
  V_segments = as.character(c(".*", "IGHV", "IGHD", "IGKV", "IGKV", "IgKINTR", "TRGV", "TRDV", "TRDD" , "TRBV")),
  J_segments = as.character(c(".*", ".*", ".*", "IGKJ", "KDE", ".*", ".*", ".*", ".*", ".*")),
  Titles = as.character(c("Total", "IGH-Vh-Jh", "IGH-Dh-Jh", "Vk-Jk", "Vk-Kde" , "Intron-Kde", "TCRG", "TCRD-Vd-Dd", "TCRD-Dd-Dd", "TCRB-Vb-Jb")),
  read_thresholds = c(0, 30, 30, 70, 40, 170, 50, 20, 30, 30),
  clones_in_locus = 0,
  clones_above_threshold = 0,
  clones_removed = 0
)

input_file = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original.txt"
input_data = read.table(file=input_file, header = T, sep = "\t", quote = "",fill = T, stringsAsFactors = F)

for(i in 1:nrow(minimum_read_threshold)){
  V_segment = minimum_read_threshold[i, "V_segments"]
  J_segment = minimum_read_threshold[i, "J_segments"]
  read_threshold = minimum_read_threshold[i, "read_thresholds"]
  title = minimum_read_threshold[i, "Titles"]
  clones_in_locus = grepl(V_segment, input_data$V_Segment_Major_Gene) & grepl(J_segment, input_data$J_Segment_Major_Gene)
  fltr = input_data$Total_Read_Count >= read_threshold
  minimum_read_threshold[i, "clones_in_locus"] = sum(clones_in_locus)
  minimum_read_threshold[i, "clones_above_threshold"] = sum(fltr & clones_in_locus)
  minimum_read_threshold[i, "clones_removed"] = (sum(clones_in_locus) - sum(fltr & clones_in_locus))
  print(paste("Clones in locus", title, ":", sum(clones_in_locus), "Clones above threshold:", sum(fltr & clones_in_locus), "Removing:", (sum(clones_in_locus) - sum(fltr & clones_in_locus))))
  input_data = input_data[fltr | !clones_in_locus,]
}

write.table(input_data, file = "D:/wd/prisca/heatmaps_vincent/Dx_R_patients_2017_from_original_filter.txt", append = F, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(minimum_read_threshold, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

#check

for(i in 1:nrow(minimum_read_threshold)){
  V_segment = minimum_read_threshold[i, "V_segments"]
  J_segment = minimum_read_threshold[i, "J_segments"]
  read_threshold = minimum_read_threshold[i, "read_thresholds"]
  title = minimum_read_threshold[i, "Titles"]
  clones_in_locus = grepl(V_segment, input_data$V_Segment_Major_Gene) & grepl(J_segment, input_data$J_Segment_Major_Gene)
  fltr = input_data$Total_Read_Count < read_threshold
  print(paste("Clones in locus", title, ":", sum(clones_in_locus), "Clones below threshold:", sum(fltr & clones_in_locus)))
}





#stable rearrengements (freq>0.01 & read>100)

complete_data_fltr = complete_data[complete_data$Frequency > 0.01 & complete_data$Normalized_Read_count > 100,]

#for document28.doc
library(reshape2)

complete_data_fltr = read.table(file = "D:/wd/prisca/heatmaps_vincent/dx_rec_from_orig_complete_data_fltr.txt", header = T, sep = "\t", quote = "", stringsAsFactors = F)

frequencies = c(0, 0.01, 0.1, 1, 5)

combined = tbl[NULL,]
combined.heatmap.data = combined.heatmap.data[NULL,]

for(freq in frequencies){
  tmp = complete_data_fltr[complete_data_fltr$Frequency_diag > freq | complete_data_fltr$type == "recid",]
  tmp$frequency = freq
  combined.heatmap.data = rbind(combined.heatmap.data, tmp)
  print(paste(freq, nrow(tmp)))
  tbl = data.frame(table(tmp$patient, tmp$type))
  names(tbl) = c("Patient", "type", "count")
  tbl$frequency = freq
  combined = rbind(combined, tbl)
}

combined$Patient.type = paste(combined$Patient, combined$type)
result.1 = dcast(combined, frequency ~ Patient.type, value.var = "count")

result.1$diag_mean = rowMeans(result.1[,names(result.1)[grepl("diag", names(result.1))]])
result.1$recid_mean = rowMeans(result.1[,names(result.1)[grepl("recid", names(result.1))]])
result.1$both_mean = rowMeans(result.1[,names(result.1)[grepl("both", names(result.1))]])

result.1$diag_positive = rowSums(result.1[,names(result.1)[grepl("diag", names(result.1))]] > 0)
result.1$recid_positive = rowSums(result.1[,names(result.1)[grepl("recid", names(result.1))]] > 0)
result.1$both_positve = rowSums(result.1[,names(result.1)[grepl("both", names(result.1))]] > 0)

result.1$diag_positve_perc = result.1$diag_positive / 43 * 100
result.1$recid_positve_perc = result.1$diag_positive / 43 * 100
result.1$both_positve_perc = result.1$both_positve / 43 * 100

#result.1[,(ncol(result.1) - 8):ncol(result.1)]



write.table(result.1, file = "clipboard", append = F, quote = F, sep = "\t", row.names = F, col.names = T)

























