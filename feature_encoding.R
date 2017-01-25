library(DNAshapeR)

# Save temporary FASTA file for DNAShapeR
prepare_fasta <- function(snps, fasta_file) {
   lines <- c()
   for(i in 1:nrow(snps)) {
      lines <- c(lines, paste(">seqA", i, sep=""))
      lines <- c(lines, snps[i, 1])
      lines <- c(lines, paste(">seqB", i, sep=""))
      lines <- c(lines, snps[i, 2])
   }
   fcon <- file(fasta_file)
   writeLines(lines, fcon)
   close(fcon)
}


# Get DNA shapes
calculate_shape_features <- function(snps) {
   fasta_file <- "tmp.fn"
   prepare_fasta(snps, fasta_file)
   shapes <- getShape(fasta_file)

   mgw <- data.frame(shapes$MGW[, 3:7]) 
   prot <- data.frame(shapes$ProT[, 3:7]) 
   helt <- data.frame(shapes$HelT[, 2:7]) 
   roll <- data.frame(shapes$Roll[, 2:7]) 

   wt_indices <- seq(1, nrow(mgw), 2)
   snp_indices <- seq(2, nrow(mgw), 2)

   wt_mgw <- mgw[wt_indices,]
   wt_prot <- prot[wt_indices,]
   wt_helt <- helt[wt_indices,]
   wt_roll <- roll[wt_indices,]

   snp_mgw <- mgw[snp_indices,]
   snp_prot <- prot[snp_indices,]
   snp_helt <- helt[snp_indices,]
   snp_roll <- roll[snp_indices,]

   colnames(wt_mgw) <- c("MGW_WT_0", "MGW_WT_1", "MGW_WT_2", "MGW_WT_3", "MGW_WT_4")
   colnames(wt_prot) <- c("PROT_WT_0", "PROT_WT_1", "PROT_WT_2", "PROT_WT_3", "PROT_WT_4")
   colnames(wt_helt) <- c("HELT_WT_0", "HELT_WT_1", "HELT_WT_2", "HELT_WT_3", "HELT_WT_4", "HELT_WT_5")
   colnames(wt_roll) <- c("ROLL_WT_0", "ROLL_WT_1", "ROLL_WT_2", "ROLL_WT_3", "ROLL_WT_4", "ROLL_WT_5")

   colnames(snp_mgw) <- c("MGW_SNP_0", "MGW_SNP_1", "MGW_SNP_2", "MGW_SNP_3", "MGW_SNP_4")
   colnames(snp_prot) <- c("PROT_SNP_0", "PROT_SNP_1", "PROT_SNP_2", "PROT_SNP_3", "PROT_SNP_4")
   colnames(snp_helt) <- c("HELT_SNP_0", "HELT_SNP_1", "HELT_SNP_2", "HELT_SNP_3", "HELT_SNP_4", "HELT_SNP_5")
   colnames(snp_roll) <- c("ROLL_SNP_0", "ROLL_SNP_1", "ROLL_SNP_2", "ROLL_SNP_3", "ROLL_SNP_4", "ROLL_SNP_5")

   shape_features <- cbind(wt_mgw, wt_prot, wt_helt, wt_roll, snp_mgw, snp_prot, snp_helt, snp_roll)
   rownames(shape_features) <- NULL

   file.remove("tmp.fn")
   file.remove("tmp.fn.HelT")
   file.remove("tmp.fn.MGW")
   file.remove("tmp.fn.ProT")
   file.remove("tmp.fn.Roll")

   return(shape_features)
}


# Calculate GC-content features
calculate_gc_content <- function(snps) {

   calc_gc <- function(seq) {
      counts = table(seq)
      if(is.na(counts["G"])) counts["G"] <- 0
      if(is.na(counts["C"])) counts["C"] <- 0
      return((counts["C"] + counts["G"])/sum(counts))
   }

   wt_gcscore_9 <- c()
   snp_gcscore_9 <- c()
   wt_gcscore_7 <- c()
   snp_gcscore_7 <- c()

   for(i in 1:nrow(snps)) {
      wt_gcscore_9 = c(wt_gcscore_9, calc_gc(strsplit(snps[i, 1], "")[[1]]))
      snp_gcscore_9 = c(snp_gcscore_9, calc_gc(strsplit(snps[i, 2], "")[[1]]))
      wt_gcscore_7 = c(wt_gcscore_7, calc_gc(strsplit(snps[i, 1], "")[[1]][2:8]))
      snp_gcscore_7 = c(snp_gcscore_7, calc_gc(strsplit(snps[i, 2], "")[[1]][2:8]))
   }

   res <- cbind(wt_gcscore_9, snp_gcscore_9, wt_gcscore_7, snp_gcscore_7)
   colnames(res) <- c("WT_GCSCORE_9", "SNP_GCSCORE_9", "WT_GCSCORE_7", "SNP_GCSCORE_7")
   return(res)
}


# Encode sequence and mutation features
encode_sequence <- function(snps) {
   res <- c()
   for(i in 1:nrow(snps)) {
      row <- c()
      seq <- strsplit(snps[i, 2], "")[[1]]
      for(j in 1:length(seq)) {
         row <- c(row, as.numeric(seq[j] == c("A", "C", "G", "T")))
      }
      snp_seq_4 <- seq[5]
      wt_seq_4 <- strsplit(snps[i, 1], "")[[1]][5]
      row <- c(row, as.numeric(wt_seq_4 == c("A", "C", "G", "T")))
      row <- c(row, as.numeric(paste0(wt_seq_4, snp_seq_4) ==
                               c("AC","AG","AT","CA","CG", "CT","GA","GC","GT","TA", "TC","TG")))
      res <- rbind(res, row)
   }
   res <- data.frame(res)
   colnames(res) <- c("SNP_SEQ_0_A", "SNP_SEQ_0_C", "SNP_SEQ_0_G", "SNP_SEQ_0_T",
                      "SNP_SEQ_1_A", "SNP_SEQ_1_C", "SNP_SEQ_1_G", "SNP_SEQ_1_T",
                      "SNP_SEQ_2_A", "SNP_SEQ_2_C", "SNP_SEQ_2_G", "SNP_SEQ_2_T",
                      "SNP_SEQ_3_A", "SNP_SEQ_3_C", "SNP_SEQ_3_G", "SNP_SEQ_3_T",
                      "SNP_SEQ_4_A", "SNP_SEQ_4_C", "SNP_SEQ_4_G", "SNP_SEQ_4_T",
                      "SNP_SEQ_5_A", "SNP_SEQ_5_C", "SNP_SEQ_5_G", "SNP_SEQ_5_T",
                      "SNP_SEQ_6_A", "SNP_SEQ_6_C", "SNP_SEQ_6_G", "SNP_SEQ_6_T",
                      "SNP_SEQ_7_A", "SNP_SEQ_7_C", "SNP_SEQ_7_G", "SNP_SEQ_7_T",
                      "SNP_SEQ_8_A", "SNP_SEQ_8_C", "SNP_SEQ_8_G", "SNP_SEQ_8_T",
                      "WT_SEQ_4_A", "WT_SEQ_4_C", "WT_SEQ_4_G", "WT_SEQ_4_T",
                      "WT_SNP_AC", "WT_SNP_AG", "WT_SNP_AT", "WT_SNP_CA", "WT_SNP_CG",
                      "WT_SNP_CT", "WT_SNP_GA", "WT_SNP_GC", "WT_SNP_GT", "WT_SNP_TA",
                      "WT_SNP_TC", "WT_SNP_TG")
   return(res)
}


# Enrich feature set with meaningful combinations
add_feature_combinations <- function(data) {
   pairs <- rbind(c("WT_GCSCORE_7", "SNP_GCSCORE_7"),
                  c("WT_GCSCORE_9", "SNP_GCSCORE_9"))
   pairs <- rbind(pairs, t(sapply(0:5, function(i) c(paste0("HELT_SNP_", i), paste0("HELT_WT_", i)))))
   pairs <- rbind(pairs, t(sapply(0:5, function(i) c(paste0("ROLL_SNP_", i), paste0("ROLL_WT_", i)))))
   pairs <- rbind(pairs, t(sapply(0:4, function(i) c(paste0("MGW_SNP_", i), paste0("MGW_WT_", i)))))
   pairs <- rbind(pairs, t(sapply(0:4, function(i) c(paste0("PROT_SNP_", i), paste0("PROT_WT_", i)))))

   for(i in 1:nrow(pairs)) {
      data[paste0(pairs[i, 1], "_", pairs[i, 2], "_diff")] <- data[pairs[i, 1]] - data[pairs[i, 2]]
      data[paste0(pairs[i, 1], "_", pairs[i, 2], "_mul")] <- data[pairs[i, 1]] * data[pairs[i, 2]]
   }

   return(data)
}


# Prepare data frame with complete features
encode_features <- function(snps) {
   colnames(snps) <- c("WT_SEQ", "SNP_SEQ")
   for(i in 1:nrow(snps)) {
      seq <- snps[i, 1]
      substring(seq, 5, 5) <- snps[i, 2]
      snps[i, 2] <- seq
   }

   data <- cbind(encode_sequence(snps), calculate_gc_content(snps), calculate_shape_features(snps))
   data <- add_feature_combinations(data)
   data <- data[, order(colnames(data))]

   return(data)
}
