library(data.table)
library(stringr)
library(dplyr)
library(tibble)

# Args command
args=commandArgs(T)
id_input_file=args[1] # id input file: file with 1 column containing ids to convert (header must be "varId" or "rsId")
type_of_id_to_convert=as.numeric(args[2]) # type of id to convert: varID to rsID = 1; rsID to varID = 2
output_file=args[3]

create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}

GRCh38_gtex_id_converter <- function(id_input_file, type_of_id_to_convert, output_file) {
  
  id_file = fread(id_input_file)
  
  if (type_of_id_to_convert == 1){
    id_file = fread("GTEx_aorta_sqtls_varIDs.txt")
    
    dummy = data.frame(str_split_fixed(id_file$varID, "_", 5))
    dummy2 = as.data.frame(table(dummy[,1]))
    colnames(dummy2) = c("Chr", "Freq")
    dummy3 = data.frame(str_split_fixed(dummy2$Chr, "r", "2"))
    dummy2$Chr = as.numeric(dummy3[,2])
    dummy2 = dummy2[order(dummy2$Chr),]
    
    ##Downloads only necessary chr data and creates a dataframe of varID and rsID per chr
    for (chr in dummy2$Chr) {
      print(paste0("Printing chr",chr," data"))
      chr_data = fread(paste0("GTEx_variants/chr",chr,"_GRCh38_gtex_ids.txt.gz"), header = F)
      var_IDs = data.frame(id_file[dummy$X1 == paste0("chr",chr),])
      colnames(var_IDs) = c("varId")
      chr_join = left_join(var_IDs, chr_data, by = c("varId" = "V1"))
      x = paste0("chr",chr,"_df")
      assign(noquote(x), data.frame("varId" = chr_join$varId, "rsID" = chr_join$V2))
    }
    
    ##Creates a table to know which chr_df files were created
    dfs = create_empty_table(nrow(dummy2), 2)
    colnames(dfs) = c("chr", "df_name")
    dfs$chr = paste0("chr",dummy2$Chr)
    dfs$df_name = paste0(dfs$chr, "_df")
    
    ##Merges all chr_dfs
    df1 = data.frame(mget(dfs$df_name[1]))
    colnames(df1) = c("varID", "rsID")
    for (i in 2:(nrow(dfs))) {
      df2 = data.frame(mget(dfs$df_name[i]))
      colnames(df2) = c("varID", "rsID")
      bind_df = rbind(df1, df2)
      df1 = bind_df
    }
  }else if (type_of_id_to_convert == 2){
    ##Downloads all the chr_data file from github
    chr1_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr1_GRCh38_gtex_ids.txt.gz", header = F)
    chr2_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr2_GRCh38_gtex_ids.txt.gz", header = F)
    chr3_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr3_GRCh38_gtex_ids.txt.gz", header = F)
    chr4_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr4_GRCh38_gtex_ids.txt.gz", header = F)
    chr5_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr5_GRCh38_gtex_ids.txt.gz", header = F)
    chr6_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr6_GRCh38_gtex_ids.txt.gz", header = F)
    chr7_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr7_GRCh38_gtex_ids.txt.gz", header = F)
    chr8_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr8_GRCh38_gtex_ids.txt.gz", header = F)
    chr9_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr9_GRCh38_gtex_ids.txt.gz", header = F)
    chr10_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr10_GRCh38_gtex_ids.txt.gz", header = F)
    chr11_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr11_GRCh38_gtex_ids.txt.gz", header = F)
    chr12_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr12_GRCh38_gtex_ids.txt.gz", header = F)
    chr13_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr13_GRCh38_gtex_ids.txt.gz", header = F)
    chr14_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr14_GRCh38_gtex_ids.txt.gz", header = F)
    chr15_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr15_GRCh38_gtex_ids.txt.gz", header = F)
    chr16_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr16_GRCh38_gtex_ids.txt.gz", header = F)
    chr17_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr17_GRCh38_gtex_ids.txt.gz", header = F)
    chr18_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr18_GRCh38_gtex_ids.txt.gz", header = F)
    chr19_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr19_GRCh38_gtex_ids.txt.gz", header = F)
    chr20_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr20_GRCh38_gtex_ids.txt.gz", header = F)
    chr21_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr21_GRCh38_gtex_ids.txt.gz", header = F)
    chr22_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chr22_GRCh38_gtex_ids.txt.gz", header = F)
    ##Downloads chrX data
    chrX_data = fread("https://raw.github.com/nbbarrientos/hg38_gtex_ids_data/master/chrX_GRCh38_gtex_ids.txt.gz", header = F)
    
    ##Merges all chr_data files
    df1 = chr1_data
    for (chrom in 2:22){
      df2 = noquote(paste0("chr",chrom,"_data"))
      temp_bind = rbind(df1, df2)
      df1 = temp_bind
    }
    chr_data = rbind(temp_bind, chrX_data)
    
    ##Joins input file with chr_data
    chr_join = left_join(id_file, chr_data, by = c("rsId" = "V2"))
    bind_df = data.frame("rsID" = id_file$rsId, "varID" = chr_join$V2)
  }else 
    print("Please specify the type of variant id you wish to convert")
  
  
  write.table(bind_df, output_file, quote = FALSE, row.names = F, col.names = TRUE, sep = "\t")
  return(bind_df)
}

test = fread("GTEx_varIDs.txt")
#unique_varids = data.frame(unique(test$varId))
#write.table(test, "test_rsID_to_varID.txt", col.names = TRUE, row.names = F)
GRCh38_gtex_id_converter("test_rsID_to_varID.txt", 2, "testforadam.txt")
