read_blast <- function(blast_file){
  read_tsv(blast_file, col_names =  c("qseqid", "qacc" ,"sseqid" ,"sacc" ,"qlen" ,"slen" ,"qstart" ,"qend" ,"sstart" ,"send" ,"length" ,"pident" ,"qcovhsp", "evalue"))
}
#read in one blast result file
read_tsv(blast_file, col_names =  c("qseqid", "qacc" ,"sseqid" ,"sacc" ,"qlen" ,"slen" ,"qstart" ,"qend" ,"sstart" ,"send" ,"length" ,"pident" ,"qcovhsp", "evalue"))

# or use the function
read_blast('blast_file')

# read in many blast result files and bind them together
all_blasts <- bind_rows(lapply(blast_files, read_blast))






