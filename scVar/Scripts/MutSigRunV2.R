library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
cosmic_ref<-read.table("/codes/COSMIC_v3.3.1_SBS_GRCh38.txt",sep = "\t",header = T)
transposed_df <- as.data.frame(t(cosmic_ref))
selected_rows <- transposed_df[-1, ]
names <- unlist(transposed_df[1,])
names(selected_rows) <-names
test<-signatures.cosmic
for (i in 1:nrow(selected_rows)) {
  for (j in 1:ncol(selected_rows)) {
    if(is.character(selected_rows[i, j])){
      test[i, j] =as.double(selected_rows[i, j])
    }
  }
}
rownames(test)<-rownames(selected_rows)
colnames(test)<-colnames(selected_rows)

###all_mutations
sample.mut.ref<-read.table("Mut.table.for.signatures_all_normal.tsv",sep = "\t",header = T)
colnames(sample.mut.ref)<-c("Sample","chr","pos","ref","alt")
sample.mut.ref$pos<-as.numeric(sample.mut.ref$pos)
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref, sample.id = "Sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt",bsg=BSgenome.Hsapiens.UCSC.hg38)
row_sums<-rowSums(sigs.input)
df<-sigs.input/row_sums
write.csv(x = df,file ="Signature_celltype_all.csv",quote = F,row.names = T)
labels_test<-list('group1','group2','group3','group4','group5','group6')
repeated_labels<-unlist(lapply(labels_test,function(x) rep(x,16)))
mysignatures <- whichSignatures(tumor.ref = sigs.input,  signatures.ref = test,  sample.id = "Sample_all",contexts.needed = T)
similarity <-as.data.frame(t(mysignatures$weights))
similarity$Signature <- rownames(similarity)
all<-rbind(similarity,c(mysignatures$unknown,'unknown'))
names(all)<- c("Similarity","Signature")
out<-all[all$Similarity!=0,]
col_index1<-which(names(out)=='Similarity')
col_index2<-which(names(out)=='Signature')
out<-out[,c(col_index2,col_index1)]
path=paste("all","cosmic.csv",sep="_")
path_final=paste("cell_type",path,sep="/")
write.csv(x = out,file =path_final,quote = F,row.names = F)
row_data<-df[1,]
h<-as.data.frame(t(row_data))
h$labels<-repeated_labels
path_sig=paste("all","signature.csv",sep="_")
path_final=paste("cell_type",path_sig,sep="/")
write.table(x=h,path_final,sep=",",col.names = F,quote = F,row.names = T)
# write the mutation signatures
#write.table(x = all,file ="SignatureSimilarityMatrix_new.txt",sep="\t",quote = F,row.names = F)
#write.table(x = sigs.input,file ="Signature.txt",sep="\t",quote = F,row.names = T)



##cell_type
sample.mut.ref.celltype<-read.table("Mut.table.for.signatures_cell_type_normal.tsv",sep = "\t",header = F)
colnames(sample.mut.ref.celltype)<-c("Sample","chr","pos","ref","alt")
sample.mut.ref.celltype$pos<-as.numeric(sample.mut.ref.celltype$pos)
sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref.celltype, sample.id = "Sample", chr = "chr", pos = "pos", ref = "ref", alt = "alt",bsg=BSgenome.Hsapiens.UCSC.hg38)
row_sums<-rowSums(sigs.input)
df<-sigs.input/row_sums

write.csv(x = df,file ="Signature_celltype.csv",quote = F,row.names = T)
#write.table(x = sigs.input,file ="Signature_celltype.txt",sep="\t",quote = F,row.names = T)
labels_test<-list('group1','group2','group3','group4','group5','group6')
repeated_labels<-unlist(lapply(labels_test,function(x) rep(x,16)))

for (i in unique(sample.mut.ref.celltype$Sample)) {
  mysignatures <- whichSignatures(tumor.ref = sigs.input,  signatures.ref = test,  sample.id =i ,contexts.needed = T)
  similarity <-as.data.frame(t(mysignatures$weights))
  similarity$Signature <- rownames(similarity)
  all<-rbind(similarity,c(mysignatures$unknown,'unknown'))
  names(all)<- c("Similarity","Signature")
  out<-all[all$Similarity!=0,]
  col_index1<-which(names(out)=='Similarity')
  col_index2<-which(names(out)=='Signature')
  out<-out[,c(col_index2,col_index1)]
  path=paste(i,"cosmic.csv",sep="_")
  path_final=paste("cell_type",path,sep="/")
  write.csv(x = out,file =path_final,quote = F,row.names = F)
  row_data<-df[i,]
  h<-as.data.frame(t(row_data))
  h$labels<-repeated_labels
  path_sig=paste(i,"signature.csv",sep="_")
  path_final_sig=paste("cell_type",path_sig,sep="/")
  write.table(x=h,path_final_sig,sep=",",col.names = F,quote = F,row.names = T)
}
#mysignatures <- whichSignatures(tumor.ref = sigs.input,  signatures.ref = signatures.cosmic,  sample.id = ,contexts.needed = T)


