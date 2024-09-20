library(AnnotationForge) #load
library(GO.db)
args <- commandArgs(T)
egg<-read.csv(args[1],header=T,sep="\t") 
egg[egg==""] <- NA

###Extract GO information from annotations###
library(dplyr)
library(stringr)
gene_info <- egg %>%dplyr::select(GID = query, GENENAME = seed_ortholog) %>% na.omit() #Extract the first two columns based on the titles of the first and second columns of the egg.
goterms <- egg %>%dplyr::select(query, GOs) %>% na.omit() %>% filter(str_detect(GOs,"GO"))
all_go_list=str_split(goterms$GOs,",") #Separate the GOs values of goterms with commas as separators
gene2go <- data.frame(GID = rep(goterms$query, times = sapply(all_go_list, length)), GO = unlist(all_go_list), EVIDENCE = "IEA") %>% filter(str_detect(GO,"GO"))

###Extract KEGG information from annotations###
koterms <- egg %>%dplyr::select(GID = query, KO=KEGG_ko)%>%na.omit()%>% filter(str_detect(KO,"ko")) #Extract KEGG annotations for genes based on the titles in the first column of the egg and the KEGG_ko column

#load("/dellfsqd2/ST_OCEAN/USER/liyanan2/03_job/single_cell_seq_analysis/make.OrgDb/kegg_info.RData")
load(args[2])
colnames(ko2pathway)=c("KO",'Pathway')
koterms$KO=str_replace_all(koterms$KO,"ko:","") #Remove the prefix 'ko:' from the KO value of 'koterm' and make it consistent with the format of 'ko2pathway'.
gene2pathway <- koterms %>% left_join(ko2pathway, by = "KO") %>%dplyr::select(GID, Pathway) %>%na.omit() #Merge koterm and ko2pathway into gene2pathway, and organize the correspondence between genes and pathways
gene2pathway_name<-left_join(gene2pathway,pathway2name,by="Pathway") #Merge gene2pathway and pathway2name
write.table(gene2pathway_name,file="gene2pathway_name.txt",sep="\t",row.names=F,quote=F) 

makeOrgPackage(gene_info=gene_info, go=gene2go, ko=koterms,  pathway=gene2pathway, version="0.0.1", maintainer='liyanan2 <liyanan2@genomics.cn>', author='liyanan2',outputDir=".", tax_id=args[3], genus=args[4], species=args[5],goTable="go")


