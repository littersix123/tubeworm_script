library(ggplot2)

rm(list=ls())

kegg <- read.table("PL-M1_GO_term.xls",stringsAsFactors=F,sep="\t",quote="")

kegg <- kegg[order(kegg$p.adjust),]

kegg <- kegg[1:10,]

top10 <- data.frame(kegg$Description,kegg$Count ,kegg$p.adjust)

colnames(top10) <- c("Description","count","padj")

p <- ggplot(data=top10,aes(x=Description,y=count,fill=padj))

p1 <- p + geom_bar(stat="identity") + coord_flip()

p2 <- p1 + theme(panel.background=element_rect(fill='transparent',color='gray'),
                 axis.text.y=element_text(color="black",size=12))

p3 <- p2 + scale_fill_gradient(low="red",high="blue")

p4 <- p3 + scale_x_discrete(limits=rev(top10[,1])) +labs(x="",y="",title="PL_M1")

pdf("PL-M1_KEGG_bar_plot.pdf",width=9)

print(p4)

dev.off()
