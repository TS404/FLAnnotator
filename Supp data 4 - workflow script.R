# Setup --------------
#> Install packages #-----------------

packagelist =       c("igraph",
                      "readxl",
                      "tidyr",
                      "stringr",
                      "stringi",
                      "ape",
                      "umap",
                      "dbscan",
                      "factoextra",
                      "rgl",
                      "gplots",
                      "ggplot2",
                      "RColorBrewer",
                      "GoodmanKruskal")
githubpackagelist = c("missuse/ragp",
                      "shaunpwilkinson/kmer",
                      "shaunpwilkinson/insect")
source("http://bioconductor.org/biocLite.R")
biocLite("Biostrings")

devtools::install_github(githubpackagelist)
install.packages(packagelist)
library(githubpackagelist,
        packagelist,
        Biostrings)


colours<-palette(c("white",           #-
                   "maroon",          #A
                   "brown",           #B
                   "red2",            #C
                   "pink2",           #D
                   "green3",          #E #5
                   "lightgreen",      #F
                   "green4",          #G
                   "yellow",          #H
                   "mediumpurple4",   #I
                   "darkcyan",        #J #10
                   "cyan",            #K
                   "grey40",          #L
                   "magenta",         #M
                   "navyblue",        #N
                   "royalblue1",      #O #15
                   "black",           #P
                   "cornflowerblue",  #Q
                   "orange",          #R
                   rep("white",20)))  #etc

#> Gather files ---------------

file.full <- "C:\\Users\\TS404\\FLAs\\Supp data 2 - fasciclin alignments subset by domain type\\umap_group.all.fa"
all_1_fas.aln <- Biostrings::readAAStringSet("C:\\Users\\TS404\\FLAs\\Supp data 2 - fasciclin alignments subset by domain type\\1-FLA(aln).fa")
all_2_fas.aln <- Biostrings::readAAStringSet("C:\\Users\\TS404\\FLAs\\Supp data 2 - fasciclin alignments subset by domain type\\2-FLA(aln).fa")

SEQF <- Biostrings::readAAStringSet(file.full)

FLAnamesfull   <- readxl::read_xlsx("C:\\Users\\TS404\\FLAs\\Supp data 1 - Domain names and labels.xlsx")
FLAnames.FLA   <- (FLAnamesfull[names(SEQF2),c(6,2,3,7)])
FLAnames.fas   <- (FLAnamesfull[labels$number,c(6,2,3,7)])
FLAnames.fas.o <- (FLAnamesfull[labels.o$number,c(6,2,3,7)])

FLAtaxafull <- read.csv("C:\\Users\\TS404\\FLAs\\species list flas.csv")
rownames(FLAtaxafull) <- rownames(labels.o)
FLAtaxafull.r <- FLAtaxafull[,ncol(FLAtaxafull):1]
FLAtaxafull.r <- t(apply(FLAtaxafull, 1, function(x) c(x[!is.na(x)], x[is.na(x)])))
FLAtaxafull.r <- t(tidyr::fill(FLAtaxafull.r.t,colnames(FLAtaxafull.r.t)))

#########.
# MAAB+ #--------------------------
#########.

#> FLA annotation -------------
SEQF2 <- lapply(SEQF,as.character)

# annot <- NULL
t     <- Sys.time()
total <- (length(SEQF2)/10)
pb    <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1:total){
  SEQF3 <- SEQF2[i*10+(1:10)]
  annot <- rbind(annot, ragp::get_hmm(sequence = SEQF3,
                                      id = names(SEQF3),
                                      verbose = FALSE,
                                      sleep = 0.1))
  print(i)
  setTxtProgressBar(pb, i)
  print(Sys.time()-t)
}

SEQN <- Biostrings::readAAStringSet(file.full)
SEQN
SEQN2 <- lapply(SEQN,as.character)

# slow HMM check of all sequences
t  <- Sys.time()
pb <- txtProgressBar(min = 0, max = length(SEQN), style = 3)
annot.notfla <- NULL
for(i in 1:length(SEQN2)){
  SEQN3 <- SEQN2[i]
  annot.notfla <- rbind(annot.notfla, ragp::get_hmm(sequence = SEQN3,
                        id = names(SEQN3),
                        verbose = FALSE,
                        sleep = 0.1))
  print(i)
  setTxtProgressBar(pb, i)
  print(Sys.time()-t)
}

annotF <- annot[annot[,2]=="Fasciclin",]
annotF <- annotF[!is.na(annotF[,2]),]

# Number of Fla domains
FLAcount<-NULL
for(i in names(SEQF2)){
  FLAcount[i] <- sum(annotF[,1]==i)
}

#> FLA extraction  ---------------
SEQF4 <- lapply(SEQF2,strsplit,"")
SEQF4 <- lapply(SEQF4,"[[",1)
SEQF4 <- lapply(SEQF4,casefold,upper=FALSE)
# uppercase FLAs
for(i in 1:nrow(annotF)){
  SEQF4[[annotF[i,1]]] [annotF[i,6]:annotF[i,7]] <- casefold(SEQF4[[annotF[i,1]]] [annotF[i,6]:annotF[i,7]],upper=TRUE)
}
as.fasta(SEQF4             ,write = "Full length (annotated).fa")
as.fasta(SEQF4[FLAcount==0],write = "Full length (annotated) 0-FLA.fa")
as.fasta(SEQF4[FLAcount==1],write = "Full length (annotated) 1-FLA.fa")
as.fasta(SEQF4[FLAcount==2],write = "Full length (annotated) 2-FLA.fa")
as.fasta(SEQF4[FLAcount>=3],write = "Full length (annotated) 3+FLA.fa")

counter<-NULL
for(i in 1:length(FLAcount)){
  if(FLAcount[i]>0){
    counter <- append(counter,1:FLAcount[i])
  }
}
FLAmax <- rep(FLAcount,FLAcount)
FLAnames <- paste0(annotF[,1],".FLA.",counter,"/",FLAmax)
SEQF5 <- SEQF4[FLAcount>0]
FLAs  <- list()
for(i in 1:nrow(annotF)){
  start <-annotF[i,6] 
  end   <-annotF[i,7] 
  if(annotF[i,8]<5){
    start <- start - annotF[i,8]          # if starts <5aa short, also collect missing end residues
  }
  if(start<0){                            # prevent negaive start numbers
    start <- 0
  }
  if((128-annotF[i,9])<5){                # if stops <5aa short, also collect missing end residues
    end <- end + (128-annotF[i,9])
  }
  FLAs[[ FLAnames[i] ]] <- SEQF5[[ annotF[i,1] ]][start:end]
}
FLAs <- lapply(FLAs, function(x) x[!is.na(x)])
FLAs <- FLAs[lapply(FLAs, length)>50]
as.fasta(FLAs,write = "ALLFLAS.fa")

FLAsX  <- list()
for(i in 1:nrow(annotF)){
  start <-annotF[i,6] - annotF[i,8]       - 5    # also collect N-ter residues not picked up by model and 5 extra
  end   <-annotF[i,7] + (128-annotF[i,8]) + 5    # also collect C-ter residues not picked up by model and 5 extra
  if(start<0){                            # prevent negaive start numbers
    start <- 0
  }
  FLAsX[[ FLAnames[i] ]] <- SEQF5[[ annotF[i,1] ]][start:end]
  
  start2 <- which(diff(unlist(FLAsX[i])==casefold(unlist(FLAsX[i]),TRUE))==1)[2]
  if(!is.na(start2)){
    FLAsX[[ FLAnames[i] ]] <- FLAsX[[ FLAnames[i] ]] [1:start2]
    }
}
FLAsX <- lapply(FLAsX, function(x) x[!is.na(x)])
FLAsX <- FLAsX[lapply(FLAsX, length)>50]
as.fasta(FLAsX,write = "ALLFLAS plus overhangs.fa")

boxplot(t(head(annotF[order(annotF[,9]-annotF[,8]),8:9],n=2658)))

#> FLA masking  -----------
SEQF6 <- SEQF4
for(i in 1:nrow(annotF)){
  SEQF6[[annotF[i,1]]] [annotF[i,6]:annotF[i,7]] <- "X"
}
as.fasta(SEQF6, write = "Full length (masked).fa") # X-masked FLA domains
SEQF6.1 <- lapply(SEQF6,paste, collapse = "")

#> AGP capitals ------
SEQF7 <- ragp::scan_ag(sequence = SEQF2,
                       id = names(SEQF2),
                       dim = 3,div = 8)$sequence

names <- paste0(">",names(SEQF2))
ord1 <- 2*(1:length(names))-1
ord2 <- 2*(1:length(SEQF7))

cat(c(names,SEQF7)[order(c(ord1,ord2))],sep = "\n",
    file =  "Full length (AGP capitals).fa") # Capitalised AGP regions

#> AG region stats -------
agregions  <- ragp::scan_ag(sequence = SEQF7,
                            id       = names(SEQF2),
                            dim      = 3,
                            div      = 6,
                            type     = "extended",
                            simplify = FALSE) $locations

agreg.seqs <- NULL
for(i in 1:length(agreg.lengths)){
  agreg.seqs[[i]] <- apply(agregions[[i]],1, function(x) substr(SEQF7[[i]], x[1]-1, x[2]+1))
}

agreg.seqs.AAbin <- ape::as.AAbin(as.list(unlist(lapply(agreg.seqs,as.list))),residues = "AA")
agreg.length     <- lapply(agregions, function(x) x[,2]-x[,1])
agreg.length.sum <- unlist(lapply(agreg.length,sum))
agreg.length.av  <- unlist(lapply(agreg.length,mean))
agreg.count      <- unlist(lapply(agregions, nrow))
agreg.numP       <- lapply(agreg.seqs,stringr::str_count,"P")
agreg.numP.sum   <- unlist(lapply(agreg.numP,sum))
overall.numP.sum <- stringr::str_count(SEQF6.1, pattern= "p")
agreg.names      <- rep(names(SEQF2),agreg.count)
agreg <- cbind(agreg.length.sum,
               agreg.length.av,
               agreg.count,
               agreg.numP.sum,
               overall.numP.sum)
rownames(agreg) <- names(SEQF2)

hist(agreg.count,breaks = -1:max(agreg.count))
hist(agreg.length.sum,breaks = -1:max(agreg.length.sum))
hist(agreg.numP.sum,breaks = -1:max(agreg.numP.sum))
hist(unlist(agreg.lengths),breaks = 0:max(unlist(agreg.lengths)),xlim=c(0,100))
hist(unlist(agreg.numPs),breaks = 0:max(unlist(agreg.numPs)))
heatmap(log10(table(agreg.count,FLAcount)+0.1),Colv=NA,Rowv=NA,scale="none",xlab="n.Fas", ylab="n.agreg")


as.matrix(head((rowsum(agreg.length.av, fasclass[,2],na.rm = 1)/
                  as.matrix(table(fasclass[,2])))
               [order(table(fasclass[,2]),decreasing = 1),],20))

as.matrix(head((rowsum(agreg.length.sum, fasclass[,2],na.rm = 1)/
                  as.matrix(table(fasclass[,2])))
               [order(table(fasclass[,2]),decreasing = 1),],20))

as.matrix(head((rowsum(agreg.count, fasclass[,2],na.rm = 1)/
                  as.matrix(table(fasclass[,2])))
               [order(table(fasclass[,2]),decreasing = 1),],20))

as.matrix(head((rowsum(agreg.numP.sum, fasclass[,2],na.rm = 1)/
                  as.matrix(table(fasclass[,2])))
               [order(table(fasclass[,2]),decreasing = 1),],20))

#>> k3 kmers ----
ape::write.FASTA(agreg.seqs.AAbin,"agreg.seqs.AAbin.fa")
agreg.seqs.AAbin.2 <- insect::readFASTA("agreg.seqs.AAbin.fa",residues = "AA")

agreg.k3 <- kmer::kcount(agreg.seqs.AAbin.2,k=3)
agreg.k3 <- do.call(rbind,by(data=agreg.k3,INDICES=agreg.names,FUN=colSums))
agreg.k3 <- rbind(agreg.k3,array(0,dim = list(sum(agreg.count==0),ncol(agreg.k3)),
                                 dimnames = list(names(SEQF2)[agreg.count==0],NULL)))
agreg.k3 <- agreg.k3[names(SEQF2),]
agreg.k3 <- agreg.k3[,order(colSums(agreg.k3),decreasing = TRUE)]
agreg.k3 <- agreg.k3[,(colSums(agreg.k3)!=0)]

plot(100*(colSums(agreg.k3)/sum(agreg.k3))[1:40],type = "b",log="y")
txt<-colnames(agreg.k3)
txt[21:400]<-""
text(100*colSums(agreg.k3)/sum(agreg.k3),txt, cex=0.6,adj=c(-0.4,0.3), srt=45)

topcols <- 1:30
agreg.k3heat <- heatmap(agreg.k3[rowSums(agreg.k3[,topcols])!=0,topcols],keep.dendro = TRUE)
barplot(100*(colSums(agreg.k3)/sum(agreg.k3))[topcols],las=2)
barplot(100*(colSums(agreg.k3[,agreg.k3heat$colInd])/sum(agreg.k3))[topcols],las=2)


# k3clust<-(kmer::cluster(SEQ,3))
# plot(k3clust)

#> non-AG region stats -----------
SEQF9 <- SEQF6
for(i in 1:length(agregions)){
  if(nrow(agregions[[i]])>0){
    for(j in 1:nrow(agregions[[i]])){
      SEQF9[[i]] [agregions[[i]][j,1]:agregions[[i]][j,2]] <- "O"
    }
  }
}

SEQF9.1           <- lapply(SEQF9,paste, collapse = "")
nagreg.seqs       <- stringi::stri_split(SEQF9.1,regex = "X+|O+")
nagreg.seqs.AAbin <- ape::as.AAbin(as.list(unlist(lapply(nagreg.seqs,as.list))),residues = "AA")
nagreg.tooshort   <- lapply(nagreg.seqs.AAbin,length)<4
nagreg.seqs.AAbin <- nagreg.seqs.AAbin[!nagreg.tooshort]
nagreg.numP       <- ( unlist(lapply(SEQF9.1,stringr::str_count,"p")) / unlist(lapply(SEQF9.1,stringr::str_length)) )
ape::write.FASTA(nagreg.seqs.AAbin,"nagreg.seqs.AAbin.fa")
nagreg.seqs.AAbin.2 <- insect::readFASTA("nagreg.seqs.AAbin.fa",residues = "AA")

nagreg.count      <- unlist(lapply(nagreg.seqs,length))
nagreg.names      <- rep(names(SEQF2),nagreg.count)

FLAsX.1 <- lapply(FLAsX,paste, collapse = "")
FLAsX.numP <- ( unlist(lapply(FLAsX.1,stringr::str_count,"p")) / unlist(lapply(FLAsX.1,stringr::str_length)) )

Pro_average <- 0.047
boxplot(FLAsX.numP,nagreg.numP, (agreg.numP.sum/agreg.length.sum))
abline(h=Pro_average, lty=3)

as.matrix(head((rowsum(agreg.length.av, fasclass[,2],na.rm = 1)/
                  as.matrix(table(fasclass[,2])))
               [order(table(fasclass[,2]),decreasing = 1),],20))

#>> k3 kmers ----
nagreg.names.filtered <- nagreg.names[!nagreg.tooshort]
nagreg.k3 <- kmer::kcount(nagreg.seqs.AAbin.2[lapply(nagreg.seqs.AAbin.2,length)>3],k=3)
nagreg.k3 <- do.call(rbind,by(data=nagreg.k3,INDICES=nagreg.names.filtered,FUN=colSums))
nagreg.k3 <- rbind(nagreg.k3,array(0,dim = list(sum(nagreg.count==0),ncol(nagreg.k3)),
                                 dimnames = list(names(SEQF2)[nagreg.count==0],NULL)))
nagreg.k3 <- nagreg.k3[names(SEQF2),]
nagreg.k3 <- nagreg.k3[,order(colSums(nagreg.k3),decreasing = TRUE)]
nagreg.k3 <- nagreg.k3[,(colSums(nagreg.k3)!=0)]

plot(100*(colSums(nagreg.k3)/sum(nagreg.k3))[1:40],type = "b",log="y")
txt<-colnames(nagreg.k3)
txt[21:400]<-""
text(100*colSums(nagreg.k3)/sum(nagreg.k3),txt, cex=0.6,adj=c(-0.4,0.3), srt=45)

topcols <- 1:30
nagreg.k3heat <- heatmap(nagreg.k3[rowSums(nagreg.k3[,topcols])!=0,topcols],keep.dendro = TRUE)
barplot(100*(colSums(nagreg.k3)/sum(nagreg.k3))[topcols], las=2)
barplot(100*(colSums(nagreg.k3[,nagreg.k3heat$colInd])/sum(nagreg.k3))[topcols],las=2)


##########################.
# Inter-proline distance #--------------------------
##########################.
hist(stringi::stri_length(SEQ),breaks = 1000,xlim = c(0,400))

interP <- lapply(stringi::stri_split(SEQF,regex = "P"),
                 stringi::stri_length)

hist(interP[[2]],breaks = -1:max(unlist(interP)),plot=FALSE)$count

interP.hist <- lapply (interP,hist,breaks = -1:max(unlist(interP)),plot=FALSE)
interP.hist <- lapply(interP.hist, `[[`, 2)
interP.hist <- t(array(unlist(interP.hist),
                       dim = c(max(unlist(interP))+1,
                               length(SEQF))))

colnames(interP.hist)<- 0:max(unlist(interP))
rownames(interP.hist) <- names(SEQF)

num = 20
interP.hist2 <- cbind(interP.hist[,1:(num+1)],
                      rowSums(interP.hist[,(num+2):max(unlist(interP))]))
colnames(interP.hist2)[(num+2)]<-paste0(">",num)
rownames(interP.hist2) <- names(SEQF)
barplot(colSums(interP.hist2))
plot(colSums(interP.hist2), type="b")

interP.hist2.norm <- interP.hist2/rowSums(interP.hist2)

interP.FLAcount <- FLAcount[names(SEQF)]
heatmap(interP.hist2.norm,
        Colv = NA,
        RowSideColors=colours[interP.FLAcount+1]) # by number of Fas domains
heatmap(interP.hist2.norm,
        Colv = NA,
        RowSideColors=colours[interP.SAPCAcluster]) # by SAPCA cluster
k=8
interP.dist.cluster <- factoextra::hcut(interP.hist2.norm,k)
hist(interP.dist.cluster$cluster,breaks = k)
h <- heatmap(interP.hist2.norm,
            Colv = NA,
            Rowv = interP.dist.cluster$order,
            RowSideColors=colours[interP.dist.cluster$cluster],
            keep.dendro=TRUE,verbose=TRUE) # by dendrogram cluster

# strip.cluster <- SAPCA$seq.space.clusters$classification
# names(strip.cluster) <- gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names)))
# interP.SAPCAcluster <- strip.cluster[names(SEQF)]
h <- gplots::heatmap.2(interP.hist2.norm,
                      Colv = NA,
                      key  = 0,
                      trace="none",
                      keep.dendro=TRUE)
#find dendrogram clusters
interP.dendro.cluster <- factoextra::hcut(ape::cophenetic.phylo(ape::as.phylo(h$rowDendrogram)),k=7)$cluster
h <- gplots::heatmap.2(interP.hist2.norm,
                      Colv = NA,
                      key  = 0,
                      trace="none",
                      RowSideColors=colours[interP.dendro.cluster[rownames(interP.hist2.norm)]],
                      keep.dendro=TRUE)



interP.dist.cluster$cluster
# main interP 'clades'
k=10
interP.dist.cluster <- factoextra::hcut(interP.hist2.norm,k)
hist(interP.dist.cluster$cluster,breaks = k)

# N-glycosylation------
#> 1-fas FLAs -----
Ngly_sites_1 <- list(X=c(1530,1532),
                     A=c(1791,1822),
                     B=c(1852,1854),
                     C=c(1868,1871),
                     Ax=c(1838,1840),
                     D=c(1905,1908))
names(Ngly_sites_1) <- (gsub("/..*","",(names(Ngly_sites_1))))

Ngly_sites_1_list <- NULL
for(i in 1:length(Ngly_sites_1)){
  Ngly_sites_1_list <- cbind(Ngly_sites_1_list,
                             grepl(ignore.case = 1,
                                   "N[^P][TS]",
                                   gsub("-","",subseq(all_1_fas.aln,
                                                      Ngly_sites_1[[i]][1],
                                                      Ngly_sites_1[[i]][2]))))
}
colnames(Ngly_sites_1_list) <- names(Ngly_sites_1)
rownames(Ngly_sites_1_list) <- names(all_1_fas.aln)
Ngly_sites_1_list_b <- cbind("X" = Ngly_sites_1_list[,"X"],
                             "A" = rowSums(Ngly_sites_1_list[,c("A","Ax")]),
                             "B" = Ngly_sites_1_list[,"B"],
                             "C" = Ngly_sites_1_list[,"C"])
barplot(colMeans(Ngly_sites_1_list_b))
heatmap((Ngly_sites_1_list_b)*1,scale="none",Colv=NA,col=RColorBrewer::brewer.pal(9,"Blues"))
plot(GoodmanKruskal::GKtauDataframe(Ngly_sites_1_list_b),corrColors = "black")

#> 2-fas FLAs -----
Ngly_sites_2 <- list(X1=c(2372,2375),
                     A1=c(2607,2609),
                     B1=c(2643,2645),
                     B1x=c(2621,2623),
                     C1=c(2672,2674),
                     C1x=c(2628,2630),
                     D1=c(2692,2695),
                     X2=c(2872,2874),
                     X2x=c(2893,2895),
                     A2=c(3024,3030),
                     B2=c(3031,3033),
                     B2x=c(3056,3058),
                     C2x=c(3069,3071))
Ngly_sites_2_list <- NULL
for(i in 1:length(Ngly_sites_2)){
  Ngly_sites_2_list <- cbind(Ngly_sites_2_list,
                             grepl(ignore.case = 1,
                                   "N[^P][TS]",
                                   gsub("-","",subseq(all_2_fas.aln,
                                                      Ngly_sites_2[[i]][1],
                                                      Ngly_sites_2[[i]][2]))))
}
colnames(Ngly_sites_2_list) <- names(Ngly_sites_2)
rownames(Ngly_sites_2_list) <- names(all_2_fas.aln)
Ngly_sites_2_list_b <- cbind("X1" = Ngly_sites_2_list[,"X1"],
                             "A1" = Ngly_sites_2_list[,"A1"],
                             "B1" = Ngly_sites_2_list[,"B1"]|Ngly_sites_2_list[,"B1x"],
                             "C1" = Ngly_sites_2_list[,"C1"]|Ngly_sites_2_list[,"C1x"],
                             "D1" = Ngly_sites_2_list[,"D1"],
                             "X2" = Ngly_sites_2_list[,"X2"]|Ngly_sites_2_list[,"X2x"],
                             "A2" = Ngly_sites_2_list[,"A2"],
                             "B2" = Ngly_sites_2_list[,"B2"]|Ngly_sites_2_list[,"B2x"],
                             "C2" = Ngly_sites_2_list[,"C2x"])
                             
barplot(colMeans(Ngly_sites_2_list_b))
heatmap((Ngly_sites_2_list_b)*1,scale="none",Colv=NA,col=RColorBrewer::brewer.pal(9,"Blues"))
plot(GoodmanKruskal::GKtauDataframe(Ngly_sites_2_list_b),corrColors = "black")

#> subset by fas domain family -----
Ngly_sites_1_list_b_fas <- cbind(Ngly_sites_1_list_b,fasclass[rownames(Ngly_sites_1_list_b),])
Ngly_sites_1_list_b_fas <- Ngly_sites_1_list_b_fas[order(Ngly_sites_1_list_b_fas[,6]),]
Ngly_sites_1_list_b_fas_tab <- NULL
for (i in 1:ncol(Ngly_sites_1_list_b)){
  Ngly_sites_1_list_b_fas_tab <- cbind(Ngly_sites_1_list_b_fas_tab,
                                       table(Ngly_sites_1_list_b_fas[,6],Ngly_sites_1_list_b_fas[,i])[,2]/table(Ngly_sites_1_list_b_fas[,6]))
}
colnames(Ngly_sites_1_list_b_fas_tab) <- colnames(Ngly_sites_1_list_b)
heatmap(Ngly_sites_1_list_b_fas_tab[nrow(Ngly_sites_1_list_b_fas_tab):1,], Rowv=NA,Colv=NA,col=blues9,scale="none")

Ngly_sites_2_list_b_fas <- cbind(Ngly_sites_2_list_b,fasclass[rownames(Ngly_sites_2_list_b),])
Ngly_sites_2_list_b_fas <- Ngly_sites_2_list_b_fas[order(Ngly_sites_2_list_b_fas[,11]),]
Ngly_sites_2_list_b_fas_tab <- NULL
for (i in 1:ncol(Ngly_sites_2_list_b)){
  Ngly_sites_2_list_b_fas_tab <- cbind(Ngly_sites_2_list_b_fas_tab,
                                       table(Ngly_sites_2_list_b_fas[,11],Ngly_sites_2_list_b_fas[,i])[,2]/table(Ngly_sites_2_list_b_fas[,11]))
}
colnames(Ngly_sites_2_list_b_fas_tab) <- colnames(Ngly_sites_2_list_b)
heatmap(Ngly_sites_2_list_b_fas_tab[nrow(Ngly_sites_2_list_b_fas_tab):1,], Rowv=NA,Colv=NA,col=blues9,scale="none")

# 2d density plots -----

temp <- cbind(as.vector(dist(interP.hist2.norm)),
              as.vector(dist(interP.hist2.norm[,1:19])))

temp <- cbind(as.vector(dist(interP.hist2.norm[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names))),])),
              as.vector(dist(SAPCA$seq.space.PCA$coordinates[,1:10])))

temp <- cbind(as.vector(dist(k3[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names))),1:10])),
              as.vector(dist(k3[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names))),1:50])))

temp <- cbind(as.vector(dist(interP.hist2.norm[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names))),])),
              as.vector(dist(k3[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names))),1:100])))

temp <- cbind(unlist(agreg.numP),unlist(agreg.length))


colnames(temp) <-c("x","y")
temp2<-temp[sample(1:nrow(temp),100000),]
ggplot2::ggplot(data.frame(temp2), ggplot2::aes(x = x, y = y)) +
  ggplot2::geom_hex(bins = 50,ggplot2::aes(colour=..count..)) +
  ggplot2::geom_smooth(method='lm',formula=y~x,color="black",se=0) +
  # ggplot2::xlim(0, 0.1) +
  # ggplot2::ylim(0, 0.1) +
  ggplot2::theme_bw()
summary(lm(temp2[,1] ~ temp2[,2]))

# UMAP-------

#> agreg.k3 --------
agreg.k3.umap.x    <- umap::umap(agreg.k3,
                           n_neighbors  = 30,
                           n_components = 3)
agreg.k3.hdbscan <- dbscan::hdbscan(agreg.k3.umap$layout,minPts = nrow(agreg.k3.umap$layout)/90)
agreg.k3.hdbscan$cluster[agreg.k3.hdbscan$cluster==0] <- 24
names(agreg.k3.hdbscan$cluster) <- rownames(agreg.k3.umap$data)
agreg.k3.umap$groupmeans <- NULL
for (i in 1:max(agreg.k3.hdbscan$cluster)){
  agreg.k3.umap$groupmeans <- rbind (agreg.k3.umap$groupmeans,
                                      colMeans(agreg.k3.umap$layout[agreg.k3.hdbscan$cluster==i,]))
}
rownames(agreg.k3.umap.x$groupmeans) <- 1:max(agreg.k3.hdbscan$cluster)
rgl::plot3d(agreg.k3.umap$layout[,1:3],
            xlab="",ylab="",zlab="",
            type     = "s",
            specular = "black",   
            col      = agreg.k3.hdbscan$cluster+1,
            radius   = (agreg.k3.hdbscan$membership_prob)/10+0.05)
rgl::axes3d(color = "black", labels = FALSE)  
rgl::text3d(agreg.k3.umap$groupmeans,
            texts = paste0("---",letters[as.numeric(rownames(agreg.k3.umap$groupmeans))]),
            adj   = -0.1)
agreg.k3.hdbscan

agreg.k3.freq <- NULL
for (i in 1:max(agreg.k3.hdbscan$cluster)){
  agreg.k3.freq <- rbind(agreg.k3.freq,
                         colSums(agreg.k3[agreg.k3.hdbscan$cluster==i,]))
}
rownames(agreg.k3.freq) <- letters[1:max(agreg.k3.hdbscan$cluster)]
agreg.k3.freq.norm <- agreg.k3.freq/as.vector(table(agreg.k3.hdbscan$cluster)[1+1:max(agreg.k3.hdbscan$cluster)])
colSums(agreg.k3)/length(agreg.k3.hdbscan$cluster)
agreg.k3.freq.norm[agreg.k3.freq.norm==0] <- NA
heatmap(agreg.k3.freq.norm[max(agreg.k3.hdbscan$cluster):1,1:40], Rowv=NA)

#> nagreg.k3 --------
nagreg.k3.umap.x    <- umap::umap(nagreg.k3,
                               n_neighbors  = 30,
                               n_components = 3,
                               min_dist = )
nagreg.k3.hdbscan <- dbscan::hdbscan(nagreg.k3.umap$layout,minPts = nrow(nagreg.k3.umap$layout)/100)
nagreg.k3.hdbscan$cluster[nagreg.k3.hdbscan$cluster==0] <- 24
names(nagreg.k3.hdbscan$cluster) <- rownames(nagreg.k3.umap$data)
nagreg.k3.umap$groupmeans <- NULL
for (i in 1:max(nagreg.k3.hdbscan$cluster)){
  nagreg.k3.umap$groupmeans <- rbind (nagreg.k3.umap$groupmeans,
                                    colMeans(nagreg.k3.umap$layout[nagreg.k3.hdbscan$cluster==i,]))
}
rownames(nagreg.k3.umap$groupmeans) <- 1:max(nagreg.k3.hdbscan$cluster)
rgl::plot3d(nagreg.k3.umap$layout[,1:3],
            xlab="",ylab="",zlab="",
            type     = "s",
            specular = "black",   
            col      = nagreg.k3.hdbscan$cluster+1,
            radius   = (nagreg.k3.hdbscan$membership_prob)/25+0.02)
rgl::axes3d(color = "black", labels = FALSE)
rgl::text3d(nagreg.k3.umap$groupmeans,
            texts = paste0("---",letters[as.numeric(rownames(nagreg.k3.umap$groupmeans))]),
            adj   = -0.1)
nagreg.k3.hdbscan

plot(GoodmanKruskal::GKtauDataframe(cbind(agreg.k3.hdbscan$cluster,
                                          nagreg.k3.hdbscan$cluster)),
     corrColors = "black")

nagreg.k3.freq <- NULL
for (i in 1:max(nagreg.k3.hdbscan$cluster)){
  nagreg.k3.freq <- rbind(nagreg.k3.freq,
                         colSums(nagreg.k3[nagreg.k3.hdbscan$cluster==i,]))
}
rownames(nagreg.k3.freq) <- letters[1:max(nagreg.k3.hdbscan$cluster)]
nagreg.k3.freq.norm <- nagreg.k3.freq/as.vector(table(nagreg.k3.hdbscan$cluster)[1+1:max(nagreg.k3.hdbscan$cluster)])
colSums(nagreg.k3)/length(nagreg.k3.hdbscan$cluster)
nagreg.k3.freq.norm[nagreg.k3.freq.norm==0] <- NA
heatmap(nagreg.k3.freq.norm[max(nagreg.k3.hdbscan$cluster):1,1:50], Rowv=NA)

#> interP ------
interP.hist2.norm.umap    <- umap::umap(interP.hist2.norm[,1:21],
                                        n_neighbors  = 50,
                                        n_components = 3)
interP.hist2.norm.hdbscan <- dbscan::hdbscan(interP.hist2.norm.umap$layout,minPts = nrow(interP.hist2.norm.umap$layout)/50)
interP.hist2.norm.hdbscan$cluster[interP.hist2.norm.hdbscan$cluster==0] <- 24
names(interP.hist2.norm.hdbscan$cluster) <- rownames(interP.hist2.norm.umap$data)
rgl::plot3d(interP.hist2.norm.umap$layout[,1:3],
            xlab="",ylab="",zlab="",
            type     = "s",
            specular = "black",   
            col      = interP.hist2.norm.hdbscan$cluster+1,
            radius   = (interP.hist2.norm.hdbscan$membership_prob)/5+0.015)
rgl::axes3d(color = "black", labels = FALSE)
interP.hist2.norm.hdbscan

interP.pca <- prcomp(interP.hist2.norm)
rgl::plot3d(interP.pca$x[,1:3],
            type     = "s",
            specular = "black",   
            col      = interP.hist2.norm.hdbscan$cluster+1,
            radius   = (interP.hist2.norm.hdbscan$membership_prob)/50+0.002)


interP.hist2.norm.freq <- NULL
for (i in 1:max(interP.hist2.norm.hdbscan$cluster)){
  interP.hist2.norm.freq <- rbind(interP.hist2.norm.freq,
                          colSums(interP.hist2.norm[interP.hist2.norm.hdbscan$cluster==i,]))
}
rownames(interP.hist2.norm.freq) <- letters[1:max(interP.hist2.norm.hdbscan$cluster)]
interP.hist2.norm.freq.norm <- interP.hist2.norm.freq/as.vector(table(interP.hist2.norm.hdbscan$cluster)[1+1:max(interP.hist2.norm.hdbscan$cluster)])
heatmap(interP.hist2.norm.freq.norm[max(interP.hist2.norm.hdbscan$cluster):1,], Rowv=NA, Colv=NA)

#> Fas sequence ------
res.prop1 <- "https://raw.githubusercontent.com/TS404/DefSpace/master/data/Amino_acid_properties.csv"

SAPCA <- PCA_MSA (MSA        = file.full,
                  res.prop   = res.prop1,
                  cys        = 0,
                  # bootstrap  = 10,
                  clusterPCs = 1:50, 
                  clusters   = 1:20,
                  model      = "VVV")

fas.msa.umap    <- umap::umap(SAPCA$numerical.alignment$MSA.scale.wide,
                              n_neighbors  = 50,
                              n_components = 3)
fas.msa.umap.x  <- umap::umap(SAPCA$numerical.alignment$MSA.scale.wide,
                              n_neighbors  = 50,
                              min_dist     = 0.8,
                              n_components = 3)
fas.msa.hdbscan <- dbscan::hdbscan(fas.msa.umap$layout,minPts = nrow(fas.msa.umap$layout)/50)
names(fas.msa.hdbscan.retain$cluster) <- rownames(fas.msa.umap$data)
fas.msa.umap$groupmeans <- NULL
for (i in 1:max(fas.msa.hdbscan.retain$cluster)){
  fas.msa.umap$groupmeans <- rbind (fas.msa.umap$groupmeans,
                                    colMeans(fas.msa.umap$layout[fas.msa.hdbscan.retain$cluster==i,]))
}
rownames(fas.msa.umap$groupmeans) <- 1:max(fas.msa.hdbscan.retain$cluster)

fas.msa.umap.x$groupmeans <- NULL
for (i in 1:max(fas.msa.hdbscan.retain$cluster)){
  fas.msa.umap.x$groupmeans <- rbind (fas.msa.umap.x$groupmeans,
                                    colMeans(fas.msa.umap.x$layout[fas.msa.hdbscan.retain$cluster==i,]))
}
rownames(fas.msa.umap.x$groupmeans) <- 1:max(fas.msa.hdbscan.retain$cluster)


#plot
plot(table(fas.msa.hdbscan.retain$cluster))

rgl::plot3d(fas.msa.umap.x$layout[,1:3],
            type     = "s",
            specular = "black",   
            col      = fas.msa.hdbscan.retain$cluster+1,
            xlab = "", ylab = "", zlab = "",
            radius   = (fas.msa.hdbscan$membership_prob/5+0.1))
rgl::axes3d(color = "black", labels = FALSE)
rgl::text3d(fas.msa.umap.x$groupmeans,
            texts = paste0("---",LETTERS[as.numeric(rownames(fas.msa.umap$groupmeans))]),
            adj   = -0.1)

plot(GoodmanKruskal::GKtauDataframe	(cbind(fas.msa.hdbscan.retain$cluster,
                                           SAPCA$seq.space.clusters$classification)),
     corrColors = "black")
umap.tree <- ape::as.phylo(hclust(dist(fas.msa.umap.retain$groupmeans[,1:3])))
umap.tree$tip.label <- LETTERS[as.numeric(umap.tree$tip.label)] 
plot(umap.tree)


#> Consensus sequences -----

consensuses <- NULL
for(i in 1:max(fas.msa.hdbscan.retain$cluster)){
  consensuses <- rbind(consensuses,
                       seqinr::consensus(subset(SAPCA$numerical.alignment$MSA,
                                                subset=fas.msa.hdbscan.retain$cluster==i)))
}
rownames(consensuses) <- 1:max(fas.msa.hdbscan.retain$cluster)
consensuses.AA <- as.AAstringSet(consensuses)

# Save msas into folders
newdir <- "fasdomains alignments subset by fas umap"
dir.create(newdir,showWarnings = 0)
temp           <- SAPCA$numerical.alignment$MSA
rownames(temp) <- paste0("G",fas.msa.hdbscan.retain$cluster,"_",rownames(temp)) 
as.fasta(temp[order(fas.msa.hdbscan.retain$cluster),],
         write = paste0(newdir,"//umap_group.all.fa"))
for (i in sort(unique(fas.msa.hdbscan.retain$cluster))){
  as.fasta(temp[fas.msa.hdbscan.retain$cluster==i,],
           write = paste0(newdir,"//umap_group.",i,".fa"))
}
as.fasta(consensuses, write=paste0(newdir,"//consensues.fa"))

conservation.all <- apply(x      = temp,
	                      MARGIN = 2,
	                      FUN    = conservation)

conservation.group <- NULL
for (i in sort(unique(fas.msa.hdbscan.retain$cluster))){
  conservation.group <- rbind(conservation.group,
  	                          apply(x      = temp[fas.msa.hdbscan.retain$cluster==i,],
	                                MARGIN = 2,
	                                FUN    = conservation))
}

# Phylo of consensuses
fasdoms.msa <- phangorn::read.phyDat(type = "AA",format = "fasta",paste0(newdir,"//consensues.fa"))
MSA.dist <- phangorn::dist.ml(fasdoms.msa)
MSA.NJ   <- phangorn::NJ(MSA.dist)
MSA.fitNJ<- phangorn::pml(MSA.NJ, fasdoms.msa,
                          model="LG", optGamma=TRUE, k=4)
MSA.fitML<- phangorn::optim.pml(MSA.fitNJ, rearrangement = "NNI",
                                optInv=TRUE, optGamma=TRUE) #calculate ML tree (fit$tree)
MSA.ML   <- phangorn::bootstrap.pml(MSA.fitML, bs=100,optNni=TRUE) #calculate bootstraps (...optNni=TRUE)
MSA.ML.BS <- phangorn::plotBS(MSA.fitML$tree,MSA.ML, type= "phylogram") #combine main tree (from fit) with bootstraps(from pml)
hist(MSA.ML.BS$node.label,xlim=c(0,100))

# nagreg painted by agreg
rgl::plot3d(nagreg.k3.umap$layout[,1:3],
            xlab="",ylab="",zlab="",
            type     = "s",
            specular = "black",   
            col      = temp+1,
            radius   = (temp2)/10+0.02)
rgl::axes3d(color = "black", labels = FALSE)
rgl::text3d(nagreg.k3.umap$groupmeans,
            texts = paste0("---",letters[as.numeric(rownames(nagreg.k3.umap$groupmeans))]),
            adj   = -0.1)

# fas by agreg
rgl::plot3d(fas.msa.umap.retain$layout[,1:3],
            type     = "s",
            specular = "black",   
            col      = agreg.k3.hdbscan$cluster[as.character(labels$number)]+1,
            xlab = "", ylab = "", zlab = "",
            radius   = (agreg.k3.hdbscan$cluster[as.character(labels$number)]>0)*0.1+0.1)
rgl::axes3d(color = "black", labels = FALSE)
rgl::text3d(fas.msa.umap.retain$groupmeans,
            texts = paste0("---",LETTERS[as.numeric(rownames(fas.msa.umap.retain$groupmeans))]),
            adj   = -0.1)

# fas by nagreg
rgl::plot3d(fas.msa.umap.retain$layout[,1:3],
            type     = "s",
            specular = "black",   
            col      = nagreg.k3.hdbscan$cluster[as.character(labels$number)]+1,
            xlab = "", ylab = "", zlab = "",
            radius   = (nagreg.k3.hdbscan$cluster[as.character(labels$number)]>0)*0.1+0.1)
rgl::axes3d(color = "black", labels = FALSE)
rgl::text3d(fas.msa.umap.retain$groupmeans,
            texts = paste0("---",LETTERS[as.numeric(rownames(fas.msa.umap.retain$groupmeans))]),
            adj   = -0.1)

# fas by inter-P
rgl::plot3d(fas.msa.umap.retain$layout[,1:3],
            type     = "s",
            specular = "black",   
            col      = interP.hist2.norm.hdbscan$cluster[as.character(labels$number)]+1,
            xlab = "", ylab = "", zlab = "",
            radius   = (interP.hist2.norm.hdbscan$cluster[as.character(labels$number)]>0)*0.1+0.1)
rgl::axes3d(color = "black", labels = FALSE)

# inter-P by agreg
rgl::plot3d(interP.hist2.norm.umap$layout[,1:3],
            xlab="",ylab="",zlab="",
            type     = "s",
            specular = "black",   
            col      = agreg.k3.hdbscan$cluster+1,
            radius   = (agreg.k3.hdbscan$membership_prob)/15+0.015)
rgl::axes3d(color = "black", labels = FALSE)

# Seqspace --------------------------
#> Analysis -----------
SAPCA$groupmeans <- NULL
for (i in 1:max(SAPCA$seq.space.clusters$classification)){
  SAPCA$groupmeans <- rbind (SAPCA$groupmeans,
                             colMeans(SAPCA$seq.space.PCA$coordinates[SAPCA$seq.space.clusters$classification==i,1:3]))
}
rownames(SAPCA$groupmeans) <- 1:max(SAPCA$seq.space.clusters$classification)

SAPCA$groupmeans.umap <- NULL
for (i in 1:max(SAPCA$seq.space.clusters$classification)){
  SAPCA$groupmeans.umap <- rbind (SAPCA$groupmeans.umap,
                             colMeans(SAPCA$seq.space.PCA$coordinates[SAPCA$seq.space.clusters$classification==i,1:3]))
}
rownames(SAPCA$groupmeans.umap) <- 1:max(SAPCA$seq.space.clusters$classification)

# labels
Flanames2       <- gsub(paste0("X",FLAnames),replacement = ".", pattern = "[/]")
counter2        <- counter
names(counter2) <- Flanames2
counter2        <- counter2[SAPCA$numerical.alignment$seq.names]

FLAmax2         <- FLAmax
FLAmax2[is.na(FLAmax2)] <- 0
names(FLAmax2)  <- Flanames2
FLAmax2         <- FLAmax2[SAPCA$numerical.alignment$seq.names]

FLAcount2         <- counter
names(FLAcount2)  <- Flanames2

labels <- data.frame(row.names  = SAPCA$numerical.alignment$seq.names,
                     number     = as.numeric(as.character(gsub("X([0-9]*)\\..*",
                                                               "\\1",
                                                               SAPCA$numerical.alignment$seq.names))),
                     fas.count  = FLAcount2[SAPCA$numerical.alignment$seq.names],
                     fas.max    = FLAmax2[SAPCA$numerical.alignment$seq.names],
                     PCA.clust  = SAPCA$seq.space.clusters$classification,
                     umap.clust = fas.msa.hdbscan.retain$cluster)
labels.o <- labels  [order(labels$fas.count),]
labels.o <- labels.o[order(labels.o$number),]

# labels.o[labels.o==0] <- 24
fasclass.FLA <- NULL
for(i in unique(labels.o$number)){
  fasclass.FLA <- rbind(names(SEQF2),
                    c(FLA.number=i,
                      family=paste0(LETTERS[labels.o[labels.o$number==i,"umap.clust"]],collapse = "-")
                    ))
}
rownames(fasclass.FLA) <- fasclass.FLA[,1]



fasclass.fas <- (fasclass.FLA)[as.character(labels.o$number),]

agregclass  <- agreg.k3.hdbscan$cluster[fasclass.FLA[,1]]
nagregclass <- nagreg.k3.hdbscan$cluster[fasclass.FLA[,1]]
interPclass <- interP.hist2.norm.hdbscan$cluster[fasclass.FLA[,1]]
agregclass[agregclass==0]   <- 24 #24="X"
nagregclass[nagregclass==0] <- 24
interPclass[interPclass==0] <- 20
fasclass.FLA <- cbind(fasclass.FLA,
                  agreg  = letters[agregclass],
                  nagreg = letters[nagregclass],
                  interP = casefold(as.roman (interPclass),0))

#FAS-agreg co-occurence
col<-3
table(fasclass.FLA[,2],fasclass.FLA[,col])[order(table(fasclass.FLA[,2]),decreasing = 1),]
table(fas.msa.hdbscan.retain$cluster,fasclass.FLA[as.character(labels$number),col])

t <- cbind(fas=fas.msa.hdbscan.retain$cluster,
           a=agreg.k3.hdbscan$cluster[SAPCA.f2fmap[,1]],
           n=nagreg.k3.hdbscan$cluster[SAPCA.f2fmap[,1]],
           p=interP.hist2.norm.hdbscan$cluster[SAPCA.f2fmap[,1]])
t[t==0] <- NA
plot(GoodmanKruskal::GKtauDataframe(t,includeNA = "no"),corrColors = "black")

#agreg kmer by FLA family
topcols <- 1:30
toprows <- 1:20
agreg.k3.by.fas.umap <- NULL
temp <- names(sort(table(fasclass.FLA[,2]),decreasing = T))[toprows]
for (i in toprows){
  subset <- agreg.k3[fasclass.FLA[,2]==temp[i],]
  agreg.k3.by.fas.umap <- rbind(agreg.k3.by.fas.umap,
                                (colSums(subset[,topcols])/sum(subset)))
}
rownames(agreg.k3.by.fas.umap) <-temp
plot(agreg.k3.by.fas.umap[[3]], type="l", col="red", lwd=2)
heatmap(agreg.k3.by.fas.umap[20:1,agreg.k3heat$colInd],
        col=colorRampPalette(colors = c("red","white","blue"))(100),
        Colv=NA,Rowv=NA,
        scale='col')

#Arabidopsis FLA 1,3,4,9.11,12
fasclass.FLA[c("1211","1115","775","1080","1184","1216"),]
#Poplar FLA6
fasclass.FLA[c("8564","8605","9804","9601","9599"),]
#Brassica FLA4
fasclass.FLA[c("1765","2069","834"),]
#Cotton FLA1
fasclass.FLA[c("4758","5194"),]
#Zea Mays
fasclass.FLA[c("12937"),]

#TAXIDs
taxids <- data.frame(monocot       = 4447,
                     rosid         = 71275,
                     asterid       = 71274,
                     eudicot       = 71240,
                     spermatophyte = 58024,
                     embryophyte   = 3193)

taxids.fas <- NULL
for (i in 1:ncol(taxids)){
  taxids.fas <- cbind(taxids.fas,
                     rowSums(FLAtaxafull==taxids[,i],na.rm = T))
}
colnames(taxids.fas) <- colnames(taxids)
rownames(taxids.fas) <- rownames(FLAtaxafull)
table(apply(taxids.fas,FUN=paste0, collapse = "", MARGIN=1)) 
taxids.bin <- apply(taxids.fas,FUN=paste0, collapse = "", MARGIN=1) #taxid row combinations
taxids.fas.table <- table(list(taxids.bin,LETTERS[labels.o$umap.clust]))
rownames(taxids.fas.table) <- c("algae","bryophyte","basal spermatophyte","basal eudicot","asterid","rosid","monocot")
# taxids.fas.table <- taxids.fas.table/rowSums(taxids.fas.table)
taxids.fas.table
taxids.fas.table[taxids.fas.table==0]<-NA
heatmap(taxids.fas.table,Colv=NA,Rowv=NA,col=RColorBrewer::brewer.pal(9,"Greens"))

sort(table(fasclass.fas[taxids.bin=="000011",2]))
sort(table(grep("A",fasclass.fas[,2],value = T)))
table(taxids.bin[fasclass.fas[,2]=="C"])

fasclass.fas.u <- unique(cbind(fasclass.fas,taxids.bin))
taxids.FLA.table <- table(list(fasclass.fas.u[,2],fasclass.fas.u[,3]))[order(table(fasclass.fas.u[,2])),]
taxids.FLA.table[taxids.FLA.table==0] <- NA
heatmap(tail(taxids.FLA.table,60),Colv=NA,Rowv=NA,scale="column",col=RColorBrewer::brewer.pal(9,"Greens")) #NOTE MULTIPLIED BY FASnumber

#commelinid addendum
commelinid = rowSums(FLAtaxafull==4734,na.rm = T)
taxids.fas2 <- cbind (taxids.fas,noncom=!commelinid & taxids.fas[,1],com=commelinid)
taxids.bin2 <- apply(taxids.fas2,FUN=paste0, collapse = "", MARGIN=1) 
taxids.fas.table2 <- table(list(taxids.bin2,LETTERS[labels.o$umap.clust]))
rownames(taxids.fas.table2) <- c("algae","bryophyte","b. spermatophyte","b. eudicot","asterid","rosid",'com','non-com')
taxids.fas.table2[taxids.fas.table2==0]<-NA
heatmap(taxids.fas.table2,Colv=NA,Rowv=NA,col=RColorBrewer::brewer.pal(9,"Greens"))

fasclass.fas.u2 <- unique(cbind(fasclass.fas,taxids.bin2))
taxids.FLA.table2 <- table(list(fasclass.fas.u2[,2],fasclass.fas.u2[,3]))[order(table(fasclass.fas.u2[,2])),]
taxids.FLA.table2[taxids.FLA.table2==0] <- NA
heatmap(tail(taxids.FLA.table2,60),Colv=NA,Rowv=NA,scale="column",col=RColorBrewer::brewer.pal(9,"Greens")) #NOTE MULTIPLIED BY FASnumber

#domain co-occurence
barplot(tail(sort(table(fasclass.FLA[,2])),40),cex.axis=0.8,horiz=TRUE,las=1)
barplot(sort(table(t[,6])),cex.axis=0.8,horiz=TRUE,las=1)

cbind(FLAnamesfull[fasclass.FLA[stringr::str_length(fasclass.FLA[,2])>3,],],
      fasclass.FLA[stringr::str_length(fasclass.FLA[,2])>3,])[1:42,c(1,8,4,7,9)]->t
t[order((t[,5]),decreasing = 1),]

#clades vs nfas
temp <- table(labels.o[c(3,5)])
colnames(temp)<-LETTERS[as.numeric(colnames(temp))]
plot(t(temp[order(nrow(temp):1),]),col=colours,main="")
legend("right", as.character(c(7,5,4,3,2,1)), fill=colours[1:6], cex=0.8, title="n.fas")
barplot((temp[order(nrow(temp):1),order(ncol(temp):1)]),col=colours,main="",horiz=TRUE,las=1)
legend("right",  as.character(c(1,2,3,4,5,7)), fill=colours[6:1], cex=0.8, title="n.fas")

#clades vs genus
# temp <- table(cbind(FLAnamesfull[labels.o$number,2],LETTERS[labels.o[,5]])
temp <- table(data.frame(x=FLAtaxafull[,5],y=LETTERS[labels.o[,5]]))
temp <- table(data.frame(x=FLAtaxafull.r[,5],y=LETTERS[labels.o[,5]]))

heatmap(temp)

plot(fas.msa.hdbscan.retain,scale=1, gradient=c("black","grey40"))

# all vs all GoodmanKruskal

t1 <- cbind(fas=fas.msa.hdbscan.retain$cluster[intersect(rownames(Ngly_sites_1_list_b),SAPCA.f2fmap[,1])],
            a=agreg.k3.hdbscan$cluster[intersect(rownames(Ngly_sites_1_list_b),SAPCA.f2fmap[,1])],
            n=nagreg.k3.hdbscan$cluster[intersect(rownames(Ngly_sites_1_list_b),SAPCA.f2fmap[,1])],
            p=interP.hist2.norm.hdbscan$cluster[intersect(rownames(Ngly_sites_1_list_b),SAPCA.f2fmap[,1])])
t1[t1==0]<-NA
t1 <- cbind(Ngly_sites_1_list_b[intersect(rownames(Ngly_sites_1_list_b),SAPCA.f2fmap[,1]),]*1,
            t1)

t2 <-  cbind(fas=fas.msa.hdbscan.retain$cluster[intersect(rownames(Ngly_sites_2_list_b),SAPCA.f2fmap[,1])],
             a=agreg.k3.hdbscan$cluster[intersect(rownames(Ngly_sites_2_list_b),SAPCA.f2fmap[,1])],
             n=nagreg.k3.hdbscan$cluster[intersect(rownames(Ngly_sites_2_list_b),SAPCA.f2fmap[,1])],
             p=interP.hist2.norm.hdbscan$cluster[intersect(rownames(Ngly_sites_2_list_b),SAPCA.f2fmap[,1])])
t2[t2==0]<-NA
t2 <- cbind(Ngly_sites_2_list_b[intersect(rownames(Ngly_sites_2_list_b),SAPCA.f2fmap[,1]),]*1,
            t2)

plot(GoodmanKruskal::GKtauDataframe(t1,includeNA="no"),corrColors = "black")
plot(GoodmanKruskal::GKtauDataframe(t2,includeNA="no"),corrColors = "black")

#> Plot pca by cluster ------------------

plot_modelfit(SAPCA,legend = FALSE)
plot(table(SAPCA$seq.space.clusters$classification))
plot_3Dclusters(SAPCA,
                plotPCs    = 1:3)

rgl::par3d()->pp
rgl::par3d(pp)

#> Plot pca by n Fas domains ------------------

plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                # radius  = FLAmax2[SAPCA$numerical.alignment$seq.names]/3,
                col     = colorRampPalette(c("white", "red"))(3)[labels$fas.max])

#> plot pca by fas umap -------
PCA.by.UMAP <- NULL
for (i in 1:max(fas.msa.hdbscan.retain$cluster)){
  PCA.by.UMAP <- rbind (PCA.by.UMAP,
                        colMeans(SAPCA$seq.space.PCA$coordinates[fas.msa.hdbscan.retain$cluster==i,]))
}
rownames(PCA.by.UMAP) <- 1:nrow(PCA.by.UMAP)

plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                col     = 1+labels$umap.clust,
                radius  = (fas.msa.hdbscan.retain$membership_prob+0.1))
rgl::text3d(PCA.by.UMAP[,1:3],
            texts = paste0("-------",LETTERS[as.numeric(rownames(fas.msa.umap$groupmeans))]),
            adj   = -0.1)

# 1 = more than 2-fas FLA = black
# 2 = 1-fas FLA = white
# 3 = second fas of 2-fas FLA = Yellow
# 4 = first fas of 2-fas FLA = Red
colours<-palette(c("black",           
                   "white", 
                   "yellow",            
                   "red"))
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                col     = 1+1*(FLAcount2[SAPCA$numerical.alignment$seq.names]==1)
                           +2*(FLAmax2[SAPCA$numerical.alignment$seq.names]==2))

rgl::par3d(pp)

#> Plot pca by P number ------------------
hist(log10(numP.AGP), breaks = 40)
Pcolours <- c(colorRampPalette(c("white", "blue"))(2.5*mean(numP)),
              rep("blue",100))
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,0.5,
                col     = Pcolours[numP])

#> Plot pca by AG region data  ------------------
#> AG count
hist((agreg.count))
Pcolours <- c(colorRampPalette(c("white", "blue"))(2.5*mean(agreg.count)),
              rep("blue",100))
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,0.5,
                col     = Pcolours[((agreg.count[SAPCA.f2fmap[,1]])+1)])

#> AG length sum
hist(log10(agreg.length.sum),breaks=40)
Pcolours <- c(colorRampPalette(c("white", "blue"))(2.5*mean(agreg.length.sum)),
              rep("blue",100))
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,0.5,
                col     = Pcolours[((agreg.length.sum[SAPCA.f2fmap[,1]])+1)])

#> AG length average
hist(log10(agreg.length.av),breaks=40)
Pcolours <- c(colorRampPalette(c("white", "blue"))(2.5*mean(agreg.length.av)),
              rep("blue",100))
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,0.5,
                col     = Pcolours[((agreg.length.av[SAPCA.f2fmap[,1]])+1)])


#> Plot by interP distance ------------------
#>?\ by number of Prolines with interP distance P \
P=3
Pcolours <- c(colorRampPalette(c("white", "blue"))(2.5*mean(interP.hist2[,P+1])),
              rep("blue",100))
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,radius=0.5,
                col     = Pcolours[1+interP.hist2[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names))),(P+1)]])

#> by interP dendrogram cluster
k=8
interP.dist.cluster <- factoextra::hcut(interP.hist2.norm,k)
hist(interP.dist.cluster$cluster,breaks = k)
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                col     = blues9[interP.dist.cluster$cluster[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names)))]])
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                col     = colours[interP.dendro.cluster[gsub("X","",gsub("\\..*","",(SAPCA$numerical.alignment$seq.names)))]])

# hct clustering of is distances
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                col     = colours[factoextra::hcut(SAPCA$seq.space.PCA$coordinates[,1:40],20)$cluster])
# hct clustering of its distance dendrogram
seqspace.dendro.cluster <- factoextra::hcut(ape::cophenetic.phylo(ape::as.phylo(gplots::heatmap.2(SAPCA$seq.space.PCA$coordinates[,1:40],trace="none")$rowDendrogram)),k=20)$cluster
plot_3Dclusters(SAPCA,
                plotPCs = 1:3,
                col     = colours[seqspace.dendro.cluster])

#> Plot by kmer cluster ------------------

# Seqspace Conversions -----
fasdoms.phylo.f2fmap <- t(data.frame(strsplit(fasdoms.phylo$tip.label,"\\.")))
rownames(fasdoms.phylo.f2fmap) <- fasdoms.phylo$tip.label
SAPCA.f2fmap <- t(data.frame(strsplit(SAPCA$numerical.alignment$seq.names,"\\.")))
rownames(SAPCA.f2fmap) <- rownames(fas.msa.umap$data)

fasdoms.phylo.f2fmap <- gsub("X","",fasdoms.phylo.f2fmap)
SAPCA.f2fmap <- gsub("X","",SAPCA.f2fmap)[,1:4]

colnames(fasdoms.phylo.f2fmap) <- c("id","FLA","fas.n","fas.max")
colnames(SAPCA.f2fmap) <- c("id","FLA","fas.n","fas.max")

# add names to k2, k3 and interP
names(k2.hdbscan$cluster) <- rownames(agreg.k2)
names(k3.hdbscan$cluster) <- rownames(agreg.k3)
names(agreg.k3.hdbscan$cluster) <- rownames(agreg.k3)
names(interP.hist2.norm.hdbscan$cluster) <- rownames(interP.hist2.norm)

k2.SAPCAorder <- agreg.k2[SAPCA.f2fmap[,1],]
k3.SAPCAorder <- agreg.k3[SAPCA.f2fmap[,1],]
interP.hist2.norm.SAPCAorder <- interP.hist2.norm[SAPCA.f2fmap[,1],]
colnames(interP.hist2.norm.SAPCAorder) <- paste0("interP_",colnames(interP.hist2.norm.SAPCAorder))
agreg.SAPCAorder <- agreg[SAPCA.f2fmap[,1],]

k2.hdbscan.SAPCAorder <- k2.hdbscan$cluster[SAPCA.f2fmap[,1]] +1
k3.hdbscan.SAPCAorder <- k3.hdbscan$cluster[SAPCA.f2fmap[,1]] +1
agreg.k3.hdbscan.SAPCAorder <- agreg.k3.hdbscan$cluster[SAPCA.f2fmap[,1]] +1
interP.hist2.norm.hdbscan.SAPCAorder <- interP.hist2.norm.hdbscan$cluster[SAPCA.f2fmap[,1]]

# Phylogeny -----------------
fasdoms.msa <- phangorn::read.phyDat(type = "AA",format = "fasta",file.full)
MSA.dist    <- phangorn::dist.ml(fasdoms.msa)
MSA.NJ      <- phangorn::NJ(MSA.dist)
MSA.fitNJ   <- phangorn::pml(MSA.NJ, fasdoms.msa,
                          model="LG", optGamma=TRUE, k=4)
MSA.fitML   <- phangorn::optim.pml(MSA.fitNJ, rearrangement = "NNI",
                                optInv=TRUE, optGamma=TRUE) #calculate ML tree (fit$tree)
MSA.ML      <- phangorn::bootstrap.pml(MSA.fitML, bs=100,optNni=TRUE) #calculate bootstraps (...optNni=TRUE)
MSA.ML.BS   <- phangorn::plotBS(MSA.fitML$tree,MSA.ML, type= "phylogram") #combine main tree (from fit) with bootstraps(from pml)
hist(MSA.ML.BS$node.label,xlim=c(0,100))

#> by cluster ----------
cluster.MSAs <- cluster.MSAs.dist <- cluster.MSAs.NJ <- cluster.MSAs.fitNJ <- cluster.MSAs.fitML <- cluster.MSAs.ML <- cluster.MSAs.ML.BS <-  NULL
for (i in 1:18){
  cluster.MSAs[[i]] <- phangorn::read.phyDat(paste0(newdir,"//umap_group.",i,".fa"),
                                             type="AA",
                                             format="fasta")
  cluster.MSAs.dist[[i]] <- phangorn::dist.ml(cluster.MSAs[[i]])
  cluster.MSAs.NJ[[i]]   <- phangorn::NJ(cluster.MSAs.dist[[i]])
  cluster.MSAs.fitNJ[[i]]<- phangorn::pml(cluster.MSAs.NJ[[i]], cluster.MSAs[[i]],
                            model="LG", optGamma=TRUE, k=4)
  cluster.MSAs.fitML[[i]]<- phangorn::optim.pml(cluster.MSAs.fitNJ[[i]], rearrangement = "NNI",
                                  optInv=TRUE, optGamma=TRUE) #calculate ML tree (fit$tree)
  cluster.MSAs.ML[[i]]   <- phangorn::bootstrap.pml(cluster.MSAs.fitML[[i]], bs=100,optNni=TRUE) #calculate bootstraps (...optNni=TRUE)
  cluster.MSAs.ML.BS[[i]]<- phangorn::plotBS(cluster.MSAs.fitML[[i]]$tree,cluster.MSAs.ML[[i]], type= "phylogram") #combine main tree (from fit) with bootstraps(from pml)
  hist(cluster.MSAs.ML.BS[[i]]$node.label,xlim=c(0,100))
  ape::write.tree(cluster.MSAs.ML.BS[[i]],
                  paste0("cluster.",LETTERS[[i]],".nwk"))
}

#> iTol annotation -----

agreg.o  <- agreg.k3.hdbscan$cluster[as.character(labels.o$number)]
nagreg.o <- nagreg.k3.hdbscan$cluster[as.character(labels.o$number)]
interP.o <- interP.hist2.norm.hdbscan$cluster[as.character(labels.o$number)]

write.csv(cbind(labels.o,
                agreg.o,
                nagreg.o,
                interP.o,
                taxids.bin2),file = file.choose())

# domain ploits ---------------
sequences <- SEQF[1:100]

gpis <- NULL
pbt  <- txtProgressBar(min = 0, max = length(sequences), style = 3)
pbw  <- winProgressBar(min = 0, max = length(sequences), title = "HMM progress")

for(i in 1:length(sequences)){
  seqsubset <- sequences[i]
  gpis <- rbind(gpis, ragp::get_big_pi(sequence = seqsubset,
                                       id = names(seqsubset),
                                       verbose = FALSE,
                                       sleep = 0))
  setTxtProgressBar(pbt, i)
  setWinProgressBar(pbw, i, title= paste("GPI scan progress:",
                                         round(i/length(sequences)*100, 0),
                                         "%      (",
                                         names(sequences[i]),
                                         ")"
  ))
}
close(pbw)
rownames(gpis)<-gpis$id
signal.full <- ragp::get_phobius_file(file.full)
rownames(signal.full)<-signal.full$Name

temp <- (table(fasclass.fas[,2],gpis.fas$is.bigpi))
temp2 <- (temp[,2]/(temp[,1]+temp[,2]))
barplot(temp2[names(tail(sort(table(fasclass.FLA[,2])),21))],cex.axis=0.8,horiz=TRUE,las=1,xlim = 0:1)

subset<- 
a2 <- head(annotations,44)
s2 <- SEQF[unique(a2$id)[!is.na(unique(a2$id))]]
g2 <- gpis.FLA[unique(a2$id)[!is.na(unique(a2$id))],]
sig2 <- signal.full[unique(a2$id)[!is.na(unique(a2$id))],]

plot_domains(sequences   = s2,
             annotations = a2,
             gpis = g2,
             predict_nglc = FALSE)     

#> Species seqs -----
seqs.phys <- SEQF[FLAnames.FLA$X__1=="Zea  mays"]
seqs.phys <- SEQF[FLAnames.FLA$X__1=="Gossypium  hirsutum"]
seqs.phys <- SEQF[FLAnames.FLA$X__1=="Populus  trichocarpa"]
seqs.phys <- SEQF[FLAnames.FLA$X__1=="Eucalyptus  grandis"]
seqs.phys <- SEQF[FLAnames.FLA$X__1=="Brassica  rapa FPsc"]
seqs.phys <- seqs.phys[order(fasclass.FLA[as.character(names(seqs.phys)),2],decreasing = TRUE),]
plot_domains(sequences   = seqs.phys,
             annotations = annot,
             clades      = fasclass.FLA,
             signal      = signal.full,
             predict_nglc = 1,
             gpis        = gpis.FLA[names(seqs.phys),],
             #labels      = TRUE,
             y_lab       = FLAnamesfull[names(seqs.phys),]$id,
             dom_colour  = colours)
cbind(letters[agreg.k3.hdbscan$cluster[names(seqs.phys)]],
      letters[nagreg.k3.hdbscan$cluster[names(seqs.phys)]],
      letters[interP.hist2.norm.hdbscan$cluster[names(seqs.phys)]])

#> arabidopsis ------
seqs.ara.names <- data.frame(number=c("1042","1046","1080","1113","1115","1121","1125","1128","1129","2392","1146","1157","1179","1184","1186","1187","1202","1203","1211","1216"),
                             accession=c("AT2G04780","AT3G52370","AT1G03870","AT2G20520","AT2G24450","AT2G35860","AT2G45470","AT3G11700","AT3G12660","AT3G46550","AT3G60900","AT4G12730","AT4G31370","AT5G03170","AT5G06390","AT5G06920","AT5G40940","AT5G44130","AT5G55730","AT5G60490"),
                             name=c("FLA7","FLA15","FLA9","FLA6","FLA3","FLA16","FLA8","FLA18","FLA14","FLA4","FLA10","FLA2","FLA5","FLA11","FLA17","FLA21","FLA20","FLA13","FLA1","FLA12"))
seqs.ara.names <- seqs.ara.names[order(fasclass.FLA[as.character(seqs.ara.names[,1]),2],decreasing = TRUE),]
seqs.ara <- SEQF[seqs.ara.names[,1]]
plot_domains(sequences   = seqs.ara,
             annotations = annot,
             clades      = fasclass.FLA,
             signal      = signal.full,
             gpis        = gpis.FLA[names(seqs.ara),],
             labels      = TRUE,
             y_lab       = seqs.ara.names[,3],
             dom_colour  = colours)
cbind(letters[agreg.k3.hdbscan$cluster[names(seqs.ara)]],
      letters[nagreg.k3.hdbscan$cluster[names(seqs.ara)]],
      letters[interP.hist2.norm.hdbscan$cluster[names(seqs.ara)]])

#> charcterised seqs--------
seqs.chara <- SEQF[c("8564", "1765", "4758", "12937")]
seqs.chara <- seqs.chara[order(fasclass.FLA[as.character(names(seqs.chara)),2],decreasing = TRUE),]
plot_domains(sequences   = seqs.chara,
             annotations = annot,
             clades      = fasclass.FLA,
             signal      = signal.full,
             gpis        = gpis.FLA[names(seqs.chara),],
             labels      = TRUE,
             dom_colour  = colours)
cbind(letters[agreg.k3.hdbscan$cluster[names(seqs.chara)]],
      letters[nagreg.k3.hdbscan$cluster[names(seqs.chara)]],
      letters[interP.hist2.norm.hdbscan$cluster[names(seqs.chara)]])

#############.
# Functions #-----------------------------
#############.

#' Darken colour shade
isDark <- function(col,threshold=123) { (sum( col2rgb(col) * c(499, 587,114))/1000 < threshold) }

#' Largest objects
tail(sort( sapply(ls(),function(x){object.size(get(x))})),30)

#' Simple scan for N-glycosylation site mofis
predict_nglc <- function(sequences){
  
  nglc <- gregexpr('N[^P][ST]', sequences, ignore.case = FALSE)
  countnglc <- unlist(lapply(nglc,length))
  
  nglc2 <- data.frame(id          = rep(names(sequences),countnglc),
                      align_start = unlist(nglc))
  
  nglc2
}

#' Sequence conservation within alignment column
conservation <- function(x){
  x <- x[!is.na(x)]
  sort(table(x), decreasing=TRUE, na.rm=TRUE)[1]/length(x)
}


#' Generate domain diagram of domains and AG regions
#'
#' 
#'@param sequences   A list with elements of class "SeqFastaAA"
#'@param annotations A data frame of annotations produced by \code{\link[ragp]{get_hmm}}
#'@param gpis A data frame of annotations produced by \code{\link[ragp]{get_phobius_file}}
#'@param signal A data frame of annotations produced by \code{\link[ragp]{predict_hyp}}
#'@param predict_nglc Predict N-glycosylation sites
#'@param labels Print name labels over domains
#'@param height height of plot rectangles 
#'@return a ggplot2 visualisation of domain structure
#' 
#'@details 
#'
#'@references 
#'
#'@seealso \code{\link[ragp]{maab} \link[ragp]{scan_ag}}
#'
#'@examples
#'annotations <- get_hmm(sequence = sequences,
#'                       id = names(sequences),
#'                       verbose = FALSE,
#'                       sleep = 0.1)
#'gpis <- ragp::get_big_pi(sequence = sequences,
#'                         id = names(sequences),
#'                         verbose = FALSE,
#'                         sleep = 0.1)   
#'plot_domains(sequences   = sequences,
#'             annotations = annotations,
#'             gpis        = gpis) 
#'      
#'@export
#'
plot_domains <- function(sequences,
                         annotations,
                         clades         = NULL,
                         gpis           = NULL,    
                         signal         = NULL,
                         predict_nglc   = FALSE,
                         labels         = FALSE,
                         labelthreshold = 110,
                         y_lab          = TRUE,
                         dom_limits     = TRUE,
                         height         = 0.8,
                         dom_colour     = NULL){
  
  # Annotations ------
  # Order to match sequences
  
  annotations2 <- NULL
  for (name in names(sequences)){
    annotations2 <- rbind(annotations2,annotations[annotations$id==name,])
  }
  annotations2 <- annotations2[!is.na(annotations2$id),]
  annotations2 <- annotations2[!duplicated.array(annotations2),]
  annotations2[is.na(annotations2)] <- 0
  annotations  <- annotations2
  
  # Clades -------------
  if(!is.null(clades)){
    clades2        <- stringi::stri_split(clades[,2],regex = '-')
    names(clades2) <- rownames(clades)
    clades3        <- NULL
    for (name in unique(annotations$id)){
      clades3 <- append(clades3,unlist(clades2[names(clades2)==name]))
    }
    annotations$name[annotations$name=="Fasciclin"]<-clades3[clades3!=""]
  }
  
  # Domains ---------
  domcount <- NULL
  for(i in names(sequences)){
    domcount[i] <- sum(annotations$id==i, na.rm = 1)
  }
  annotations <- data.frame(y = rep(1:(length(sequences)),domcount),
                            annotations)
  annotations <- annotations[!is.na(annotations$name),]
  
  # Domains data to plot
  y1 <- annotations$y - height/2
  y2 <- annotations$y + height/2
  x1 <- annotations$align_start
  x2 <- annotations$align_end
  Domain  <- annotations$name
  
  d.domains <- data.frame(x1,x2,y1,y2,Domain)
  
  # Domains (full-length) -------
  dommax <- NULL
  for(i in unique(annotations$desc)){
    dommax[i] <- max(annotations$model_end[annotations$desc==i],na.rm = 1)
  }
  annotations$full_seq_length <- unlist(lapply(sequences[annotations$id],length))
  annotations$max_start <- annotations$align_start - annotations$model_start
  annotations$max_end   <- annotations$align_end   + (dommax[annotations$desc] - annotations$model_end)
  annotations$max_start[annotations$max_start <=0] <- 0
  annotations$max_end[annotations$max_end >= annotations$full_seq_length] <-
    annotations$full_seq_length[annotations$max_end >= annotations$full_seq_length]
  
  # Full-length domains data to plot
  y1 <- annotations$y - height/2
  y2 <- annotations$y + height/2
  x1 <- annotations$max_start
  x2 <- annotations$max_end
  Domain  <- annotations$name
  
  d.domains.max <- data.frame(x1,x2,y1,y2,Domain)
  
  # AGP regions  -----
  agregions  <- ragp::scan_ag(sequence = sequences,
                              id       = names(sequences),
                              dim      = 3,
                              div      = 6,
                              type     = "extended",
                              simplify = FALSE) $locations
  
  agregions2 <- matrix(unlist(lapply(agregions,t)), ncol = 2, byrow = TRUE)
  agregcount <- unlist(lapply(agregions, nrow))
  agregions3 <- data.frame(y           = rep(1:(length(sequences)),agregcount),
                           id          = rep(names(sequences),agregcount),
                           align_start = agregions2[,1],
                           align_end   = agregions2[,2])
  
  # AGP regions data to plot
  y1 <- agregions3$y - height/2
  y2 <- agregions3$y + height/2
  x1 <- agregions3$align_start
  x2 <- agregions3$align_end
  Domain <- rep("AG region",nrow(agregions3))
  
  d.agregions <- data.frame(x1,x2,y1,y2,Domain)
  
  ## Hyp predition ---------------------
  # stringi::stri_locate_all_regex(pattern = '[PAST]P|P[PAST]',
  #                                str=sequences$`133` [agregions$`133`[1,1]:agregions$`133`[1,2]])[[1]][,1] + 
  #   agregions$`133`[1,1]
  
  # hyp <- ragp::predict_hyp(sequence  = sequences,
  #                          id        = names(sequences),
  #                          tprob     = 0.3) $prediction
  
  # hyp <- hyp[hyp$HYP=="Yes",]
  # hyp <- hyp[!is.na(hyp$HYP),]
  # hypcount <- NULL
  # for(i in names(sequences)){
  #   hypcount[i] <- sum(hyp$id==i, na.rm = 1)
  # }
  # hyp <- data.frame(y = rep(1:(length(sequences)),hypcount),
  #                   hyp)
  #
  # # Hyp data to plot
  # y1 <- hyp$y - height/2
  # y2 <- hyp$y + height/2
  # x1 <- hyp$P_pos
  # x2 <- hyp$P_pos
  
  # y <- c(y1,y2)
  # x <- c(x1,x2)
  # P <- c(1:length(x1),1:length(x1))
  #
  # d.hyp=data.frame(x,y,P)
  
  ## GPI anchors ------
  if(!is.null(gpis)){
    gpi <- data.frame(y=1:nrow(gpis), gpis)
    gpi <- gpi[gpi$is.bigpi,]
    
    # GPI data to plot
    y <- gpi$y
    x <- as.numeric(gpi$omega_site)
    d.gpi=data.frame(x,y)
  }
  # N-Glycosylation ---------
  if (predict_nglc){
    nglc <- predict_nglc (sequences)
    nglccount <- NULL
    for(i in names(sequences)){
      nglccount[i] <- sum(nglc$id==i, na.rm = 1)
    }
    nglc <- data.frame(y = rep(1:(length(sequences)),nglccount),
                       nglc)
    
    # N-glc data to plot
    y <- nglc$y
    x <- as.numeric(nglc$align_start)
    d.nglc <- data.frame(x,y)
  }
  
  # Signal peptide -----------------
  
  if(!is.null(signal)){
    signal <- signal[names(sequences),]
    signal <- data.frame(y = 1:nrow(signal),
                         signal)
    
    # rownames(signal) <- names(sequences)
    signal <- signal[!is.na(signal$cut_site),]
    
    # Signal data to plot
    y1 <- signal$y - height/2
    y2 <- signal$y + height/2
    x1 <- rep(0,nrow(signal))
    x2 <- as.numeric(signal$cut_site)
    Domain  <- rep("Signal",nrow(signal))
    
    d.signal <- data.frame(x1,x2,y1,y2,Domain)
    rownames(d.signal) <- rownames(signal)
  }
  
  # Plot data ---------------------------
  # y labels
  if(all(y_lab==TRUE)){
    y_lab <- names(sequences)
  }
  if(all(y_lab==FALSE)){
    y_lab <- rep("",length(sequences))
  }
  
  # Backbones data to plot
  y1 <- 1:(length(sequences))
  y2 <- 1:(length(sequences))
  x1 <- rep(0,length(sequences))
  x2 <- nchar(as.character(sequences))
  
  y <- c(y1,y2)
  x <- c(x1,x2)
  
  d.backbone <- data.frame(x,y)
  
  # Generate plot with ggplot2
  p <- ggplot2::ggplot() 
  p <- p + ggplot2::scale_x_continuous(name="Length", breaks = seq(0,max(nchar(as.character(sequences))), by = 100))
  p <- p + ggplot2::scale_y_continuous(name="Protein", breaks = 1:length(sequences), labels=y_lab)                 
  p <- p + ggplot2::geom_line (d.backbone,  mapping=ggplot2::aes(x=x, y=y, group=y))                                               
  
  if(dom_limits){
    p <- p + ggplot2::geom_rect(d.domains.max, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Domain),
                                color="black", alpha = 0.3)    
  }
  
  p <- p + ggplot2::geom_rect (d.agregions, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2),fill="White", color="black")
  p <- p + ggplot2::geom_rect (d.domains,   mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=Domain), color="black")
  if(!is.null(dom_colour)){
    names(dom_colour) <- c("X",LETTERS[1:(length(dom_colour)-1)])
    p <- p + ggplot2::scale_fill_manual(values = dom_colour)
  }
  # p <- p + ggplot2::geom_line (d.hyp,       mapping=ggplot2::aes(x=x, y=y, group=P), color="#333333")
  p <- p + ggplot2::geom_rect (d.agregions, mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", alpha=0)
  if(!is.null(gpis)){
    p <- p + ggplot2::geom_point(d.gpi,            mapping=ggplot2::aes(x, y),shape = 21, colour = "black", fill = "yellow", size = height*2.5)
  }
  if(predict_nglc){
    p <- p + ggplot2::geom_point(d.nglc,           mapping=ggplot2::aes(x, y),shape = 21, colour = "black",  size = height*1.5)
  }
  
  if(!is.null(signal)){
    p <- p + ggplot2::geom_rect (d.signal,   mapping=ggplot2::aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="darkblue", color="black")
  }
  
  p <- p + ggplot2::theme_bw()                                                                                                
  p <- p + ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                          panel.grid.minor.y = ggplot2::element_blank(),
                          panel.grid.minor.x = ggplot2::element_blank())
  
  if(labels){
    darktextsize  <- rep(3,length(d.domains$Domain))
    darktextsize[stringi::stri_length(d.domains$Domain)>1] <- 0
    if(!is.null(dom_colour)){
      lighttextsize <- darktextsize
      lighttextsize[is.na(match(d.domains$Domain, names(dom_colour)[vapply(dom_colour,isDark,T,threshold=labelthreshold)]))] <- 0
      darktextsize[!is.na(match(d.domains$Domain, names(dom_colour)[vapply(dom_colour,isDark,T,threshold=labelthreshold)]))] <- 0
      p <- p + ggplot2::geom_text(data=d.domains, ggplot2::aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=d.domains$Domain,
                                                               fontface = "bold"), colour="white", size=lighttextsize)
    }
    p <- p + ggplot2::geom_text(data=d.domains, ggplot2::aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=d.domains$Domain,
                                                             fontface = "bold"), size=darktextsize)
  }
  
  p
} 
