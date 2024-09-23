#Alicia Petrany
#Last Edit - 9-14-2024
#Bioinformatics Club Workshop data simulation

#This file simulates our fake genotype and rna-seq data!
#there are 20 total mice. 10 of them have wisker wobble syndrome (wws), and 10
#are normal mice (wt)

#first, I'm going to make the genotyping data. This dataset will ultimately have
#the following columns:
#pos - the position of a snp
#chrom - chromosomal a snp falls on
#gene - the gene the snp affects
#mice 1-20 - each mouse will have its own column, with each entry containing the
#            nucleotide each mouse has at a chromosomal location. The first 10
#            mice have wws, while the second 10 are WT.
#p.val - the probability that the snps at a given genomic coordinate are
#        are signficantly different between the wws and wt mice


#basic genotype parameters:
nmarkers = 10000 #number of surveyed sites
multi.allele.prob = 0.1 #probability for binomial dist to simulate number of allelles
n.mice.sd = 4 #for position with alternate alleles, this is the SD for the normal distribution simulating # of affected mice
assignment.bias.prob = 0.6 #see section where I simulate this for explanation
possible.alleles = c("A", "T", "C", "G", "H", "N")
n.wws = 10 #the number of wws mice
n.wt = 10 #the number of wt mice

#simulate number of markers with multiple alleles
multi.allele.dist <- rbinom(n = nmarkers,
                            size = length(possible.alleles),
                            prob = multi.allele.prob)

#simulate the number of mice with a certain snp along a normal dist
n.affected.mice <- floor(abs(rnorm(n = nmarkers,
                                   mean = 0,
                                   sd = n.mice.sd)))

#get number of mice in each condition group to have each allele
assignment.bias <- rbinom(n = nmarkers,
                         size = 1:20,
                         prob = assignment.bias.prob)/20
n.alt.affected <- round(n.affected.mice*assignment.bias)
n.main.affected <- n.affected.mice - n.alt.affected

n.alt.affected[n.alt.affected > 10] <- 10
n.main.affected[n.main.affected > 10] <- 10

#time to actually start assembling the dataframe. start with all mice having
#all the same allele at all positions
alleles <- sample(possible.alleles,
                  size = nmarkers,
                  replace = T,
                  prob = c(0.3, 0.3, 0.15, 0.15, 0.07, 0.03))

gdf <- matrix(data = rep(alleles, n.wws + n.wt),
               ncol = n.wws + n.wt,
               nrow = nmarkers)

scolnames <- paste0("wws.mouse.", c(1:n.wws))
nscolnames <- paste0("wt.mouse.", c(1:n.wt))
colnames(gdf) <- c(scolnames, nscolnames)
gdf <- as.data.frame(gdf)

#introduce alt alleles
for(i in 1:nrow(gdf)){
  n.alts <- multi.allele.dist[i]
  if(n.alts > 0){
    possible.alts <- possible.alleles[!(possible.alleles %in% gdf[i,])]
    alts <-sample(possible.alts, multi.allele.dist[i])
    alts.prob <- runif(length(alts))
    alts.prob <- alts.prob/sum(alts.prob)

    main <- sample(c("wws", "s"), 1)
    if(main == "s"){
      main.indices = c(1:n.wws)
      alt.indices = c((n.wws+1):(n.wws + n.wt))
    }else{
      main.indices = c((n.wws+1):(n.wws + n.wt))
      alt.indices = c(1:n.wws)
    }

    main.indices <- sample(main.indices, n.main.affected[i])
    alt.indices<- sample(alt.indices, n.alt.affected[i])

    if(length(main.indices) == 1){
      gdf[,main.indices][i] <- sample(alts,
                                       length(main.indices),
                                       prob = alts.prob,
                                       replace = T)
    }
    if(length(main.indices) > 1){
      gdf[,main.indices][i,] <- sample(alts,
                                       length(main.indices),
                                       prob = alts.prob,
                                       replace = T)
    }
    if(length(alt.indices) == 1){
      gdf[,alt.indices][i] <- sample(alts,
                                      length(alt.indices),
                                      prob = alts.prob,
                                      replace = T)
    }
    if(length(alt.indices) > 1){
      gdf[,alt.indices][i,] <- sample(alts,
                                       length(alt.indices),
                                       prob = alts.prob,
                                       replace = T)
    }
  }
}


#simulate other information about each position
genenames <- make.unique(replicate(nmarkers,
                                    paste0(sample(c(letters, 0:9),
                                                  sample(c(3:7), 1)),
                                           collapse = "")))

g.metadata <- data.frame(
  chrom = sample(c(1:19), nrow(gdf), replace = T),
  pos = sample(30000:100000, nrow(gdf)),
  gene = genenames)

#add in our gene of interest, plus genes it interacts with randomly
g.metadata$gene[1] <- "Wbbl1"
wbbl.interactors <- c("Fzlz3", "Sqk1", "Snf5",
                      "Jmpy2", "Nzl1", "Wvr",
                      "Twst", "Prpp", "Scth",
                      "Bll1", "Flkr2", "Snzz",
                      "Tgl", "Bnk5", "Sprt")
gdf[1,][,1:n.wws] <- "N"
gdf[1,][,(n.wws + 1):(n.wws + n.wt)] <- "A"
g.metadata$gene[2:(length(wbbl.interactors)+1)] <- wbbl.interactors

#peform fishers testing
fisher.testing <- function(x, n.wws, n.wt){
  letters <- as.character(x)
  REF <- names(which.max(table(letters)))
  bools <- letters == REF
  wws.ref <- length(which(bools[1:n.wws]))
  wt.ref <- length(which(bools[(n.wws + 1):length(x)]))
  cont_matrix <- matrix(c(wws.ref, n.wws-wws.ref,
                          wt.ref, n.wt-wt.ref),
                        nrow = 2, byrow = TRUE)
  return(fisher.test(cont_matrix)$p)
}

g.metadata$p.val <- as.numeric(apply(gdf,
                                     MARGIN = 1,
                                     FUN = fisher.testing,
                                     n.wws = n.wws,
                                     n.wt = n.wt))
gdf <- cbind(g.metadata, gdf)


#yay, thats done! now lets move onto the rna seq
#I'm using a simulation framework from a package I recently developed called DEGage
#there are three samples for each mouse
library(DEGage)

#generate counts
ndegs = 1000
nEEs = nmarkers - ndegs
lfcs <- runif(ndegs, -5, -1.5)
lfcs[1:(length(wbbl.interactors)+1)] <- runif((length(wbbl.interactors)+1), 5,7.5)
lfcs[1] <- 10
means = runif(nEEs + ndegs, 10, 100)
disps = runif(nEEs + ndegs, 0.1, 10)

rna.seq <- DEGage_Simulation(ngenes = nEEs, ndegs = ndegs,
                  ngroup1 = n.wws*3,
                  ngroup2 = n.wt*3,
                  lfc = lfcs,
                  min.prop.zeros = 0,
                  max.prop.zeros = 0,
                  means = means,
                  dispersions = disps)

rownames(rna.seq) <- make.unique(gdf$gene)

rna.colnames <- c()
for(lab in c(scolnames, nscolnames)){
  rna.colnames <- c(rna.colnames, paste0(lab, ".sample.", 1:3))
}
colnames(rna.seq) <- rna.colnames
#perform DE analysis
library(DESeq2)
condition <- factor(c(rep("wws", n.wws*3), rep("wt", n.wt*3)))
colData <- data.frame(row.names = colnames(rna.seq), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = rna.seq, colData = colData, design = ~ condition)
res <- DESeq(dds)
res <- results(res)

rna.met <- data.frame(res[,c(2, 6)])


#get some pca coords

pca <- prcomp(rna.seq[1:ndegs,])
percent.variance <-( pca$sdev^2 / sum(pca$sdev^2) * 100 )[1:2]
coords <- pca$rotation[,c(1,2)]
percent.variance

genders <- c()
for(i in 1:(n.wws + n.wt)){
  genders <- c(genders, rep(sample(c('m', "f"), 1), 3))
}

mouse.met <- data.frame(
  mouse = substr(rownames(coords), 1, nchar(rownames(coords))-9 ),
  sample = rownames(coords),
  pca1.859 = as.numeric(coords[,1]),
  pca2.520 = as.numeric(coords[,2]),
  gender = genders
)


#final organization and writing out files

gdf <- gdf[order(gdf$pos),]
gdf <- gdf[order(gdf$chrom), ]

rna.seq <- rna.seq[match(gdf$gene, rownames(rna.seq)),]
rna.met <- rna.met[match(gdf$gene, rownames(rna.met)),]

write.csv(gdf, "genotype_data.csv")
write.csv(rna.seq, "rnaseq_counts.csv")
write.csv(rna.met, "rnaseq_DE_results.csv")
write.csv(mouse.met, "mouse_metadata.csv")
