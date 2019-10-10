#################################################################


#####################Extract_all#################################
Extract_all <- function(chromosome,id_all, all_vcf,gap,filter=T){
  counter = 0
  for( i in 1:length(chromosome)){

    if(counter==0){
      out <-  Extract_by_chr(id_all = id_all,all_vcf = all_vcf,chr = chromosome[i],filter=T,gap=gap)

    }else{
      out <-  cbind.data.frame(out,Extract_by_chr(id_all = id_all,all_vcf = all_vcf,chr = chromosome[i],gap=gap,filter=T))
    }
    counter <- counter +1
    cat(i,"\n")
  }
  return(out)
}
#################################################################



#####################Extract_by_chr##############################
Extract_by_chr <- function(id_all, all_vcf,chr,filter=T,gap=NULL){
  if(!require(data.table)){
    require(data.table)
  }
  start <- c()
  end <- c()
  for ( i in 1:length(id_all)){
    path_co <- paste0(all_vcf[i],"/",id_all[i],".genotype.",chr,".rough_COs.refined.breaks.txt")
    if(file.exists(path_co)){
      co <- data.frame(fread(path_co))
      if(filter){
        co <- filter_co(co = co,gap=gap)}
      colnames(co) <- c("V1","V2","V3","V4","V5")
      #check chr
      start <- unique(c(co$V3,start))
      end <- unique(c(co$V4,end))
    }
  }

  boarder <- sort(unique(c(start,end)))
  pos.start <- boarder[1:c(length(boarder) -1)]
  pos.end <- boarder[2:length(boarder)]
  pos <- (pos.start + pos.end)/2
  chromosome  <- rep(chr, length(pos))
  bin_name <- paste(chromosome,pos,sep = "-")
  # generate holder
  holder <- data.frame(array(data = NA,dim = c(length(id_all)+2,length(pos))))
  colnames(holder) <- bin_name
  rownames(holder) <- c(c("chr","pos"),id_all)
  holder[1,] <- chromosome
  holder[2,] <- pos

  for( i in 3:(length(id_all)+2)){
    # read in co file
    path_co <- paste0(all_vcf[i-2],"/",id_all[i-2],".genotype.",chr,".rough_COs.refined.breaks.txt")
    if(file.exists(path_co)){
      co <- data.frame(fread(path_co))
      if(filter)
        co <- filter_co(co = co,gap=gap)
      colnames(co) <- c("V1","V2","V3","V4","V5")

      for( j in 1:nrow(co)){
        id <- findInterval(holder[2,],c(co$V3[j],co$V4[j])) == 1
        holder[i,id] <- rep(co$V5[j],sum(id))
      }
    }
    #cat(i,"-",nrow(co),rep(co$V5[j],sum(id))[1],"-",length(rep(co$V5[j],sum(id))),"\n")

  }
  return(holder)
}
#####################################################################


############################ get.info ###############################
get.info <- function(chroms,chroms.len,bin.size=1e6){
  num.bin <- rep(NA,length(chroms))
  index <- data.frame(array(NA,dim = c(length(chroms),2)))
  colnames(index) <- c("start","end")
  rownames(index) <- chroms
  num.bin <- ceiling(chroms.len/bin.size)
  names(num.bin) <- chroms

  if(any(num.bin==1))
    warnings("a few chr only have one bin, start and end are set to the same","\n")

  index$start[1] <- 1
  index$end[1] <- num.bin[1]
  loca <- seq(from = 0.5,to =num.bin[1])
  loca.chr <- rep(chroms[1] , num.bin[1])
  for( i in 2:length(chroms)){
    index$start[i] <- sum(num.bin[1:c(i-1)])+1
    index$end[i] <- sum(num.bin[1:i])
    loca <- c(loca,seq(from = 0.5,to =num.bin[i]))
    loca.chr <- c(loca.chr,rep(chroms[i],num.bin[i]))
  }
  chr.loca <- rep(NA,max(index$end))
  for( i in 1:length(loca)){
    chr.loca[i] <- paste(loca.chr[i],"-",loca[i],sep = "")
  }
  return(list("num.bin" = num.bin,"index"=index,"loca"=loca,"loca.chr"=loca.chr,"chr.loca"=chr.loca))
}

#############################################################################


############################# wrap_get_density ##############################
wrap_get_density<- function(chr.ranks, binsize ,test){
  if(is.character(chr.ranks)){
    chr.ranks <- fread(chr.ranks, sep="\t")
    setnames(chr.ranks, old=c("V1","0","1","rank"), new=c("idx", "name","size","rank" ))
    chr.ranks <- chr.ranks[order(chr.ranks$rank)]
    chr.ranks <- data.frame(chr.ranks)
    }
  chroms <- chr.ranks$name
  chroms.len <- chr.ranks$size
  out.info <- get.info(chroms = chroms, chroms.len = chroms.len, bin.size = binsize)
  output.count <- data.frame(array(0, dim = c(1, max(out.info$index$end))))
  colnames(output.count) <- 1:ncol(output.count)
  #output.gap <-
  for( i in 1:length(names(test))){
    # chromsome in my test
    chr.now <- names(test)[i]
    count <- test[[i]]$density
    names(count) <- as.character(format(as.numeric(names(count))/binsize,digit=0))

    #output.count[,]

    loca.now <- out.info$loca[out.info$loca.chr == chr.now]
    index.now <- out.info$index[chr.now,]
    idx.now <- c(index.now$start:index.now$end)
    idx.now2 <- c(1:length(loca.now))
    mced.id <- match(as.numeric(names(count)),idx.now2)
    output.count[1,idx.now[mced.id]] <- count
  }
  #output.count[,which(output.count > 50)] <- 50
  return(output.count)
}

################################################################################

########################## Get_density_input ###################################

get_density_input <- function(inputfile,binsize,chr.ranks,cutoff=10, minsize=2000000){
  # compute the marker density on the input for Tiger
  # NB : This function is only implemented for contigs larger that 2 Mb!
  # takes:
  #          inputfile = Tiger input file (e.g. path/to/sample.genotype)
  #          binsize   = integer, size of the bins in basepairs for which the marker density will be computed.
  #          chr.ranks = .tsv or dataframe with index, name and size in basepairs and sizerank. here we use the one created by stripes_2/scripts/make_contig_names.py
  #          cutoff    = number of informative markers per bin that is considered acceptable.
  # returns:
  #          output.count =  ##TODO

  ## load requirements
  if(!require(data.table))
    require(data.table)
  input <- fread(inputfile)

  ## check if chr.rank is filename or data.frame, if filename, load as dataframe.
if(is.character(chr.ranks)){
  chr.ranks <- fread(chr.ranks, sep="\t")
  setnames(chr.ranks, old=c("V1","0","1","rank"), new=c("idx", "name","size","rank" ))
  chr.ranks <- chr.ranks[order(chr.ranks$rank)]
  chr.ranks <- data.frame(chr.ranks)
  }
  ## format chr.ranks
  ## substitute chromosome names in input file with numeric rank
  input$V1 <- chr.ranks$rank[match(input$V1,chr.ranks$name)]
  ## filter input for chromosome minimum size
  chr.all <- chr.ranks$rank[chr.ranks$size>minsize]
  input <- subset(input,subset = input$V1 %in% chr.all)


  ## calculate the number of mrker along the genome
  chr.num <- sort(unique(input$V1))
  out.put <- list()

  for( i in 1:length(chr.num)){  ##for each chromosome
    chrom <- chr.num[i]
    input.chr <- data.frame(subset(input,input$V1 %in% chrom))  # get all information on that chromosome from the input file / sample
    bins <- seq(1,chr.ranks$size[as.numeric(chr.ranks$rank) == chrom],by=binsize)  # create bins from 1 to chromosome size, with size binsize
    index <- findInterval(x = input.chr$V2,vec = bins)  # assign all sites to a bin
    index.pos <- seq(0,length(bins)+1,by=1)*binsize  #
    names(index.pos) <- seq(0,length(index.pos)-1,by=1)
    out <- aggregate(x = index,by=list(index),FUN=length)
    density <- out$x
    names(density) <- index.pos[as.character(out$Group.1)]
    if(length(density) <3){
      warning("too few windows\n")
      dis <- ""
    }else{
      dis <- as.numeric(names(density)[seq(2,length(density),by = 1)]) - as.numeric(names(density)[seq(1,length(density)-1,by = 1)] )

    }
    ratio <- sum(density > cutoff)/length(bins)
    out.put[[i]] <- list("density"=density,"gap"=dis,"coverage.chrom"=ratio)
  }
  names(out.put) <- chr.num

  names(out.put) <- chr.ranks$name[match(chr.num,chr.ranks$rank)]
  return(out.put)
}

##################################################################


################## filter_co #####################################

filter_co <- function(co,gap=10e6){
  if(is.null(gap))
    stop("gap must be provided to filter the genotype")

  colnames(co) <- c("sample","chr","start","end","genotype")
  co <- co[order(co$chr,co$start,decreasing = F),]
  chr <- unique(co$chr)
  #xo <- 0
  dis <- c()
  co.chr <- co
  line.filter <- c()
  if(nrow(co) <3){
    return(co)
  }else{
    for (j in 1:(nrow(co.chr)-2)){
      #Hom <- co.chr$genotype[j]==co.chr$genotype[j+2]  & co.chr$genotype[j] %in% c("CC" ,"LL") & co.chr$genotype[j+1] =="CL"
      #het <- co.chr$genotype[j]==co.chr$genotype[j+2]  & co.chr$genotype[j]  == "CL" & co.chr$genotype[j+1] %in% c("CC" ,"LL")
      dis <- co.chr$start[j+2] - co.chr$end[j]
      #if(Hom | het){
      if(dis < gap)
        line.filter <- c(line.filter,j+1)
      #}
    }

    if(length(line.filter)>=1){
      co.now <- co[-line.filter,]
      return(co.now)
    }else{
      return(co)
    }
  }
}

################################################################################

##################### Mat2geno #################################################

Mat2geno <- function(out=out_fam,id=idv.reliable,binsize=1e6,length.chr){

  chr  <- as.vector(t(out[1,]))
  pos <- as.vector(t(out[2,]))
  chr.lev <- as.character(sort(as.numeric(unique(chr))))
  out.idv <- out[id,]
  #Create out put holder
  #dime <- c(length(id),sum(length.chr[as.numeric(chr.lev)])*1e6/binsize)
  #output <- data.frame(array(NA,dime))
  for( i in 1:length(chr.lev)){
    index <- chr==chr.lev[i]
    out.idv.now <- out.idv[,index]
    pos.now <- pos[index]
    vec <- seq(0,length.chr[chr.lev[i]]*binsize,by=binsize)
    index.bin <- findInterval(pos.now,vec)
    # create holder
    if(i==1){
      #if(i==1 ){
      out.put <- data.frame(array(NA,c(length(id),length(unique(index.bin)))))
      colnames(out.put) <- paste0(chr.lev[i],"-",unique(index.bin))
      rownames(out.put) <- id
      for( j in 1:nrow(out.idv.now)){
        # fill in genotype
        # set to numeric
        geno <- gsub(pattern = "CC",replacement = 1,out.idv.now[j,])
        geno <- gsub(pattern = "LL",replacement = -1,geno)
        geno <- gsub(pattern = "CL",replacement = 0,geno)
        geno <- aggregate(as.numeric(geno),by=list("index"=index.bin),FUN= function(x) return(mean(x,na.rm=T)))
        rownames(geno) <- geno[,1]
        out.put[j,] <- geno[as.character(unique(index.bin)),2]
        #cat(i,"_",j,"\n")
      }

    }else{
      out.put1 <- data.frame(array(NA,c(length(id),length(unique(index.bin)))))
      colnames(out.put1) <- paste0(chr.lev[i],"-",unique(index.bin))
      rownames(out.put1) <- id
      for( j in 1:nrow(out.idv.now)){
        geno <- gsub(pattern = "CC",replacement = 1,out.idv.now[j,])
        geno <- gsub(pattern = "LL",replacement = -1,geno)
        geno <- gsub(pattern = "CL",replacement = 0,geno)
        geno <- aggregate(as.numeric(geno),by=list("index"=index.bin),FUN= function(x) return(mean(x,na.rm=T)))
        rownames(geno) <- geno[,1]
        out.put1[j,] <- geno[as.character(unique(index.bin)),2]
        #cat(i,"_",j,"\n")
      }
      out.put <- cbind(out.put,out.put1)
    }
  }
  return(out.put)
}
########################################################################################################

################################## Remove_het_mat ######################################################

Remove_het_mat <- function(mat.fam){
  mat.fam.re <- mat.fam
  for( i in 1:nrow(mat.fam)){
    index.fam  <- which(mat.fam.re[i,] < 1 & mat.fam.re[i,] >0)
    mat.fam.re[i,index.fam] <- rep(NA,length(index.fam))
    index.fam2  <- which(mat.fam.re[i,] < 0 & mat.fam.re[i,] >-1)
    mat.fam.re[i,index.fam2] <- rep(NA,length(index.fam2))
    index.fam3 <- which(is.nan(unlist(mat.fam.re[i,])))
    mat.fam.re[i,index.fam3] <- rep(NA,length(index.fam3))
    # # prob
    # index.prob <- which(mat.prob.re[i,] >= 0.8 )
    # mat.prob.re[i,index.prob] <- rep(1,length(index.prob))
    #
    # index.prob2 <- which(mat.prob.re[i,] <= -0.8 )
    # mat.prob.re[i,index.prob2] <- rep(-1,length(index.prob2))
    #
    # index.prob3<- which(mat.prob.re[i,] >= -0.4 & mat.prob.re[i,] <= 0.4 )
    # mat.prob.re[i,index.prob3] <- rep(0,length(index.prob3))
    #
    # index.prob4<- which(mat.prob.re[i,] > 0.4 & mat.prob.re[i,] < 0.8 )
    # mat.prob.re[i,index.prob4] <- rep(NA,length(index.prob4))
    #
    # index.prob5<- which(mat.prob.re[i,] > -0.8 & mat.prob.re[i,] < -0.4 )
    # mat.prob.re[i,index.prob5] <- rep(NA,length(index.prob5))
    # cat(i,"\n")
  }
  return(mat.fam.re)
}

#########################################################################################################


###################################### Export2qtl #######################################################
Export2rqtl <- function(genoFile,phenoFile,output.name=paste0("/Users/yanjunzan/Documents/impute/results/F2_833.within.fam._Tiger_1Mb_bins_cut_", 5,  "_",  Sys.Date())){
  # check if the id are the same
  if(!identical(rownames(genoFile),as.character(phenoFile$id)))
    stop("individual id coded in rownames(genoFile) and phenoFile$id must match")

  # recode autosome name
  #all.chr  <-matchingNames
  ## CAREFUL: Remove a.b from the end (after the -), carreful that the name still have this form
  chr <- as.numeric(sapply(strsplit(colnames(genoFile), split='-', fixed=TRUE), function(x) (x[1])))
  pos <- as.numeric(sapply(strsplit(colnames(genoFile), split='-', fixed=TRUE), function(x) (x[2])))

  #
  line1 <- c(colnames(phenoFile), colnames(genoFile))
  line2 <- c(rep("", length(colnames(phenoFile))), chr)
  test2 <- data.frame(array(NA, dim = c(nrow(genoFile)+2, length(line1))))
  test2[1, ] <- line1
  test2[2, ] <- line2
  test2[3:nrow(test2), 5:length(line1)] <- genoFile
  test2[3:nrow(test2), 1] <- rownames(genoFile)
  test2[3:nrow(test2),2:4] <- phenoFile[,2:4]
  write.table(test2, file = paste0(output.name,  ".csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  phe.name <- colnames(phenoFile)[!(colnames(phenoFile) %in% "id")]
  geno <- test2[,-c(2,3,4)] #tilman
  #print(geno[1,1]) # <- "id" # tilman
  geno[1,1] <- "id"
  #print(geno)
  write.table(geno, file = paste0(output.name,  "geno.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE) # tilman

  #write.table(test2[,-c(2,3,4)], file = paste0(output.name,  "geno.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = FALSE)
  #phenoFile[,1] <- NULL  ### addition by Tilman! somehow it generates an extra column named "V1"
  write.table(phenoFile, file = paste0(output.name,  "pheno.csv"), quote = FALSE, sep = ",", row.names = FALSE, col.names = T)

  #cat("file", paste0(output.name,  ".csv"),"generated","\n" )
}
