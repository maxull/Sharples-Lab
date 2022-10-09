################################################################
###                                                          ###
### ebGSEA functions (9.10.22)                               ###
###                                                          ###
################################################################


# https://github.com/aet21/ebGSEA/blob/master/ebGSEA/R/convertIDs.R


convertIDs <- function( ids, from, to, db, ifMultiple=c("putNA", "useFirst")) {
        stopifnot( inherits( db, "AnnotationDb" ) )
        ifMultiple <- match.arg( ifMultiple )
        suppressWarnings( selRes <- AnnotationDbi::select(
                db, keys=ids, keytype=from, columns=c(from,to) ) )
        
        if ( ifMultiple == "putNA" ) {
                duplicatedIds <- selRes[ duplicated( selRes[,1] ), 1 ]
                selRes <- selRes[ ! selRes[,1] %in% duplicatedIds, ]
        }
        
        return( selRes[ match( ids, selRes[,1] ), 2 ] )
}



#####################################################################################################################

# https://github.com/aet21/ebGSEA/blob/master/ebGSEA/R/doGSEAft.R

doGSEAft <- function(selEID.v,ptw.ls,allEID.v,ncores=4,minN=5,adjPVth=0.05){
        
        message(" Running Fisher's Exact Test... ")
        gseaFT.m <- matrix(unlist(mclapply(ptw.ls,gseaFTfn,selEID.v,allEID.v,minN,mc.cores=ncores)),ncol=5,byrow=TRUE)
        message(" Done ")
        colnames(gseaFT.m) <- c("nREP","nOVL","OR","P","Genes");
        rownames(gseaFT.m) <- names(ptw.ls);
        
        tmp.s <- sort(as.numeric(gseaFT.m[,4]),decreasing=FALSE,index.return=TRUE);
        sgseaFT.m <- gseaFT.m[tmp.s$ix,];
        padj.v <- p.adjust(as.numeric(sgseaFT.m[,4]),method="BH");
        sel.idx <- which(padj.v <= adjPVth);
        
        topGSEAft.m <- cbind(sgseaFT.m[sel.idx,1:4],padj.v[sel.idx],sgseaFT.m[sel.idx,5]);
        try(colnames(topGSEAft.m) <- c("nREP","nOVL","OR","P","adjP","Genes"))
        
        topGSEAft.lm <- list();
        topGSEAft.lm[[1]] <- topGSEAft.m;
        tmp.s <- sort(as.numeric(topGSEAft.m[,3]),decreasing=TRUE,index.return=TRUE);
        topGSEAft.lm[[2]] <- topGSEAft.m[tmp.s$ix,];
        names(topGSEAft.lm) <- c("Rank(P)","Rank(OR)");
        
        return(topGSEAft.lm);
}

gseaFTfn <- function(termEID.v,selEID.v,allEID.v,minN=5){
        
        enrEID.v <- as.character(intersect(termEID.v,selEID.v));
        presEID.v <- intersect(termEID.v,allEID.v);
        nrep <- length(enrEID.v);
        if(nrep >= minN){
                #print(enrEID.v);
                tmp.m <- matrix(nrow=2,ncol=2);
                tmp.m[1,1] <- nrep;
                tmp.m[1,2] <- length(selEID.v)-nrep;
                tmp.m[2,1] <- length(intersect(allEID.v,termEID.v))-nrep;
                tmp.m[2,2] <- length(allEID.v) - tmp.m[1,1] - tmp.m[1,2] - tmp.m[2,1];
                ft.o <- fisher.test(tmp.m,alt="gr");
                
                mapSYM.v <- convertIDs(enrEID.v, "ENTREZID", "SYMBOL", org.Hs.eg.db,ifMultiple="useFirst");
                #      print(mapSYM.v);
                out.v <- c(sum(tmp.m[,1]),nrep,ft.o$est,ft.o$p.value,PasteVector(mapSYM.v));
        }else {
                out.v <- c(length(presEID.v),nrep,1,1,"NoGenes");
        }
        return(out.v);
}

PasteVector <- function(v){
        
        vt <- v[1];
        if(length(v) > 1){
                for(g in 2:length(v)){
                        vt <- paste(vt,v[g],sep=" ")
                        
                }
        }
        vt <- paste(vt," EnD",sep="");
        out.v <- sub(" EnD","",vt);
        out.v <- sub("NA , ","",out.v);
        out.v <- sub(" , NA","",out.v);
        out.v <- sub(" , NA , "," , ",out.v);
        return(out.v);
}


#####################################################################################################################


### https://github.com/aet21/ebGSEA/blob/master/ebGSEA/R/doGSEAwt.R

doGSEAwt <- function(rankEID.m,ptw.ls,ncores=4,minN=5,adjPVth=0.05){
        rankEID.v<-rownames(rankEID.m)
        message(" Running Wilcox Test and Known Population Median Test... ")
        gseaWT.m <- matrix(unlist(mclapply(ptw.ls,gseaWTfn,rankEID.v,mc.cores=ncores,minN=minN)),ncol=4,byrow=TRUE)
        message(" Done ")
        colnames(gseaWT.m) <- c("nREP","AUC","P(WT)","P(KPMT)");
        rownames(gseaWT.m) <- names(ptw.ls);
        
        tmp.s <- sort(gseaWT.m[,3],decreasing=FALSE,index.return=TRUE);
        sgseaWT.m <- gseaWT.m[tmp.s$ix,];
        padj.v <- p.adjust(sgseaWT.m[,3],method="BH");
        
        sel.idx <- which(padj.v <= adjPVth);
        topGSEAwt.lm <- list();
        sym.v<-convertIDs(rankEID.v, 'ENTREZID', 'SYMBOL', org.Hs.eg.db, ifMultiple="useFirst")
        del.idx<-which(is.na(sym.v))
        sym.v<-sym.v[-del.idx]
        rankEID.m<-rankEID.m[-del.idx,]
        rankEID.v<-rankEID.v[-del.idx]
        if(length(sel.idx)>1){
                
                topGSEAwt.m <- cbind(sgseaWT.m[sel.idx,],padj.v[sel.idx]);
                colnames(topGSEAwt.m) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");
                
                topGSEAwt.lm[[1]] <- topGSEAwt.m;
                tmp.s <- sort(topGSEAwt.m[,2],decreasing=TRUE,index.return=TRUE);
                topGSEAwt.lm[[2]] <- topGSEAwt.m[tmp.s$ix,];
                topGSEAwt.lm[[3]] <- list()
                for (i in 1:nrow(topGSEAwt.m)){
                        EID.v<-intersect(ptw.ls[[match(rownames(topGSEAwt.m)[i],names(ptw.ls))]],rankEID.v)
                        pathgene.m<-matrix(NA,nrow=length(EID.v),ncol=2)
                        rownames(pathgene.m)<-sym.v[match(EID.v,rankEID.v)]
                        colnames(pathgene.m)<-c('Pvalue','Statistic')
                        pathgene.m[,1:2]<-rankEID.m[match(EID.v,rankEID.v),1:2]
                        topGSEAwt.lm[[3]][[i]]<-pathgene.m
                }
                names(topGSEAwt.lm[[3]])<-rownames(topGSEAwt.m)
                names(topGSEAwt.lm) <- c("Rank(P)","Rank(AUC)","Genestat");
        }
        if (length(sel.idx)==1) {
                topGSEAwt.v <- as.vector(c(sgseaWT.m[sel.idx,],padj.v[sel.idx]));
                names(topGSEAwt.v) <- c("nREP","AUC","P(WT)","P(KPMT)","adjP");
                EID.v<-intersect(ptw.ls[[match(rownames(sgseaWT.m)[sel.idx],names(ptw.ls))]],rankEID.v)
                pathgene.m<-matrix(NA,nrow=length(EID.v),ncol=2)
                rownames(pathgene.m)<-sym.v[match(EID.v,rankEID.v)]
                colnames(pathgene.m)<-c('Pvalue','Statistic')
                pathgene.m[,1:2]<-rankEID.m[match(EID.v,rankEID.v),1:2]
                topGSEAwt.lm <- list("Rank(P)"=topGSEAwt.v,"Rank(AUC)"=topGSEAwt.v,"POI"=rownames(sgseaWT.m)[sel.idx],
                                     "Genestat"=pathgene.m);
        }
        
        return(topGSEAwt.lm);
}



gseaWTfn <- function(termEID.v,rankEID.v,minN=5){
        commonEID.v <- intersect(termEID.v,rankEID.v);
        nrep <- length(commonEID.v);
        if(length(commonEID.v)>=minN){
                otherEID.v <- setdiff(rankEID.v,termEID.v);
                match(commonEID.v,rankEID.v) -> rank1.idx;
                match(otherEID.v,rankEID.v) -> rank2.idx;
                wt.o <- wilcox.test(rank1.idx,rank2.idx,alt="less");
                pv <- wt.o$p.value;
                n1 <- length(rank1.idx);
                n2 <- length(rank2.idx);
                auc <- 1 - wt.o$stat/(n1*n2);
                ### now do kpmt
                pop.v <- 1:length(rankEID.v);
                names(pop.v) <- rankEID.v;
                obs.v <- commonEID.v;
                pvKPMT <- kpmt(pop=pop.v,obs=obs.v,tail="lower")[[4]];
                out.v <- c(nrep,auc,pv,pvKPMT);
        }else {
                out.v <- c(nrep,0,1,1);
        }
        return(out.v);
}


#####################################################################################################################

# ://github.com/aet21/ebGSEA/blob/master/ebGSEA/R/doGT.R


doGT <- function(pheno.v,data.m,model=c("linear"),array=c("450k","850k"),ncores=4){
        if(array=="450k"){
                message(" Mapping 450k probes to genes... ")
                data("dualmap450kEID");
                subsets <- mclapply(mapEIDto450k.lv,intersect,rownames(data.m),mc.cores = ncores);
                message(" Done ")
        }else {
                message(" Mapping EPIC probes to genes... ")
                data("dualmap850kEID");
                subsets <- mclapply(mapEIDto850k.lv,intersect,rownames(data.m),mc.cores = ncores);
                message(" Done ")
        }
        nrep.v <- unlist(lapply(subsets,length));
        selG.idx <- which(nrep.v>0);
        message(" Running Global Test... ")
        gt.o <- gt(response=pheno.v,alternative=t(data.m),model=model,directional = FALSE, standardize = FALSE, permutations = 0, subsets=subsets[selG.idx],trace=F);
        message(" Done ")
        resGT.m <- as.matrix(result(gt.o));
        tmp.s <- sort(resGT.m[,2],decreasing=TRUE,index.return=TRUE);
        sresGT.m <- resGT.m[tmp.s$ix,];
        return(sresGT.m);
}


#####################################################################################################################

### https://github.com/aet21/ebGSEA/blob/master/ebGSEA/R/selEIDfromSelCpG.R

selEIDfromSelCpG <- function(selCpG.v,allCpG.v,pvth=0.3/length(selCpG.v),array=c("450k","850k")){
        
        if(array=="450k"){
                message(" Mapping 450k probes to genes... ")
                data("dualmap450kEID");
                unqEID.v <- unique(unlist(map450ktoEID.lv[match(selCpG.v,names(map450ktoEID.lv))]));
                tmp.idx <- match(unqEID.v,names(mapEIDto450k.lv));
                tmp.l <- mapEIDto450k.lv[tmp.idx];
                message(" Done ")
        }else {
                message(" Mapping EPIC probes to genes... ")
                data("dualmap850kEID");
                unqEID.v <- unique(unlist(map850ktoEID.lv[match(selCpG.v,names(map850ktoEID.lv))]));
                tmp.idx <- match(unqEID.v,names(mapEIDto850k.lv));
                tmp.l <- mapEIDto850k.lv[tmp.idx];
                message(" Done ")
        }
        
        ncpg.v <- unlist(lapply(tmp.l,length));
        obsN.v <- unlist(lapply(lapply(tmp.l,intersect,selCpG.v),length));
        prob <- length(selCpG.v)/length(allCpG.v);
        
        pv.v <- pbinom(obsN.v,size=ncpg.v,prob,lower.tail=FALSE);
        names(pv.v) <- unqEID.v;
        sel.idx <- which(pv.v < pvth);
        selEID.v <- unqEID.v[sel.idx];
        return(list(selEID=selEID.v,selPV=pv.v[sel.idx],allPV=pv.v));
}









