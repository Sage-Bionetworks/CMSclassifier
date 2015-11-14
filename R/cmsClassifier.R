# Rd
# description >> predict the Consensus Molecular Subtype (CMS) of colorectal cancer samples, based on log2_scaled GEP   
# argument
# item >> Exp >> a dataframe with log2_scaled Gene Expression Profiles (GEP) data values, samples in columns, genes in rows, rownames corresponding to Entrez IDs  
# item >> method >> a character vector, accepted values are: RF = random forest  predictor (here data will be automatically row-centered) ; SSP = single sample predictor , based on correlation to centroids (the data won't be row-centered)
# value >> a dataframe, samples in rows, columns : output of the predictor(s)
# author >> Aurelien de Reynies
# keyword >> internal
# end
classifyCMS <- function(Exp,method=c("RF","SSP")){
      res1<-res2<-NULL
      if("RF"  %in% method & ncol(Exp)>1)res1 <- classifyCMS.RF(Exp)  
      if("SSP" %in% method)res2 <- classifyCMS.SSP(Exp)
      L <- list("predictedCMS"=NULL,"nearestCMS"=NULL,"RF.details"=NULL,"SSP.details"=NULL)
      if(!is.null(res1)){
           L[[1]] <- res1[,"RF.predictedCMS",drop=F]
           L[[2]] <- res1[,"RF.nearestCMS",drop=F]
           L[[3]] <- res1
      }
      if(!is.null(res2)){
           if(is.null(L[[1]])){
               L[[1]] <- res2[,"SSP.predictedCMS",drop=F]
           }else{
               L[[1]] <- cbind(L[[1]],res2[,"SSP.predictedCMS",drop=F])
           }
           if(is.null(L[[2]])){
               L[[2]] <- res2[,"SSP.nearestCMS",drop=F]
           }else{
               L[[2]] <- cbind(L[[2]],res2[,"SSP.nearestCMS",drop=F])
           }
           L[[4]] <- res2
      }
      names(L[[1]]) <- gsub(".predictedCMS","",names(L[[1]]))
      names(L[[2]]) <- gsub(".nearestCMS","",names(L[[2]]))
      L <- L[which(!sapply(L,is.null))]
      L  
}


# Rd
# description >> Random Forest (RF) predictor of the Consensus Molecular Subtype (CMS) of colorectal cancer samples, based on log2_scaled GEP   
# argument
# item >> Exp >> a dataframe with log2_scaled Gene Expression Profiles (GEP) data values, samples in columns (n>1 mandatory), genes in rows, rownames corresponding to Entrez IDs  
# item >> center >>  boolean : should GEP be row-centered (default : TRUE)
# item >> minPosterior >>  numeric : minimal posterior probablity to classify a sample in a CMS
# value >> a dataframe, samples in rows, columns = posterior probability to be classified in each of the four CMS centroids, nearest CMS (ie CMS with highest posterior prob), predicted CMS 
# author >> Justin Guinney, Aurelien de Reynies
# keyword >> internal
# end
classifyCMS.RF <- function(Exp,center=TRUE,minPosterior=.5){
        data(model)
        Exp <- naImpute(Exp,G=listModelGenes("RF"))
        if(center) Exp <- Exp - rowMeans(Exp) 
        res <- as.data.frame(predict(finalModel, t(Exp),type="prob"))
        names(res) <- paste(names(res),".posteriorProb",sep="")
        res$predictedCMS <- res$nearestCMS <- apply(res,1,function(z){y<-which(z==max(z));paste(c("CMS1","CMS2","CMS3","CMS4")[y],collapse=",")})
        res$predictedCMS[which(apply(res[,1:4],1,max)<minPosterior)]<-NA
        names(res) <- paste("RF.",names(res),sep="")
        res
}



# Rd
# description >> nearest-centroid Single Sample Predictor (SSP) of the Consensus Molecular Subtype (CMS) of colorectal cancer samples, based on log2_scaled GEP   
# argument
# item >> Exp >> a dataframe with log2_scaled Gene Expression Profiles (GEP) data values, samples in columns, genes in rows, rownames corresponding to Entrez IDs  
# item >> minCor >> minimal correlation (sample x CMS centroid) for a sample to be classified in a CMS (NB : for a correlation below this threshold the sample will remain unclassified)
# item >> minDelta >> minimal difference between the correlation of a sample to its two nearest centroids, to be classified in a CMS (NB: for a delta of correlation below this threshold the sample will remain unclassified)
# value >> a dataframe, samples in rows, columns = correlation of the GEP to the four CMS centroids, nearest CMS, correlation to the nearest CMS, delta of correlation between the 2 nearest CMS, predicted CMS (NA value in case of uncertainty)
# author >> Aurelien de Reynies
# keyword >> internal
# end
classifyCMS.SSP <- function(Exp,minCor=.15,minDelta=.06){

              data(centroids)
              
              n <- ncol(centroids)
              G <- intersect(rownames(centroids),rownames(Exp))
              centroids <- centroids[G,]
              Exp <- Exp[G,,drop=F]

              nam <- sapply(strsplit(names(centroids),".",fixed=T),function(z)z[length(z)])

              tmp <-  as.data.frame(cor(Exp,centroids,use="pair"))
              names(tmp) <- paste("corTo",names(centroids),sep="")

              fun <- function(tmp,agreg.fun = min,minCor=.15,minDelta=.06,n=24){

                  tmpa <-as.data.frame(sapply(1:4,function(i)apply(tmp[,seq(i,n,by=4)],1,agreg.fun)))
                  names(tmpa) <- paste("corToCMS",1:4,sep="")
            
                  tmpa$nearestCMS <-  paste("CMS",1:4,sep="")[apply(tmpa,1,which.max)]
                  tmpa$corToNearest <- apply(tmpa[,1:4],1,max)
                  tmpa$deltaSecondNearestCMS <- sapply(1:nrow(tmpa),function(i) {
                                   w <- setdiff(paste("CMS",1:4,sep=""),tmpa[i,"nearestCMS"])
                                   tmpa[i,"corToNearest"] - max(tmpa[i,paste("corTo",w,sep="")]) })
                  tmpa$predictedCMS <- tmpa$nearestCMS
                  try(tmpa[ which(tmpa$corToNearest< minCor | tmpa$deltaSecondNearestCMS < minDelta),"predictedCMS"] <-  NA)
            
                  
                  tmpa
              }



             res = cbind ("min"= fun(tmp,agreg.fun =min,minCor=minCor,minDelta=minDelta,n=n),
                          "median"=fun(tmp,agreg.fun =median,minCor=minCor,minDelta=minDelta,n=n),
                          "max"=fun(tmp,agreg.fun =max,minCor=minCor,minDelta=minDelta,n=n))
              
             res <- as.data.frame(res)
              
             w1 <- which(is.na(res$min.predictedCMS))
             w2 <- which(is.na(res$min.predictedCMS) & !is.na(res$median.predictedCMS) & res$median.predictedCMS==res$min.nearestCMS & res$median.predictedCMS==res$max.nearestCMS)
             
             res$predictedCMS <- res$nearestCMS <- res$min.nearestCMS
             res$predictedCMS[w1] <- NA
             res$predictedCMS[w2] <- res$median.predictedCMS[w2]
             res <- res[,-grep("ToNearest",names(res))]
             res <- res[,-grep("delta",names(res))]
             res <- res[,-grep(".nearest",names(res))]
             res <- res[,-grep(".predicted",names(res))]
           
             names(res) <- paste("SSP.",names(res),sep="")
             res
}



# Rd
# description >> aggregate probe level GEP to gene level (entrez id)  GEP 
# argument
# item >> Exp >> a dataframe with log2_scaled GEP data values, samples in columns, genes in rows, rownames corresponding to Entrez IDs  
# item >> Gpl >> a dataframe giving the mapping between probes and Entrez gene Ids, probes correspond to rownames, Entrez Gene Ids are given in the column entrez
# item >> entrez >> a character string, the name of the column of the Entrez Gene Ids
# value >>  a dataframe with GEP data - 1 row per entrez ID 
# author >> Aurelien de Reynies
# keyword >> internal
# end
probesToEntrez <- function(Exp,Gpl,entrez="Entrez"){
        rn <- setdiff(intersect(rownames(Exp),rownames(Gpl)),rownames(Gpl)[which(Gpl[,entrez] %in% c("","---") | is.na(Gpl[,entrez]))])
        N <- length(rn)
        p <- ncol(Exp)
        data <- Exp[rn,,drop=F]
        partition <- Gpl[rn,entrez]
        nna <- function(z){as.numeric(!is.na(z))}
        res <- rowsum(data, partition,na.rm=T)  / rowsum(apply(data,2,nna), partition)
        res <- as.data.frame(apply(res,2,function(z){try(z[which(is.nan(z))]<-NA);z}))
        res
}


# Rd
# description >> genes used in the RF or SSP CMS predictor
# argument
# item >> method >> a character string : either RF (Random Forest predictor) or SSP (Single Sample Predictor)
# value >> a vector of character strings : entrez gene ids of the genes used in the RF or SSP predictor (depending on the chosen method)
# author >> Justin Guinney, Rodrigo Dienstmann, Aurelien de Reynies
# keyword >> internal
# end
listModelGenes <- function(method=c("RF","SSP")[1]){
  if(method=="RF"){
        data(model)
        return(mgenes)
  }
  if(method=="SSP"){
        data(centroids)
        return(rownames(centroids))
  }
  NULL
}

# Rd
# description >> impute NA values in a dataframe
# argument
# item >> Exp >>  a dataframe of GEP measures, samples in columns, genes in rows, rownames corresponding to Entrez IDs
# item >> G >> a vector of character strings (optionnal) : entrez gene ids that must be found in the rownames of Exp (will be added if absent, imputed value = global mean)  
# value >> the imputed datataframe
# author >> Aurelien de Reynies
# keyword >> internal
# end
naImpute <- function(Exp,G=NULL){ 
        m <- rowMeans(Exp,na.rm=T)
        mm <- mean(m,na.rm=T)
        m[which(is.na(m))] <- mm
        mna <- rowMeans(Exp)
        w <- which(is.na(mna))
        for(i in w){
            try(Exp[i,which(is.na(Exp[i,]))] <- m[i])
        }
        G <- setdiff(G,rownames(Exp))
        if(!is.null(G)) {
            tmp <- as.data.frame(array(mm,dim=c(length(G),ncol(Exp))))
            rownames(tmp) <- G
            names(tmp) <- names(Exp)
            colnames(tmp)<-colnames(Exp)##added by SJCG: rbind fails if Data Frame has sample names
            Exp <- rbind(Exp,tmp)
        }       
        Exp
}


