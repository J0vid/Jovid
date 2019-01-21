#' read TPS files
#'
#' @param data A .TPS file
#' @return A matrix of the landmarks for each observation
#' @export
read.tps <- function(data) {
  # Reads the .tps file format produced by TPSDIG
  # (http://life.bio.sunysb.edu/morph/ into a single data frame
  # USAGE: R> read.tps("filename.tps")
  a = readLines(data) # so we can do some searching and indexing
  LM = grep("LM", a) # find the line numbers for LM
  ID.ind = grep("ID", a) # find the line numbers for ID
  # and the ID values, SCALE values, and image names
  ID = gsub("(ID=)(.*)", "\\2", grep("ID", a, value=T))
  SCALE = gsub("(SCALE=)(.*)", "\\2", grep("SCALE", a, value=T))
  images = basename(gsub("(IMAGE=)(.*)", "\\2", a[ID.ind - 1]))
  # FOR EACH LOOP
  skip = LM # set how many lines to skip
  # and how many rows to read
  nrows = as.numeric(gsub("(LM=)(.*)", "\\2", grep("LM", a, value=T)))
  l = length(LM) # number of loops we want

  landmarks = vector("list", l) # create an empty list

  for (i in 1:l) {
    landmarks[i] = list(data.frame(
      read.table(file=data, header=F, skip=LM[i],
                 nrows=nrows[i], col.names=c("X", "Y")),
      IMAGE = as.character(images[i]),
      ID = ID[i],
      SCALE = SCALE[i]))
  }
  do.call(rbind, landmarks) # rbind the list items into a data.frame
}


#' Calculates procrustes distance for row format
#'
#' @param data A set of data you want to get distances for
#' @param target some reference you want every observation to be measured against
#' @return Procrustes distances
#' @export
procdist= function(data, target){

  pdist= rep(0, length(data[,1]))

  for(i in 1:length(data[,1])){

    pdist[i]=	sqrt(sum((as.vector(as.numeric(data[i,])) - as.vector(as.numeric(target)))^2))

  }
  return(pdist)
}


#' Calculates procrustes distance for array format
#'
#' @param data A set of data you want to get distances for
#' @param reference some reference you want every observation to be measured against
#' @return Procrustes distances
#' @export
procdist.array= function(data, reference){

  pdist= rep(0, length(data[1,1,]))

  for(i in 1:length(data[1,1,])){

    pdist[i]=  sqrt(sum((as.matrix(as.numeric(data[,,i])) - as.matrix(as.numeric(reference)))^2))

  }
  return(pdist)

}

#' Converts tps imports to array format
#'
#' @return An array
#' @export

#####tps2array
tps2array= function(tpsfile, arrayname= NA){

  Nlandmarks= sum(as.numeric(tpsfile$ID)== 1)
  ID.nums = rep(0, length(tpsfile$ID))
  coord.array = array(dim= c(Nlandmarks,2, length(ID.nums)/Nlandmarks), dimnames= arrayname)
  #doesn't work with repeats ID.nums = unique(as.vector(tpsfile$ID))

  for(i in 0:(length(tpsfile$ID)/Nlandmarks)){
    ID.nums[((i*Nlandmarks) + 1):(((i+1)*Nlandmarks) + 1)]= rep(i, length(Nlandmarks))
  }


  for(ind in 1:length(coord.array[1,1,])) {
    coord.array[,,ind] = as.matrix(tpsfile[ID.nums== unique(ID.nums)[ind] ,1:2])
  }
  return(coord.array)
}


#' Converts row formats to array format for 2D data
#'
#' @return An array
#' @export
row2array= function(data, Nlandmarks= length(data[1,])/2){
  if(is.matrix(data)| is.data.frame(data)== T){
    xseq= seq(1,Nlandmarks*2, 2)
    yseq= seq(2,Nlandmarks*2, 2)
    a= array(0, dim= c(Nlandmarks, 2, length(data[,1])))

    for(i in 1:length(data[,1])){

      a[,,i]= as.matrix(cbind(data[i, xseq], data[i, yseq]))

    } }

  else {
    if(is.matrix(data)== F){
      xseq= seq(1,Nlandmarks*2, 2)
      yseq= seq(2,Nlandmarks*2, 2)

      a= as.matrix(cbind(as.matrix(data)[1, xseq], as.matrix(data)[1, yseq]))

    }
  }
  return(a)
}#end function

#' Converts row formats to array format for 3D data
#'
#' @return An array
#' @export
row2array3d= function(data, Nlandmarks= length(data[1,])/3){
  if(is.matrix(data)| is.data.frame(data)== T){
    xseq= seq(1,Nlandmarks*3, 3)
    yseq= seq(2,Nlandmarks*3, 3)
    zseq = seq(3,Nlandmarks*3, 3)
    a= array(0, dim= c(Nlandmarks, 3, length(data[,1])))

    for(i in 1:length(data[,1])){

      a[,,i]= as.matrix(cbind(data[i, xseq], data[i, yseq], data[i, zseq]))

    } }

  else {
    if(is.matrix(data)== F){
      xseq= seq(1,Nlandmarks*3, 3)
      yseq= seq(2,Nlandmarks*3, 3)
      zseq = seq(3,Nlandmarks*3, 3)

      a= as.matrix(cbind(data[xseq], data[yseq], data[zseq]))

    }
  }
  return(a)
}#end function


#' Converts array mean into a matrix
#'
#' @return
#' @export
array.mean= function(data){

  mean.matrix= matrix(0, nrow= length(data[,1,1]), ncol= length(data[1,,1]))

  for(i in 1:length(data[,1,1])){
    for(j in 1:length(data[1,,1])){

      mean.matrix[i,j]= mean(data[i,j,])

    }
  }
  return(mean.matrix)
}#end function


#' Calculates centroid size
#'
#' @return
#' @export
csize= function(data){


  csize= rep(0, length(data[1,1,]))

  for(i in 1:length(data[1,1,])){
    centroid = colMeans(data[,,i])
    csize[i]=  sqrt(sum((as.matrix(data[,,i]) - centroid)^2))

  }
  return(csize)

}


#' Calculates mean shape
#'
#' @return
#' @export
mean.shape= function(data){
  m= matrix(nrow=dim(data)[1], ncol=dim(data)[2])
  m= apply(data, c(1,2), mean)
}


#' Converts tps files to row format
#'
#' @return
#' @export

tps2row= function(coord){
  ##only for 2d TPS data

  Nlandmarks= sum(as.numeric( coord$ID)== 1)
  x.y= as.matrix(coord[,1:2])
  tps.matrix= matrix(0, length(x.y[,1])/Nlandmarks, 2)
  xseq= seq(1, Nlandmarks, 2)
  yseq= seq(2, Nlandmarks, 2)

  for(ind in 1:dim(tps.matrix)[1]){

    tps.matrix[ind, xseq]= coord[ind:(ind+Nlandmarks),1]


  }

}




#' Converts arrays to row formats for 2D data
#'
#' @return
#' @export
array2row= function(data){

  Nlandmarks= dim(data)[2] * dim(data)[1]
  x.y= matrix(0, dim(data)[3], dim(data)[1]*dim(data)[2])

  xseq= seq(1, Nlandmarks, 2)
  yseq= seq(2, Nlandmarks, 2)

  for(ind in 1:dim(x.y)[1]){

    x.y[ind, xseq]= data[,1,ind]
    x.y[ind, yseq]= data[,2,ind]


  }
  return(x.y)
}




#' Converts arrays to row formats for 3D data
#'
#' @return
#' @export
array2row3d <- function(data){
  Nlandmarks= dim(data)[2] * dim(data)[1]

  if(length(dim(data)) > 2) {
    x.y= matrix(0, dim(data)[3], dim(data)[1]*dim(data)[2])
  } else if(length(dim(data)) == 2) x.y= matrix(0, 1, dim(data)[1]*dim(data)[2])

  xseq= seq(1, Nlandmarks, 3)
  yseq= seq(2, Nlandmarks, 3)
  zseq= seq(3, Nlandmarks, 3)

  for(ind in 1:dim(x.y)[1]){

    x.y[ind, xseq]= data[,1,ind]
    x.y[ind, yseq]= data[,2,ind]
    x.y[ind, zseq]= data[,3,ind]

  }
  return(x.y)
}

#' Calculates angle between 2 vectors
#'
#' @return
#' @export
angle = function( a, b){
  #acos(cor(a,b)) * 180/pi
  acos(sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) ) *180/pi
}#end function

#' CVA from julian claude's book
#'
#' @return
#' @export
can.var = function(yourarray, gnames){ require(MASS)
  row.aligned= array2row3d(yourarray)
  canvar= lda(row.aligned, gnames)
  LD = canvar$scaling

  ldascore = predict(canvar)$x

  n= dim(row.aligned)[1]
  regr.resid = lm(row.aligned~ gnames)
  dfw= n- length(levels(gnames))
  SSw = var(regr.resid$residuals) * (n-1)
  VCVw = SSw/dfw

  LDs = VCVw %*% LD

  rgbmat= rbind(c(1,0,0), c(0,0,1))
  hist(ldascore[c(which(gnames== levels(gnames)[1])),1 ], xlim= range(ldascore[,1]), main= "Score distribution", xlab= "Discriminant score", col= rgb(rgbmat[,1][1],rgbmat[,2][1], rgbmat[,3][1], alpha= .8))

  for(i in 2:length(unique(gnames))){

    hist(ldascore[c(which(gnames== levels(gnames)[i])), 1 ], col= rgb(rgbmat[,1][i],rgbmat[,2][i], rgbmat[,3][i], alpha= .8), add= T)

  }

  ###more than 3 groups will cause a crash here


  summary(aov(ldascore[,1]~ gnames))

  row.m= colMeans(row.aligned)

  lda.estimate.neg = row2array(t(min(ldascore[,1]) * LDs + row.m))
  lda.estimate.pos = row2array(t(max(ldascore[,1]) * LDs + row.m))
  lda.estimate.mean = row2array(t( 0* LDs + row.m))

  plot(lda.estimate.pos[,,1], typ= "n", axes= F, xlab= "", ylab= "", ylim= range(lda.estimate.neg[,2,1]))

  points(lda.estimate.mean[,,1], lwd=2, lty= 2, col= "red", pch=19)
  points(lda.estimate.neg[,,1], lwd= 3, col="blue", pch=19)
  points(lda.estimate.pos[,,1], lwd= 3, col="green", pch=19)

  legend("bottomleft", legend= c("mean", "negative", "positive"), pch= 19, col= c("red", "blue", "green"))


  ###estimates of lda shape
  return(list(scores= ldascore, neg.est= lda.estimate.neg, mean.est= lda.estimate.mean, pos.est= lda.estimate.pos, LDs= LDs))
}

##need to solve the color histogram problem more generally


#' Calculates the mode
#'
#' @return
#' @export
###mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' PLS from julian claude's book
#'
#' @return
#' @export
pls<-function(M1, M2){
  p1<-dim(M1)[2]
  p2<-dim(M2)[2]
  n<-dim(M1)[1]
  sM12<-svd(var(cbind(M1,M2))[1:p1,(p1+1):(p1+p2)])
  vM12<-var(cbind(M1,M2))[1:p1,(p1+1):(p1+p2)]
  vM21<-var(cbind(M1,M2))[(p1+1):(p1+p2), 1:p1]
  v11<-var(M1)
  v22<-var(M2)
  D<-sM12$d
  F1<-sM12$u
  F2<-sM12$v
  Rv<-sum(diag(vM12%*%vM21))/sqrt(sum(diag(v11%*%v11))* sum(diag(v22%*%v22)))
  x1= matrix(0,nrow= dim(M1)[1], ncol= dim(M1)[2])
  for(i in 1:dim(M1)[2]){x1[,i]= apply(X = M1, MARGIN = 1, function(X) {sum(X*F1[,i])})}
  x2= matrix(0,nrow= dim(M2)[1], ncol= dim(M1)[2])
  for(i in 1:dim(M1)[2]){x2[,i]= apply(X = M2, MARGIN = 1, function(X) {sum(X*F2[,i])})}
  list(Rv=Rv, F1=F1, F2=F2, D=D, x1=x1,x2=x2)
}#end function


#' Checks if a number is even
#'
#' @return
#' @export
is.even <- function(x) {x %% 2 == 0 }



#' Calculates where LMs end up for a given rotation
#'
#' @return
#' @export
rotate.coords = function(X1,Y1, Px, Py, D1){ ##X1=Xcoords, Y1=Ycoords, Px=Rotation pt X, Py = Rotation pt Y, D1=clockwise degrees rotation

  d1= -D1*pi/180
  #xRot = Px + cos(d1) * (X1 – Px) – sin(d1) * (Y1 – Py)
  #yRot = Py + sin(d1) * (X1 – Px) + cos(d1) * (Y1 – Py)
  rot=cbind(Px + cos(d1) * (X1 - Px) - sin(d1) * (Y1 - Py), Py + sin(d1) * (X1 - Px) + cos(d1) * (Y1 - Py))
}

#' Wrapper for Delaunay triangulation function
#'
#' @return
#' @export
tri.surf = function(tri.object, num.passes){
  require(tripack)# triangulations
  require(sp)#point in polygon
  tri.object.prime=tri.object
  for(i in 1:num.passes){
    gen.tri= tri.mesh(tri.object[,1], tri.object[,2])
    tri.cent = matrix(0, nrow= length(triangles(gen.tri)[,1]), ncol= 2) #get centroids from triangulation
    for(i in 1:length(tri.cent[,1])){
      tri.cent[i,]= round(colMeans(tri.object[triangles(gen.tri)[i,1:3],]))
    }

    are.you.in = point.in.polygon(tri.cent[,1], tri.cent[,2], tri.object.prime[,1], tri.object.prime[,2]) #index for out of boundary triangles caused by concavities

    tri.cent= tri.cent[are.you.in==1,] #keep only triangles in original boundaries
    tri.object = rbind(as.matrix(tri.object), tri.cent) #for each iteration of num.passes, bind the new centroid coordinates with the starting coords
  }
  return(tri.object)
}#end function

#' Calculates weighted inverse from DOQTL
#'
#' @return
#' @export
# weighted inverse from DOQTL####
W.inv<- function(W, symmetric=TRUE,inverse=TRUE){
  eW <- eigen(W, symmetric=symmetric)
  d <- eW$values
  if (min(d) <0  && abs(min(d))>sqrt(.Machine$double.eps))
    stop("'W' is not positive definite")
  else d[d<=0]<- ifelse(inverse, Inf, 0)
  A <- diag(d^ifelse(inverse, -0.5, 0.5)) %*% t(eW$vector)
  A # t(A)%*%A = W^{-1}
}


#' Calculates least squares from DOQTL
#'
#' @return
#' @export
lmGls<- function (formula, data, A, ...) {
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$A <- NULL
  m[[1L]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  yy <- model.response(m)
  y<- A%*%yy
  xx <- model.matrix(Terms, m, contrasts)
  x<- A%*%xx; colnames(x)[1]<- "Intercept"
  dtf<- data.frame(y=y,x)
  fit<- lm(y~.-1, data=dtf, ...)

  fit
}

#' rQTL in the DO
#'
#' @return
#' @export
DO.rQTL <- function(pheno, probs, K, addcovar, intcovar, snps) {

  #pheno = as.matrix(LDs.dataframe[,1])
  #K = as.matrix(K.rqtl)
  retval = NULL
  #snps = MM_snps
  #probs = all.probs.rqtl
  #addcovar = covar.rqtl
  #intcovar = DO.weights$Weight


  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]


  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = NULL, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(length(pheno)))


  vTmp$AA = 2 * K


  # This tells QTLRel to fit the additive model.
  class(prdat) = c(class(prdat), "addEff")
  vc = estVC(y = pheno, x = addcovar, v = vTmp)

  # Additive and interactive covariates.
  res = scanOne.rQTL(y = pheno, x = addcovar, prdat = prdat, vc = vc, intcovar = intcovar)

  # Convert the model coefficients to a matrix.
  coef = matrix(unlist(res$parameters), length(res$parameters), length(res$parameters[[1]]), dimnames = list(res$snp, names(res$parameters[[1]])), byrow = TRUE)

  # Return the LRS, LOD, p-value and -log10(p-value).
  p = pchisq(q = res$p, df = dim(probs)[[2]] - 1, lower.tail = FALSE)
  return(list(lod = cbind(snps[,1:4], perc.var = res$v, lrs = res$p, lod = res$p / (2 * log(10)), p = p, neg.log10.p = -log(p, 10)), coef = coef, ss = res$ss))

} # qtl.qtlrel()


#' scanone rqtl
#'
#' @return
#' @export
scanOne.rQTL <- function(y,x,prdat,vc,intcovar = NULL){
  #taken from QTLREL's scanone function
  # prdat$pr: n by ? by ? matrix, allele probabilities
  # vc: object from estVC or aicVC
  # test: “Chisq”, “F” or “Cp”
  nb<- length(vc$par) - sum(vc$nnl)
  nr<- nrow(vc$y)
  cov<- matrix(0,nrow=nr,ncol=nr)
  for(i in 1:vc$nv)if(vc$nnl[i]) cov<- cov + vc$v[[i]]*vc$par[nb+vc$nn[i]]

  diag.cov<- diag(cov)
  if( max( abs( cov-diag(diag.cov) ) ) < min(1e-5,1e-5*max(diag.cov)) ){
    if( max(diag.cov-min(diag.cov)) < min(1e-5,1e-5*max(diag.cov)) ){
      weights<- NULL
    }else weights<- 1/diag.cov
  }else weights<- NA
  gcv<- W.inv(cov)


  nsnp<- dim(prdat$pr)[3]
  if(!is.null(intcovar)) nint<- ncol(as.matrix(intcovar))
  model.par<- vector("list",nsnp)
  names(model.par)<- prdat$snp
  P<- rep(Inf,nsnp)
  names(P)<- prdat$snp
  V<- P
  V2<-P

  for(k in 1:nsnp){
    #null model
    oTmp<- data.frame(y=y,x,intcovar,prdat$pr[,-1,k])
    g0<- lmGls(y~.,data=oTmp,A=gcv)
    P0<- logLik(g0)
    #full model
    oTmp<- data.frame(y=y,x,intcovar,prdat$pr[,-1,k])
    nc<- ncol(oTmp); nq<- ncol(prdat$pr[,-1,k])-1
    str<- paste(paste("(",paste(colnames(oTmp)[nc-nq-(nint:1)],collapse="+"),")",sep=""),
                paste("(",paste(colnames(oTmp)[(nc-nq):nc],collapse="+"),")",sep=""),
                sep=":")
    str<- paste("y~.+",str,sep="")


    g<- lmGls(formula(str),data=oTmp,A=gcv)

    model.par[[k]]<- g$coef
    P[k]<- logLik(g)
    #V[k]<- sum(g$res^2)
    V2[k]<- sum(g$res^2)
    P[k] <- 2*(P[k]-P0)
    V[k] <- sum(g0$res^2) - V2[k]
    V[k]<- V[k]/sum(anova(g0)[,"Sum Sq"])
  }


  list(snp=prdat$snp,
       chr=prdat$chr,
       dist=prdat$dist,
       p=P,
       v=V*100,
       parameters=model.par)
}


#' manhattan plots
#'
#' @return
#' @export
manhattan.plot<-function(chr, pos, pvalue,
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2,
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("grey39","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {

  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")

  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }

  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;

  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }

  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=NULL, pos=3, offset=0.5,
                      col=NULL, fontface=NULL, fontsize=NULL, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }

  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    }
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }

  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)

  if (length(ann.settings)>1) {
    lcols<-trellis.par.get("superpose.symbol")$col
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch,
                              col=lcols[(i-2) %% length(lcols) +1 ],
                              fill=lfills[(i-2) %% length(lfills) +1 ],
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label,
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label,
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]],
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)

  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places),
      pos=round(genpos,thin.pos.places),
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()

  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr),
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }

  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) {
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+.5;
    A$ylim=c(0,maxy);
    A;
  }

  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings,
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab,
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}



#' @param X A landmark set
#' @return A matrix of symmetric landmarks for each observation
#' @export
OSymm <- function(X, midline, right, left) {
  ncl <- ncol(X)
  Xr <- cbind(X[,-ncl], -X[,ncl])
  Xow <- Xo <- rbind(X[c(midline, right, left),])
  Xrw <- Xr <- rbind(Xr[c(midline, left, right),])
  rownames(Xrw) <- rownames(Xr) <- rownames(X)
  Xo[which(is.na(Xr))] <- NA
  Xr[which(is.na(Xo))] <- NA
  mo <- matrix(apply(Xo, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xow), nc=ncol(Xow))
  mr <- matrix(apply(Xr, 2, mean, na.rm=TRUE), byrow=TRUE, nr=nrow(Xrw), nc=ncol(Xrw))
  Xrwc <- Xrw-mr
  SVD <- svd(t(na.omit(Xr-mr)) %*% na.omit(Xo-mo))
  L <- diag(SVD$d)
  S <- ifelse(L<0, -1, L)
  S <- ifelse(L>0, 1, L)
  RM <- SVD$v %*% S %*% t(SVD$u)
  Xrot <- (Xow-mo) %*% RM
  SC <- apply(array(c(Xrwc,Xrot), dim=c(nrow(Xrot),ncol(Xrot),2), dimnames=list(rownames(Xrot),colnames(Xrot))), 1:2, mean, na.rm=TRUE)
  Xrot[which(is.na(Xrot))] <- Xrwc[which(is.na(Xrot))]
  Xrwc[which(is.na(Xrwc))] <- Xrot[which(is.na(Xrwc))]
  list(rec.orig=Xrot, symmconf=SC, rec.ref=Xrwc)}



#' probablistic rank based sensitivity plots
#'
#' @param model A model fit in base/caret/H2O (Must be probablistic)
#' @return Sensitivity and Balanced Accuracy plots for the model
#' @export
rank.sens <- function(model, validation_set = testing){

  if(class(model) == 'H2OMultinomialModel'){
    model.stack <- h2o.predict(model, validation_set)
    model.probs <- as.matrix(model.stack[,-1])
    top1sens <- confusionMatrix(as.factor(as.matrix(model.stack[,1])), testing$Syndrome)$byClass[,1]
    top1acc <- confusionMatrix(as.factor(as.matrix(model.stack[,1])), testing$Syndrome)$byClass[,11]
  } else {model.stack <- predict(model, validation_set)
  model.probs <- predict(model, validation_set, type = "prob")
  top1sens <- confusionMatrix(model.stack, testing$Syndrome)$byClass[,1]
  top1acc <- confusionMatrix(model.stack, testing$Syndrome)$byClass[,11]}


  top5.check <- rep(NA, length(testing$Syndrome))

  for(i in 1: length(testing$Syndrome)){
    if(sum(names(sort(model.probs[i,], decreasing = T)[1:5]) == testing$Syndrome[i]) > 0){
      top5.check[i] <- as.character(testing$Syndrome[i])
    } else{top5.check[i] <- names(sort(model.probs[i,], decreasing = T)[1])}
  }

  #top 10
  top10.check <- rep(NA, length(testing$Syndrome))

  for(i in 1: length(testing$Syndrome)){
    if(sum(names(sort(model.probs[i,], decreasing = T)[1:10]) == testing$Syndrome[i]) > 0){
      top10.check[i] <- as.character(testing$Syndrome[i])
    } else{top10.check[i] <- names(sort(model.probs[i,], decreasing = T)[1])}
  }



  top5sens <- confusionMatrix(as.factor(top5.check), testing$Syndrome)$byClass[,1]
  top10sens <- confusionMatrix(as.factor(top10.check), testing$Syndrome)$byClass[,1]

  model.rank.sens <- rbind(cbind(top1sens, 1),
                           cbind(top5sens - top1sens, 5),
                           cbind(top10sens - top5sens, 10)
  )

  model.rank.sens <- data.frame(model.rank.sens, syndrome = sort(unique(testing$Syndrome)))
  model.rank.sens[,2] <- as.factor(model.rank.sens[,2])
  colnames(model.rank.sens)[1:2] <- c("sensitivity", "rank")

  fill <- c("black", "darkgrey", "lightgrey")

  p1 <- ggplot() +
    geom_bar(aes(y = sensitivity, x = syndrome, fill = rank), data = model.rank.sens, stat="identity",  position = position_stack(reverse = T)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_hline(yintercept = .7, color = 2) +
    xlab("Syndrome") +
    ylab("Sensitivity") +
    scale_fill_manual(values=fill) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))



  top5acc <- confusionMatrix(as.factor(top5.check), testing$Syndrome)$byClass[,11]
  top10acc <- confusionMatrix(as.factor(top10.check), testing$Syndrome)$byClass[,11]

  model.rank.acc <- rbind(cbind(top1acc, 1),
                          cbind(top5acc - top1acc, 5),
                          cbind(top10acc - top5acc, 10)
  )

  model.rank.acc <- data.frame(model.rank.acc, syndrome = sort(unique(testing$Syndrome)))
  model.rank.acc[,2] <- as.factor(model.rank.acc[,2])
  colnames(model.rank.acc)[1:2] <- c("accuracy", "rank")

  p2 <- ggplot() +
    geom_bar(aes(y = accuracy, x = syndrome, fill = rank), data = model.rank.acc, stat="identity",  position = position_stack(reverse = T))  +
    geom_hline(yintercept = .7, color = 2) +
    xlab("Syndrome") +
    ylab("Balanced Accuracy") +
    scale_fill_manual(values=fill) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

  return(list(p1, p2))

}





