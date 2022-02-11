#====================================================================
# Plot the top 2 PCs of the K matrix showing tst and trn points (user-level)
#====================================================================
# Z = NULL; tst = U = d = group = group.shape = set.color = set.size = df = NULL
# axis.labels = TRUE; curve = FALSE; bg.color = "gray20"; unified = TRUE; ntst = 36;
# line.color = "gray90"; line.tick = 0.3; legend.pos="right";
# point.color = "gray20"; sets = c("Testing","Supporting","Non-active")

net.plot <- function(object, B, Z = NULL, K, tst = NULL,
           U = NULL, d = NULL, group = NULL, group.shape = NULL,
           set.color = NULL, set.size = NULL, df = NULL, main, axis.labels = TRUE,
           curve = FALSE, bg.color = "gray20", unified = TRUE, ntst = 36,
           line.color = "gray90", line.tick = 0.3, legend.pos="right",
           point.color = "gray20", sets = c("Testing","Supporting","Non-active"))
{
  set <- PC1 <- PC2 <- PC1_TRN <- PC1_TST <- PC2_TRN <- PC2_TST <- NULL
  legend.pos <- match.arg(legend.pos,
    choices=c("right","bottomright","bottomleft","topleft","topright","none"))

  if(!inherits(object, "SSI")) stop("The input object is not of the class 'SSI'")

  if(is.null(U) & is.null(d))
  {
    if(is.character(K)){
      K <- readBinary(K)
    }
    if(is.null(K))
      stop("Matrix 'K' must be a positive semi definite matrix\n")
    if(!is.null(Z)) {
      if(length(dim(Z))!=2) stop("Object 'Z' must be a matrix with ncol(Z)=nrow(K)\n")
      K <- float::tcrossprod(Z,float::tcrossprod(Z,K))   # Z%*%K%*%t(Z)
    }
    tmp <- float::svd(K,nu=2,nv=0)
    d <- tmp$d
    U <- tmp$u
    expvarPC <- 100*d/sum(d)  # 100*d/sum(float::diag(K))
  }else{
    if(is.null(U)){
      stop("You are providing the eigevalues, but not the eigenvectors")
    }else{
      if(is.null(d)){
        message("You are providing the eigenvectors, but not the eigenvalues\n",
                "No variance explained can be calculated")
        expvarPC <- NULL
      }else{
        if(nrow(U) == length(d)){
          expvarPC <- 100*d/sum(d)
        }else expvarPC <- NULL
      }
    }
  }
  tmp <- paste0(" (",sprintf('%.1f',expvarPC),"%)")
  if(length(tmp)<2) tmp <- NULL
  labelsPC <- paste0("PC ",1:2,tmp[1:2])

  if(!is.null(tst)){
      if(any(!tst %in% object$tst))
          stop("Some elements in 'tst' are not contained in set 'object$tst'")
  }else tst <- object$tst
  if(!unified & length(tst) >= ntst){
   cat("Large number of testing individuals. Only the first",ntst,"are shown\n")
   tst <- tst[1:ntst]
  }

  justx <- ifelse(length(grep("left",legend.pos))>0,0,1)
  justy <- ifelse(length(grep("bottom",legend.pos))>0,0,1)
  if(!legend.pos %in% c("none","right")) legend.pos <- c(abs(justx-0.01),abs(justy-0.01))

  theme0 <- ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    legend.box.spacing = ggplot2::unit(0.4, "lines"),
    legend.background = ggplot2::element_rect(fill = "gray95"),
    panel.background = ggplot2::element_rect(fill = bg.color),
    legend.justification = c(justx,justy),
    legend.position=legend.pos,
    legend.key.height = ggplot2::unit(0.9,"line"),
    legend.key.width = ggplot2::unit(0.9, "lines"),
    legend.title = ggplot2::element_blank(),
    legend.margin = ggplot2::margin(t=0,b=0.25,l=0.25,r=0.25,unit='line'),
    strip.text = ggplot2::element_blank(), panel.spacing = ggplot2::unit(0.1,"lines")
  )

  if(missing(B)){
    if(is.null(df)) df <- summary.SSI(object)$optCOR['df']
    if(0 > df | df > range(object$df)[2])
      stop("Parameter 'df' must be greater than zero and no greater than trn size")
    B <- as.matrix(coef.SSI(object, df=df))
  }else{
    stopifnot(length(dim(B))==2L)
    df <- mean(do.call(c, lapply(1:nrow(B), function(i) sum(abs(B[i,]) > 0))))
  }

  flagGp <- !is.null(group)
  if(is.null(group)) group <- data.frame(group=rep(1,nrow(U)))
  gpName <- colnames(group)

  if(!(class(sets) == "character" & length(sets) == 3))
   stop("Parameter 'sets' must be a triplet of 'character' type")

  dat <- data.frame(id=1:nrow(U),set=sets[3],group=group,float::dbl(U[,1:2]))
  dat$set <- as.character(dat$set)

  # Testing and training (active) set
  dat$set[tst] <- sets[1]
  index <- do.call(c, lapply(1:ncol(B), function(j) any(abs(B[object$tst %in% tst,,drop=FALSE][,j]) > 0)))
  dat$set[object$trn[index]] <- sets[2]
  dat$set[object$trn[!index]] <- sets[3]

  colnames(dat) <- c("id","set","group","PC1","PC2")

  dat$group <- factor(as.character(dat$group))
  dat$set <- factor(dat$set,levels=c(sets))

  # Shape and color for the levels of group
  if(!flagGp) dat$group <- dat$set
  levelsGp <- levels(dat$group)
  if(length(levelsGp) > 5)
   stop("Number of levels of 'group' must be at most 5")

  if(is.null(group.shape)){
    if(flagGp){
      group.shape <- c(21,22,23,24,25)
    }else group.shape <- c(21,21,21)
  }
  group.shape <- group.shape[1:length(levelsGp)]

  if(is.null(set.color)){
    set.color <- c("#E69F00","#56B4E9","#999999")
  }
  set.color <- set.color[1:length(sets)]

  if(is.null(set.size)){
    set.size <- c(2.5,1.5,1)
  }
  set.size <- set.size[1:length(sets)]

  if(any(is.na(group.shape)))
    stop("The number of elements in 'group.shape' must be of length ",length(levelsGp))

  if(any(is.na(set.size)) | any(is.na(set.color)))
    stop("The number of elements in 'set.size' and 'set.color' must be of length ",length(sets))

  if(missing(main)){
     main0 <- bquote(.(object$name)*". Support set size="*.(round(df)))
     theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }else{
    main0 <- main
    if(is.null(main)){
      theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_blank())
    }else{
      theme0 <- theme0 + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }
  }

  if(!axis.labels){
    theme0 <- theme0 + ggplot2::theme(axis.text=ggplot2::element_blank(),
                        axis.ticks=ggplot2::element_blank())
  }

  names(group.shape) <- levelsGp
  names(set.color) <- names(set.size) <- sets

  # If unified plot
  if(unified)
  {
    pt <- ggplot2::ggplot(dat,ggplot2::aes(x=PC1,y=PC2)) +
           ggplot2::geom_point(data=dat[dat$set==sets[3],],ggplot2::aes(shape=group,fill=set),
                    color=point.color,size=set.size[3])

    for(i in 1:length(tst))
    {
      indexTRN <- which(abs(B[which(object$tst == tst[i]), ]) > 0)
      if(length(indexTRN)>0)
      {
        dat1 <- dat[object$trn, c("PC1","PC2")][indexTRN,]
        dat2 <- dat[tst, c("PC1","PC2")][i,]
        colnames(dat1) <- paste0(colnames(dat1),"_TRN")
        colnames(dat2) <- paste0(colnames(dat2),"_TST")
        dat1 <- data.frame(dat2[rep(1,nrow(dat1)),],dat1)
        if(curve){
          pt <- pt + ggplot2::geom_curve(ggplot2::aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                        data=dat1,alpha=0.4,size=line.tick,color=line.color,curvature=0.4)
        }else{
          pt <- pt + ggplot2::geom_segment(ggplot2::aes(x=PC1_TST,y=PC2_TST,xend=PC1_TRN,yend=PC2_TRN),
                        data=dat1,alpha=0.4,size=line.tick,color=line.color)
        }
      }
    }

    pt <- pt  +
      ggplot2::geom_point(data=dat[dat$set==sets[1],],ggplot2::aes(shape=group,fill=set),color=point.color,size=set.size[1]) +
      ggplot2::geom_point(data=dat[dat$set==sets[2],],ggplot2::aes(shape=group,fill=set),color=point.color,size=set.size[2]) +
      ggplot2::theme_bw() + theme0
  }else{
      dat2 <- c()
      for(i in 1:length(tst))
      {
        indexTRN <- which(abs(B[which(object$tst == tst[i]), ]) > 0)
        if(length(indexTRN) > 0)
        {
          tmp <- dat[object$trn, ][-indexTRN,]
          tmp$set <- sets[3]
          tmp <- rbind(dat[object$trn, ][indexTRN,], tmp, dat[tst, ][i,])
          dat2 <- rbind(dat2,data.frame(tmp, ind = i))
        }
      }

      pt <- ggplot2::ggplot(dat2,ggplot2::aes(x=PC1,y=PC2)) + ggplot2::facet_wrap(~ind) +
             ggplot2::geom_point(data=dat2[dat2$set==sets[3],],ggplot2::aes(fill=set,shape=group),color=point.color,size=set.size[3]) +
             ggplot2::geom_point(data=dat2[dat2$set==sets[2],],ggplot2::aes(fill=set,shape=group),color=point.color,size=set.size[2]) +
             ggplot2::geom_point(data=dat2[dat2$set==sets[1],],ggplot2::aes(fill=set,shape=group),color=point.color,size=set.size[1]) +
             ggplot2::theme_bw() + theme0

  }

  pt <- pt + ggplot2::labs(title=main0, x=labelsPC[1],y=labelsPC[2]) +
    ggplot2::scale_shape_manual(values = group.shape,
              guide=ggplot2::guide_legend(override.aes=list(size=2,fill="white"))) +
    ggplot2::scale_fill_manual(values = set.color,
              guide=ggplot2::guide_legend(override.aes=list(shape=21,size=2)))

 if(!flagGp) pt <- pt + ggplot2::guides(shape="none")

  pt
}
