ContTable=function(tab, title, chisqtest=F,xlabL="X label",ylabL="Y label", scaleV=T){
  library("RColorBrewer")
  if (chisqtest==T){
    a1=chisq.test(tab)
    tit2=sprintf("Chisq = %g", round(a1$p.value*100)/100)
  } else {
    tit2=" "
  }
  nr=nrow(tab)
  nc=ncol(tab)
  if (chisqtest==T){
    l1=(a1$observed-a1$expected)
    l1[which(l1<0, arr.ind=T)]=0
    image(l1, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
          xlab=xlabL, ylab=ylabL,
          main=sprintf("%s %s", title, tit2))
  }else{
    if (scaleV==T){
      image(scale(tab), col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
            xlab=xlabL, ylab=ylabL,
            main=sprintf("%s %s", title, tit2))
    }else{
      image(tab, col=brewer.pal(9, "BuGn"), xaxt="n", yaxt="n",
            xlab=xlabL, ylab=ylabL,
            main=sprintf("%s %s", title, tit2))   
    }}
  xval=seq(0, 1, length=nr)
  yval=seq(0, 1, length=nc)
  axis(1, at=xval, rownames(tab))
  axis(2, at=yval, colnames(tab))
  for(i in 1:nr){
    for (j in 1:nc){
      text(xval[i], yval[j], tab[i,j], cex=1.5)
    }
  }
}

## Create variant key
variant_key <- function(df) {
  paste(df$CHROM, df$POS, df$Tumor_Sample_Barcode, sep = ":")
}

## Count overlap per sample
count_overlap <- function(key1, samp1, key2, samp2, samples) {
  sapply(samples, function(s) {
    sum(key1[samp1 == s] %in% key2[samp2 == s])
  })
}

## Safe numeric filter
af_filter <- function(df, cutoff) {
  mgrb_af_num <- as.numeric(df$mgrb_af)
  mgrb_af_num2 <- ifelse(is.na(mgrb_af_num), 0, mgrb_af_num)
  gn_af_num   <- as.numeric(df$gn_af)
  gn_af_num2 <- ifelse(is.na(gn_af_num),   0, gn_af_num)
  
  keep <- ifelse(is.na(mgrb_af_num), 0, mgrb_af_num) <= cutoff &
    ifelse(is.na(gn_af_num),   0, gn_af_num)   <= cutoff
  
  keep2 <- which(mgrb_af_num2<=cutoff & gn_af_num2<=cutoff)
  test=df[keep, ]
}

