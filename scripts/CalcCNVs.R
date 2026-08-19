## Calc CIN

CalcNBP=function(DM, chrid, thresh=0.05){
  # inputs: DM: chromosome CNS (vector)
  # genomic ranges object which contains chromosome location for each gene
  # thresh: minimum threshold to define a CN change
  # chrid is a genomic ranges object with locations of genes includ chr no
  # calculate the number of breakpoints
  DiffId=DM[2:length(DM)]-DM[1:(length(DM)-1)]
  # figure out where chromosome ends are and remove these samples
  t2=seqlevels(chrid)
  mx=match(t2, seqnames(chrid))
  nb=which(DiffId[-na.omit(mx-1)]>0.05)
  length(nb)
}


CalcGGI=function(DM, thresh=0.5,Glength=NULL){
  ## INPUTS
  ## DM: vector of CNgains/losses sorted by location in genome
  ## thresh: threshold to determine a gain or loss
  ## Glength: length of each gene (weighted CN according to length)
  l1=which(abs(DM)>=(thresh))
  if(is.null(Glength)==T){
    l2=length(l1)/length(DM)
  }else{
    l2=sum(Glength[l1])/sum(Glength)
  }
  l2
}