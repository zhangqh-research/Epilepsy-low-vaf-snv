args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1)
module <- normalizePath(args[1],winslash='/')
d <- read.delim(file.path(module,'tables','patient_group_metric_medians.tsv'),check.names=FALSE,na.strings='NA')
only <- c('Mutect2-only','Strelka2-only','MuSE-only','FreeBayes-only')
consensus <- 'Multi-caller consensus'
results <- list()
for (metric in unique(d$metric)) {
  a <- d[d$metric==metric,]
  reference <- a[a$group==consensus,c('pair_id','patient_median')]
  names(reference)[2] <- 'consensus_median'
  batch <- list()
  for (group in only) {
    x <- a[a$group==group,c('pair_id','patient_median')]
    z <- merge(x,reference,by='pair_id')
    z <- z[is.finite(z$patient_median) & is.finite(z$consensus_median),]
    delta <- z$patient_median-z$consensus_median
    nz <- delta[delta!=0]
    test <- NULL;status <- 'OK';warnings <- character()
    if (length(delta)<10) status <- 'descriptive_only_fewer_than_10_patient_pairs'
    else if (length(nz)<5) status <- 'descriptive_only_fewer_than_5_nonzero_differences'
    else test <- withCallingHandlers(wilcox.test(z$patient_median,z$consensus_median,paired=TRUE,alternative='two.sided',exact=FALSE,correct=TRUE),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    qq <- if (length(delta)) quantile(delta,c(.25,.5,.75),names=FALSE,type=7) else rep(NA_real_,3)
    skew <- if (is.finite(qq[3]-qq[1]) && qq[3]>qq[1]) (qq[3]+qq[1]-2*qq[2])/(qq[3]-qq[1]) else NA_real_
    ranks <- rank(abs(nz),ties.method='average')
    rbc <- if(length(nz)) (sum(ranks[nz>0])-sum(ranks[nz<0]))/sum(ranks) else NA_real_
    batch[[group]] <- data.frame(metric=metric,only_group=group,reference_group=consensus,
      n_evaluable_patient_pairs=nrow(z),n_nonzero_differences=length(nz),
      n_positive=sum(delta>0),n_negative=sum(delta<0),n_zero=sum(delta==0),
      paired_delta_q25=qq[1],paired_delta_median=qq[2],paired_delta_q75=qq[3],
      signed_rank_biserial=rbc,wilcoxon_V=if(is.null(test))NA_real_ else unname(test$statistic),
      p_value=if(is.null(test))NA_real_ else test$p.value,p_holm=NA_real_,holm_family_size=4,status=status,
      delta_quartile_skewness=skew,asymmetry_caution=is.finite(skew) && abs(skew)>.5,
      patient_ids=paste(z$pair_id,collapse=';'),warnings=paste(unique(warnings),collapse=';'),stringsAsFactors=FALSE)
  }
  out <- do.call(rbind,batch)
  out$p_holm <- p.adjust(out$p_value,method='holm',n=4)
  results[[metric]] <- out
}
out <- do.call(rbind,results)
options(digits=17)
target <- file.path(module,'tables','patient_paired_comparisons.tsv')
stopifnot(!file.exists(target))
con <- file(target,open='wx',encoding='UTF-8')
write.table(out,con,sep='\t',row.names=FALSE,quote=FALSE,na='NA');close(con)
cat('Comparisons:',nrow(out),'\n')
