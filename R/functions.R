#### Common functions used in the analysis

renameMarker <- function(oldname){
  oldname <- ifelse(oldname == 'HLA_DR','HLA-DR', oldname)
  oldname <- ifelse(oldname == 'Alpha_SMA','SMA', oldname)
  oldname <- ifelse(oldname == 'E_Cadherin','E-cadherin', oldname)
  oldname <- ifelse(oldname == 'CA_IX','CAIX', oldname)
  oldname <- ifelse(oldname == 'Collagen_I','Collagen-I', oldname)
  oldname <- ifelse(oldname == 'Ki67','Ki-67', oldname)
  oldname <- ifelse(oldname == 'CD274_PDL1','PDL1', oldname)
  oldname <- ifelse(oldname == 'CD134_OX40','OX40', oldname)
  oldname <- ifelse(oldname == 'CD279_PD1','PD1', oldname)
  oldname <- ifelse(oldname == 'CD223_LAG3','LAG3', oldname)
  oldname <- ifelse(oldname == 'CD366_TIM3','TIM3', oldname)
  oldname <- ifelse(oldname == 'CD278_ICOS','ICOS', oldname)
  oldname <- ifelse(oldname == 'CD278','ICOS', oldname)
  return(oldname)
}

renamePID <- function(oldname,IDtype = 'ROI'){
  newname <- c()
  for (each in oldname){
    if (stringr::str_detect(each, 'F')){
      a = paste0(each)
    }else{
      if (IDtype == 'patient'){
        a = paste0(each,'F')
      }else{
        splitID = str_split(each,'_')[[1]]
        a = paste0(splitID[1], 'F_', splitID[2])
      }
    }
    newname = c(newname,a)
  }
  return(newname)
}


renameCluster <- function(clustern){
  clustern <- ifelse(clustern=='Double-positive T cells','DPT',clustern)
  clustern <- ifelse(clustern=='MC1 (HLADR-CD14+CD16-CD11C+CD11b-)','MC1',clustern)
  clustern <- ifelse(clustern=='MC2 (HLADR-CD14+CD16-CD11C-CD11b-)','MC2',clustern)
  clustern <- ifelse(clustern=='MC3 (HLADR-CD14-CD16-CD11C-CD11b^{hi})','MC3',clustern)
  clustern <- ifelse(clustern=='MC4 (HLADR^{hi}CD14^{hi}CD16+CD11c+CD11b+)','MC4',clustern)
  clustern <- ifelse(clustern=='MC5 (HLADR+CD14^{hi}CD16+CD11c-CD11b+)','MC5',clustern)
  clustern <- ifelse(clustern=='MC6 (HLADR+CD14+CD16-CD11c-CD11b-)','MC6',clustern)
  clustern <- ifelse(clustern=='Regulatory T cells','Treg',clustern)
  clustern <- ifelse(clustern=='B cells','B',clustern)
  clustern <- ifelse(clustern=='CD4 T cells','CD4T',clustern)
  clustern <- ifelse(clustern=='CD8 T cells','CD8T',clustern)
  clustern <- ifelse(clustern=='Lymphocyte','Lym.',clustern)
  return(clustern)
}

renameCluster2.forcox <- function(clustern){
  clustern <- ifelse(clustern=='Stroma (Collagen+)','S1.Collagen',clustern)
  clustern <- ifelse(clustern=='Stroma (FAP+)','S2.FAP',clustern)
  clustern <- ifelse(clustern=='Stroma (PDGFRb+)','S3.PDGFRb',clustern)
  clustern <- ifelse(clustern=='Stroma (SMA+)','S4.SMA',clustern)
  clustern <- ifelse(clustern=='Stroma (Vimentin+)','S5.Vimentin',clustern)
  clustern <- ifelse(clustern=='Tumor (CA9+)','T1.CA9',clustern)
  clustern <- ifelse(clustern=='Tumor (Ki67+)','T2.Ki67',clustern)
  clustern <- ifelse(clustern=='Tumor (VEGF+)','T3.VEGF',clustern)
  clustern <- ifelse(clustern=='Tumor n.c.','T4.n.c', clustern)
  clustern <- ifelse(clustern=='Double-positive T cells','DPT',clustern)
  clustern <- ifelse(clustern=='MC1 (HLADR-CD14+CD16-CD11C+CD11b-)','MC1',clustern)
  clustern <- ifelse(clustern=='MC2 (HLADR-CD14+CD16-CD11C-CD11b-)','MC2',clustern)
  clustern <- ifelse(clustern=='MC3 (HLADR-CD14-CD16-CD11C-CD11b^{hi})','MC3',clustern)
  clustern <- ifelse(clustern=='MC4 (HLADR^{hi}CD14^{hi}CD16+CD11c+CD11b+)','MC4',clustern)
  clustern <- ifelse(clustern=='MC5 (HLADR+CD14^{hi}CD16+CD11c-CD11b+)','MC5',clustern)
  clustern <- ifelse(clustern=='MC6 (HLADR+CD14+CD16-CD11c-CD11b-)','MC6',clustern)
  clustern <- ifelse(clustern=='Regulatory T cells','Treg',clustern)
  clustern <- ifelse(clustern=='B cells','B',clustern)
  clustern <- ifelse(clustern=='CD4 T cells','CD4T',clustern)
  clustern <- ifelse(clustern=='CD8 T cells','CD8T',clustern)
  return(clustern)
}


volcanoPlot2 <- function(pval_df, p_val = "padj", clustern = 'cluster', fcn='logfc', p.threshold = 0.05, fc.threshold=2, change=NULL, c.color=NULL, region=FALSE, countn='count'){
  #pval_df <- pvalue.4.inv
  plotdf <- data.frame(logfc = pval_df[,fcn], pvalue = pval_df[,p_val], cluster = pval_df[,clustern], count= pval_df[,countn])
  plotdf <- plotdf[!is.infinite(plotdf$logfc), ]
  
  if (is.null(change)){
    plotdf$change = as.factor(ifelse(plotdf$pvalue < p.threshold & abs(round(plotdf$logfc,1)) >= fc.threshold,
                                     ifelse(round(plotdf$logfc,1) >= fc.threshold ,'UP','DOWN'),'NOT'))
    change.color <- c('blue',"grey", "red2")
  }else{
    plotdf$change <- pval_df$change
    change.color <- c.color
  }
  
  if (region){
    plotdf$region <- pval_df$region
    gg_logFC = ggplot(plotdf, aes(x = logfc, y = -log10(pvalue))) + 
      geom_point(aes(color = change, size = count)) + theme_bw() + 
      theme(legend.position = 'right', text = element_text(size = 25)) +
      scale_color_manual(values = change.color) + 
      facet_wrap(~region, nrow = 1, scales = 'free_x') + 
      geom_label_repel(data = plotdf, 
                       aes(x = logfc, y = -log10(pvalue), 
                           label = ifelse(((-log10(pvalue) >= -log10(p.threshold)) & (abs(round(logfc,1)) >= fc.threshold)), as.character(cluster), "")), 
                       point.padding = 1, force = 10) +
      xlab("log2(FoldChange)") + ylab("-log10(FDR)") +
      geom_vline(xintercept=c(-fc.threshold,fc.threshold),lty=3,col="black",lwd=0.5) + #add line,|FoldChange|>2
      geom_hline(yintercept = -log10(p.threshold),lty=3,col="black",lwd=0.5) #add line, padj<0.05
    #+ theme(legend.position = "none")
  }else{
    gg_logFC = ggplot(plotdf, aes(x = logfc, y = -log10(pvalue))) + 
      geom_point(aes(color = change, size = 1)) + theme_bw() + 
      scale_color_manual(values = change.color) + 
      geom_label_repel(data = plotdf, 
                       aes(x = logfc, y = -log10(pvalue), 
                           label = ifelse(((-log10(pvalue) >= -log10(p.threshold)) & (abs(round(logfc,1)) >= fc.threshold)), as.character(cluster), "")), 
                       point.padding = 1, force = 10) + 
      xlab("log2(FoldChange)") + ylab("-log10(FDR)") +
      geom_vline(xintercept=c(-fc.threshold,fc.threshold),lty=3,col="black",lwd=0.5) + #add line, |FoldChange|>2
      geom_hline(yintercept = -log10(p.threshold),lty=3,col="black",lwd=0.5) #add line, padj<0.05
    #+ theme(legend.position = "none")
  }
  return(gg_logFC)
}


getP <- function(result, count.table=NULL, reg = 'invasive', percentcol = 'count',
                 clustercol = "cluster", treatcol = "response"){
  #result <- tmp
  #reg = 'invasive'
  p.all <- data.frame(marker = c(), prop_R = c(), prop_NR = c(), logfc = c(), p=c())
  colnames(result)[which(colnames(result)==percentcol)] <- 'percent' 
  colnames(result)[which(colnames(result)==clustercol)] <- 'cluster'
  colnames(result)[which(colnames(result)==treatcol)] <- 'response'
  
  for (each in unique(result$cluster)){
    tmp <- subset(result, cluster == each & region==reg)
    p1 = wilcox.test(percent ~ response, data = tmp)$p.value
    
    prop_R <- mean(subset(tmp, cluster == each & response=='R')$percent)
    prop_NR <- mean(subset(tmp, cluster == each & response=='NR')$percent)
    logFC <- log2(prop_R/prop_NR)
    
    p.all = rbind(p.all, data.frame(marker = each, prop_R = prop_R, prop_NR = prop_NR, logfc = logFC, p = p1))
  }
  
  if (!is.null(count.table)){
    p.all <- merge(p.all, count.table, by.x='marker',by.y='X')
    colnames(p.all)[which(colnames(p.all)=='X0')] <- 'count'
  }
  
  p.all$p.adj.BH <- p.adjust(p.all$p, method = 'BH')
  
  return(p.all)
}

##### for survival analysis:
setGroup <- function(survival_df, variables=c('UP','DOWN','UP.DOWN'), statuscol='status', timecol='time'){
  colnames(survival_df)[which(colnames(survival_df)==statuscol)] = 'status'
  colnames(survival_df)[which(colnames(survival_df)==timecol)] = 'time'
  
  res.cut <- surv_cutpoint(survival_df, time = "time", event = "status",
                           variables = variables)
  res.cat <- surv_categorize(res.cut)
  for (each in variables){
    res.cat[,each] <- as.character(res.cat[,each])
  }
  res.cat[res.cat=='high'] = 'High'
  res.cat[res.cat=='low'] = 'Low'
  survival_df <- res.cat
  return(survival_df)
}
###################################

# for neighborhood analysis
getP1 <- function(result, reg='cluster1'){
  #result <- tmp
  #reg = 'invasive'
  #result <- result[result$FirstLabel != 'Other PDGFRb+' & result$SecondLabel != 'Other PDGFRb+',]
  #result <- result[result$FirstLabel != 'Other FAP+' & result$SecondLabel != 'Other FAP+',]
  
  p.all <- data.frame(cluster1 = c(), cluster2=c(), ct_base = c(), ct_perm = c(), logfc = c(), p=c())
  #colnames(result)[which(colnames(result)=='log2')] <- 'percent' 
  #colnames(result)[which(colnames(result)=='treatment')] <- 'response'
  
  #result <- data_sub
  
  for (each1 in unique(result$FirstLabel)){
    for (each2 in unique(result$SecondLabel)){
      
      tmp1 <- subset(result, FirstLabel == each1 & SecondLabel == each2)
      
      if (dim(tmp1)[1] != 0){
        tmp <- melt(tmp1, id.vars = c("ID","group","FirstLabel","SecondLabel","p_gt","p_lt",
                                      "direction","p","sig","sigval","log2"))
        #colnames(tmp)
        p1 = wilcox.test(value ~ variable, data = tmp)$p.value
        
        ct_perm <- mean(subset(tmp, variable=='ct_perm_mean')$value)
        ct_base <- mean(subset(tmp, variable=='ct_baseline')$value)
        logFC <- log2(ct_base/ct_perm)
        
        p.all = rbind(p.all, data.frame(cluster1 = each1, cluster2= each2, ct_base=ct_base, ct_perm = ct_perm, logfc = logFC, p = p1))
      }
    }
  }
  p.all$p.adj.BH <- p.adjust(p.all$p, method = 'BH')
  p.all$region <- reg
  #p.all$tumor.p.adj.BH <- p.adjust(p.all$tumor.p, method = 'BH')
  return(p.all)
}

plotDatp <- function(dat_p){
  pmat = dcast(dat_p, 'FirstLabel ~ SecondLabel', value.var = 'sigval', fun.aggregate = sum,
               fill=0, drop=F)
  
  rname = pmat$FirstLabel
  
  pmat = pmat %>%
    select(-c('FirstLabel')) %>%
    as.matrix()
  row.names(pmat) <- rname
  return(pmat)
}

###################################

# for DEGs 
runedgeR <- function(exprSet,group_list,
                     g1="control",g2="case",
                     pro='test'){
  print(table(group_list))
  colnames(exprSet)
  cat(paste0('Now process the project : ',pro))
  
  ### ---------------
  ###
  ### Then run edgeR 
  ###
  ### ---------------
  library(edgeR)
  g=factor(group_list)
  g=relevel(g,g1)
  d <- DGEList(counts=exprSet,group=g)
  keep <- rowSums(cpm(d)>1) >= 2
  table(keep)
  d <- d[keep, , keep.lib.sizes=FALSE]
  d$samples$lib.size <- colSums(d$counts)
  d <- calcNormFactors(d)
  d$samples
  
  dge=d
  design <- model.matrix(~0+factor(group_list))
  rownames(design)<-colnames(dge)
  colnames(design)<-levels(factor(group_list))
  
  dge <- estimateGLMCommonDisp(dge,design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  
  fit <- glmFit(dge, design)
  # https://www.biostars.org/p/110861/
  lrt <- glmLRT(fit,  contrast=c(-1,1)) 
  nrDEG=topTags(lrt, n=nrow(dge))
  nrDEG=as.data.frame(nrDEG)
  head(nrDEG)
  DEG_edgeR =nrDEG 
  nrDEG=DEG_edgeR[,c(1,5)]
  colnames(nrDEG)=c('log2FoldChange','pvalue')
  #draw_h_v(exprSet,nrDEG,paste0(pro,'_edgeR')) 
  
  save(DEG_edgeR,
       exprSet,group_list,
       file = paste0(pro,'_DEG_results_edgeR.Rdata')) 
  
  write.csv(nrDEG,file =  paste0(pro,'_DEG_results_edgeR.csv'))
  return(nrDEG)
}

getDEGs <- function(DEG_result,DEG_type = 'Deseq2',thre_logFC=1,thre_p=0.05, usePadj = TRUE){
  if (DEG_type == 'Deseq2'){
    colnames(DEG_result)[which(colnames(DEG_result)=='log2FoldChange')] = 'logFC'
    colnames(DEG_result)[which(colnames(DEG_result)=='pvalue')] = 'p.value'
    colnames(DEG_result)[which(colnames(DEG_result)=='padj')] = 'p.adj'
  }
  
  if (DEG_type == 'edgeR'){
    #colnames(DEG_result)[which(colnames(DEG_result)=='log2FoldChange')] = 'logFC'
    colnames(DEG_result)[which(colnames(DEG_result)=='PValue')] = 'p.value'
    colnames(DEG_result)[which(colnames(DEG_result)=='FDR')] = 'p.adj'
  }
  
  if (DEG_type == 'limma'){
    #colnames(DEG_result)[which(colnames(DEG_result)=='log2FoldChange')] = 'logFC'
    colnames(DEG_result)[which(colnames(DEG_result)=='P.Value')] = 'p.value'
    colnames(DEG_result)[which(colnames(DEG_result)=='adj.P.Val')] = 'p.adj'
  }
  
  if (is.null(thre_logFC)){
    thre_logFC <- with(DEG_result,mean(abs(logFC)) + 2*sd(abs(logFC)))
  }
  if (usePadj){
    #u1=rownames(DEG_result[with(DEG_result,logFC>thre_logFC & p.adj<thre_p),])
    #d1=rownames(DEG_result[with(DEG_result,logFC<-thre_logFC & p.adj<thre_p),])
    a = subset(DEG_result, logFC >= (thre_logFC) &  p.adj<thre_p)
    b = subset(DEG_result, logFC < (-thre_logFC) &  p.adj<thre_p)
    
  }else{
    #u1=rownames(DEG_result[with(DEG_result,logFC>thre_logFC & p.value<thre_p),])
    #d1=rownames(DEG_result[with(DEG_result,logFC< -thre_logFC & p.value<thre_p),])
    a = subset(DEG_result, logFC >= (thre_logFC) &  p.value<thre_p)
    b = subset(DEG_result, logFC < (-thre_logFC) &  p.value<thre_p)
  }
  
  a <- a[,c('logFC','p.value','p.adj')]
  a$direction <- 'UP'
  a <- a[order(a$logFC,decreasing = TRUE),]
  
  b <- b[,c('logFC','p.value','p.adj')]
  b$direction <- 'DOWN'
  b <- b[order(b$logFC,decreasing = FALSE),]
  
  DEGdf <- rbind(a,b)
  
  return(DEGdf)
}

volcanoPlot <- function(pval_df, p_val = "padj", clustern = 'cluster', fcn='logfc', p.threshold = 0.05, text_table = NULL, 
                        fc.threshold=2, change=NULL, c.color=NULL){
  #pval_df <- DEG_edgeR.sub
  plotdf <- data.frame(logfc = pval_df[,fcn], pvalue = pval_df[,p_val], cluster = pval_df[,clustern])
  plotdf <- plotdf[!is.infinite(plotdf$logfc), ]
  
  if (is.null(change)){
    plotdf$change = as.factor(ifelse(plotdf$pvalue < p.threshold & abs(round(plotdf$logfc,1)) >= fc.threshold,
                                     ifelse(round(plotdf$logfc,1) >= fc.threshold ,'UP','DOWN'),'NOT'))
    change.color <- c('blue',"black", "red")
  }else{
    plotdf$change <- pval_df$change
    change.color <- c.color
  }
  
  gg_logFC = ggplot(plotdf, aes(x = logfc, y = -log10(pvalue))) + 
    geom_point(aes(color = change),alpha = 0.5, size=1.75) + theme_bw() + 
    theme(legend.position = 'right', text = element_text(size = 25)) +
    geom_text_repel(data = text_table, aes(x = logFC, y = -log10(FDR),label=gene), force = 10, size=5) +
    #geom_label_repel(data = text_table, aes(x = logFC, y = -log(p.adj),label=X),color="black")+
    scale_color_manual(values = change.color) +
    geom_vline(xintercept=c(-fc.threshold,fc.threshold),lty=3,col="black",lwd=1) + #添加横线|FoldChange|>2
    geom_hline(yintercept = -log10(p.threshold),lty=3,col="black",lwd=1) #添加竖线padj<0.05
  return(gg_logFC)
}

## for correlation between bulk RNA
getCorrandP <- function(sig.result, cell.result, corrm = 'pearson'){
  #sig.result, cell.result: row - sample, column - signaturescore/celltypefraction
  sig.result.sub <- sig.result[rownames(cell.result),]
  
  sigs <- colnames(sig.result)
  cells <- colnames(cell.result)
  
  scatter.df <- data.frame(sig=NULL, celltype=NULL, corr = NULL, rsq=NULL, pvalue=NULL)
  for (sig in sigs) {
    #i=12
    for (cell in cells){
      x = sig.result.sub[,sig]
      y = cell.result[,cell]
      
      cor_value <- cor(x, y, method = corrm)
      model<-lm(y~x)
      r2 <- summary(model)[["r.squared"]]
      f <- summary(model)$fstatistic
      pvalue <- pf(f[1], f[2], f[3], lower.tail=F)
      tmp <- data.frame(sig=sig, celltype=cell, corr = cor_value ,rsq = r2, pvalue = pvalue)
      
      scatter.df <- rbind(scatter.df, tmp)
    }
  }
  return(scatter.df)
}
