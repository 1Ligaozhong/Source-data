
  dir.create("PDFs")
  dir.create("PDFs/图片")
  dir.create("files")
  dir.create("files/文件")
  dir.create("origin_datas")
  dir.create("origin_datas/TCGA")
  dir.create("origin_datas/GEO")
}

library(ggsci)
scales::show_col(pal_nejm()(9)[1:9])
scales::show_col(pal_lancet()(9)[1:9])
scales::show_col(pal_npg()(9)[1:9])
scales::show_col(pal_d3(palette = 'category10')(8)[1:8])

###############
subtype.color=pal_nejm(alpha = 1)(9)[c(3,4,2)]
names(subtype.color)=c('C1','C2','C3')
risk.group.color=pal_npg()(8)[c(5,2)]
names(risk.group.color)=c('High','Low')

heatmap.color=(pal_jama(palette = c('default'), alpha =1)(9))[c(1,2)]
mycolor <- c(pal_npg()(9)[c(2,3,5,1,9)])
scales::show_col(mycolor)

corrmat.color=c("#009E73","white","#E69F00")

PPARs.color=pal_nejm()(9)[c(7,4,8)]
############################################################################
########## TCGA-STAD 数据处理
## 无虚线
ggplotKMCox=function(dat,title='Groups',labs=NULL,add_text=NULL,pal=NULL){
  library(ggplot2)
  library(ggsci)
  colnames(dat)=c('time','status','groups')
  #sdf<-survdiff(Surv(time,status) ~ groups,data=dat)
  #print((sdf))
  #summary(sdf)
  #p<-pchisq(sdf$chisq,length(sdf$n)-1,lower.tail=FALSE)
  sf<-survfit(Surv(time,status) ~ groups,data=dat)
  if(is.null(pal)){
    surv=survminer::ggsurvplot(sf, data = dat
                               , palette =pal_lancet()(9)[c(2,4,3,1,5:6)]
                               ,pval = TRUE, surv.median.line='hv'
                               # ,conf.int = T
                               , xlab="Time(years)"
                               ,linetype = "strata"
                               ,conf.int.style ='step'
                               , pval.coord=c(0, 0.2) #Add p-value 
                               , risk.table.title=""
                               , risk.table = TRUE
                               , legend.title = title
                               , legend.labs = labs
    )
    
  }else{
    if(pal=='npg'){
      surv=survminer::ggsurvplot(sf, data = dat
                                 , palette = pal_npg()(9)[c(1,2,3,4:9)] #jco palette 
                                 ,pval = TRUE, surv.median.line='hv'
                                 # ,conf.int = T
                                 , xlab="Time(years)"
                                 ,linetype = "strata"
                                 ,conf.int.style ='step'
                                 , pval.coord=c(0, 0.2), #Add p-value 
                                 risk.table = TRUE, 
                                 legend.title = title
                                 ,legend.labs = labs
      )
    }else{
      surv=survminer::ggsurvplot(sf, data = dat
                                 , palette = pal #jco palette 
                                 ,pval = TRUE, surv.median.line='hv'
                                 # ,conf.int = T
                                 ,linetype = "strata"
                                 , xlab="Time(years)"
                                 ,conf.int.style ='step'
                                 , pval.coord=c(0, 0.2), #Add p-value 
                                 risk.table = TRUE, 
                                 legend.title = title
                                 ,legend.labs = labs
      )
    }
    
  }
  
  p1=surv$plot+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                ,axis.text.x=element_blank()
                                ,axis.title.x=element_blank()
                                ,plot.margin=unit(c(0.2, 0.2, 0, 0.1), "inches")
                                #,axis.title.y=element_blank()
                                ,legend.position=c(1,1), legend.justification=c(1,1)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  #p1=p1+text()
  #tms=data.frame(Group=tms.gp,value=tms.tps,Attribute=rep(data_m[1,1],length(tms.gp))
  #               ,ymax=rep(max(ylim),length(tms.gp)))
  #p4=p4+geom_text(data=tms,aes(x=Group, y=ymax, label=value),color="yellow")
  if(!is.null(add_text)){
    text.tb=surv$data.survplot[1,]
    text.tb[1,1]=0
    text.tb[1,5]=0
    text.tb$Text=add_text
    p1=p1+geom_text(data=text.tb,aes(x=time, y=surv, label=Text),color="yellow",hjust =0)
  }
  
  p2=surv$table+theme_bw()+theme(axis.text.y=element_text(family="Times",face="plain")
                                 #,axis.text.x=element_blank()
                                 #,axis.title.x=element_blank()
                                 #,axis.title.y=element_blank()
                                 ,plot.margin=unit(c(0, 0.2, 0.2, 0.1), "inches")
                                 ,plot.title=element_blank()
                                 ,legend.position=c(1,1), legend.justification=c(1,1)
                                 #,legend.background = element_rect(fill = NA, colour = NA)
                                 ,legend.title = element_text(family="Times",face="plain")
                                 ,legend.text = element_text(family="Times",face="plain"))
  p1=surv$plot
  p2=surv$table
  g2=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(0.9,0.45),align = "v")
  return(g2)
}

############### risk score
plotRiskScoreModel_use=function(riskScore,dat,time,event,cutoff,pal=NULL,hetTitle='z-score of expression',hetColor=c('green','black','red')){
  srt.inds=order(riskScore)
  dat=dat[srt.inds,]
  time=time[srt.inds]
  event=event[srt.inds]
  riskScore=riskScore[srt.inds]
  library(ggplot2)
  dt1=data.frame(V1=1:length(riskScore),V2=riskScore,RiskType=ifelse(riskScore>cutoff,'High','Low')) 
  # p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_bar(stat = 'identity', position = 'dodge')+ggsci::scale_fill_npg()+theme_bw()
  # p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_point(stat = 'identity', position = 'identity')+ggsci::scale_fill_npg()+theme_bw()
  p1=ggplot(dt1, aes(x = V1, y = V2, colour = RiskType,fill=RiskType)) +geom_point(stat = 'identity', position = 'identity')+scale_colour_manual(values = pal,aesthetics = c("colour", "fill"))+theme_bw()
  p1=p1+geom_hline(aes(yintercept=cutoff),colour='black',linetype="dashed")
  p1=p1+ylab('RiskScore')+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
                                ,axis.title.x=element_blank(),legend.position=c(1,0), legend.justification=c(1,0)
                                ,legend.background = element_rect(fill = NA, colour = NA)
                                ,plot.margin=unit(c(0.1, 0.1, 0, 0.1), "inches")
                                ,legend.title = element_text(family="Times",face="plain")
                                ,legend.text = element_text(family="Times",face="plain"))
  
  dt2=data.frame(V1=1:length(riskScore),V2=time,Status=ifelse(event==1,'Dead','Alive'))  
  p2=ggplot(dt2, aes(x = V1, y = V2, colour = Status,shape =Status)) +geom_point()+ggsci::scale_fill_npg()+theme_bw()
  p2=p2+ylab('Time')+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_blank()
                           ,axis.title.x=element_blank(),legend.position=c(1,1), legend.justification=c(1,1)
                           ,legend.background = element_rect(fill = NA, colour = NA)
                           ,plot.margin=unit(c(0, 0.1, 0, 0.1), "inches")
                           ,legend.title = element_text(family="Times",face="plain")
                           ,legend.text = element_text(family="Times",face="plain"))
  
  data=as.data.frame(scale(dat))
  hc.r = hclust(dist(t(data)))
  data=data[,hc.r$order]
  data$ID <- 1:nrow(dat)
  #colnames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  colnames(data_m)=c('ID','V1','V2')
  data_m$V2[which(data_m$V2>mean(data_m$V2)+3*sd(data_m$V2))]=mean(data_m$V2)+3*sd(data_m$V2)
  data_m$V2[which(data_m$V2<mean(data_m$V2)-3*sd(data_m$V2))]=mean(data_m$V2)-3*sd(data_m$V2)
  
  data_m$V1=mg_str_outline(data_m$V1,isCut = T,n=50)
  #print(data_m$V1)
  #print(head(data_m))
  #data_m[1:20,]
  p3 <- ggplot(data_m, aes(x=ID,y=V1)) 
  p3 <- p3 + geom_tile(aes(fill=V2))
  p3=p3+scale_fill_gradient2(low = hetColor[1],mid=hetColor[2], high = hetColor[3])
  p3=p3+theme_bw()
  p3=p3+ labs(fill=hetTitle) 
  #p3=p3+guides(fill=guide_legend(title="New Legend Title"))
  p3=p3+xlab('Samples')
  #p3 <- p3 + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
  p3=p3+theme(axis.text.y=element_text(family="Times",face="plain")
              ,axis.text.x=element_blank()
              #,axis.title.x=element_blank()
              ,axis.title.y=element_blank()
              ,legend.position='bottom'
              #,legend.justification=c(1,1)
              #,legend.background = element_rect(fill = NA, colour = NA)
              ,plot.margin=unit(c(0, 0.1, 0.1, 0.1), "inches")
              ,legend.title = element_text(family="Times",face="plain")
              ,legend.text = element_text(family="Times",face="plain"))
  
  g1=ggpubr::ggarrange(p1,p2,p3, ncol = 1, nrow = 3,heights = c(0.5,0.5,1),align = "v")
  return(g1)
}

################ timeROC
ggplotTimeROC_use=function(time,status,score,mks=c(1,3,5)){
  #time=g.os
  #status=g.ev
  #score=as.numeric(cpm.score)
  #cx=coxRun(data.frame(time,status,score))
  #if(cx[1]<=1){
  #  score=-1*score
  #}
  roc.tm=mg_surv_pROC(time,status,score,mks)
  print('roc.tm')
  print((roc.tm))
  library(survival)
  library(ggplot2)
  mks=mg_predict_time_ymd(time,mks)
  print(mks)  
  ROC.DSST=timeROC::timeROC(T=time,
                            delta=status
                            ,marker=score,
                            cause=1,weighting="marginal",
                            times=mks,
                            iid=TRUE)
  print(ROC.DSST)
  mks=mks[which(!is.na(ROC.DSST$AUC)&ROC.DSST$AUC>0)]
  print(mks)
  if(length(mks)>0){
    if(max(ROC.DSST$AUC)<0.5){
      score=-1*score
    }
    ROC.DSST=timeROC::timeROC(T=time,
                              delta=status
                              ,marker=score,
                              cause=1,weighting="marginal",
                              times=mks,
                              iid=TRUE)
    print(ROC.DSST$times)
    if(max(ROC.DSST$times)<20){
      lb=paste0(ROC.DSST$times,'-Years')
    }else if(max(ROC.DSST$times)<365){
      lb=paste0(round(ROC.DSST$times/12,0),'-Years')
    }else{
      lb=paste0(round(ROC.DSST$times/365,0),'-Years')
    }
    
    lbs=paste0(lb,',AUC=',round(ROC.DSST$AUC,2),',95%CI(',paste0(round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,1]/100,2),'-',
                                                                 round(confint(ROC.DSST,level = 0.95,na.rm=T)$CI_AUC[,2]/100,2)),')')
    #roc.tm=ROC.DSST$times[which(ROC.DSST$times>0)]
    
    #p.dat=rbind()
    #for(i in which(ROC.DSST$times>0)){
    #los=lowess(ROC.DSST$FP[,i], y=ROC.DSST$TP[,i], f = 1/3, iter = 100)
    #los$x=c(0,los$x,1)
    #los$y=c(0,los$y,1)
    # p.dat=rbind(p.dat,data.frame(los$x, y=los$y,rep(lbs[i],length(los$y)),stringsAsFactors = F))
    #}
    
    p.dat=rbind()
    print(length(roc.tm))
    for(i in 1:length(roc.tm)){
      #print(i)
      r1=roc.tm[[i]]
      x1=1-r1$specificities
      y1=r1$sensitivities
      #print(cbind(1-r1$specificities,r1$sensitivities))
      nx1=unique(x1)
      ny1=c()
      for(x in unique(x1)){
        x.inds=which(x1==x)
        if(length(x.inds)>0&x<0.5){
          ny1=c(ny1,min(y1[x.inds]))
        }else if(length(x.inds)>0){
          ny1=c(ny1,max(y1[x.inds]))
        }else{
          ny1=c(ny1,y1[x.inds][1])
        }
      }
      #print(cbind(nx1,ny1))
      p.dat=rbind(p.dat,data.frame(x=nx1, y=ny1,rep(lbs[i],length(nx1)),stringsAsFactors = F))
    }
    colnames(p.dat)=c('V1','V2','Type')
    p.dat=as.data.frame(p.dat)
    
    p1=ggplot(p.dat, aes(x=V1,y=V2, fill=Type))
    p1=p1+geom_line(aes(colour=Type),lwd=1.1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    # p1=p1+stat_smooth(aes(colour=Type),se = FALSE, size = 1)+theme_bw()+xlab('False positive fraction')+ylab('True positive fraction')
    
    p1=p1+scale_colour_manual(values = c(pal_npg(alpha =0.8)(9)[c(8,4,3,5)]))
    
    p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),axis.text.x=element_text(family="Times",face="plain")
                ,axis.title.x=element_text(family="Times",face="plain"),axis.title.y=element_text(family="Times",face="plain")
                ,plot.title=element_blank()
                ,plot.margin=unit(c(0.1, 0.1, 0.1, 0.1), "inches")
                ,legend.position=c(1,0)
                ,legend.justification=c(1,0)
                ,legend.background = element_rect(fill = NA, colour = NA)
                ,legend.title = element_text(family="Times",face="plain")
                ,legend.text = element_text(family="Times",face="plain"))
    return(p1)
  }else{
    return(mg_getplot_bank('No data plot by ROC!'))
  }
}

################ 
plotCoxModel_Batch_use=function(riskScore,dat,time,event,cutoff,pal=NULL,title='Groups',hetTitle='z-score of expression',hetColor=c('green','black','red'),mks=c(1,3,5),labs=NULL){
  g1=plotRiskScoreModel_use(riskScore,dat,time,event,cutoff,pal,hetTitle,hetColor)
  dat=data.frame(time,event,ifelse(riskScore>cutoff,'High','Low'),stringsAsFactors = F)
  dat=dat[order(dat[,3]),]
  cx=coxRun(dat = data.frame(time,event,riskScore))
  txt=paste0('HR=',round(cx[2],2),' 95CI%(',round(cx[3],2),'-',round(cx[4],2),')')
  dat[,3]=factor(dat[,3],levels = c('High','Low'))
  print(unique(dat[,3]))
  g3=ggplotKMCox(dat,title=title,add_text = txt,pal,labs = labs)
  g2=ggplotTimeROC_use(time,event,riskScore,mks)
  g23=ggpubr::ggarrange(g2,g3, ncol = 1, nrow = 2,heights = c(1,1)
                        #,align = "hv"
                        ,labels = toupper(letters)[2:3])
  g23.row=ggpubr::ggarrange(g3,g2,ncol = 2, nrow = 1,heights = c(1,1)
                            #,align = "hv"
                            ,labels = toupper(letters)[2:3])
  
  gal=ggpubr::ggarrange(g1,g23, ncol = 2, nrow = 1,widths =  c(1.5,1)
                        #,align = "hv"
                        ,labels = c(toupper(letters)[1],''))
  return(list(gal,g23,g23.row))
}

##################
createCoxModel_use=function(dat,time,event,isStep=F,direction=c("both", "backward", "forward")[1],check=T){
  cls=colnames(dat)
  dat1=cbind(dat,time,event)
  colnames(dat1)=c(paste0('g',1:ncol(dat)),'time','status')
  dat1=as.data.frame(dat1)
  if(ncol(dat)>nrow(dat)&check){
    print('gene count > sample count')
    return(NULL)
  }
  #nas=apply(dat1, 1, function(x){
  #  return(sum(is.na(x)))
  #})
  #dat1=dat1[which(nas==0),]
  
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(dat1)[1:ncol(dat)],collapse = '+')))
  library(survival)
  cox <- coxph(fmla, data = dat1)
  
  if(isStep){
    tryCatch({
      cox=step(cox,direction =direction,steps = 10000)
    },error = function(e) {
      print(conditionMessage(e))
      return(NULL)
    })
  }
  #score=predict(cox,data=dat1)
  sig.genes=cls[as.numeric(gsub('g','',names(cox$coefficients)))]
  fls=c('RiskScore=')
  for(i in 1:length(sig.genes)){
    if(cox$coefficients[i]>0){
      fls=c(fls,'+',round(cox$coefficients[i],3),'*',sig.genes[i])
    }else{
      fls=c(fls,round(cox$coefficients[i],3),'*',sig.genes[i])
    }
  }
  score=predictRiskScore(cox$coefficients,sig.genes,dat)
  return(list(Cox=cox,Score=score,Genes=sig.genes,Coef=cox$coefficients,fmla=paste0(fls,collapse = '')))
}

mg_plot_lasso_use <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
  if(is.null(lambda)){
    lmda=cv_fit$lambda.min
  }else{
    lmda=lambda
  }
  fit.coef=fit$beta[(apply(fit$beta,1,function(x){
    return(sum(x!=0))
  })>0),]
  
  fit.coef=as.matrix(fit.coef)
  colnames(fit.coef)=fit$lambda
  #fit$lambda==cv_fit$lambda
  library(ggplot2)
  dat=data.table::melt(t(as.matrix(fit.coef)))
  dat_z=dat[which(dat$value==0),]
  dat=dat[which(dat$value!=0),]
  dat.sv=rbind()
  for (u in unique(dat_z[,2])) {
    t.z=dat_z[which(dat_z[,2]==u),1]
    t.zx=max(t.z)
    dat.sv=rbind(dat.sv,c(t.zx,u,0))
    t.zn=min(t.z)
    if(t.zx!=t.zn){
      dat.sv=rbind(dat.sv,c(t.zn,u,0))
    }
  }
  colnames(dat.sv)=colnames(dat_z)
  #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
  dat=crbind2DataFrame(rbind(dat,dat.sv))
  mn=min(-log(dat$Var1))
  mx=max(-log(dat$Var1))
  if(show_text){
    mx=(mx-mn)*0.1+mx
  }
  p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
  p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
  if(show_text){
    fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
    for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
    p=p+ggrepel::geom_label_repel(
      aes(label = Var2,color=Var2),
      data = for_label,hjust = 0
    )
  }
  p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
  p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
  tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                 ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
  p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Parial Likelihood Deviance')+
    geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
    geom_point(aes(colour=col))
  p1=p1+theme_bw()+theme(legend.position = "none")
  # gal=p+p1
  gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                        ,align = "hv"
                        ,labels = figLabels)
  return(gal)
}

mg_lasso_cox_use=function(dat,time,event,nfolds=3,lambda.min=T,show_text=T,figLabels=c('A','B')){
  library("glmnet") 
  library('survival')
  t.inds=which(!is.na(time)&!is.na(event)&time>0)
  dat=dat[t.inds,]
  time=as.numeric(time[t.inds])
  event=as.numeric(event[t.inds])
  y=Surv(time,event)
  set.seed(333333333.3333)
  fit1_cv = cv.glmnet(as.matrix(dat), y, family = "cox", nfolds=nfolds)
  fit<-glmnet(dat, y, family = "cox")
  if(lambda.min){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=fit1_cv$lambda.1se
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  genes=row.names(coefficients)[Active.Index]
  Active.coefficients<-coefficients[Active.Index]  
  g=mg_plot_lasso(fit,fit1_cv,lambda = lambda,show_text=show_text,figLabels=figLabels)
  return(list(Mode1=fit,Model2=fit1_cv,Genes=genes,Coef=Active.coefficients,lambda=lambda,plot=g))
}

plot_ggviolin=function(data,melt=F,group.col=NULL,ylab='',xlab='',leg.title='Group'){
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  #s设置比较组
  group=levels(factor(data_m$Group))
  data_m$Group=factor(data_m$Group, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  library(ggpubr)
  p1=ggviolin(data_m, x="Group", y="value", fill = "Group", 
              xlab=xlab, ylab=ylab,
              legend.title=leg.title,
              add = "boxplot", add.params = list(fill="white"))+ 
    stat_compare_means(comparisons = my_comparisons,label='p.signif',method='wilcox.test')
  if(!is.null(group.col)){
    p1=p1+scale_fill_manual(values = group.col)
  }else{
    p1=p1+scale_fill_manual(values = c(pal_lancet('lanonc',alpha =0.6)(9)[c(2,1,3,4:9)]))
  }
  return(p1)
}

######################################################################################################
######################################################################################################

ann.pcg=ann[which(ann$gene_type=='protein_coding'),]
dim(ann.pcg)

######################### pancancer

# pancan.exp=c()
# for(f in files){
#   df=data.table::fread(file = f,data.table=T,sep='\t',check.names = F)
#   df=data.frame(df,check.names = F,stringsAsFactors = F)
#   rownames(df)=df$Tags
#   df=as.matrix(df[,-1])
#   pancan.exp=cbind(pancan.exp,df)
# }
################# TCGA-STAD

# stad.tcga.exp.tpm=data.frame(stad.tcga.exp.tpm,check.names = F)
# rownames(stad.tcga.exp.tpm)=stad.tcga.exp.tpm$Tags
# stad.tcga.exp.tpm=stad.tcga.exp.tpm[,-1]
# range(stad.tcga.exp.tpm)
# stad.tcga.exp.tpm=log2(stad.tcga.exp.tpm+1)
# stad.tcga.exp.tpm=exp_ensg2symbol(stad.tcga.exp.tpm)
# save(stad.tcga.exp.tpm,file = 'origin_datas/TCGA/stad.tcga.exp.tpm.RData')
load('origin_datas/TCGA/stad.tcga.exp.tpm.RData')
stad.tcga.exp=stad.tcga.exp.tpm
stad.tcga.exp[1:4,1:4]
write.table(stad.tcga.exp,file = 'files/文件/stad.tcga.exp.txt',sep = '\t',row.names = T,col.names = T,quote = F)

#######################################
########### GDC 下载最新的临床病理数据
stad.tcga.cli=readMatrix('origin_datas/TCGA/TCGA-STAD_clinical.txt',row = F)
dim(stad.tcga.cli)
write.table(stad.tcga.cli,file = 'files/文件/stad.tcga.cli.txt',sep = '\t',row.names = T,col.names = T,quote = F)

stad.tcga.cli$Tumor_Sample_Barcode=stad.tcga.cli$A0_Samples
stad.tcga.cli$os=stad.tcga.cli$A1_OS
stad.tcga.cli$os_status=stad.tcga.cli$A2_Event
dim(stad.tcga.cli)
stad.tcga.cli=stad.tcga.cli[,c(172:174,1:171)]
names(table(stad.tcga.cli$os_status))

stad.tcga.cli$os_status[which(stad.tcga.cli$os_status=='Dead')]=1
stad.tcga.cli$os_status[which(stad.tcga.cli$os_status=='Alive')]=0

stad.tcga.cli=data.frame(SampleID=paste0(stad.tcga.cli$Tumor_Sample_Barcode,"-01"),stad.tcga.cli,stringsAsFactors = F,check.names = F)
rownames(stad.tcga.cli)=stad.tcga.cli$SampleID

############## TCGA-STAD 表达谱数据
dim(stad.tcga.exp)
range(stad.tcga.exp)
table(substr(colnames(stad.tcga.exp),14,15))

######### -01 原发肿瘤样本
stad.tcga.exp.t.inds=grep('-01[A-Z]?',colnames(stad.tcga.exp))
length(stad.tcga.exp.t.inds) ### 414

stad.tcga.t.exp=stad.tcga.exp[,stad.tcga.exp.t.inds]
dim(stad.tcga.t.exp)

##############################
stad.tcga.t.cli=stad.tcga.cli[grep("-01[AD]?",stad.tcga.cli$SampleID),]
dim(stad.tcga.t.cli)

stad.tcga.t.cli=dplyr::distinct(stad.tcga.t.cli,SampleID,.keep_all=TRUE)
rownames(stad.tcga.t.cli)=stad.tcga.t.cli$SampleID
dim(stad.tcga.t.cli)

##########
setdiff(colnames(stad.tcga.t.exp),stad.tcga.t.cli$SampleID)

stad.tcga.t.cli=stad.tcga.t.cli[colnames(stad.tcga.t.exp),]
dim(stad.tcga.t.cli)

stad.tcga.t.cli=dplyr::rename(stad.tcga.t.cli,c('OS.time'='A1_OS'
                                                ,'OS'='A2_Event'))
table(stad.tcga.t.cli$OS)

stad.tcga.t.cli$OS[stad.tcga.t.cli$OS=='Alive']=0
stad.tcga.t.cli$OS[stad.tcga.t.cli$OS=='Dead']=1
stad.tcga.t.cli$OS[stad.tcga.t.cli$OS=='']=NA
table(stad.tcga.t.cli$OS)

stad.tcga.t.cli=crbind2DataFrame(stad.tcga.t.cli)

table(stad.tcga.t.cli$OS)
table(is.na(stad.tcga.t.cli$OS))
table(is.na(stad.tcga.t.cli$OS.time))

dim(stad.tcga.t.cli)

range(stad.tcga.t.cli$OS.time)

table(stad.tcga.t.cli$OS.time<30)
table(stad.tcga.t.cli$OS.time>365*10)
table(stad.tcga.t.cli$OS.time>365*15)

###
stad.tcga.t.cli=stad.tcga.t.cli[intersect(stad.tcga.t.cli$SampleID,colnames(stad.tcga.t.exp)), ]
dim(stad.tcga.t.cli)
stad.tcga.t.exp=stad.tcga.t.exp[,stad.tcga.t.cli$SampleID]
dim(stad.tcga.t.exp)

## TCGA 生存数据
stad.tcga.t.exp.os=stad.tcga.t.cli[,c('OS','OS.time')]
table(stad.tcga.t.exp.os$OS.time>0)
table(stad.tcga.t.exp.os$OS.time>30)
table(stad.tcga.t.exp.os$OS.time<365*10)
dim(stad.tcga.t.exp.os)

stad.tcga.t.cli=stad.tcga.t.cli[rownames(stad.tcga.t.exp.os),]
stad.tcga.t.exp=stad.tcga.t.exp[,rownames(stad.tcga.t.exp.os)]
dim(stad.tcga.t.cli)
dim(stad.tcga.t.exp)

## TCGA 临床数据
clin.selected=c('A17_Age','A18_Sex','A3_T','A4_N','A5_M','A6_Stage','A7_Grade')
setdiff(clin.selected,colnames(stad.tcga.t.cli))

stad.tcga.t.exp.cli=stad.tcga.t.cli[,clin.selected]

stad.tcga.t.exp.cli.tnm=clean_TNMStage(st=stad.tcga.t.exp.cli$A3_T
                                       , sn = stad.tcga.t.exp.cli$A4_N
                                       , sm = stad.tcga.t.exp.cli$A5_M
                                       , ss = stad.tcga.t.exp.cli$A6_Stage
                                       , gd = stad.tcga.t.exp.cli$A7_Grade
                                       , sex = stad.tcga.t.exp.cli$A18_Sex
                                       , age = stad.tcga.t.exp.cli$A17_Age
                                       , age_cut = NULL)
row.names(stad.tcga.t.exp.cli.tnm)=rownames(stad.tcga.t.cli)
dim(stad.tcga.t.exp.cli.tnm)

colnames(stad.tcga.t.exp.cli.tnm)

table(stad.tcga.t.exp.cli.tnm$Clinical_T)
table(stad.tcga.t.exp.cli.tnm$Clinical_N)
table(stad.tcga.t.exp.cli.tnm$Clinical_M)
table(stad.tcga.t.exp.cli.tnm$Clinical_Stage)
table(stad.tcga.t.exp.cli.tnm$Clinical_Grade)

#############
tcga.t.cli_use=cbind(SampleCode=rownames(stad.tcga.t.exp.os),stad.tcga.t.exp.os,stad.tcga.t.exp.cli.tnm)
colnames(tcga.t.cli_use)
median(tcga.t.cli_use$Age,na.rm = T)

tcga.t.cli_use$Age1=ifelse(tcga.t.cli_use$Age>60,'>60','<=60')
tcga.t.cli_use$Status=ifelse(tcga.t.cli_use$OS==0,'Alive','Dead')

tcga.t.cli_use=dplyr::rename(tcga.t.cli_use,c('T Stage'='Clinical_T'
                                              ,'N Stage'='Clinical_N'
                                              ,'M Stage'='Clinical_M'
                                              ,'Stage'='Clinical_Stage'
                                              ,'Grade'='Clinical_Grade'
))
dim(tcga.t.cli_use)
tcga.t.cli_use$Stage=gsub("Stage ","",tcga.t.cli_use$Stage)

stad.tcga.t.cli=stad.tcga.t.cli[rownames(tcga.t.cli_use),]
stad.tcga.t.exp=stad.tcga.t.exp[,rownames(tcga.t.cli_use)]
stad.tcga.t.exp.os=stad.tcga.t.exp.os[rownames(tcga.t.cli_use),]
dim(stad.tcga.t.exp.os)

writeMatrix(stad.tcga.t.cli,outpath = 'files/文件/stad.tcga.t.cli.txt')
writeMatrix(stad.tcga.t.exp,outpath = 'files/文件/stad.tcga.t.exp.txt')

################################################################
################################################################
tcga.t.exp=stad.tcga.t.exp
################################################################
################################
dir.create('01_PPARs_Pan-cancer')
gene.selected=c('PPARA','PPARD','PPARG')

PPARs.exp=mg_mysql_getGenes_Exp_by_Symbol(geneSymbols = gene.selected
                                          ,db = c('TCGA'))
PPARs.cli=getTCGAClinicalBySamples(PPARs.exp$SampleCode)

table(PPARs.exp$Site)
PPARs.exp=PPARs.exp[which(PPARs.exp$Site %in% c('Normal','Tumor')),]
PPARs.exp[,gene.selected]=log2(PPARs.exp[,gene.selected]+1)
range(PPARs.exp[,c(4:5)])
head(PPARs.exp)
write.table(PPARs.exp,file = 'files/文件/Pancancer.PPARs.exp.txt',sep = '\t',row.names = F,col.names = T,quote = F)
write.table(PPARs.exp,file = '01_PPARs_Pan-cancer/Pancancer.PPARs.exp.txt',sep = '\t',row.names = F,col.names = T,quote = F)

###############
library(ggsci)
library(ggpubr)
sampleType.color=pal_lancet('lanonc')(9)[c(3,7)]

colnames(PPARs.exp)[4]
groupViolin(PPARs.exp[,c(3,2,4)]
                  ,melt = T
                  # ,ylim = c(0,150)
                  ,ylab = 'PPARA Expression'
                  ,group_col=sampleType.color
            )

p.all=list()
for(gene in gene.selected){
  gene=sym(gene)
  p=ggplot(PPARs.exp,
         aes(x = TCGA_CODE,
             y = !!gene, 
             fill = Site))+
    geom_boxplot()+
    # geom_jitter()+
    scale_fill_manual(values = sampleType.color) +
    xlab('') + #labs(fill = 'Cancers') +
    stat_compare_means(aes(group=Site), label = "p.signif",method='wilcox.test') +
    theme_bw() +
    labs(fill = 'Group') + ylab(paste('The expression of gene',gene,"\nlog2(TPM+1)")) +
    theme_classic()+
    theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times")
          ,axis.ticks = element_line(color = "black"))+labs(fill = "SampleType")
  p.all=c(p.all,list(p))
}

fig1a=ggarrange(plotlist = p.all
                , labels = c('A', 'B','C')
                , nrow = 3,common.legend=T)
fig1a

PPARs.t.exp=PPARs.exp[which(PPARs.exp$Site=='Tumor'),]
head(PPARs.t.exp)

PPARs.t.exp_melt=reshape2::melt(PPARs.t.exp[,c('TCGA_CODE',gene.selected)], id.vars=c("TCGA_CODE"))
head(PPARs.t.exp_melt)

fig1b=ggplot(PPARs.t.exp_melt,
             aes(x = TCGA_CODE,
                 y = value, 
                 fill = variable))+
  geom_boxplot()+
  stat_compare_means(aes(group=variable), label = "p.signif") +
  scale_fill_manual(values = PPARs.color) +
  xlab('') + 
  theme_bw() +
  labs(fill = 'Group') + ylab("log2(TPM+1)") +
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times")
        ,axis.ticks = element_line(color = "black")
        ,legend.position = 'top')+labs(fill = "PPARs")
fig1b

########## 肿瘤细胞系 CCLE数据库

table(ccle.smp$primary_disease)
table(ccle.smp$lineage)

ccle.smp=ccle.smp[which(ccle.smp$lineage!='unknown'),]
head(ccle.smp)

###### RNASeq表达谱


colnames(ccle.tpm.exp)=gsub(" \\(.*\\)","",colnames(ccle.tpm.exp))
ccle.tpm.exp=as.matrix(ccle.tpm.exp)
ccle.tpm.exp=t(ccle.tpm.exp)

dim(ccle.tpm.exp)
range(ccle.tpm.exp)
head(ccle.tpm.exp[,1:10])

ccle.tpm.exp[rownames(ccle.tpm.exp) %in% gene.selected,]

##########
setdiff(colnames(ccle.tpm.exp),rownames(ccle.smp))
ccle.smp=ccle.smp[intersect(colnames(ccle.tpm.exp),rownames(ccle.smp)),]
ccle.tpm.exp=ccle.tpm.exp[,rownames(ccle.smp)]
dim(ccle.tpm.exp)

##########################
PPARs.ccle.exp=ccle.tpm.exp[rownames(ccle.tpm.exp) %in% gene.selected,rownames(ccle.smp)]
PPARs.ccle.exp=PPARs.ccle.exp[gene.selected,]
head(PPARs.ccle.exp)

######### 
PPARs.ccle.exp_use=crbind2DataFrame(cbind(cell.line=ccle.smp$lineage
                                          ,t(PPARs.ccle.exp)))
dim(PPARs.ccle.exp_use)
head(PPARs.ccle.exp_use)
write.table(PPARs.ccle.exp_use,file = 'files/文件/ccle.PPARs.exp.txt',sep = '\t',quote = F,row.names = T,col.names = T)
write.table(PPARs.ccle.exp_use,file = '01_PPARs_Pan-cancer/ccle.PPARs.exp.txt',sep = '\t',quote = F,row.names = T,col.names = T)

##### 排序
library(dplyr)
head(PPARs.ccle.exp_use)
med=PPARs.ccle.exp_use %>% group_by(cell.line) %>% 
  dplyr::mutate(sum = rowSums(across(PPARA:PPARG), na.rm = T)) %>% 
  summarise_at(vars(sum),median,na.rm=T)
med=data.frame(med)
head(med)
dim(med)

PPARs.ccle.exp_use$cell.line=factor(PPARs.ccle.exp_use$cell.line
                              , levels=med[order(med[,"sum"],decreasing = T),"cell.line"])
head(PPARs.ccle.exp_use)

PPARs.ccle.exp_melt=reshape2::melt(PPARs.ccle.exp_use, id.vars=c("cell.line"))
head(PPARs.ccle.exp_melt)

fig1c=ggplot(PPARs.ccle.exp_melt,
         aes(x = cell.line,
             y = value, 
             fill = variable))+
  geom_boxplot()+
  stat_compare_means(aes(group=variable), label = "p.signif") +
  scale_fill_manual(values = PPARs.color) +
  xlab('') + 
  theme_bw() +
  labs(fill = 'Group') + ylab("log2(TPM+1)") +
  theme_classic()+
  theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times")
        ,axis.ticks = element_line(color = "black")
        ,legend.position = 'top')+labs(fill = "PPARs")
fig1c

fig1=mg_merge_plot(fig1a,fig1b,fig1c
                   ,nrow = 3,ncol = 1,heights = c(3,1,1.2)
                   ,labels = c("",'D','E'))
savePDF('PDFs/Fig1.pdf',fig1,height = 16,width = 12)
savePDF('01_PPARs_Pan-cancer/Fig1.pdf',fig1,height = 16,width = 12)

######################
dir.create('02_PPARs_Clinical')
########### 配对样本中的表达
get_tcga_exp_pair=function(genes,CODE){
  # genes=gene.selected
  # CODE='STAD'
  df_expr=mg_mysql_getGenes_Exp_by_Symbol(geneSymbols = genes,db = 'TCGA')
  df_expr=df_expr[which(df_expr$TCGA_CODE==CODE),]
  dim(df_expr)
  head(df_expr)
  normal.smp=df_expr[grep('-11$',df_expr[,'SampleCode']),'SampleCode']
  tumor.smp=gsub('-11$','-01',normal.smp)
  
  inds=which(!is.na(match(tumor.smp,df_expr[,'SampleCode'])))
  normal.smp=normal.smp[inds]
  tumor.smp=tumor.smp[inds]
  all(substr(normal.smp,1,12)==substr(tumor.smp,1,12))
  
  pair.smp=c(normal.smp,tumor.smp)
  
  df_expr=df_expr[match(pair.smp,df_expr$SampleCode),]
  df_expr$Patient=substr(df_expr$SampleCode,1,12)
  return(df_expr)
}

PPARs.exprs.pair=get_tcga_exp_pair(genes = gene.selected,CODE = 'STAD')
PPARs.exprs.pair[,gene.selected]=log2(PPARs.exprs.pair[,gene.selected]+1)
head(PPARs.exprs.pair)

PPARs.exprs.pair_m=PPARs.exprs.pair[,-c(1:2)] %>% group_by(Patient,Site) %>%
  tidyr::gather(PPARA:PPARG, key = gene, value = exp) %>%
  reshape2::dcast(gene+Patient~Site,value.var='exp',mean)
head(PPARs.exprs.pair_m)

p.all=list()
for(gene in gene.selected){
  # gene='PPARA'
  df=PPARs.exprs.pair_m[which(PPARs.exprs.pair_m$gene==gene),]
  p=ggpaired(df[,c(3,4)], cond1 = 'Normal', cond2 = 'Tumor', fill = "condition", palette = sampleType.color,
             legend.title="Group",xlab="",ylab = paste('The expression of gene',gene,"\nlog2(TPM+1)"))+
    stat_compare_means(method = 'wilcox.test',paired = TRUE, label = "p.format",)+labs(fill='SampleType')
  p.all=c(p.all,list(p))
}

mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = LETTERS[1:3],common.legend = T)

####################################################
################# STAD
####################################################
stad.tcga.PPARs.exp=PPARs.exp[which(PPARs.exp$TCGA_CODE=='STAD'),]
rownames(stad.tcga.PPARs.exp)=stad.tcga.PPARs.exp$SampleCode
all(rownames(stad.tcga.PPARs.exp)==colnames(stad.tcga.exp))

stad.tcga.PPARs.exp=stad.tcga.PPARs.exp[colnames(stad.tcga.exp),]

ind1=which(stad.tcga.PPARs.exp$Site=='Tumor')
stad.tcga.t.PPARs.exp=stad.tcga.PPARs.exp[ind1,]
dim(stad.tcga.t.PPARs.exp)

#################
setdiff(stad.tcga.t.PPARs.exp$SampleCode,rownames(stad.tcga.t.cli))
mg_venn_plot(list(Exprs=stad.tcga.t.PPARs.exp$SampleCode
             ,Clin=rownames(stad.tcga.t.cli))
             ,fill=mg_colors[1:2])

###########
dim(tcga.t.cli_use)
dim(stad.tcga.t.PPARs.exp)
smp=intersect(rownames(tcga.t.cli_use),stad.tcga.t.PPARs.exp$SampleCode)
length(smp) ######### 375例样本

stad.tcga.t.cli=stad.tcga.t.cli[smp,]
stad.tcga.t.exp=stad.tcga.t.exp[,smp]

tcga.t.cli_use=tcga.t.cli_use[smp,]
stad.tcga.t.PPARs.exp=stad.tcga.t.PPARs.exp[match(smp,stad.tcga.t.PPARs.exp$SampleCode),]
dim(stad.tcga.t.PPARs.exp)
dim(tcga.t.cli_use)

table(tcga.t.cli_use$Stage)
table(tcga.t.cli_use$Grade)

tcga.t.cli_2=tcga.t.cli_use
tcga.t.cli_2=data.frame(tcga.t.cli_2)

tcga.t.cli_2$Age1=factor(tcga.t.cli_2$Age1,levels = c('<=60','>60'))
tcga.t.cli_2$Gender=as.factor(tcga.t.cli_2$Gender)

tcga.t.cli_2$Status=ifelse(tcga.t.cli_2$OS==1,'Dead','Alive')

tcga.t.cli_2$T.Stage=ifelse(is.na(tcga.t.cli_2$T.Stage),NA,ifelse(tcga.t.cli_2$T.Stage %in% c('T1','T2'),'T1+T2','T3+T4'))
tcga.t.cli_2$N.Stage=ifelse(is.na(tcga.t.cli_2$N.Stage),NA,ifelse(tcga.t.cli_2$N.Stage %in% c('N0'),'N0','NX'))
tcga.t.cli_2$M.Stage=ifelse(is.na(tcga.t.cli_2$M.Stage),NA,ifelse(tcga.t.cli_2$M.Stage %in% c('M0'),'M0','M1'))
tcga.t.cli_2$Stage=ifelse(!is.na(tcga.t.cli_2$Stage),ifelse(tcga.t.cli_2$Stage %in% c('I','II'),'I+II','III+IV'),NA)
tcga.t.cli_2$Grade=ifelse(!is.na(tcga.t.cli_2$Grade),ifelse(tcga.t.cli_2$Grade %in% c('G1','G2'),'G1+G2','G3+G4'),NA)

colnames(tcga.t.cli_2)

p.all=list()
for(gene in gene.selected){
  p=mg_violin(data.frame(tcga.t.cli_2$Stage,stad.tcga.t.PPARs.exp[,gene])
              , melt = T
              , xlab = 'Stage'
              , ylab = paste('The expression of gene', gene, "\nlog2(TPM+1)")
              , legend.pos = 'top'
              ,test_method='kruskal.test')
  p=p+theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    # plot.subtitle=element_text(size=8),
    plot.caption=element_blank(),
    #plot.title = element_text(color = "red", size = 18, face = "bold", vjust = 0, hjust = 0.5),
    axis.ticks = element_blank(),
    axis.line = element_line(colour = "grey50"),
    panel.grid = element_line(color = "#b4aea9"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    # panel.grid.major.y = element_line(linetype = "dashed"),
    panel.grid.major.y=element_blank(),
    panel.background = element_rect(fill = "white", color = "white"),
    plot.background = element_rect(fill = "white", color = "white")
  )
  p.all=c(p.all,list(p))
}
length(p.all)
mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'B')

########
stad.tcga.t.PPARs.exp_melt=reshape2::melt(stad.tcga.t.PPARs.exp[,c('SampleCode',gene.selected)], id.vars=c("SampleCode"))
head(stad.tcga.t.PPARs.exp_melt)

stad.tcga.t.PPARs.exp_melt=right_join(stad.tcga.t.PPARs.exp_melt,tcga.t.cli_2)
range(stad.tcga.t.PPARs.exp_melt$Exprs)
head(stad.tcga.t.PPARs.exp_melt)
colnames(stad.tcga.t.PPARs.exp_melt)[2:3]=c('PPARs','Exprs')
head(stad.tcga.t.PPARs.exp_melt)

library(ggstatsplot)
colnames(tcga.t.cli_2)
colnames(tcga.t.cli_2)[c(7,9,11,10)]

p.all=list()
for(i in c(7,9,11,10)){
  x=colnames(tcga.t.cli_2)[i]
  comparisons=names(table(tcga.t.cli_2[i]))
  for(gene in gene.selected){
    p=stad.tcga.t.PPARs.exp_melt %>%
      dplyr::filter(PPARs %in% c(gene)) %>%
      ggbetweenstats(
        x = !!sym(x),
        y = Exprs,
        plot.type = "boxviolin",
        type = "parametric",
        pairwise.comparisons = TRUE,
        p.adjust.method = "fdr",
        ylab="log2(TPM+1)",
        title=gene,
        # ggplot.component = ggplot2::scale_y_continuous(
        #   breaks = seq(0, 9, 1),
        #   limits = (c(0, 9))
        # ),
        ggtheme = theme(
          axis.text = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12),
          # plot.subtitle=element_text(size=8),
          plot.caption=element_blank(),
          #plot.title = element_text(color = "red", size = 18, face = "bold", vjust = 0, hjust = 0.5),
          # axis.ticks = element_blank(),
          axis.line = element_line(colour = "grey50"),
          panel.grid = element_line(color = "#b4aea9"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          # panel.grid.major.y = element_line(linetype = "dashed"),
          panel.grid.major.y=element_blank(),
          panel.background = element_rect(fill = "white", color = "white"),
          plot.background = element_rect(fill = "white", color = "white")
        )
      )+ggsignif::geom_signif(comparisons = list(c(comparisons)),map_signif_level  = F)
    p.all=c(p.all,list(p))
  }
}
length(p.all)
fig2a=mg_merge_plot(p.all,nrow = 4,ncol = 3,labels = LETTERS[1:4])
savePDF('PDFs/Fig2A.pdf',fig2a,height = 16,width = 12)
savePDF('02_PPARs_Clinical/Fig2A.pdf',fig2a,height = 16,width = 12)

####### PPARs基因表达与临床病理特征之间的关系，三线表展示
stad.tcga.t.PPARs.exp_use=stad.tcga.t.PPARs.exp[,c('SampleCode',gene.selected)]
head(stad.tcga.t.PPARs.exp_use)
stad.tcga.t.PPARs.exp_use$PPARA=ifelse(stad.tcga.t.PPARs.exp_use$PPARA>median(stad.tcga.t.PPARs.exp_use$PPARA),'High','Low')
stad.tcga.t.PPARs.exp_use$PPARD=ifelse(stad.tcga.t.PPARs.exp_use$PPARD>median(stad.tcga.t.PPARs.exp_use$PPARD),'High','Low')
stad.tcga.t.PPARs.exp_use$PPARG=ifelse(stad.tcga.t.PPARs.exp_use$PPARG>median(stad.tcga.t.PPARs.exp_use$PPARG),'High','Low')

stad.tcga.t.PPARs.exp_use=reshape2::melt(stad.tcga.t.PPARs.exp_use[,c('SampleCode',gene.selected)], id.vars=c("SampleCode"))
head(stad.tcga.t.PPARs.exp_use)
colnames(stad.tcga.t.PPARs.exp_use)[2:3]=c('PPARs','Type')

stad.tcga.t.PPARs.exp_use=right_join(stad.tcga.t.PPARs.exp_use,tcga.t.cli_2)
range(stad.tcga.t.PPARs.exp_use$Exprs)
save(stad.tcga.t.PPARs.exp_use,file = 'origin_datas/TCGA/stad.tcga.t.PPARs.exp_use.RData')
head(stad.tcga.t.PPARs.exp_use)

write.table(stad.tcga.t.PPARs.exp_use
            ,file = '02_PPARs_Clinical/stad.tcga.t.PPARs.clinical.txt'
            ,sep = '\t',quote = F,row.names = F,col.names = T)

# setwd('Z:/users/lishuang/Work1/Work/20221209_STAD_PPARs')
# load('origin_datas/TCGA/stad.tcga.t.PPARs.exp_use.RData')
# head(stad.tcga.t.PPARs.exp_use)
# library(gtsummary)
# stad.tcga.t.PPARs.exp_use=stad.tcga.t.PPARs.exp_use[,-10]
# colnames(stad.tcga.t.PPARs.exp_use)[12]='Age'
# 
# t1=stad.tcga.t.PPARs.exp_use %>%
#   filter(PPARs=='PPARA') %>%
#   select(Type,T.Stage,N.Stage,M.Stage,Stage,Grade,Age,Gender)%>%
#   tbl_summary(by = Type) %>%
#   add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
#   add_overall() %>%
#   add_n() %>%
#   modify_header(label ~ "**Variable**") %>%
#   modify_spanning_header(c("stat_1", "stat_2") ~ "**PPARA**") %>%
#   modify_footnote(
#     all_stat_cols() ~ "Median (IQR) or Frequency (%)"
#   ) %>%
#   modify_caption("**Table 1. Relationship Between PPARs family member expression and clinicopathological features in the TCGA-STAD cohort.**") %>%
#   bold_labels()
# 
# t2=stad.tcga.t.PPARs.exp_use %>%
#   filter(PPARs=='PPARD') %>%
#   select(Type,T.Stage,N.Stage,M.Stage,Stage,Grade,Age,Gender)%>%
#   tbl_summary(by = Type) %>%
#   add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
#   # add_overall() %>%
#   # add_n() %>%
#   modify_header(label ~ "**Variable**") %>%
#   modify_spanning_header(c("stat_1", "stat_2") ~ "**PPARD**") %>%
#   modify_footnote(
#     all_stat_cols() ~ "Median (IQR) or Frequency (%)"
#   ) %>%
#   modify_caption("**Table 1. Relationship Between PPARs family member expression and clinicopathological features in the TCGA-STAD cohort.**") %>%
#   bold_labels()
# 
# t3=stad.tcga.t.PPARs.exp_use %>%
#   filter(PPARs=='PPARG') %>%
#   select(Type,T.Stage,N.Stage,M.Stage,Stage,Grade,Age,Gender)%>%
#   tbl_summary(by = Type) %>%
#   add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
#   # add_overall() %>%
#   # add_n() %>%
#   modify_header(label ~ "**Variable**") %>%
#   modify_spanning_header(c("stat_1", "stat_2") ~ "**PPARG**") %>%
#   modify_footnote(
#     all_stat_cols() ~ "Median (IQR) or Frequency (%)"
#   ) %>%
#   modify_caption("**Table 1. Relationship Between PPARs family member expression and clinicopathological features in the TCGA-STAD cohort.**") %>%
#   bold_labels()
# 
# 
# tbl_merge <-
#   tbl_merge(tbls = list(t1, t2, t3),
#             tab_spanner = c("**PPARA**", "**PPARD**","**PPARG**"))
# 
# tbl_merge
# 
# tbl_merge %>%    # build gtsummary table
#   as_gt() %>%             # convert to gt table
#   gt::gtsave(             # save table as image
#     filename = "Z:/users/lishuang/Work1/Work/20221209_STAD_PPARs/Table1.html",
#   )
# 
# 
# library(gdtools)
# library(flextable)
# library(officer)
# library(webshot2)
# library(writexl)
# 
# f='Z:/users/lishuang/Work1/Work/20221209_STAD_PPARs/files/Table1.docx'
# sect_properties <- prop_section(
#   page_size = page_size(
#     orient = "landscape",
#     width = 20, height = 16
#   )
#   # type = "continuous"
#   # page_margins = page_mar()
# )
# 
# tbl_merge %>%
#   as_flex_table() %>%
#   flextable::save_as_docx(path = f,pr_section = sect_properties)
# 
# ###############
# f='Z:/users/lishuang/Work1/Work/20221209_STAD_PPARs/files/Table1.xlsx'
# tbl_merge %>%
#   gtsummary::as_tibble() %>%
#   writexl::write_xlsx(., path = f)
# 
# tbl_merge %>%
#   as_gt() %>%             # convert to gt table
#   gt::gtsave(             # save table as image
#     filename = "Table1.png",
#     path='Z:/users/lishuang/Work1/Work/20221209_STAD_PPARs/PDFs/'
#   )


###################### PPARs基因表达与预后的关系
dir.create('03_PPARs_Prognosis')
library(survival)
dim(stad.tcga.t.PPARs.exp)
head(stad.tcga.t.PPARs.exp_use)

p.all=list()
for(gene in gene.selected){
  df_m=stad.tcga.t.PPARs.exp[,c('SampleCode',gene)]
  all(tcga.t.cli_2$SampleCode==stad.tcga.t.PPARs.exp$SampleCode)
  cutoff=survminer::surv_cutpoint(data.frame(time=tcga.t.cli_2$OS.time/365,
                                                     event=tcga.t.cli_2$OS,
                                                     risk=stad.tcga.t.PPARs.exp[,gene]), time = "time", event = "event",
                                          variables = c("risk"))
  cutoff=cutoff$cutpoint$cutpoint
  df_m$Type=ifelse(df_m[,gene]>cutoff,'High','Low')
  
  # df_forCox=cbind(df_m,tcga.t.cli_2)
  # fmla <- as.formula("Surv(OS.time, OS) ~ PPARG")
  # cox <- coxph(fmla, data =as.data.frame(df_forCox))
  # survminer::ggforest(cox,data=df_forCox,noDigits = 3)
  
  p=ggplotKMCox(data.frame(time = tcga.t.cli_2$OS.time/365
                           , event = tcga.t.cli_2$OS
                           , groups=df_m$Type)
                , title=gene
                , pal = risk.group.color
                , labs = c('High','Low')
                , add_text = '')
  p.all=c(p.all,list(p))
}

fig3a=mg_merge_plot(p.all,nrow = 1,ncol = 3)
fig3a

p.all=list()
for(gene in gene.selected){
  df_m=stad.tcga.t.PPARs.exp[,c('SampleCode',gene)]
  all(tcga.t.cli_2$SampleCode==stad.tcga.t.PPARs.exp$SampleCode)
  
  cutoff=survminer::surv_cutpoint(data.frame(time=stad.tcga.t.cli$A8_New_Event_Time/365,
                                             event=stad.tcga.t.cli$A8_New_Event,
                                             risk=stad.tcga.t.PPARs.exp[,gene]), time = "time", event = "event",
                                  variables = c("risk"))
  cutoff=cutoff$cutpoint$cutpoint
  df_m$Type=ifelse(df_m[,gene]>cutoff,'High','Low')
  
  p=ggplotKMCox(data.frame(time = stad.tcga.t.cli$A8_New_Event_Time/365
                           , event = stad.tcga.t.cli$A8_New_Event
                           , groups=df_m$Type)
                , title=gene
                , pal = risk.group.color
                , labs = c('High','Low')
                , add_text = '')
  p.all=c(p.all,list(p))
}

fig3b=mg_merge_plot(p.all,nrow = 1,ncol = 3)
fig3b

fig3ab=mg_merge_plot(fig3a,fig3b,nrow = 2,ncol = 1,labels = LETTERS[1:2])
fig3ab

savePDF('PDFs/Fig3AB.pdf',fig3ab,height = 8,width = 12)
savePDF('03_PPARs_Prognosis/Fig3AB.pdf',fig3ab,height = 8,width = 12)

################# 基因共表达分析
dir.create('04_PPARs_Co-expression')
library(parallel)
detectCores()
cal_correlation_parallel <- function(x,y,threads=8,method='spearman'){
  library(foreach)
  library(doParallel)
  library(abind)
  if(method=='spearman'){
    otu_table1 <- apply(x,2,rank)
    otu_table2 <- apply(y,2,rank)
  }else{
    otu_table1 <- x
    otu_table2 <- y
  }
  
  r <- function(rx,ry){
    n <- length(rx)
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
    return(c(r,p))
  }
  arraybind <- function(...){
    abind(...,along = 3,force.array=TRUE)
  }
  nc1 <- ncol(otu_table1)
  nc2 <- ncol(otu_table2)
  registerDoParallel(cores = threads)
  corr <- foreach (i = 1:nc1,.combine = "arraybind") %dopar%{
    corr1 <- matrix(rep(0,2*nc1),nrow = 2,ncol=nc2)
    for(j in 1:nc2) {
      # if(j > i) 
      corr1[,j] <- r(otu_table1[,i],otu_table2[,j])
    }
    corr <- corr1
  }
  if(is.matrix(corr)){
    rr=corr[1,]
    pp <- corr[2,]
    rr=matrix(rr)
    pp=matrix(pp)
  }else if(is.array(corr)){
    rr <- corr[1,,]
    # rr <- rr+t(rr)
    # diag(rr) <- 1
    pp <- corr[2,,]
  }
  lp <- lower.tri(pp)
  pa <- pp[lp]
  # pa <- p.adjust(pa, "fdr")
  pp[lower.tri(pp, diag = FALSE)] <- pa
  # pp <- pp+t(pp)
  rownames(pp) <- colnames(otu_table2)
  colnames(pp) <- colnames(otu_table1)
  rownames(rr) <- colnames(otu_table2)
  colnames(rr) <- colnames(otu_table1)
  # return(list(r = rr,p = pp))
  
  df_cor=rr
  df_pval=pp
  ### Pivot data from wide to long
  library(tidyverse)
  g = pivot_longer(data=rownames_to_column(as.data.frame(df_cor),var = "to"),
                   cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                   names_to = "from",
                   values_to = "cor")
  gp = pivot_longer(data=rownames_to_column(as.data.frame(df_pval)),
                    cols = 2:(ncol(df_pval)+1),
                    names_to = "gene",
                    values_to = "p")
  g$p.adj = p.adjust(gp$p, "fdr")
  g=g[,c(2,1,3:4)]
  g=g[order(g$from),]
  return(g)
}
stad.tcga.t.pcg.exp=stad.tcga.t.exp[rownames(stad.tcga.t.exp) %in% ann.pcg$gene_name,]
dim(stad.tcga.t.pcg.exp)
dim(stad.tcga.t.PPARs.exp)

rownames(stad.tcga.t.PPARs.exp)=stad.tcga.t.PPARs.exp$SampleCode
all(rownames(stad.tcga.t.PPARs.exp)==colnames(stad.tcga.t.pcg.exp))

tcga.PPARs.PCG.cor=cal_correlation_parallel(x=stad.tcga.t.PPARs.exp[,gene.selected],y=t(stad.tcga.t.pcg.exp),method = 'pearson')
tcga.PPARs.PCG.cor=data.frame(tcga.PPARs.PCG.cor)
head(tcga.PPARs.PCG.cor)
######## PPARs 与 PPARs相关性
tcga.PPARs.PCG.cor=tcga.PPARs.PCG.cor[-which(tcga.PPARs.PCG.cor$from==tcga.PPARs.PCG.cor$to),]

cor.cutoff=0.4
p.cutoff=0.05
tcga.PPARs.PCG.cor$type=case_when(tcga.PPARs.PCG.cor$cor>cor.cutoff & tcga.PPARs.PCG.cor$p.adj<p.cutoff ~ "positive",
                                  tcga.PPARs.PCG.cor$cor< (-cor.cutoff) & tcga.PPARs.PCG.cor$p.adj<p.cutoff ~ "negative",
                                  TRUE~"not" )
tcga.PPARs.PCG.cor=tcga.PPARs.PCG.cor[order(-tcga.PPARs.PCG.cor$cor),]
head(tcga.PPARs.PCG.cor)
write.table(tcga.PPARs.PCG.cor
            ,file = '04_PPARs_Co-expression/tcga.PPARs.PCG.cor.txt'
            ,row.names = F,col.names = T,sep = '\t',quote = F)

write.table(tcga.PPARs.PCG.cor
            ,file = 'files/文件/tcga.PPARs.PCG.cor.txt'
            ,row.names = F,col.names = T,sep = '\t',quote = F)

######## 相关性显著的
tcga.PPARs.PCG.cor.sig=tcga.PPARs.PCG.cor[which(abs(tcga.PPARs.PCG.cor$cor)>cor.cutoff & tcga.PPARs.PCG.cor$p.adj<p.cutoff),]
dim(tcga.PPARs.PCG.cor.sig)
write.table(tcga.PPARs.PCG.cor.sig
            ,file = '04_PPARs_Co-expression/tcga.PPARs.PCG.cor.sig.txt'
            ,row.names = F,col.names = T,sep = '\t',quote = F)
write.table(tcga.PPARs.PCG.cor.sig
            ,file = 'files/文件/tcga.PPARs.PCG.cor.sig.txt'
            ,row.names = F,col.names = T,sep = '\t',quote = F)

table(tcga.PPARs.PCG.cor.sig$from)
table(tcga.PPARs.PCG.cor.sig$type)
table(tcga.PPARs.PCG.cor.sig$from,tcga.PPARs.PCG.cor.sig$type)

p.all=list()
for(gene in gene.selected){
  df=tcga.PPARs.PCG.cor[which(tcga.PPARs.PCG.cor$from==gene),]
  ind1=head(which(df$cor>cor.cutoff),n=10)
  ind2=tail(which(df$cor<(-cor.cutoff)),n=10)
  for_label <-df[c(ind1,ind2),]
  df[ind1,]
  
  dim(for_label)
  
  options(ggrepel.max.overlaps = Inf)
  p <- ggplot(data = df,
               aes(x = cor,
                   y = -log10(p.adj)))
  p=p+geom_point(alpha=0.4, size=2, aes(color=type))
  p=p+scale_color_manual(values=c(mg_colors[1],'black',mg_colors[2]),limits = c("positive",'not', "negative"),name='State')
  p=p+geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)
  p=p+ylab('-log10(FDR)')+xlab("Pearson's correlation")+labs(title = gene)
  p=p+geom_point(size =1, shape = 2, data = for_label)+
    ggrepel::geom_label_repel(aes(label = to,label.size =0.5),
                              data = for_label,
                              color="black")
  p=p+theme_bw()
  p=p+theme(
    axis.text.y=element_text(family="Times",face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
    axis.title.y=element_text(family="Times",face="plain"), #设置y轴标题的字体属性
    legend.text=element_text(face="plain", family="Times", colour="black"), #设置图例的子标题的字体属性
    legend.title=element_text(face="plain", family="Times", colour="black"), #设置图例的总标题的字体属性
    legend.justification=c(1,1), legend.position=c(1,1),
    legend.background = element_rect(fill = NA, colour = NA)
  )
  p
  p.all=c(p.all,list(p))
}
fig4a=mg_merge_plot(p.all,nrow = 1,ncol = 3,common.legend = T)
fig4a

dotplot_batch_use=function(enrichmentORA,dbs=c('pathway_KEGG','pathway_Wikipathway','pathway_Reactome'
                                               ,'geneontology_Biological_Process','geneontology_Cellular_Component'
                                               ,'geneontology_Molecular_Function','disease_Disgenet','disease_OMIM'
                                               ,'drug_DrugBank','phenotype_Human_Phenotype_Ontology'),top=20,FDR=T
                           ,high_col='red',low_col='blue'){
  #library('ggpubr')
  library(ggplot2)
  #library(ggsci)
  #db='pathway_KEGG'
  plts=list()
  p.ct=0
  all.enrich=enrichmentORA
  dbs=intersect(dbs,enrichmentORA$DB)
  #windowsFonts(mgFont = windowsFont("Times New Roman"))
  for(db in dbs){
    tl=paste0('All ',db)
    pathway_data=all.enrich[all.enrich$DB==db,]
    if(nrow(pathway_data)>0){
      if(nrow(pathway_data)>top){
        pathway_data=pathway_data[order(pathway_data$FDR)[1:top],]
        tl=paste0('Top',top,' ',db)
      }
      desc=pathway_data$description
      ndesc=c()
      for(de in desc){
        if(nchar(de)>50&length(grep(' ',de))>0){
          de1=unlist(strsplit(de,' '))
          d2=paste0(de1[(ceiling(length(de1)/2)+1):length(de1)],collapse = ' ')
          if(nchar(d2)>50){
            d2=paste0(substr(d2,0,47),'...')
          }
          de2=paste0(paste0(de1[1:ceiling(length(de1)/2)],collapse = ' '),'\n'
                     ,d2)
          ndesc=c(ndesc,de2)
        }else{
          ndesc=c(ndesc,de)
        }
      }
      pathway_data$description=ndesc      
      # pathway_data$description <- factor(pathway_data$description
      #                                    ,levels=pathway_data$description[order(pathway_data$size,decreasing = T)], ordered=TRUE)
      pathway_data$description <- factor(pathway_data$description
                                         ,levels=pathway_data$description[order(pathway_data$enrichmentRatio,decreasing = F)], ordered=TRUE)
      
      pathway_data$pValue[pathway_data$pValue<1e-16]=1e-16
      pathway_data$FDR[pathway_data$FDR<1e-16]=1e-16
      
      bubble=ggplot(data = pathway_data, aes(x = enrichmentRatio, y = description)) 
      
      bubble=bubble+xlab('Enrichment Ratio')
      if(FDR){
        
        bubble=bubble+geom_point(aes(size = size,color = -log10(FDR))) +scale_color_gradient(low = low_col, high = high_col)
      }else{
        bubble=bubble+geom_point(aes(size = size,color = -log10(pValue))) +scale_color_gradient(low = low_col, high = high_col)
      }
      
      bubble=bubble+ggtitle(tl)
      bubble=bubble+ggsci::scale_fill_npg()+ggplot2::theme_bw()+theme(axis.title.y=element_blank()
                                                                      ,axis.text.y=element_text(family="Times",face="plain")
                                                                      ,axis.text.x=element_text(family="Times",face="plain")
                                                                      ,plot.title = element_text(hjust = 0.5,family="Times",face="plain")
                                                                      ,axis.title.x=element_text(family="Times",face="plain")
                                                                      ,legend.title = element_text(family="Times",face="plain")
                                                                      ,legend.text = element_text(family="Times",face="plain"))
      plts=c(plts,list(bubble))
      p.ct=p.ct+1
    }
  }
  if(p.ct>1){
    nc=2
    nr=ceiling(p.ct/2)
  }else{
    nc=p.ct
    nr=1
  }
  g1=ggpubr::ggarrange(plotlist=plts, ncol = nc, nrow = nr,labels = toupper(letters)[1:p.ct],align = "hv")
  #detach('package:ggpubr')
  #detach('package:ggplot2')
  #detach('package:ggsci')
  
  return(g1)
}

length(unique(tcga.PPARs.PCG.cor.sig$to))

tcga.PPARs.cor.sig.enrich.res=mg_clusterProfiler(genes =unique(tcga.PPARs.PCG.cor.sig$to))
tcga.PPARs.cor.sig.enrich.res$Enrich_tab[1:4,1:5]

write.table(tcga.PPARs.cor.sig.enrich.res$Enrich_tab
            ,file = '04_PPARs_Co-expression/tcga.PPARs.cor.sig.enrich.res.txt'
            ,row.names = T,col.names = T,sep = '\t',quote = F)

write.table(tcga.PPARs.cor.sig.enrich.res$Enrich_tab
            ,file = 'files/文件/tcga.PPARs.cor.sig.enrich.res.txt'
            ,row.names = T,col.names = T,sep = '\t',quote = F)

fig4b=dotplot_batch_use(tcga.PPARs.cor.sig.enrich.res$Enrich_tab
                       , high_col = 'red'
                       , low_col = 'purple'
                       , top = 20)
fig4b
############ 与PPARs均存在显著相关的基因
head(tcga.PPARs.PCG.cor.sig)
dim(tcga.PPARs.PCG.cor.sig)
length(unique(tcga.PPARs.PCG.cor.sig$to))

plot_PointHeatmap=function(df_m,univar,cutoff=NULL,df_cor=NULL,pal=NULL,hetTitle='z-score',hetColor=c('#3B4992FF', 'white', '#EE0000FF')){
  dat=df_m
  srt.inds=order(dat[,univar])
  dat=dat[srt.inds,]
  dim(dat)
  ############
  if(is.null(cutoff)){
    cutoff=median(dat[,univar])
  }else{
    cutoff=cutoff
  }
  library(ggplot2)
  dt1=data.frame(V1=1:nrow(dat),V2=dat[,univar],Group=ifelse(dat[,univar]>cutoff,'High','Low')) 
  p1=ggplot(dt1, aes(x = V1, y = V2, colour = Group,fill=Group)) +
    geom_point(stat = 'identity', position = 'identity')+
    scale_colour_manual(values = risk.group.color,aesthetics = c("colour", "fill"))+
    theme_classic()
  p1=p1+coord_cartesian(expand = FALSE)
  p1=p1+geom_hline(aes(yintercept=cutoff),colour='black',linetype="dashed")
  p1=p1+ylab(univar)
  p1=p1+theme(axis.text.y=element_text(family="Times",face="plain"),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.line.x = element_blank(),
              axis.ticks.x = element_blank(),
              plot.margin=unit(c(0.1, 0.1, 0, 0.1), "inches"),
              legend.position = c(1, 0),
              legend.justification = c(1, 0),
              legend.background = element_rect(fill = NA, colour = NA),
              legend.title = element_text(family="Times",face="plain",colour = 'black'),
              legend.text = element_text(family="Times",face="plain",colour = 'black'))
  p1
  #######################
  ### scale
  data=as.data.frame(scale(dat[,-which(colnames(dat)==univar)]))
  ### row cluster
  hc.r = hclust(dist(t(data)))
  data=data[,hc.r$order]
  ################
  data$ID <- 1:nrow(dat)
  #colnames(data)
  data_m <- reshape2::melt(data, id.vars=c("ID"))
  colnames(data_m)=c('ID','V1','V2')
  data_m$V2[which(data_m$V2>mean(data_m$V2)+3*sd(data_m$V2))]=mean(data_m$V2)+3*sd(data_m$V2)
  data_m$V2[which(data_m$V2<mean(data_m$V2)-3*sd(data_m$V2))]=mean(data_m$V2)-3*sd(data_m$V2)
  
  data_m$V1=mg_str_outline(data_m$V1,isCut = T,n=50)
  
  if(!is.null(df_cor)){
    df_cor=df_cor[which(df_cor$from==univar),]
    df_cor=df_cor[match(colnames(data),df_cor$to),]
    df_cor=na.omit(df_cor)
    rownames(df_cor)=NULL
  }
  ######
  p2 <- ggplot(data_m, aes(x=ID,y=V1))
  p2 <- p2 + geom_tile(aes(fill=V2))
  p2=p2+coord_cartesian(expand = FALSE)
  p2=p2+scale_fill_gradient2(low = hetColor[1],mid=hetColor[2], high = hetColor[3])
  p2=p2+theme(axis.text.y=element_text(family="Times",face="plain",size = 12,hjust=1,vjust=0.5),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              plot.margin=unit(c(0.1, 0.1, 0, 0.1), "inches"),
              legend.position = 'bottom',
              legend.margin = margin(t=-25),
              legend.key.height= unit(0.1, 'inches'),
              legend.key.width = unit(0.2, 'inches'),
              legend.title.align=0.5,
              panel.background = element_rect(fill = "white", color = "white"),
              plot.background = element_rect(fill = "white", color = "white")
  )+labs(fill=hetTitle) + xlab('')+ylab('')
  p2
  
  g1=ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2,heights = c(1,1),align = "v",vjust=)
  g1
  return(g1)
}

df=tcga.PPARs.PCG.cor.sig[tcga.PPARs.PCG.cor.sig$to %in% names(which(table(tcga.PPARs.PCG.cor.sig$to)>2)),]
dim(df)
p.all=list()
for(gene in gene.selected){
  df_m=cbind(t(stad.tcga.t.pcg.exp[unique(df$to),]),subset(stad.tcga.t.PPARs.exp,select=gene))
  dim(df_m)
  head(df_m)
  p=plot_PointHeatmap(df_m = df_m,univar = gene)
  p
  p.all=c(p.all,list(p))
}

fig4c=mg_merge_plot(p.all,nrow = 1,ncol = 3)
fig4c

fig4=mg_merge_plot(fig4a,fig4c,fig4b,nrow = 3,ncol = 1,heights = c(1,1,2),labels = c('A','B',''))
savePDF('PDFs/Fig4.pdf',fig4,height = 16,width = 12)
savePDF('04_PPARs_Co-expression/Fig4.pdf',fig4,height = 16,width = 12)

######## 通路分析
dir.create('05_PPARs_Pathway')
# library(GSVA)
# library(GSEABase)
# stad.tcga.exp=stad.tcga.exp.tpm
# gmtFile='origin_datas/h.all.v7.5.1.symbols.gmt'
# c2KEGG <- getGmt(gmtFile,
#                  collectionType=BroadCollection(category="c2"),
#                  geneIdType=SymbolIdentifier())
# tcga.h.all.ssgsea <- gsva(as.matrix(stad.tcga.exp),
#                           c2KEGG,
#                           method = 'ssgsea',
#                           min.sz = 10,
#                           max.sz = 500,
#                           verbose = TRUE)
# save(tcga.h.all.ssgsea,file='origin_datas/immune/tcga.h.all.ssgsea.RData')
# gmtFile='origin_datas/c2.cp.kegg.v7.5.1.symbols.gmt'
# c2KEGG <- getGmt(gmtFile,
#                  collectionType=BroadCollection(category="c2"),
#                  geneIdType=SymbolIdentifier())
# tcga.kegg.ssgsea <- gsva(as.matrix(stad.tcga.exp),
#                          c2KEGG,
#                          method = 'ssgsea',
#                          min.sz = 10,
#                          max.sz = 500,
#                          verbose = TRUE)
# save(tcga.kegg.ssgsea,file='origin_datas/immune/tcga.kegg.ssgsea.RData')
load('origin_datas/immune/tcga.kegg.ssgsea.RData')
load('origin_datas/immune/tcga.h.all.ssgsea.RData')
dim(tcga.h.all.ssgsea)
dim(tcga.kegg.ssgsea)

all(colnames(tcga.h.all.ssgsea)==colnames(stad.tcga.exp.tpm))
all(colnames(stad.tcga.exp.tpm)==rownames(stad.tcga.PPARs.exp))

rownames(tcga.h.all.ssgsea)=gsub("HALLMARK_","",rownames(tcga.h.all.ssgsea))
rownames(tcga.kegg.ssgsea)=gsub("KEGG_","",rownames(tcga.kegg.ssgsea))

tcga.t.h.all.ssgsea=tcga.h.all.ssgsea[,rownames(tcga.t.cli_use)]
tcga.t.kegg.ssgsea=tcga.kegg.ssgsea[,rownames(tcga.t.cli_use)]

write.table(tcga.h.all.ssgsea,file = '05_PPARs_Pathway/tcga.h.all.ssgsea.txt',sep = '\t',quote = F)
write.table(tcga.h.all.ssgsea,file = 'files/文件/tcga.h.all.ssgsea.txt',sep = '\t',quote = F)

########
table(substr(colnames(tcga.h.all.ssgsea),14,15))
all(colnames(tcga.h.all.ssgsea)==stad.PPARs.exp$SampleCode)
#########
tcga.h.all.ssgsea_use=data.frame(t(tcga.h.all.ssgsea))
tcga.h.all.ssgsea_use$SampleName=rownames(tcga.h.all.ssgsea_use)
tcga.h.all.ssgsea_use[1:4,1:5]

tcga.h.all.ssgsea_m=reshape2::melt(tcga.h.all.ssgsea_use,id.vars='SampleName')
head(tcga.h.all.ssgsea_m)
colnames(tcga.h.all.ssgsea_m)=c('SampleName','pathway','ssGSEA')
head(tcga.h.all.ssgsea_m)

all(tcga.h.all.ssgsea_use$SampleName==stad.tcga.PPARs.exp$SampleCode)

tcga.h.all.ssgsea_m=merge(tcga.h.all.ssgsea_m,stad.tcga.PPARs.exp[,c('SampleCode','Site')],by.x='SampleName',by.y='SampleCode')
head(tcga.h.all.ssgsea_m)

##########
tcga.h.all.cmp=compare_means(ssGSEA ~ Site,
                             group.by='pathway',
                             data=tcga.h.all.ssgsea_m)
tcga.h.all.cmp=data.frame(tcga.h.all.cmp,check.names = F,stringsAsFactors = F)
tcga.h.all.cmp=crbind2DataFrame(tcga.h.all.cmp)
head(tcga.h.all.cmp)

write.table(tcga.h.all.cmp
            ,file = '05_PPARs_Pathway/tcga.h.all.cmp.txt'
            ,row.names = T,col.names = T,sep = '\t',quote = F)

table(tcga.h.all.cmp$p.adj<0.05)
table(tcga.h.all.cmp$p.adj<0.01)

all(tcga.h.all.cmp$pathway==rownames(tcga.h.all.ssgsea))
tcga.h.all.cmp=tcga.h.all.cmp[match(rownames(tcga.h.all.ssgsea),tcga.h.all.cmp$pathway),]
all(tcga.h.all.cmp$pathway==rownames(tcga.h.all.ssgsea))

###############
df=stad.tcga.PPARs.exp[,gene.selected]

ind1=which(stad.tcga.PPARs.exp$Site=='Tumor')
ind2=which(stad.tcga.PPARs.exp$Site=='Normal')

# hc.r = hclust(dist(scale(df[ind1,])))
# t.smp.order=rownames(df[ind1,])[hc.r$order]
t.smp.order=names(sort(rowMeans(df[ind1,])))

# hc.r = hclust(dist(scale(df[ind2,])))
# n.smp.order=rownames(df[ind2,])[hc.r$order]
n.smp.order=names(sort(rowMeans(df[ind2,])))

smp.order=c(t.smp.order,n.smp.order)

min_max=apply(df,2,function(x){range(x)})
min_max=floor(min_max)
min_max

library(circlize)
library(ComplexHeatmap)

column_ha <- HeatmapAnnotation(df = data.frame(stad.tcga.PPARs.exp[smp.order,c('Site',gene.selected)])
                               , col = list(PPARA=colorRamp2(min_max[,1], c('white', PPARs.color[1]))
                                            ,PPARD=colorRamp2(min_max[,2], c('white', PPARs.color[2]))
                                            ,PPARG=colorRamp2(min_max[,3], c('white', PPARs.color[3]))
                                            ,Site=c('Tumor'=sampleType.color[2],'Normal'=sampleType.color[1])))

row_ha <- rowAnnotation(p.signif = anno_text(tcga.h.all.cmp$p.signif, 
                                             just = "left",
                                             # location = unit(0, "npc"),
                                             gpar(fontsize = 8),
                                             show_name = TRUE), 
                        annotation_name_rot = 0)

fig5a=Heatmap(as.matrix(t(scale(t(tcga.h.all.ssgsea[,smp.order]))))
              , name = "ssGSEA"
              , cluster_rows = T
              , cluster_row_slices = T
              , row_names_gp = gpar(fontsize = 8)
              , show_row_names = TRUE
              , row_names_side ='left'
              , show_row_dend = F
              , cluster_columns = F
              , cluster_column_slices = F
              , show_column_dend = F
              , show_column_names = F
              , col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))
              , top_annotation = column_ha
              , right_annotation =row_ha
              , border = TRUE)
fig5a

######## 相关性分析
df=cal_correlation_parallel(x=stad.tcga.t.PPARs.exp[,gene.selected],y=t(tcga.t.h.all.ssgsea)
                            ,method = 'pearson')
df <- df %>%
  mutate(col = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*","ns"),
                        right = FALSE, include.lowest = TRUE))
head(df)

#### 去掉均不显著的通路
df=df[df$to %in% (names(which(table(df$to,df$p.signif)[,4]<3))),]

########
df1=reshape2::dcast(df[,1:3],to ~ from)
rownames(df1)=df1$to
df1=df1[,-1]

hc.r=hclust(dist(df1))
hc.r$order
plot(hc.r)

pw.order=rownames(df1)[hc.r$order]
df$to=factor(df$to,pw.order)

fig5b=ggplot(df,aes(y = to,x = from))+
  geom_point(aes(size=abs(cor),color=col))+
  geom_text(aes(y = to,x = from,label=p.signif))+
  theme_bw()+scale_color_manual(values = ggsci::pal_lancet(alpha = 0.5)(9))+
  ylab('')+xlab('')+labs(color="pearson's cor",size="pearson's cor")
fig5b

######## 相关性分析-正常样本
all(colnames(tcga.h.all.ssgsea)==stad.tcga.PPARs.exp$SampleCode)
inds=which(stad.tcga.PPARs.exp$Site=='Normal')

df=cal_correlation_parallel(x=stad.tcga.PPARs.exp[inds,gene.selected],y=t(tcga.h.all.ssgsea[,inds])
                            ,method = 'pearson')
df <- df %>%
  mutate(col = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*","ns"),
                        right = FALSE, include.lowest = TRUE))
head(df)

#### 去掉均不显著的通路
df=df[df$to %in% (names(which(table(df$to,df$p.signif)[,4]<3))),]
#########
df1=reshape2::dcast(df[,1:3],to ~ from)
rownames(df1)=df1$to
df1=df1[,-1]

hc.r=hclust(dist(df1))
hc.r$order

pw.order=rownames(df1)[hc.r$order]

df$to=factor(df$to,pw.order)

fig5c=ggplot(df,aes(y = to,x = from))+
  geom_point(aes(size=abs(cor),color=col))+
  geom_text(aes(y = to,x = from,label=p.signif))+
  theme_bw()+scale_color_manual(values = ggsci::pal_lancet(alpha = 0.5)(9))+
  ylab('')+xlab('')+labs(color="pearson's cor",size="pearson's cor")
fig5c

fig5bc=mg_merge_plot(fig5b,fig5c,nrow = 1,ncol = 2,labels = LETTERS[2:3])
savePDF('PDFs/Fig5BC.pdf',fig5bc,height = 6,width = 12)
savePDF('05_PPARs_Pathway/Fig5BC.pdf',fig5bc,height = 6,width = 12)

pdf('PDFs/Fig5A.pdf',height = 6,width = 12)
fig5a
dev.off()

pdf('05_PPARs_Pathway/Fig5A.pdf',height = 6,width = 12)
fig5a
dev.off()

###################
################# 与免疫的相关性分析
# immunomodulators 基因 
dir.create('06_PPARs_Immune')

                                       header = F, check.names = F,stringsAsFactors = F)
head(immunomodulators.geneset)
immunomodulators.geneset <- subset(immunomodulators.geneset,
                                   V2 != "Immunoinhibitor" &
                                     V2 != "TIL")
colnames(immunomodulators.geneset) <- c('Genes', 'Type', 'Alisa')
head(immunomodulators.geneset)
table(immunomodulators.geneset$Type)
immunomodulators.geneset <- immunomodulators.geneset[!duplicated(immunomodulators.geneset$Genes), ]
rownames(immunomodulators.geneset) <- immunomodulators.geneset$Genes
immunomodulators.geneset$Alisa[which(immunomodulators.geneset$Alisa == "C10orf54")] <- "VSIR"
head(immunomodulators.geneset)
dim(immunomodulators.geneset)

dim(stad.tcga.t.exp)
tcga.t.immu.exp=stad.tcga.t.exp[rownames(stad.tcga.t.exp) %in% immunomodulators.geneset$Alisa,]
tcga.t.immu.exp=tcga.t.immu.exp[match(immunomodulators.geneset$Alisa,rownames(tcga.t.immu.exp)),]
tcga.t.immu.exp=na.omit(tcga.t.immu.exp)

dim(tcga.t.immu.exp)
length(unique(immunomodulators.geneset$Alisa))
########
immunomodulators.geneset=immunomodulators.geneset[match(rownames(tcga.t.immu.exp),immunomodulators.geneset$Alisa),]
dim(immunomodulators.geneset)

intersect(tcga.PPARs.PCG.cor.sig$to,immunomodulators.geneset$Alisa)

############## 热图
library(ComplexHeatmap)
row.anno = rowAnnotation(Type = immunomodulators.geneset$Type
                         , col=list(Type=c('chemokine'=mycolor[1],
                                           'receptor' = mycolor[2],
                                           'MHC' = mycolor[3],
                                           'Immunostimulator' = mycolor[4]))
                         , annotation_width = unit(c(1,2), 'cm')
                         , annotation_height = unit(0.2, "cm")
                         , gap = unit(1, 'mm'))

df=(stad.tcga.t.PPARs.exp[,gene.selected])
smp.order=names(sort(rowMeans(df)))

# hc.r = hclust(dist(scale(stad.tcga.t.PPARs.exp[,gene.selected])))
# smp.order=rownames(stad.tcga.t.PPARs.exp)[hc.r$order]

ha <- HeatmapAnnotation(df = scale(stad.tcga.t.PPARs.exp[,gene.selected])
                        , col = list(PPARA=circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))
                                     ,PPARD=circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))
                                     ,PPARG=circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))))

Heatmap(as.matrix(t(scale(t(tcga.t.immu.exp))))
        , name = "log2(TPM+1)"
        , row_split = immunomodulators.geneset$Type
        , cluster_rows = T
        , cluster_row_slices = T
        , row_title_gp = gpar(fill = mycolor)
        , show_row_dend = F
        # , column_split = tcga.group$group
        , cluster_columns = F
        # , cluster_column_slices = T
        , show_column_dend = F
        , show_column_names = F
        , col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))
        , left_annotation = row.anno
        , top_annotation = ha
        , border = TRUE)


intersect(tcga.PPARs.PCG.cor$to,immunomodulators.geneset$Alisa)
df=tcga.PPARs.PCG.cor[tcga.PPARs.PCG.cor$to %in% immunomodulators.geneset$Alisa,]
df=data.frame(df)

df <- df %>%
  mutate(lty = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.adj, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
df=df[which(df$p.adj<0.05),]
unique(df$to)

library(ggcor)
immunomodulators.geneset$color = case_when(immunomodulators.geneset$Type=='chemokine' ~ mycolor[1],
                                           immunomodulators.geneset$Type=='receptor' ~ mycolor[2],
                                           immunomodulators.geneset$Type=='MHC' ~ mycolor[3],
                                           immunomodulators.geneset$Type=='Immunostimulator' ~ mycolor[4])
all(rownames(tcga.t.immu.exp)==immunomodulators.geneset$Alisa)

quickcor(t(tcga.t.immu.exp),cor.test = TRUE,type = "lower",show.diag = T) +
  geom_square(data = get_data(p.value < 0.05,type = "lower", show.diag = T)) +
  # geom_number(data = get_data(p.value < 0.05, type = "lower"),aes(num = r),size=1)+
  ggcor::anno_link(df, mapping = aes(colour = col,
                                     size = abs(cor),
                                     linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = rev(red_blue())) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  theme(
    axis.text.y=element_text(family="Times",face="plain",size = 8,colour =immunomodulators.geneset$color), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
    axis.title.y=element_text(family="Times",face="bold",size = 12), #设置y轴标题的字体属性
    legend.text=element_text(face="plain", family="Times", colour="black"), #设置图例的子标题的字体属性
    legend.title=element_text(face="plain", family="Times", colour="black"), #设置图例的总标题的字体属性
    legend.justification=c(1,1), legend.position=c(1,1),
    legend.background = element_rect(fill = NA, colour = NA)
  )+ylab('')+
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 3),
    linetype = "none"
  )
############
######### 免疫特征
get.IOBR.immu.format=function(tcga.t.exp.cibersort){
  tcga.t.exp.cibersort = data.frame(tcga.t.exp.cibersort)
  rownames(tcga.t.exp.cibersort) = tcga.t.exp.cibersort$ID
  tcga.t.exp.cibersort = tcga.t.exp.cibersort[, -1]
  colnames(tcga.t.exp.cibersort) = gsub('(.*)_.*', "\\1", colnames(tcga.t.exp.cibersort))
  return(tcga.t.exp.cibersort)
}
load('origin_datas/immune/tcga.exp.cibersort.RData')
load('origin_datas/immune/tcga.exp.estimate.RData')

tcga.exp.cibersort=get.IOBR.immu.format(tcga.exp.cibersort)
tcga.exp.estimate=get.IOBR.immu.format(tcga.exp.estimate)
dim(tcga.exp.cibersort)

tcga.t.cibersort=tcga.exp.cibersort[rownames(tcga.t.cli_use),1:22]
tcga.t.estimate=tcga.exp.estimate[rownames(tcga.t.cli_use),1:4]
dim(tcga.t.cibersort)
tcga.t.cibersort[1:4,1:5]

write.table(tcga.t.cibersort
            ,file = '06_PPARs_Immune/tcga.t.cibersort.txt'
            ,sep = '\t',row.names = T,col.names = F,quote = F)

df=cal_correlation_parallel(x=stad.tcga.t.PPARs.exp[,gene.selected],y=tcga.t.cibersort
                            ,method = 'pearson')
df=data.frame(df)
head(df)

df <- df %>%
  mutate(lty = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         col = cut(p.adj, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)

########
line.color=c(pal_nejm()(8)[c(1,2)],'grey60')
########
fig6a=quickcor(tcga.t.cibersort, type = "lower") +
  geom_square() + 
  ggcor::anno_link(df, mapping = aes(colour = col,
                                     size = abs(cor),
                                     linetype = lty),diag.label = TRUE) +
  scale_colour_manual(values = line.color)+
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_size_area(max_size = 1) +
  scale_fill_gradient2n(colours = rev(red_blue())) +
  remove_x_axis()+
  theme(
    axis.text.y=element_text(family="Times",face="plain",size = 10), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
    axis.title.y=element_text(family="Times",face="plain",size = 10), #设置y轴标题的字体属性
    legend.text=element_text(family="Times",face="plain", colour="black"), #设置图例的子标题的字体属性
    legend.title=element_text(family="Times",face="plain", colour="black"), #设置图例的总标题的字体属性
    # legend.justification=c(1,1), legend.position=c(1,1),
    legend.key.size =unit(c(0.1), "inches"),
    legend.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = 'red'),
  )+ylab('')+
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "pearson's p", order = 2),
    size = guide_legend(title = "pearson's r", order = 3),
    linetype = guide_legend(title = "pearson's r", order = 2)
  )
fig6a

############################
tme.genesets=readxl::read_excel('origin_datas/TME.geneSets.PMID34019806.xlsx')
tme.genesets=data.frame(tme.genesets)
head(tme.genesets)
#############
tme.type=readxl::read_excel('origin_datas/TME.geneSets.classification.PMID34019806.xlsx',sheet = "Sheet2")
tme.type=data.frame(tme.type,check.names = F,stringsAsFactors = F)
unique(tme.type$Group)
tme.type$color = case_when(tme.type$Group=='Angiogenesis Fibroblasts' ~ mycolor[4],
                           tme.type$Group=='Pro-tumor Immune infiltrate' ~ mycolor[3],
                           tme.type$Group=='Anti-tumor Immune infiltrate' ~ mycolor[1],
                           tme.type$Group=='EMT signature Proliferation rate' ~ mycolor[5])
unique(tme.type$Group)

####
load('origin_datas/immune/tcga.tme.ssgsea.RData')
dim(tcga.tme.ssgsea)
tcga.t.tme.ssgsea=tcga.tme.ssgsea[tme.type$`Process/Signature`,rownames(tcga.t.cli_use)]
dim(tcga.t.tme.ssgsea)

df=cal_correlation_parallel(x=stad.tcga.t.PPARs.exp[,gene.selected],y=t(tcga.t.tme.ssgsea)
                            ,method = 'pearson')
df=data.frame(df)
head(df)

df <- df %>%
  mutate(lty = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         col = cut(p.adj, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
dim(tcga.t.tme.ssgsea)
head(df)

write.table(df,file = '06_PPARs_Immune/tcga.PPARs.tme.cor.res.txt'
            ,sep = '\t',row.names = F,col.names = T,quote = F)

all(rownames(tcga.t.tme.ssgsea)==tme.type$`Process/Signature`)
fig6b=quickcor(t(tcga.t.tme.ssgsea),cor.test = TRUE,type = "lower",show.diag = T) +
  geom_square(data = get_data(p.value < 0.05,type = "lower", show.diag = T)) +
  # geom_number(data = get_data(p.value < 0.05, type = "lower"),aes(num = r),size=1)+
  ggcor::anno_link(df, mapping = aes(colour = col,
                                     size = abs(cor),
                                     linetype = lty),diag.label = TRUE) +
  scale_colour_manual(values = line.color)+
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_size_area(max_size = 1) +
  scale_fill_gradient2n(colours = rev(red_blue())) +
  remove_x_axis()+
  theme(
    axis.text.y=element_text(family="Times",face="plain",size = 10,colour =rev(tme.type$color)), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
    axis.title.y=element_text(family="Times",face="plain",size = 10), #设置y轴标题的字体属性
    legend.text=element_text(family="Times",face="plain", colour="black"), #设置图例的子标题的字体属性
    legend.title=element_text(family="Times",face="plain", colour="black"), #设置图例的总标题的字体属性
    # legend.justification=c(1,1), legend.position=c(1,1),
    legend.key.size =unit(c(0.1), "inches"),
    legend.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = 'red'),
  )+ylab('')+
  guides(
    fill = guide_colourbar(title = "corr", order = 1),
    colour = guide_legend(title = "pearson's p", order = 2),
    size = guide_legend(title = "pearson's r", order = 3),
    linetype = guide_legend(title = "pearson's r", order = 2)
  )
fig6b

######### 与免疫检查点基因的相关性分析

rownames(icg.genes)=icg.genes$`Checkpoint gene`
icg.genes$`Activator(A)/Inhibitor(I)`[icg.genes$`Activator(A)/Inhibitor(I)`=='A']='Immune Activation'
icg.genes$`Activator(A)/Inhibitor(I)`[icg.genes$`Activator(A)/Inhibitor(I)`=='I']='Immune Inhibition'

######### TCGA
setdiff(icg.genes$`Checkpoint gene`,rownames(stad.tcga.t.exp))
setdiff(icg.genes$`Checkpoint gene`,rownames(stad.tcga.t.exp))
rownames(stad.tcga.t.exp)[which(rownames(stad.tcga.t.exp)=='VSIR')]='VISTA'

icg.genes=c('CD274','CTLA4','HAVCR2','LAG3','PDCD1LG2','TIGIT','SIGLEC15')
setdiff(icg.genes,rownames(stad.tcga.t.exp))

tcga.t.exp.icg=stad.tcga.t.exp[rownames(stad.tcga.t.exp) %in% icg.genes,]
dim(tcga.t.exp.icg)

df=cal_correlation_parallel(x=stad.tcga.t.PPARs.exp[,gene.selected],y=t(tcga.t.exp.icg)
                            ,method = 'pearson')
df=data.frame(df)
df2=reshape2::dcast(data = df[,1:3],from ~ to)
library(ggradar)
ggradar(df2,grid.min=-1,grid.mid=0,grid.max = 1,
        values.radar=c(-1,0,1),
        group.line.width = 1, 
        group.point.size = 2,
        group.colours = c("#00AFBB", "#E7B800", "#FC4E07"),
        # Background and grid lines
        background.circle.colour = "white",
        gridline.mid.colour = "grey",
        legend.position = "bottom",
        legend.title = 'PPARs')
###########
##########
library(fmsb)
create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, #线条颜色
    # pfcol = scales::alpha(color, 0.5), 
    plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "black", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
df=cal_correlation_parallel(x=stad.tcga.t.PPARs.exp[,gene.selected],y=t(tcga.t.exp.icg)
                            ,method = 'pearson')
head(df)
write.table(df,file = 'files/文件/tcga.PPARs.immunecheckpoints.txt',row.names = F,col.names = T,sep = '\t',quote = F)
write.table(df,file = '06_PPARs_Immune/tcga.PPARs.immunecheckpoints.txt',row.names = F,col.names = T,sep = '\t',quote = F)

df=data.frame(df)
df2=reshape2::dcast(data = df[,1:3],from ~ to)
rownames(df2)=df2$from
df2=df2[,-1]
head(df2)

maxValue=ceiling(max(abs(df2))*10)/10  #设置刻度
df2=rbind(rep(maxValue,ncol(df2)),rep(-maxValue,ncol(df2)),df2)
head(df2)

pdf(file='PDFs/Fig6C.pdf',height=4.5,width=4.5)
op <- par(mar = c(1, 1, 1, 1))
# Create the radar charts
create_beautiful_radarchart(
  data = df2, caxislabels = seq(-maxValue,maxValue,by = 0.1),
  seg=length(seq(-maxValue,maxValue,by = 0.1))-1,
  color = PPARs.color,
)
# Add an horizontal legend
legend(
  x = "bottom", legend = rownames(df2[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)
dev.off()

pdf(file='06_PPARs_Immune/Fig6C.pdf',height=4.5,width=4.5)
op <- par(mar = c(1, 1, 1, 1))
# Create the radar charts
create_beautiful_radarchart(
  data = df2, caxislabels = seq(-maxValue,maxValue,length.out = 7),
  color = PPARs.color,
)
# Add an horizontal legend
legend(
  x = "bottom", legend = rownames(df2[-c(1,2),]), horiz = TRUE,
  bty = "n", pch = 20 , col = c("#00AFBB", "#E7B800", "#FC4E07"),
  text.col = "black", cex = 1, pt.cex = 1.5
)
par(op)
dev.off()

######### 与PPARs最相关的免疫检查点基因的展示
head(df2)
pairs_use=apply(df2[-c(1,2),], 1,function(x){c(rownames(x),names(which.max(abs(x))))})
pairs_use[1]

library(ggplot2)
library(ggpubr)
library(ggExtra)

PPARs.color
icgs.color=mycolor[c(1,4,3)]

p.all=list()
for(i in 1:length(pairs_use)){
  # i=1
  gene1=names(pairs_use[i])
  gene2=pairs_use[i]
  ######
  x=as.numeric(stad.tcga.t.exp[gene1,])
  y=as.numeric(stad.tcga.t.exp[gene2,])
  #相关性分析
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="pearson")
  cor=corT$estimate
  pValue=corT$p.value
  p1=ggplot(df1, aes(x, y)) + 
    xlab(gene1)+ylab(gene2)+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x =x, y =y))
  p2=ggMarginal(p1, type = "density", xparams = list(fill = PPARs.color[i]),yparams = list(fill = icgs.color[i]))
  p2
  p.all=c(p.all,list(p2))
}
fig6d=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'D')
fig6d
###########
p.all=list()
for(gene in gene.selected){
  df=data.frame(group=stad.tcga.t.PPARs.exp[,gene],t(tcga.t.exp.icg))
  # cutoff=survminer::surv_cutpoint(data.frame(time=tcga.t.cli_2$OS.time/365,
  #                                            event=tcga.t.cli_2$OS,
  #                                            risk=stad.tcga.t.PPARs.exp[,gene]), time = "time", event = "event",
  #                                 variables = c("risk"))
  # cutoff=cutoff$cutpoint$cutpoint
  # df$group=ifelse(df$group>cutoff,'High','Low')
  df$group=ifelse(df$group>median(df$group),'High','Low')
  
  head(df)
  df_melt=reshape2::melt(df, id.vars=c("group"))
  head(df_melt)
  p=ggplot(df_melt,
           aes(x = variable,
               y = value, 
               fill = group))+
    geom_boxplot()+
    stat_compare_means(aes(group=group), label = "p.signif",method='wilcox.test') +
    scale_fill_manual(values = risk.group.color) +
    xlab('') + 
    theme_bw() +
    labs(fill = 'Group') + ylab("log2(TPM+1)") +
    theme_classic()
  if(gene=='PPARG'){
    p=p+theme(axis.text.x=element_text(angle=0,hjust = 0.5,colour="black",family="Times")
              ,axis.ticks = element_line(color = "black")
              ,legend.position = 'right')+labs(fill = gene)
  }else{
    p=p+theme(axis.text.x=element_blank()
              ,axis.ticks = element_line(color = "black")
              ,legend.position = 'right')+labs(fill = gene)
    
  }
  p.all=c(p.all,list(p))
}
fig6e=mg_merge_plot(p.all,nrow = 3,ncol = 1,labels = 'E')
fig6e

p.all=list()
for(gene in gene.selected){
  df=data.frame(group=stad.tcga.t.PPARs.exp[,gene],(tcga.t.cibersort))
  df$group=ifelse(df$group>median(df$group),'High','Low')
  df_melt=reshape2::melt(df, id.vars=c("group"))
  head(df_melt)
  p=ggplot(df_melt,
           aes(x = variable,
               y = value, 
               fill = group))+
    geom_boxplot()+
    stat_compare_means(aes(group=group), label = "p.signif",method='wilcox.test') +
    scale_fill_manual(values = risk.group.color) +
    xlab('') + 
    theme_bw() +
    labs(fill = 'Group') + ylab("Estimated Proportion") +
    theme_classic()
  if(gene=='PPARG'){
    p=p+theme(axis.text.x=element_text(angle=45,hjust = 1,colour="black",family="Times")
              ,axis.ticks = element_line(color = "black")
              ,legend.position = 'right')+labs(fill = gene)
  }else{
    p=p+theme(axis.text.x=element_blank()
              ,axis.ticks = element_line(color = "black")
              ,legend.position = 'right')+labs(fill = gene)
    
  }
  
  p.all=c(p.all,list(p))
}

fig6f=mg_merge_plot(p.all,nrow = 3,ncol = 1,heights = c(1,1,1.5),labels = 'F')
fig6ef=mg_merge_plot(fig6e,fig6f,widths = c(8,10),labels = c('E','F'))

savePDF('PDFs/Fig6A.pdf',fig6a,height = 8,width = 8)
savePDF('PDFs/Fig6B.pdf',fig6b,height = 8,width = 10)
savePDF('PDFs/Fig6D.pdf',fig6d,height = 4.5,width = 13.5)
savePDF('PDFs/Fig6EF.pdf',fig6ef,height = 9,width = 18)

savePDF('06_PPARs_Immune/Fig6A.pdf',fig6a,height = 8,width = 8)
savePDF('06_PPARs_Immune/Fig6B.pdf',fig6b,height = 8,width = 10)
savePDF('06_PPARs_Immune/Fig6D.pdf',fig6d,height = 4.5,width = 13.5)
savePDF('06_PPARs_Immune/Fig6EF.pdf',fig6ef,height = 9,width = 18)

#################### 分子亚型与突变之间的关系
dir.create('07_PPARs_Mutation')
############################
STAD=getTCGAMAFByCode('STAD')
############ Estimates Tumor Mutation Burden in terms of per megabases
tcga.tmb=tmb(STAD,logScale = F)
tcga.tmb=data.frame(tcga.tmb)
write.table(tcga.tmb,file = '07_PPARs_Mutation/tcga.tmb.txt'
            ,sep = '\t',row.names = T,col.names = T,quote = F)
tcga.tmb$SampleCode=paste0(tcga.tmb$Tumor_Sample_Barcode,"-01")
rownames(tcga.tmb)=tcga.tmb$SampleCode
head(tcga.tmb)

df=right_join(tcga.tmb,stad.tcga.t.PPARs.exp,by='SampleCode')

p.all=list()
for(i in 1:length(gene.selected)){
  ######
  gene=gene.selected[i]
  x=as.numeric(df[,gene])
  y=as.numeric(df[,'total_perMB_log'])
  #相关性分析
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="pearson")
  cor=corT$estimate
  pValue=corT$p.value
  p1=ggplot(df1, aes(x, y)) + 
    xlab(gene)+ylab('log10(TMB)')+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x =x, y =y))
  p2=ggMarginal(p1, type = "density", xparams = list(fill = PPARs.color[i]),yparams = list(fill = mycolor[5]))
  p2
  p.all=c(p.all,list(p2))
}
fig7a=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'A')
fig7a

df=right_join(tcga.tmb,stad.tcga.t.PPARs.exp,by='SampleCode')

p.all=list()
for(gene in gene.selected){
  group=df[,gene]
  group=ifelse(group>median(group),'High','Low')
  
  p=plot_ggviolin(data.frame(group
                             ,df[,'total_perMB_log'])
                  ,melt = T,group.col = risk.group.color,ylab='log10(TMB)',leg.title = gene)
  p.all=c(p.all,list(p))
}
fig7b=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'B')

##########
fig7ab=mg_merge_plot(fig7a,fig7b,nrow = 2,ncol = 1,labels = LETTERS[1:2])

savePDF('PDFs/Fig7AB.pdf',fig7ab,height = 8,width = 12)

####### 在线绘制
library(dplyr)
df=right_join(tcga.tmb,stad.tcga.t.PPARs.exp,by='SampleCode')
df[,gene.selected]=apply(df[,gene.selected], 2, function(x){ifelse(x>median(x),'High','Low')})
df=na.omit(df)
head(df)

writeMatrix(df[,c('Tumor_Sample_Barcode',gene.selected)],'files/文件/tcga.PPARs.forMut.txt',row = F)
writeMatrix(df[,c('Tumor_Sample_Barcode','PPARA')],'files/tcga.PPARA.forMut.txt',row = F)
writeMatrix(df[,c('Tumor_Sample_Barcode','PPARD')],'files/tcga.PPARD.forMut.txt',row = F)
writeMatrix(df[,c('Tumor_Sample_Barcode','PPARG')],'files/tcga.PPARG.forMut.txt',row = F)

writeMatrix(df[,c('Tumor_Sample_Barcode',gene.selected)],'07_PPARs_Mutation/tcga.PPARs.forMut.txt',row = F)
writeMatrix(df[,c('Tumor_Sample_Barcode','PPARA')],'07_PPARs_Mutation/tcga.PPARA.forMut.txt',row = F)
writeMatrix(df[,c('Tumor_Sample_Barcode','PPARD')],'07_PPARs_Mutation/tcga.PPARD.forMut.txt',row = F)
writeMatrix(df[,c('Tumor_Sample_Barcode','PPARG')],'07_PPARs_Mutation/tcga.PPARG.forMut.txt',row = F)

table(df$PPARA,df$PPARD)
table(df$PPARA,df$PPARG)
table(df$PPARD,df$PPARG)

########
tcga.MATH=math.score(maf = STAD)
tcga.MATH=data.frame(tcga.MATH)
tcga.MATH=data.frame(tcga.MATH)
tcga.MATH$SampleCode=paste0(tcga.MATH$Tumor_Sample_Barcode,"-01")
rownames(tcga.MATH)=tcga.MATH$SampleCode
head(tcga.MATH)

df=right_join(tcga.MATH,stad.tcga.t.PPARs.exp,by='SampleCode')

p.all=list()
for(i in 1:length(gene.selected)){
  ######
  gene=gene.selected[i]
  x=as.numeric(df[,gene])
  y=as.numeric(df[,'MATH'])
  #相关性分析
  df1=as.data.frame(cbind(x,y))
  corT=cor.test(x,y,method="pearson")
  cor=corT$estimate
  pValue=corT$p.value
  p1=ggplot(df1, aes(x, y)) + 
    xlab(gene)+ylab('Intra-tumor genetic heterogeneity')+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x =x, y =y))
  p2=ggMarginal(p1, type = "density", xparams = list(fill = PPARs.color[i]),yparams = list(fill = mycolor[5]))
  p2
  p.all=c(p.all,list(p2))
}
mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'B')

####################
tcga.immu.lands.p1=readMatrix(paste0(MG_Grobal_baseFolder,'/source/PMC5982584_supplement_2.txt'),row = F)
tcga.immu.lands.p1=data.frame(array=paste0(tcga.immu.lands.p1$`TCGA Participant Barcode`,"-01"),tcga.immu.lands.p1,stringsAsFactors = F,check.names = F)
dim(tcga.immu.lands.p1)
#############

colnames(tcga.immu.lands.p2)[1]=c('array')

dim(tcga.immu.lands.p2)
dim(tcga.immu.lands.p3)

tcga.immu.lands.p2p3=merge(tcga.immu.lands.p2,tcga.immu.lands.p3,by='array',all.x=TRUE,all.y=TRUE,sort = F)
tcga.immu.lands.p2p3=tcga.immu.lands.p2p3[,c('array','LOH_frac_altered','purity','ploidy')]
dim(tcga.immu.lands.p2p3)

setdiff(tcga.immu.lands.p2$array,tcga.immu.lands.p2p3$array)
setdiff(tcga.immu.lands.p3$array,tcga.immu.lands.p2p3$array)

tcga.immu.lands=merge(tcga.immu.lands.p1,tcga.immu.lands.p2p3,by='array',all.x=TRUE,all.y=TRUE,sort = F)
dim(tcga.immu.lands)

setdiff(tcga.immu.lands.p2p3$array,tcga.immu.lands$array)
setdiff(tcga.immu.lands.p1$array,tcga.immu.lands$array)
table(tcga.immu.lands$`TCGA Study`)

tcga.immu.lands<-tcga.immu.lands[which(tcga.immu.lands$`TCGA Study`=='STAD'),]
rownames(tcga.immu.lands)=tcga.immu.lands$array
dim(tcga.immu.lands)
head(tcga.immu.lands)
tcga.immu.lands$SampleCode=tcga.immu.lands$array

setdiff(stad.tcga.t.PPARs.exp$SampleCode,tcga.immu.lands$SampleCode)

col.selected=c('Aneuploidy Score','Homologous Recombination Defects','Intratumor Heterogeneity','LOH_frac_altered','purity','ploidy')

df=full_join(stad.tcga.t.PPARs.exp,tcga.immu.lands[,c('SampleCode',col.selected)],by='SampleCode')
head(df)

p.all=list()
for(i in 1:length(gene.selected)){
  for(j in 1:length(col.selected)){
    ######
    gene=gene.selected[i]
    x=as.numeric(df[,gene])
    y=as.numeric(df[,col.selected[j]])
    #相关性分析
    df1=as.data.frame(cbind(x,y))
    corT=cor.test(x,y,method="pearson")
    cor=corT$estimate
    pValue=corT$p.value
    p1=ggplot(df1, aes(x, y)) + 
      xlab(gene)+ylab(col.selected[j])+
      geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
      stat_cor(method = 'spearman', aes(x =x, y =y))
    p2=ggMarginal(p1, type = "density", xparams = list(fill = PPARs.color[i]),yparams = list(fill = mycolor[5]))
    p2
    p.all=c(p.all,list(p2))
  }
  
}
length(p.all)
mg_merge_plot(p.all,nrow = 6,ncol = 3,labels = 'B')

##################### 药物IC50
############# GDSC数据库细胞系表达谱
dir.create('08_Drug_sensitivity')
GDSC_exp<-read.table('origin_datas/GDSC/Cell_line_RMA_proc_basalExp/Cell_line_RMA_proc_basalExp.txt',header = T,check.names = F,sep = '\t')
GDSC_exp[1:4,1:5]
setdiff(gene.selected,GDSC_exp$GENE_SYMBOLS)
GDSC_exp$GENE_SYMBOLS[grep("PPARgamma",GDSC_exp$GENE_SYMBOLS)]

table(duplicated(GDSC_exp$GENE_SYMBOLS))
######## 去掉重复的GENE_SYMBOLS
GDSC_exp<-GDSC_exp[!duplicated(GDSC_exp$GENE_SYMBOLS),]
rownames(GDSC_exp)<-GDSC_exp$GENE_SYMBOLS
GDSC_exp$GENE_title<-NULL
GDSC_exp$GENE_SYMBOLS<-NULL
GDSC_exp[1:4,1:5]
range(GDSC_exp)
dim(GDSC_exp)
############# 药物IC50
gdsc_drug_ic<-read.csv('origin_datas/GDSC/PANCANCER_IC_Fri Jun 10 02_20_40 2022.csv',header = T)
dim(gdsc_drug_ic)
head(gdsc_drug_ic)
colnames(gdsc_drug_ic)

table(gdsc_drug_ic$Tissue.sub.type)
gdsc_drug_ic=gdsc_drug_ic[which(gdsc_drug_ic$Tissue.sub.type=='stomach'),]
head(gdsc_drug_ic)

gdsc_drug_ic[,10:13]<-NULL
gdsc_drug_ic[,c(5:8)]<-NULL
head(gdsc_drug_ic)
# colnames(gdsc_drug_ic)
# gdsc_drug_ic=gdsc_drug_ic[,c('Drug.name','Drug.Id','Cell.line.name','Cosmic.sample.Id','AUC')]
colnames(gdsc_drug_ic)[4]<-"sample"
head(gdsc_drug_ic)

#################################
dim(GDSC_exp)
GDSC.PPARs.exp=t(GDSC_exp[rownames(GDSC_exp) %in% gene.selected,])
GDSC.PPARs.exp=data.frame(GDSC.PPARs.exp)
GDSC.PPARs.exp$sample=gsub("DATA.","",rownames(GDSC.PPARs.exp))
head(GDSC.PPARs.exp)

GDSC_drug<-merge(gdsc_drug_ic,GDSC.PPARs.exp,by='sample')
head(GDSC_drug)

write.table(GDSC_drug,file = '08_Drug_sensitivity/GDSC.PPARs.drug.txt'
            ,sep = '\t',row.names = F,col.names = T,quote = T)

drug_name<-as.data.frame(table(GDSC_drug$Drug.name))
colnames(drug_name)<-c("drug","counts")

library(psych)
GDSC.PPARs.drug.cor=c()
for(gene in gene.selected){
  if(!gene %in% colnames(GDSC_drug)){
    next
  }
  df_cor=GDSC_drug %>%
    group_by(Drug.name) %>%
    summarise(cor=corr.test(AUC,!!sym(gene),method = "pearson",adjust = "fdr")$r
              ,pvalue=corr.test(AUC,!!sym(gene),method = "pearson",adjust = "fdr")$p
              ,FDR=corr.test(AUC,!!sym(gene),method = "pearson",adjust = "fdr")$p.adj)
  df_cor=as.data.frame(df_cor)
  df_cor$PPARs=gene
  GDSC.PPARs.drug.cor=rbind(GDSC.PPARs.drug.cor,df_cor)
}
head(GDSC.PPARs.drug.cor)

write.table(GDSC.PPARs.drug.cor,file = '08_Drug_sensitivity/GDSC.PPARs.drug.cor.txt'
            ,sep = '\t',row.names = F,col.names = T,quote = T)

########
cor.cutoff=0
cor.p.cutoff=0.05
GDSC.PPARs.drug.cor$Class=case_when(GDSC.PPARs.drug.cor$cor>cor.cutoff & GDSC.PPARs.drug.cor$pvalue<cor.p.cutoff ~ "positive",
                                   GDSC.PPARs.drug.cor$cor < (-cor.cutoff) & GDSC.PPARs.drug.cor$pvalue < cor.p.cutoff ~ "negative",
                                   TRUE ~ "ns" )
head(GDSC.PPARs.drug.cor)

p.all=list()
for(gene in unique(GDSC.PPARs.drug.cor$PPARs)){
  # gene='PPARA'
  df=GDSC.PPARs.drug.cor[which(GDSC.PPARs.drug.cor$PPARs==gene),]
  df<-df[which(df$Class != 'ns'),]
  df<-df[order(df$cor),]
  df$`-log10(pvalue)`<-(-log10(df$pvalue))
  range(df$`-log10(pvalue)`)
  
  library(tidyverse)
  p=df %>% 
    ggplot(aes(reorder(Drug.name, cor), cor)) + 
    geom_col(aes(fill = `-log10(pvalue)`)) + 
    coord_cartesian(ylim = c(-1, 1))+
    scale_fill_gradient2(low = "blue", 
                         mid = 'white',
                         high = "red", 
                         midpoint = 1.3) + 
    labs(x = "")+
    labs(y = paste("Correlation of drug sensitivity and",gene))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1,family ='Times',colour ='black',size =10),
          axis.title.x = element_text(family ='Times',colour ='black',size =10))
  p
  p.all=c(p.all,list(p))
}
fig8a=mg_merge_plot(p.all,nrow = 1,ncol = 2,widths = c(2,1),labels = LETTERS[1:2])
fig8a

########################### 药物靶向的通路
drug_pathway<-read.csv('origin_datas/GDSC/PANCANCER_ANOVA_Fri Jun 10 02_18_35 2022.csv',header = T,stringsAsFactors =F)
dim(drug_pathway)
colnames(drug_pathway)
drug_pathway[,c(5:22)]<-NULL
drug_pathway[,c(2:3)]<-NULL
head(drug_pathway)
drug_pathway$target_pathway
#########
inds=which(GDSC.PPARs.drug.cor$Class != 'ns' & GDSC.PPARs.drug.cor$PPARs=='PPARA')
GDSC.PPARs.drug.cor.sig=GDSC.PPARs.drug.cor[inds,]
head(GDSC.PPARs.drug.cor.sig)
#########
drug_pathway<-drug_pathway[which(drug_pathway$drug_name %in% GDSC.PPARs.drug.cor.sig$Drug.name),]
unique(drug_pathway$drug_name)

head(drug_pathway)
drug_pathway_use<-data.frame(table(drug_pathway),stringsAsFactors = F)
drug_pathway_use$target_pathway=as.character(drug_pathway_use$target_pathway)
drug_pathway_use$Drugs=as.character(drug_pathway_use$Drugs)

head(drug_pathway_use)
drug_pathway_use<-drug_pathway_use[which(drug_pathway_use$Freq != 0),]
colnames(drug_pathway_use)<-c("Drugs","pathway","freq")
head(drug_pathway_use)
dim(drug_pathway_use) ### 15个通路

drug_pathway_use_final<-merge(drug_pathway_use,
                              GDSC.PPARs.drug.cor.sig,
                              by.x='Drugs',by.y='Drug.name')
head(drug_pathway_use_final)

drug_pathway_use_final[,3]<-NULL
drug_pathway_use_final[,6:8]<-NULL
head(drug_pathway_use_final)

length(drug_pathway_use_final$pathway)
#####
table(drug_pathway_use_final$cor<0)
##### 敏感
drug_pathway_use_final=drug_pathway_use_final[which(drug_pathway_use_final$cor<0),]

pathway<-as.data.frame(table(drug_pathway_use_final$pathway))
pathway_use<-as.character(pathway$Var1)
pathway_use
length(pathway_use)

########## 相关性
head(drug_pathway_use_final)
table(drug_pathway_use_final$Drugs,drug_pathway_use_final$pathway)

down=reshape2::dcast(drug_pathway_use_final,Drugs ~ pathway,value.var="cor",fill = 0)
rownames(down)=down$Drugs
down=down[,-1]
down=t(down)
down=as.matrix(down)
down

######### 相关性显著性p值
head(drug_pathway_use_final)

up=reshape2::dcast(drug_pathway_use_final,Drugs ~ pathway,value.var='pvalue',fill=0)
rownames(up)=up$Drugs
up=up[,-1]
up=t(up)
up=(-log10(up))
up[which(is.infinite(up),arr.ind = )]=0
up=as.matrix(up)

############ 可视化
library("ComplexHeatmap")
library("circlize")
plotMutiHeatmap=function(up,down,up_break,up_colors,down_break,down_colors,title){
  UpColor <- colorRamp2(breaks = up_break, colors = up_colors)
  DnColor <- colorRamp2(breaks = down_break, colors = down_colors)
  
  DiagFunc <- function(up, down){
    function(j, i, x, y, width, height, fill){
      grid.polygon(unit.c(x - 0.5*width, x - 0.5*width, x + 0.5*width), 
                   unit.c(y - 0.5*height, y + 0.5*height, y - 0.5*height),
                   gp = gpar(fill = DnColor(down[i, j]), col = "grey")) 
      
      grid.polygon(unit.c(x + 0.5*width, x + 0.5*width, x - 0.5*width), 
                   unit.c(y + 0.5*height, y - 0.5*height, y + 0.5*height),
                   gp = gpar(fill = UpColor(up[i, j]), col = "grey"))
      ##### 显著性
      if(up[i, j]>=1.3){
        txt="***"
        if(up[i, j]>=1.3&up[i, j]<2){
          txt='*'
        }else if(up[i, j]>=2&up[i, j]<3){
          txt='**'
        }else if(up[i, j]>=3&up[i, j]<4){
          txt='***'
        }
        grid.text(label=txt,x=(x + 0.5*width),
                  y=(y+ 0.5*height),just = c('right','top'))
      }
      if(down[i, j]>0.1){
      }
    }
  }
  ######## 热图
  p1 <- Heatmap(up, column_title = title
                , rect_gp = gpar(type = "none")
                , show_heatmap_legend = F
                , cluster_rows = T
                , show_row_dend = F
                , row_names_gp = gpar(fontsize = 8)
                , column_names_gp = gpar(fontsize = 8)
                , cluster_columns = T
                , show_column_dend = F
                , cell_fun = DiagFunc(up = up, down = down) 
  ) 
  ############## legend correlation
  col_fun = colorRamp2(down_break, down_colors) 
  lgd <- Legend(title = "cor", 
                col_fun = col_fun, 
                grid_width = unit(1, "mm"),
                # at = c(-0.3,0,0.3), 
                # labels = c("-0.3","0","0.3"),  
                direction = "vertical" 
  )
  ############## legend correlation p value
  col_fun2 = colorRamp2(up_break, up_colors) 
  lgd2 <- Legend(title = "p.val", 
                 col_fun = col_fun2, 
                 grid_width = unit(1, "mm"),
                 at = c(0,1,2,3,4,5), 
                 labels = c('0',"1","2","3","4",">5"),  
                 direction = "vertical"
  )
  
  draw(p1, annotation_legend_list = list(lgd,lgd2), annotation_legend_side = "right"
       ,heatmap_legend_side = "right", merge_legend = TRUE)
}

range(up)
range(down)

up_break=c(0, 2)
down_break=c(-1,0,1)
up_colors=c("#FFFFFF","#6f9a8d")
down_colors=c("blue",'white',"red")
range(down)
range(up)

pdf('PDFs/Fig8C.pdf',height = 4,width = 6)
plotMutiHeatmap(up=up,down=down
                ,up_break=up_break
                ,up_colors=up_colors
                ,down_break=down_break
                ,down_colors=down_colors,title='')
dev.off()

pdf('08_Drug_sensitivity/Fig8C.pdf',height = 4,width = 6)
plotMutiHeatmap(up=up,down=down
                ,up_break=up_break
                ,up_colors=up_colors
                ,down_break=down_break
                ,down_colors=down_colors,title='')
dev.off()

########### 药物IC50预测
get_drug_IC50=function(exp.dat,drug){
  library(pRRophetic)
  ic50.res=pRRopheticPredict(as.matrix(exp.dat)
                             , drug
                             , "urogenital_system"
                             , selection = 1
                             , dataset = "cgp2016")
  return(ic50.res)
}

library(pRRophetic)
data(PANCANCER_IC_Tue_Aug_9_15_28_57_2016)

gdsc.drugs=unique(as.character(drug_pathway_use_final$Drugs))
gdsc.drugs=intersect(gdsc.drugs,drugData2016$Drug.name)

tcga.drug.ic50=c()
for(drug in gdsc.drugs){
  ic50=get_drug_IC50(exp.dat = stad.tcga.t.exp,drug = drug)
  ic50=as.matrix(ic50)
  colnames(ic50)=drug
  tcga.drug.ic50=cbind(tcga.drug.ic50,ic50)
}
head(tcga.drug.ic50)

#相关性分析
all(stad.tcga.t.PPARs.exp$SampleCode==rownames(tcga.drug.ic50))

plot_cor_point=function(x,y,method='Pearson',top_col='#D55E00',right_col='#009E73'
                      ,ylab='y expression',xlab='x expression',title=NULL
                      ,marginal.type=c("histogram", "boxplot", "density", "violin", "densigram")[1]){
  library(ggplot2)
  library(ggpubr)
  library(ggExtra)
  corT=cor.test(x,y,method="pearson")
  cor=corT$estimate
  pValue=corT$p.value
  df=data.frame(x,y)
  p=ggplot(df, aes(x, y)) + 
    xlab(xlab)+ylab(ylab)+
    geom_point()+ geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
    stat_cor(method = 'pearson', aes(x =x, y =y))
  p=ggMarginal(p, type = marginal.type, xparams = list(fill = top_col),yparams = list(fill = right_col))
  return(p)
}

p.all=list()
for(i in 1:length(gene.selected)){
  gene=gene.selected[i]
  df=data.frame(cbind(stad.tcga.t.PPARs.exp[,gene],tcga.drug.ic50),check.names = F)
  colnames(df)
  p=plot_cor_point(x = df[,2],y = df[,1],
                   ylab = gene,xlab = paste(colnames(df)[2],'IC50'),
                   right_col = PPARs.color[i],top_col = pal_npg()(9)[2],
                   marginal.type='density')
  p.all=c(p.all,list(p))
}
fig8d=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'D')
fig8d

p.all=list()
df=as.data.frame(cbind(tcga.drug.ic50,stad.tcga.t.PPARs.exp))
for(gene in gene.selected){
  group=df[,gene]
  group=ifelse(group>median(group),'High','Low')
  
  p=plot_ggviolin(data.frame(group
                             ,df[,"OSI-027"])
                  ,melt = T,group.col = risk.group.color,ylab='Estimated IC50 of OSI-027',leg.title = gene)
  p
  p.all=c(p.all,list(p))
}
fig8e=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = 'B')
fig8e

###### 常用的化疗药物
drug.selected_use=c('Docetaxel','Vinorelbine','Paclitaxel','Cisplatin')

for(drug in drug.selected_use){
  ic50=get_drug_IC50(exp.dat = stad.tcga.t.exp,drug = drug)
  ic50=as.matrix(ic50)
  colnames(ic50)=drug
  tcga.drug.ic50=cbind(tcga.drug.ic50,ic50)
}
head(tcga.drug.ic50)
tcga.drug.ic50=data.frame(tcga.drug.ic50)

write.table(tcga.drug.ic50,file = 'files/文件/tcga.drug.ic50.txt'
            ,sep = '\t',row.names = T,col.names = T,quote = F)
write.table(tcga.drug.ic50,file = '08_Drug_sensitivity/tcga.drug.ic50.txt'
            ,sep = '\t',row.names = T,col.names = T,quote = F)

df=cbind(tcga.drug.ic50,stad.tcga.t.PPARs.exp[,gene.selected])
p.all=list()
for(drug in drug.selected_use){
  df1=df[,c(drug,gene.selected)]
  df2=pivot_longer(data=rownames_to_column(as.data.frame(df1)),
               cols = 3:(ncol(df1)+1),
               names_to = "PPARs",
               values_to = "exprs")
  p=ggplot(df2, aes(x = !!sym(drug), y = exprs,color=PPARs))+
    geom_point(size=1)+
    geom_smooth(method = "lm")+
    stat_cor(method = "pearson",aes(x = !!sym(drug), y = exprs, color = PPARs))+
    scale_color_manual(values = PPARs.color)+
    theme_bw()+
    xlab(paste(drug,'IC50'))+ylab('log2(TPM+1)')
  p
  p.all=c(p.all,list(p))
}
length(p.all)

fig8f=mg_merge_plot(p.all,nrow = 1,ncol = 4,labels = 'F',common.legend = T)
fig8f

fig8def=mg_merge_plot(fig8d,fig8e,fig8f,nrow = 3,ncol = 1)
savePDF('PDFs/Fig8DEF.pdf',fig8def,height = 12,width = 12)
savePDF('PDFs/Fig8A.pdf',fig8a,height = 4,width = 8)

savePDF('08_Drug_sensitivity/Fig8DEF.pdf',fig8def,height = 12,width = 12)
savePDF('08_Drug_sensitivity/Fig8A.pdf',fig8a,height = 4,width = 8)

############## 甲基化分析
load(paste0(MG_Grobal_baseFolder,'/source/genome_info/cpg.450.27.cpg2enstbytss.Rdata'))
load(paste0(MG_Grobal_baseFolder,'/source/gencode.v33.id.tab.RData'))
dat=cbind(cpg2enstbytss,gencode.v33.id.tab[match(cpg2enstbytss[,2],gencode.v33.id.tab[,1]),-1])
colnames(dat)=c('CpG','ENST','TSS','ENSG','Symbol','ENST_Type','LncRNA')
#########
cpg2symbol=dat[match(gene.selected,dat$Symbol),]
cpg2symbol=cpg2symbol[,c(1,5)]
head(cpg2symbol)

######
stad.methy=getTCGAMethyCpGByCode('STAD')
stad.methy450=stad.methy$M450k
dim(stad.methy450)
setdiff(cpg2symbol$CpG,rownames(stad.methy450))
head(cpg2symbol)

stad.methy450=stad.methy450[match(cpg2symbol$CpG,rownames(stad.methy450)),]
stad.methy450=na.omit(stad.methy450)
stad.methy450=data.frame(t(stad.methy450))
stad.methy450$SampleCode=rownames(stad.methy450)
head(stad.methy450)

df=merge(stad.tcga.t.PPARs.exp[,c('SampleCode','PPARD')],stad.methy450,by='SampleCode')
head(df)
plot_cor_point(x = df[,2],y=df[,3])

save.image(file = '20221209_STAD_PPARs_final.RData')
