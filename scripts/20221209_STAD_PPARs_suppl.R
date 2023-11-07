
head(stad.tcga.t.PPARs.exp_use)
library(gtsummary)
stad.tcga.t.PPARs.exp_use=stad.tcga.t.PPARs.exp_use[,-10]
colnames(stad.tcga.t.PPARs.exp_use)[12]='Age'

t1=stad.tcga.t.PPARs.exp_use %>%
  filter(PPARs=='PPARA') %>%
  select(Type,T.Stage,N.Stage,M.Stage,Stage,Grade,Age,Gender)%>%
  tbl_summary(by = Type) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  add_overall() %>%
  add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**PPARA**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  modify_caption("**Table 1. Relationship Between PPARs family member expression and clinicopathological features in the TCGA-STAD cohort.**") %>%
  bold_labels()

t2=stad.tcga.t.PPARs.exp_use %>%
  filter(PPARs=='PPARD') %>%
  select(Type,T.Stage,N.Stage,M.Stage,Stage,Grade,Age,Gender)%>%
  tbl_summary(by = Type) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  # add_overall() %>%
  # add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**PPARD**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  modify_caption("**Table 1. Relationship Between PPARs family member expression and clinicopathological features in the TCGA-STAD cohort.**") %>%
  bold_labels()

t3=stad.tcga.t.PPARs.exp_use %>%
  filter(PPARs=='PPARG') %>%
  select(Type,T.Stage,N.Stage,M.Stage,Stage,Grade,Age,Gender)%>%
  tbl_summary(by = Type) %>%
  add_p(pvalue_fun = ~style_pvalue(.x, digits = 2)) %>%
  # add_overall() %>%
  # add_n() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**PPARG**") %>%
  modify_footnote(
    all_stat_cols() ~ "Median (IQR) or Frequency (%)"
  ) %>%
  modify_caption("**Table 1. Relationship Between PPARs family member expression and clinicopathological features in the TCGA-STAD cohort.**") %>%
  bold_labels()


tbl_merge <-
  tbl_merge(tbls = list(t1, t2, t3),
            tab_spanner = c("**PPARA**", "**PPARD**","**PPARG**"))

tbl_merge

tbl_merge %>%    # build gtsummary table
  as_gt() %>%             # convert to gt table
  gt::gtsave(             # save table as image
  


library(gdtools)
library(flextable)
library(officer)
library(webshot2)
library(writexl)


sect_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 20, height = 16
  )
  # type = "continuous"
  # page_margins = page_mar()
)

tbl_merge %>%
  as_flex_table() %>%
  flextable::save_as_docx(path = f,pr_section = sect_properties)

###############

tbl_merge %>%
  gtsummary::as_tibble() %>%
  writexl::write_xlsx(., path = f)

tbl_merge %>%
  as_gt() %>%             # convert to gt table
  gt::gtsave(             # save table as image
    filename = "Table1.png",
    path='Z:/users/lishuang/Work1/Work/20221209_STAD_PPARs/PDFs/'
  )


################################


load('origin_datas/TCGA/stad.tcga.exp.tpm.RData')
stad.tcga.exp=stad.tcga.exp.tpm
dim(stad.tcga.exp)

library(GSVA)
library(GSEABase)

gmtFile='origin_datas/h.all.v7.5.1.symbols.gmt'

c2KEGG <- getGmt(gmtFile,
                 collectionType=BroadCollection(category="c2"),
                 geneIdType=SymbolIdentifier())
tcga.h.all.ssgsea <- gsva(as.matrix(stad.tcga.exp),
                         c2KEGG,
                         method = 'ssgsea',
                         min.sz = 10,
                         max.sz = 500,
                         verbose = TRUE)
dir.create('origin_datas/immune')
save(tcga.h.all.ssgsea,file='origin_datas/immune/tcga.h.all.ssgsea.RData')

gmtFile='origin_datas/c2.cp.kegg.v7.5.1.symbols.gmt'
c2KEGG <- getGmt(gmtFile,
                 collectionType=BroadCollection(category="c2"),
                 geneIdType=SymbolIdentifier())
tcga.kegg.ssgsea <- gsva(as.matrix(stad.tcga.exp),
                         c2KEGG,
                         method = 'ssgsea',
                         min.sz = 10,
                         max.sz = 500,
                         verbose = TRUE)
save(tcga.kegg.ssgsea,file='origin_datas/immune/tcga.kegg.ssgsea.RData')

#####################
tme.genesets=readxl::read_excel('origin_datas/TME.geneSets.PMID34019806.xlsx')
tme.genesets=data.frame(tme.genesets)
head(tme.genesets)

pathway.genesets=readxl::read_excel('origin_datas/TME.geneSets.PMID34079012.xlsx')
pathway.genesets=data.frame(pathway.genesets)
head(pathway.genesets)

tme.genesets.list=split(x=tme.genesets,f=tme.genesets$Gene.signature)
tme.genesets.list=sapply(tme.genesets.list, function(x){subset(x,select='Gene',drop=TRUE)})

pathway.genesets.list=split(x=pathway.genesets,f=pathway.genesets$Pathway)
pathway.genesets.list=sapply(pathway.genesets.list, function(x){subset(x,select='Gene',drop=TRUE)})

tcga.tme.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = stad.tcga.exp
                                                ,genelist = tme.genesets.list)

tcga.pathway.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = stad.tcga.exp
                                                    ,genelist = pathway.genesets.list)


save(tcga.tme.ssgsea,file = 'origin_datas/immune/tcga.tme.ssgsea.RData')
save(tcga.pathway.ssgsea,file = 'origin_datas/immune/tcga.pathway.ssgsea.RData')

cancer.immu.genesets=readxl::read_excel('origin_datas/cancer-immunity cycle.geneSets.PMID30154154.xlsx')
cancer.immu.genesets=data.frame(cancer.immu.genesets)
head(cancer.immu.genesets)

cancer.immu.genesets.list=split(x=cancer.immu.genesets,f=cancer.immu.genesets$Signature.set)
cancer.immu.genesets.list=sapply(cancer.immu.genesets.list, function(x){subset(x,select='GeneSymbol',drop=TRUE)})

tcga.cancer.immu.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = stad.tcga.exp
                                                        ,genelist = cancer.immu.genesets.list)

save(tcga.cancer.immu.ssgsea,file = 'origin_datas/immune/tcga.cancer.immu.ssgsea.RData')
##################
library(IOBR)

#### ESTIMATE
tcga.exp.estimate<-deconvo_estimate(eset=stad.tcga.exp)
save(tcga.exp.estimate,file='origin_datas/immune/tcga.exp.estimate.RData')

### CIBERSORT
tcga.exp.cibersort<-deconvo_cibersort(eset=stad.tcga.exp,arrays=T)
save(tcga.exp.cibersort,file='origin_datas/immune/tcga.exp.cibersort.RData')

