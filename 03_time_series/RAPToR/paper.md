
# Real age prediction from the transcriptome with RAPToR

Romain Bulteauand Mirko Francesconi✉

Transcriptomic data is often affected by uncontrolled variation among samples that can obscure and confound the effects of interest. This variation is frequently due to unintended differences in developmental stages between samples. The transcriptome itself can be used to estimate developmental progression, but existing methods require many samples and do not estimate a specimen's real age. Here we present real-age prediction from transcriptome staging on reference (RAPToR), a computational method that precisely estimates the real age of a sample from its transcriptome, exploiting existing time-series data as reference. RAPToR works with whole animal, dissected tissue and single-cell data for the most common animal models, humans and even for non-model organisms lacking reference data. We show that RAPToR can be used to remove age as a confounding factor and allow recovery of a signal of interest in differential expression analysis. RAPToR will be especially useful in large-scale single-organism profiling because it eliminates the need for accurate staging or synchronisation before profiling.

G enome-wide profiling of gene expression is a powerful technique that provides a global and unbiased view of the transcriptional state of a biological system. However, the analysis of  gene-expression  data  can  be  complicated  by  uncontrolled  and unknown sources of variance-which may be technical or biological in nature 1 -that can mask or confound the effects of variables of interest.

To  tackle  this  problem,  several  methods  have  been  developed to  learn  and  remove  hidden  covariates  (or  surrogate  variables) from the data, such as remove unwanted variance 2 , surrogate variable analysis 3 ,  or  probabilistic  estimation  of  expression  residuals 4 . However, a drawback of these methods is that the sources of variance usually remain obscure; therefore, potentially interesting biological variance might also be removed.

Unintended  differences  in  developmental  progression  across biological replicates or experimental conditions are a major source of variance in gene-expression data of developing systems (Fig. 1a), which can confound (Fig. 1b) or mask (Fig. 1c) the effect of the variable of interest. This is especially true in organisms with rapid life  cycles  and  highly  variable  growth  speed  such  as  worms,  fruit fly, or zebrafish, where numerous factors like genetic background, temperature,  diet,  crowding 5-9 ,  or  even  the  physiological  state  of the previous generation 9 substantially impact developmental speed. Carefully controlling for all conditions influencing development is therefore particularly challenging, but failing to do so can strongly impact  gene  expression.  For  example,  in Caenorhabditis  elegans even a few hours of development may result in 10,000 differentially expressed genes 10 . Hence, it is not surprising that around 50% of variance in gene expression in the profiling of a large panel of C. elegans recombinant  inbred  lines 11 is  due  to  unintended  developmental variation 12 and that almost 38% of the datasets that did not intend to include development in a C. elegans gene-expression database 13 show substantial developmental variation in gene expression 10 .

Estimating the real physiological age of the samples and identifying hidden developmental variation between them is important first to quantify the impact of the perturbation of interest on developmental speed; second, to distinguish perturbation-specific from unspecific changes in gene expression caused by development; third, to uncover time-specific effects of the perturbations under study 12 by including inferred age as a covariate in expression data analyses (such as differential expression analysis). In yeast, analogous ideas were successfully implemented to identify genetic and environmental  perturbations  impacting specific phases of the cell cycle 14 and direct and specific effects of 700 gene deletions on gene expression after removing the main source of variance (25%): a shared expression signature of cell cycle and growth rate 15 .

Extracting developmental progression from transcriptomes has recently  become  a  topic  of  intense  research,  especially  after  the advent of single-cell RNA sequencing. Many algorithms have been developed  that  learn  developmental  progression  from  large-scale bulk,  single-cell,  or  whole-organism  transcriptomic  data  and  sort samples  along  those  trajectories  (for  example  Slingshot 16 ,  DPT 17 , Monocle 18 ,  and  BLIND 19 ).  However,  a  major  drawback  of  these trajectory-learning algorithms is that they require large numbers of samples to learn the trajectory of developmental expression changes directly from the data. Moreover, they only output dataset-specific ranks  or  arbitrary  values  usually  referred  to  as  'pseudo-times' rather  than  real-age  predictions,  making  it  difficult  to  compare results across datasets or conditions.

To  overcome  these  limitations,  we  developed  a  computational framework that exploits available time-series gene-expression data as  reference  to  determine  the  absolute  age  of  even  a  single  sample  from  its  transcriptome  with  high  precision.  We  implemented RAPToR  in  R  (available  at  https://github.com/LBMC/RAPToR), providing references to stage C.  elegans , Drosophila melanogaster, Danio rerio , and Mus musculus development from gene expression.

We show that RAPToR successfully estimates age during development and ageing, estimates tissue-specific age from whole-organism data,  works  in  dissected  tissue  and  single-cell  profiling,  and  can also estimate age of one species using another species as reference. Finally, we show how to use estimated ages to quantify a perturbation effect on developmental or ageing speed, and recover the specific effects of the variables of interest on gene expression even when completely confounded by age.

Fig. 1 | estimating age from the transcriptome using RAPToR. a , Cartoon of individuals sampled at identical chronological times in two conditions, resulting in different developmental age between groups owing to the condition impacting developmental speed. b -c , Cartoons of differential expression analysis situations where hidden developmental variation is either misinterpreted as an effect of the condition due to development and condition being confounded ( b ), or masking an effect of the condition owing to developmental spread ( c ). d -f , RAPtoR staging exploits existing reference time-series expression data ( d ). this data is first decomposed into principal/independent components, which are interpolated with respect to time ( e ). Comp., component. Interpolated reference is then reconstructed at gene level with interpolated components and gene loadings ( f ). g , h , For each sample, a correlation profile is built by computing genome-wide Spearman correlation with every time point of the reference ( g ). Corr., correlation. the reference time with maximal correlation becomes the estimate, and bootstrapping on random gene subsets defines a confidence interval (see Methods) ( h ).

## Results

RAPToR design. We set out to develop a strategy to continuously stage  development  from  gene  expression  that  would  be  effective even for experiments with a limited number of samples, in which trajectory learning methods are not applicable. We reasoned that we could exploit existing developmental time-series data as reference to stage samples individually by taking the time point of the reference that has maximum correlation with a sample transcriptome as the age estimate. In this way, the age of each sample is inferred independently from others and outliers only influence their own staging (as opposed to trajectory-based approaches). Furthermore, age estimates acquired on the same reference should be comparable even across different experiments, conditions, and genetic backgrounds.

However,  using  the  time  of  maximum  correlation  with  the reference  as  the  age  estimate  limits  its  precision  to  the  temporal resolution of the reference. To overcome this limitation, we interpolate reference gene expression (Fig. 1d) with respect to time in a dimensionally-reduced space (Fig. 1e and Supplementary Note 1), generating interpolated expression profiles potentially at any time between the original reference time points (Fig. 1f).

The sample age estimate is simply the time point of maximum Spearman  correlation  between  the  interpolated  reference  and the  sample  gene  expression  (Fig.  1g).  We  then  compute  a  confidence interval of the estimate by bootstrapping on genes (Fig. 1h; Methods).

We implemented this strategy in RAPToR, an R package where we provide functions to interpolate references and stage samples.

## RAPToR accurately  infers  developmental  age  of  model  organ-

isms. To test RAPToR in the most commonly used animal model organisms,  we  built  interpolated  references  exploiting  existing time-series  data  on C.  elegans roundworm  embryonic  and  larval development 20-22 ,  zebrafi s h  embryonic  and  larval  development 23 , mouse 24 ,  and  fruit  fly 25 embryonic  development  (Supplementary Table 1) and then staged independent time-series experiments of C. elegans late-larval development 26 and zebrafish 27,28 , mouse 29 , and fruit fly 27  embryonic development.

We found RAPToR age estimates accurately match chronological  age  for  both C.  elegans and  zebrafish  ( R 2 &gt; 0.99;  Figs.  2a,b), and morphological staging (somite number) for mice ( R 2 = 0.95; Fig. 2c). Age estimates of fly single embryos 27 less accurately match chronological  age,  especially  at  later  stages  ( R 2 = 0.74;  Fig.  2d). However,  this  is  likely  due  to  the  single-individual  nature  of data  as  any  inter-individual  variation  in  developmental  speed would not be averaged out as in bulk data. Indeed, the authors used  BLIND 19 -a  trajectory-learning  method-to  re-rank  their samples 27 similarly  to  RAPToR ( ρ &gt; 0.99;  Extended Data Fig. 1). RAPToR estimates do in fact enhance expression dynamics captured by principal components (Figs. 2e,f and Extended Data Fig. 1) and for the majority of genes (Extended Data Fig. 1) in comparison to chronological age (Methods). Thus, RAPToR estimates the true physiological age of individuals and reveals the heterogeneity of their developmental speeds.

Reference interpolation greatly improves staging accuracy. Crucially, reference interpolation allows staging with an accuracy far beyond the original sampling resolution of reference time series. Indeed, RAPToR accurately stages a dense zebrafish time course 28 with over 40 times the temporal resolution of the reference before interpolation  (Extended  Data  Fig.  2  and  Supplementary  Note  1; Methods).  RAPToR  estimates  also  stay  remarkably  accurate  and precise even when staging samples with a few hundred genes or with noisy data (Supplementary Note 1 and Supplementary Figs. 1-4), and are robust to reference interpolation parameters (Supplementary Note 1, Supplementary Figs. 5 and 6, and Supplementary Table 2).

## RAPToR  correctly  infers  developmental  speed  scaling  factors.

RAPToR estimates  are  relative  to  the  reference  chronological  age. Thus, one can use RAPToR to stage samples with known chronological  age  to  estimate  developmental  speed  differences  or  scaling factors with a reference. For example, staging a C. elegans developmental time series grown at 25 °C 26 on the reference grown at 20 °C 20 recapitulates the expected 1.5-fold increase in developmental speed owing to temperature increase 20 (Fig. 2a and Supplementary Note 1).

Fig. 2 | RAPToR precisely stages development and ageing, and works from whole-organism to single-cell data. a , b , Chronological age versus RAPtoR estimates of C. elegans late-larval samples 26 (linear model is y /thinspace = /thinspace -4.7/thinspace + /thinspace1.6 x ; a ) and D. rerio embryo samples 27  (linear model is y /thinspace = /thinspace0.7 x ; b ). c , Somite number versus RAPtoR estimates of M. musculus embryo samples 29  (linear model is y /thinspace = /thinspace9 .2/thinspace + /thinspace0.05 x ). d , Chronological age versus RAPtoR estimates of D. melanogaster embryo samples 27  (linear model is y /thinspace = /thinspace1.3/thinspace + /thinspace0.77 x ). e , f , Selected principal components of the data staged in d , plotted in black along chronological age ( e ) and in red along RAPtoR estimates ( f ). g -j , Chronological age versus RAPtoR estimates of adult C. elegans bulk samples 33 ( g ) and single-worms 34  ( h ), of adult D. melanogaster 35 ( i ), and of human brain tissue 36 ( j ). k , Chronological age versus RAPtoR estimates of dissected samples of upper jaw first molars from M. musculus embryos staged using the lower jaw samples as reference 37,38 . l , Inferred pseudo-time versus RAPtoR estimates of H. sapiens embryo single cells 39 . a -d , g , h , Staged samples 26,27,29,33,34 and references 20,23-25  are from independent time-series experiments. Original time points of the references within the plot area are shown to the right (blue), but the references can span much longer coverage.

RAPToR performs well on ageing. While RAPToR works very well with  robust  expression  changes  during  development,  ageing  and ageing-related changes in gene expression are widely known to be heterogeneous, stochastic, and strongly influenced by environmental factors 30-32 . This can potentially limit the applicability of RAPToR to ageing. In fact, RAPToR performs poorly between independent ageing time series (but works within experiments; Supplementary Note 1 and  Supplementary  Fig.  7)  with  references  built  using  the  whole transcriptome.  We  reasoned  RAPToR  performance  could  increase by strengthening the ageing signal in the reference. Indeed, by building RAPToR references restricted to genes with robust monotonous trends along ageing (Methods), we could successfully estimate ageing in C. elegans bulk 33 ( R ² = 0.94; Fig. 2g) and single-worm 34 samples ( R ² = 0.67,  Fig.  2h), Drosophila 35 ( R ² = 0.96;  Fig.  2i),  and  humans 36 ( R ² = 0.74; Fig.2j, Supplementary Note 1 and Supplementary Fig. 8). Importantly, single worms staged older than their chronological age behaved like older individuals 34 and vice versa (Fig. 2h), moreover age estimates of flies under caloric restriction are consistent with the expected lifespan extension 35 (Fig. 2i). This shows that RAPToR age estimates recapitulate true differences in biological age across individuals or environmental conditions. We conclude that RAPToR reliably infers ageing from transcriptomic data.

RAPToR  accurately  stages  dissected  tissue  samples. We  tested RAPToR  on  expression  profiling  from  dissected  tissues-where variation  in  cell-type  composition  and  relative  amount  might potentially  confound  staging-using  time  series  of M.  musculus upper- and lower-jaw first molar development 37,38 . Since these two organs have very similar development 37 , we built a lower-jaw reference to stage the upper-jaw samples (Methods). RAPToR not only accurately estimates age ( R 2 &gt; 0.99; Fig. 2k), but also correctly infers the known developmental delay of upper molars in comparison to lower molars 37,38 .  Thus, despite potential confounders, RAPToR is effective and precise on dissected tissue samples.

Fig. 3 | Tissue-specific staging. a , b , Selected independent components from ICA on joint C. elegans RILs 11 (dots) and reference data 21  they were staged on (orange line). c -f , As in ( a , b ), with RILs plotted in red along soma age ( c , d ), and in blue along germline age ( e , f ). g , Root mean square error between RILs and reference for independent components 2-8 when using soma, global or germline age estimates.

RAPToR  accurately  stages  single-cell  data. Single-cell  expression data is usually much sparser and noisier than bulk data, which might potentially limit RAPToR performance. We therefore tested whether RAPToR could stage single cells from human early embryo development 39 by building a reference with a random subset of the data and staging the remaining cells (Methods). Age estimates not only match the chronological time ( R ² = 0.87; Supplementary Fig. 9), but also strongly correlate with the pseudo-times computed by the authors 39 ( R ² = 0.97;  Fig.  2l).  Therefore, despite the sparsity of the data, RAPToR ranks cells as well as pseudo-time methods specifically designed for single-cell data but at the same time provides a real biological time for each cell.

RAPToR is robust to genetic variation in gene expression. Variable genetic  background  is  another  potential  confounder,  so  we  tested RAPToR performance  on  expression  data  for  over  200 C.  elegans recombinant inbred lines (RILs) showing extensive genetic variation in gene expression 11 .  RAPToR closely matches previous estimates from a  trajectory-learning  approach 13 ( R ² = 0.94;  Extended  Data  Fig.  3), thus  confirming  the  RILs  span  mid-larval  to  young-adult  stage,  a period with vast expression changes both in the soma (molting) and the germline (spermatogenesis and oogenesis) of worms.

We  noticed  that  some  gene-expression  dynamics  in  the  RILs are  advanced  and  others  delayed  in  comparison  to  the  reference (Figs. 3a,b, Extended Data Fig. 4 and Supplementary Note 1). Shifts between soma and germline development (soma-germline heterochrony)  are  easily  induced  by  environmental  and  physiological changes in C. elegans 9,40 .  Indeed, a consistent enrichment of soma and germline genes in advanced and delayed dynamics respectively suggests  soma-germline heterochrony between the reference and the RILs (Extended Data Fig. 4).

Quantifying heterochrony with tissue-specific staging. To  confirm this we used germline- and soma-specific gene sets 22,26 to separately stage the germline and soma of the RILs (Methods; Extended Data Fig. 3). We find germline- and soma-specific dynamics align better on the reference when staged with the corresponding gene set (Figs. 3d,e,g) while they are otherwise shifted (Figs. 3c,f), confirming heterochrony  between  reference  and  RILs.  Thus  tissue-specific staging outperforms global staging in case of heterochrony between the reference and the samples to stage.

Beyond differences between the reference and RILs, we noticed that tissue-specific staging also decreases variance among the RILs. Indeed, germline genes are better fit by germline than soma age and vice versa, suggesting soma-germline heterochrony among the RILs (Extended Data Fig. 5). However, when searching for the genetic basis of this heterochrony with a multivariate quantitative trait loci (QTL) analysis, we found no significant genetic locus (even at a false discovery rate (FDR) of 0.5) and overall no significant amount of genetic variance in heterochrony (Supplementary Note 1), which is therefore likely due to unknown and uncontrolled environmental differences or to a very complex genetic architecture not captured by the model.

In summary, by using tissue-specific gene sets RAPToR provides accurate tissue-specific age estimates from whole-organism expression despite varying genetic background.

Staging  on  references  of  a  different  species. Developmental time-series  data  are  often  unavailable  for  non-model  organisms. However, gene-expression dynamics during development are often well-conserved across related species, especially during the phylotypic stage 41 . Seeing the robustness of RAPToR to genetic variation within species, we decided to test how well RAPToR can stage one species on a related species.

Staging time series of embryo development across six Drosophila species 41 on a D. melanogaster reference using orthologs  indeed  results  in  accurate  age  estimates  ( R 2 &gt; 0.99;  Fig.  4a) despite  decreasing  overall  correlation  with  increasing  phylogenetic  distance  (Fig.  4b).  Moreover,  we  infer  between-species growth  speed  factors  matching  those  found  by  the  authors (Supplementary  Table  3)  and  account  for  small  age  differences between replicates  of  each  time  point,  which  refines  expression dynamics (Extended Data Fig. 6) and reduces unexplained variance in the data (Supplementary Fig. 10).

Encouraged  by  this,  we  probed  RAPToR  limits  by  staging samples on more distant reference species. We were able to stage

Fig. 4 | staging samples cross-species. a , Chronological age versus RAPtoR estimates for time series of embryo development of six Drosophila species 41 staged on a D. melanogaster reference 25 (Extended Data Fig. 6). b , Spearman correlation between samples from a and the reference at age estimate, along RAPtoR estimates. c , Chronological age versus RAPtoR estimates for single cells of M. musculus embryos 42 , staged on a H. sapiens single-cell reference 39 using orthologs (Extended Data Fig. 7). d , Chronological age versus RAPtoR estimates for C. elegans embryo samples 27 , staged on a D. melanogaster reference 25 using orthologs.

early embryo mouse single cells 42 on a human reference 39 ,  matching  both  chronological  age  ( R ² = 0.86;  Fig.  4c)  and  pseudo-times ( ρ = 0.95; Extended Data Fig. 7), as well as human 43 on cow 44 whole embryos  ( R ² = 0.83;  Supplementary  Fig.  11).  To  our  surprise, we could even successfully stage C.  elegans embryogenesis 27 on  a D. melanogaster reference ( R ² = 0.96; Fig. 4d, Extended Data Fig. 8 and Supplementary Note 1), two species separated by 600 million years of evolution 45 .

Which  biological  processes  with  highly  conserved  dynamics during embryogenesis could account for this accurate staging? We found that gene-expression signatures of decreasing cell proliferation  shared  across  phyla 27 and  signatures  of  muscle  development or cell differentiation are necessary and almost sufficient for accurate  staging  (Supplementary  Note  1,  Extended  Data  Fig.  7,8  and Supplementary Tables 4-11; Methods).

Thus, RAPToR can stage non-model organisms using available close species data and perform well even in extremely distant species, when applied to developmental stages with highly conserved developmental dynamics.

To  summarize,  RAPToR  performs  well  across  the  organisms, sample  types,  and  diverging  genetic  backgrounds  and  species  we tested, yielding estimates that are precise, accurate thanks to interpolation, and robust to gene-set size changes.

## RAPToR finds  hidden  drug  effects  on  germline  development.

RAPToR  absolute  age  estimates  are  useful  in  many  ways.  First, rather than just getting a list of differentially expressed genes from profiling  data,  RAPToR  precisely  quantifies  the  effect  of  perturbations  on  developmental  timing,  including  in  a  tissue-specific way. For example, tissue-specific staging of C.  elegans exposed to three concentrations of mefloquine, dichlorvos, and fenamiphos 46 found that all three drugs induce a similar germline-specific and dose-dependent developmental delay (Fig. 5a, Supplementary Note 2 and Supplementary Fig. 12).

RAPToR improves differential expression analyses. Even when known chronological age is included as a model covariate in differential expression analyses, replacing it by RAPToR age estimates increases statistical power. For example, using RAPToR estimates instead  of  chronological  age  when  analyzing  expression  changes in C. elegans pash-1 versus wild type (WT) 47 (Fig. 5b), detects up to  60%  more  differentially  expressed  genes  in pash-1 and  10% more  differentially  expressed  genes  across  development  thanks to  overall  better  model  fits  (Fig.  5c,  Supplementary  Fig.  13  and Supplementary Note 2).

Quantifying developmentally driven changes in gene expression. If an experimental condition strongly impacts developmental speed but perturbed and control samples are collected at the same chronological-and therefore different physiological-time (Fig. 1a), the variable  of  interest  will  be  completely  confounded  with  development. Thus, purely developmental expression changes are wrongly attributed to the perturbation of interest (Fig. 1b). As an example, we reanalyze a dataset comparing young-adult C. elegans that developed through dauer state  (post-dauer)  to  controls  that  did  not 48 . The authors found a downregulation of spermatogenesis-associated genes  and  an  upregulation  of  oogenesis-associated  genes  from which  they  concluded  that  post-dauer  animals  have  reduced spermatogenesis and increased oogenesis.  However,  as C.  elegans switch  from  sperm  to  egg  production  during  development,  this pattern  could  simply  be  explained  by  post-dauer  samples  being physiologically  older  than  controls.  This  is  indeed  what  RAPToR found (Fig. 5d, Supplementary Fig. 14 and Supplementary Note 2). Furthermore,  strong  correlation  ( r &gt; 0.8)  between  the  observed expression  changes  in  germline  genes  and  the  expected  developmental expression changes calculated from matching time points in the reference (Fig. 5e, Extended Data Fig. 9, Supplementary Fig. 14 and Supplementary Note 2) suggests that, despite synchronization efforts, most of the initially observed differential expression is due to uncontrolled differences in developmental progression.

Recovering  confounded  perturbation-specific  effects. We  reasoned  that  integrating  RAPToR  age  estimates  and  developmental gene  expression  from  the  reference  in  the  differential  expression analysis should allow us to extract perturbation-specific expression changes even when the variable of interest is completely confounded with  development  (Supplementary  Note  2  and  Extended  Data Fig. 10). We tested this using a C. elegans larval development time series  of xrn-2 mutant  and  relative  WT  control  sampled  every 1.5 h 49 . We defined a gold standard of truly differentially expressed genes  in  the  mutant,  which  allowed  us  to  vary  the  age  difference between mutant and WT and quantify first the amount, intensity, and variance of expression changes owing to development (Fig. 5f,g, Extended Data Fig. 10 and Supplementary Note 2); second, the deleterious impact of these developmental expression changes on the performance of a standard analysis in detecting truly differentially expressed genes; and third, the improvement obtained by integrating RAPToR estimates and reference expression data in the model. As expected, with increasing age differences between mutant and WT, the performance of a standard test of differential expression sharply decreases (Fig. 5h). However, performance is greatly recovered by

(hours past L4, 25 °C)

Fig. 5 | Quantifying and correcting for developmental effects using RAPToR age estimates. a , Effect of increasing drug dose exposure 46  on RAPtoR estimates of C. elegans germline age (Methods; Supplementary Fig. 12). P values are derived from two-sided t -tests on linear model coefficients for each drug. From top to bottom, respectively: mefloquine, P /thinspace = /thinspace6.32/thinspace × /thinspace10 -5 , P /thinspace = /thinspace0.0252, P /thinspace = /thinspace0.8248; dichlorvos, P /thinspace = /thinspace1.54/thinspace × /thinspace10 -4 , P /thinspace = /thinspace0.0067 , P /thinspace = /thinspace0.0620; fenamiphos, P /thinspace = /thinspace0.0040, P /thinspace = /thinspace0.0946, P /thinspace = /thinspace0.7702. n.s., not significant. b , RAPtoR estimates versus reported chronological age highlight large developmental spread within time points of Wt C. elegans and pash-1 ts time series 47  (Supplementary Note 2 and Supplementary Fig. 13). c , R 2 per gene of identical models with chronological age, or RAPtoR estimates. Genes and gene counts above and below dashed line ( x /thinspace = /thinspace y ) are indicated in red and black, respectively. d , Germline age estimates of control and post-dauer (PD) C. elegans adults 48 , P /thinspace = /thinspace0.025, derived from a two-tailed t -test. e , Germline gene logFCs between control and PD from d in comparison to logFCs expected from developmental time difference only (Extended Data Fig. 10). f , Chronological age versus RAPtoR estimates of a time series of Wt C. elegans and xrn-2 late larval development 49 . Sample subsets defining a gold standard of truly DE genes and shifted Wt sets used in subsequent panels are color-coded. g , Correlation of observed logFCs and expected developmental logFCs computed from the interpolated reference between the xrn-2 subset and increasingly shifted Wt sets from f (Supplementary Note 2). h , PR curves showing the performance of a standard differential-expression model P value for each shifted Wt subset in detecting gold-standard DE genes. i , AUPRC of standard differential expression model P value from h , or of the age-corrected classifier for each shifted Wt subset in detecting gold-standard DE genes (Supplementary Note 2). DE, differentially expressed. In a , b , d , central bar denotes group mean.

integrating reference data in the model, especially for large age differences when the mutant effect would be fully confounded by development (Fig. 5i, Extended Data Fig. 10 and Supplementary Note 2).

In summary, we showed that using RAPToR and reference data it is possible to measure the impact of development in gene-expression analyses and recover the specific effect of a perturbation even when completely confounded with development.

## discussion

We present RAPToR, a computational strategy to accurately estimate the age of samples from their gene-expression profile. Unlike trajectory-based  methods,  RAPToR  exploits  existing  reference time-series  data  to  continuously  stage  each  sample  separately, providing several advantages: first, it eliminates the need for large datasets  to  infer  developmental  trajectories;  second,  it  provides absolute developmental times that are comparable across datasets, conditions, genetic backgrounds, profiling technologies and other covariates; and third, outliers have no impact on the staging of other samples as each sample is staged independently.

While RAPToR is limited by the existence of reference time-series data, interpolation allows precise staging well beyond the resolution of the original reference data, enabling the use of sparse time series as references. RAPToR estimates age both during development or ageing in most common animal models and humans, from bulk, single-individual, dissected-tissue, or single-cell expression profiles, and can also infer tissue-specific age from whole-organism profiles. Importantly, RAPToR can stage one species using a close species as  reference,  which  dramatically  expands  the  scope  of  RAPToR, including  to  non-model  organisms.  We  showed  how  RAPToR absolute  estimates  can  be  exploited  in  many  ways:  to  detect  the effect  of  a  perturbation  or  treatment  on  developmental or ageing speed;  as  model  covariates  to  increase  statistical  power  to  detect differential expression; finally, we showed that even in the extreme scenario  when  the  perturbation  of  interest  is  completely  confounded  with  development,  it  is  still  possible  to  recover  genuine perturbation-specific  expression  changes  by  integrating  reference data in differential expression analysis.

RAPToR can currently only stage on one developmental or ageing trajectory, so a future improvement will be to provide RAPToR with the ability to stage on multiple branching trajectories. Another avenue for improvement is to adapt this approach to data other than genome-wide gene expression, such as genome-wide binding data.

We  anticipate our strategy of staging post-profiling with RAPToR  will  be  especially  useful  in  large-scale  single-organism profiling experiments since it eliminates the need for synchronization or for tedious and potentially difficult steps of accurate staging before profiling.

To conclude, we remark that our approach is not restricted to development or ageing but can in principle be applied to any process  with  robust  underlying  reference  gene-expression  dynamics

## NATuRe MeThods

(for example, cell differentiation, cell cycle, disease progression, and drug response) and its scope will only expand with the increasing availability of time-series profiling data.

## online content

Any  methods,  additional  references,  Nature  Research  reporting  summaries, source data, extended data, supplementary information,  acknowledgments,  peer  review  information;  details  of author  contributions  and  competing  interests;  and  statements  of data and code availability are available at https://doi.org/10.1038/

s41592-022-01540-0.

Received: 8 September 2021; Accepted: 25 May 2022; Published online: 11 July 2022

## ARTICLes

## Methods

Analyses were all performed using the R statistical software (v.4.1.2)

Data pre-processing. Probe or gene IDs of datasets were converted to standard IDs (WBGene IDs for C. elegans , FBgn IDs for D. melanogaster , Ensembl IDs for D. rerio , M. musculus , H. sapiens and B. taurus ). When multiple probes or IDs matched a single standard ID, they were mean-aggregated for microarray, sum-aggregated for RNA-seq counts. IDs with no standard ID match were dropped.

For RNA-seq datasets, gene-level transcripts per kilobase million (TPM) data was used when available, or computed from raw counts using transcript lengths from the Ensembl biomart (v.99). No remapping of the transcriptomes was done, aside from the M. musculus tooth data (see below). No background correction was applied to microarray data.

Samples were considered of poor quality and discarded when the 99 th percentile of the distribution of their Spearman correlation coefficients with other samples fell below a threshold defined below for each dataset.

Expression values for all datasets were quantile normalized using the normalizeBetweenArrays function from limma 50  (v.3.50.1) on log( X + 1) transformed values unless otherwise specified.

RAPToR implementation. Our method is implemented in an R package, RAPToR (v.1.1.6), which can be downloaded and installed from https://github.com/LBMC/ RAPToR.

Functions for staging samples, plotting results, interpolation, and building references are included in the package. Detailed vignettes on general usage, reference building, and showcases are also provided with the package.

Auxiliary R data packages include references for C. elegans (embryonic, larval, and young-adult to adult development, https://github.com/LBMC/wormRef), D. melanogaster (embryonic development, https://github.com/LBMC/drosoRef), D. rerio , (embryonic and larval development, https://github.com/LBMC/zebraRef) and M. musculus (embryonic development, https://github.com/LBMC/mouseRef).

Reference interpolation. Let X ( m × n ) be the gene-expression matrix of m genes by n samples. The matrix is first gene-centered such that X0 = X -rowMeans( X ). We then use independent component analysis (ICA) ('ica' function, 'icafast' library v.1.0.2) or principal component analysis (PCA) ('prcomp' base R function) to decompose the data into a component space of dimension c such that X0 = G S T , with G ( m × c ) the gene loadings, and S ( n × c ) the sample scores. Columns of S are interpolated on with respect to time (and other potential variables of interest, for example, batch), forming a new matrix T ( l × c ) of l new time points in component space. The full interpolated expression matrix Y ( m × l ) is then reconstructed by multiplying the gene loadings matrix by the transposed T and by adding the gene centers Y = G T T + rowMeans( X ).

To interpolate the components, we fit generalized additive models to handle non-linear dynamics through splines with the 'gam' function in the 'mgcv' package (v.1.8.39) using a single model formula for all components selected by cross validation (CV) as follows. CV training sets are built with 80% of samples, with proportional representation of any covariate group (for example, batch). The model is evaluated using the average relative error, mean squared error (MSE), and average root MSE 51 . We compared generalized additive models fitted with different splines (cubic, thin plate, and duchon), and chose the model with minimal CV and prediction errors. Automatic spline parameter estimation from 'gam' function was used, unless the model was clearly performing poorly with automatic parameter estimation (overfitting, predictions not matching the component dynamics), in which case we performed further CV on reasonable spline parameter spaces to tweak the model (defining a number of knots). We further verified that RAPToR age estimates match chronological age of the original reference data and of independent time series when staged on the interpolated reference, using the R ² of linear models (Supplementary Note 1).

The number of components to fit was selected by setting a cutoff on cumulative explained variance (for example, 99%). The cutoff was adjusted according to the number of components with intelligible dynamics with respect to time (Supplementary Note 1 and Supplementary Fig. 6). Interpolation (and subsequent staging) is robust to variation in the number of components used (Supplementary Note 1).

We implemented reference interpolation with the 'ge\_im' function in the 'RAPToR' package. Model formulas and parameters for building all the references used in this study are displayed in Supplementary Table 1.

Age estimation. To perform age estimation, we implemented the 'ae' function that takes the gene-expression matrix to stage (genes as rows, samples as columns), the reference matrix (genes as rows and time points as columns), and the reference times (time values associated with the columns of the reference matrix) as inputs. The 'ae' function then finds common genes between sample and reference and computes the Spearman correlation between each sample and each reference time point. The age estimate for each sample is simply the reference time point with the highest correlation.

When an age estimate lands within 5% of the reference edges, RAPToR suggests to stage the samples on another appropriate reference if possible.

Confidence intervals on age estimates are computed by repeated staging on bootstrap gene sets of default size of one third of the total. Unless stated otherwise, the number of bootstraps is 30. A confidence interval is given by the median absolute deviation (MAD) of bootstrap estimates (est boot ) from the global estimate (est), and the resolution of the interpolation (res, time interval between two points of the interpolated reference): [est -(median(|est -est boot |) + res/2); est + (median(|est -est boot |) + res/2)].

Staging using a prior probability. We implemented the possibility of providing a prior probability in the form of parameters for a Gaussian distribution per sample (mean, standard deviation) which must be given in the time scale of the reference. A Gaussian density function over the reference time is defined per sample from these parameters. During staging, all correlation peaks of the profile are determined and ranked by averaging their scaled correlation score (height of the peak in the correlation profile scaled to [0, 1]) and prior score (value of the Gaussian density function scaled to [0, 1], at the peak time point). The first peak of the ranking is then kept as the estimate. Since the ranking is determined by averaging normalized priors and correlation scores, changing the prior standard deviation parameter results in scaling the importance of the prior with respect to the correlation information.

No priors were used for staging unless explicitly stated.

Evaluating RAPToR performance. Staging C. elegans larval development . We built the reference from a time series of WT larval development at 20 °C sampled at 26 time points from L1 feeding to 48 h 20 (Supplementary Table 1), we set the number of interpolated time points to 500.

Staged samples are WT C. elegans collected during mid to late larval development at 25 °C from 22 to 37 h after L1 feeding 26 . Only samples aged below 32 h (corresponding to about 48 h at 20 °C) were staged, to stay within the reference boundaries.

Staging D. melanogaster embryonic development . We staged a Drosophila developmental time series 27 on an interpolated reference from another embryo developmental time series 25 (Dme\_embryo reference of the drosoRef package; Supplementary Table 1). Samples were discarded when the 99 th percentile of the distribution of their Spearman correlation coefficients with others samples fell below 0.6, leaving 90 samples to stage. The number of interpolated time points in the reference was set to 500.

We compared our rankings with the BLIND 19  rankings provided in the supplementary data 27 (restricting to 77 samples as the authors used a more stringent quality filter).

To test if our age estimates capture physiological development better than chronological age, we fit identical linear models using the 'lmFit' function of 'limma' with either chronological age or RAPToR estimates as the predictor. Age is modeled using a natural cubic spline with two to eight degrees of freedom (built with the ns function of the splines package). For each gene, we use R ² to compare the goodness of fit of the models with chronological age or RAPToR age estimates.

Staging D. rerio embryonic development . We used the interpolated reference we built from embryo and larval development data 23 (Dre\_emb\_larv reference of the zebraRef package; Supplementary Table 1) to stage a zebrafish time series of embryonic development from fertilization to 72 h post-fertilization 27 . Samples were discarded when the 99 th percentile of the distribution of their Spearman correlation coefficients with others samples fell below 0.6, leaving 93 samples. The number of interpolated time points in the reference was set to 1,000.

We then used the same reference, increasing the interpolation resolution between 0 and 15 h to 800 time points (resulting in a reference time density of around one time point per minute instead of the previous one time point per hour) to stage an additional dense embryonic time series of 180 zebrafish embryos around gastrula 28 . We compare RAPToR staging to rankings (Extended Data Fig. 2a) previously determined 28 as following: the ten youngest and oldest embryos (determined through the morphological criterion of epiboly coverage) are used to select the genes with the largest decrease in expression from start to end of the time series. The average expression of these genes then determines the ranking. To show the benefit of reference interpolation, we also staged the embryos on the non-interpolated reference time series (Extended Data Fig. 2c,d).

Staging M. musculus embryonic development . We used the interpolated reference we built from mouse embryonic development time series data 24  (Mmu\_embryo reference of the mouseRef package; Supplementary Table 1) to stage an independent mouse somite-staged developmental time course 29 . The number of interpolated time points was set to 500. We compare RAPToR staging with the provided embryo somite number as no chronological age is given 29 .

Staging M. musculus first-molar embryonic development . First and second data replicates for mouse first molar embryonic development are from Pantalacci et al. 37 , and Sémon et al. 38 , respectively. Reads from both replicates were processed together, trimmed with trimmomatic 52 (v.0.39) to remove adapters, and mapped using salmon 53 (v.0.14.1) and the Ensembl 98 version of the mouse transcriptome to obtain TPM values.

## NATURE METHODS

Genes with a median expression of log(TPM + 1) &lt; 0.5 across all samples were filtered out, leaving 15,362 genes. A reference was built from both replicates of the lower jaw samples (Supplementary Table 1) and used to stage all 32 samples.

Estimating developmental speed factors and resolution increase factors. Developmental speed factors and R ² between chronological and estimated age of samples are estimated with linear models.

We call 'resolution increase factor' the factor between sampling frequencies of a reference before interpolation and of a successfully staged independent time series.

C. elegans larval development is sampled every 2 h at 20 °C (0.5/h) in the reference 20 and every hour at 25 °C (1/h, 1.5 development speed factor) in the staged time series 26 resulting in a resolution increase factor rf = (1.5 × 1)/0.5 = 3.

Drosophila embryo development is sampled every 2 h (0.5/h) in the reference 25 and every 15 min (4/h) in the staged time series 27 , resulting in a resolution increase factor rf = 4/0.5 = 8.

Mouse embryo development is sampled every 1.5 d (0.66/d) in the reference 24 and somite-staged in the target time series 29 . Since the first 30 somites of M. musculus grow in ~2.5 d 54 , the somite-staged times series has a resolution of 12 time points per day (12/d) determining a resolution increase factor rf = 12/0.66 = 18.2.

Zebrafish embryo development is sampled every hour (1/h) in the reference 23 and at a rate equivalent to 47 per hour (47/h) in the staged samples 28 (180 samples are roughly evenly staged between 5.7 and 9.5 h post-fertilization: 180/(9.5 -5.7) = 47/h), resulting in a resolution increase factor rf = 47.

Building ageing references. To build RAPToR references capable of staging adults across independent studies, we select genes with monotonous expression along chronological age. Monotonous genes are defined as those with Spearman correlation with chronological age above a threshold given for each dataset. A threshold of √ 0 . 33 selects genes where approximately a third of expression variance is explained by ageing progression. When less than 200 genes were kept by the filter, we used a more lenient threshold of √ (0.25). We then interpolate as described above, using only the first component (which is monotonous given the gene selection).

Staging C. elegans ageing. We built a C. elegans ageing reference with an unpublished time series (GSE93826) using monotonous genes, defined as those with spearman correlation with chronological age above √ 0 . 33 (see Supplementary Table 1). The number of interpolated time points was set to 500.

We then staged an RNA-seq bulk time series 33 , and a single-worm microarray profiling 34 . Single-worm behavior was provided by the authors 34 .

Staging D. melanogaster ageing . Pletcher et al. 35 profiled ageing Drosophila in food ad libitum and caloric restriction conditions. We built an ageing reference with the ad libitum food time series using monotonous genes, defined as those with Spearman correlation with chronological age above √ 0 . 33 (see Supplementary Table 1). The number of interpolated time points was set to 500. We then staged all flies (control and caloric restriction).

Staging human ageing from brain tissue. We removed outliers from Chen et al. 36 expression data when the 99 th percentile of the distribution of their Spearman correlation coefficients with other samples fell below median + 2 s.d., and when the reported RNA integrity score was below seven. The remaining 360 samples were split per tissue (BA47 or BA11), and genes with Spearman correlation with chronological age above √ 0 . 25 were defined as monotonous. For each tissue, half of the samples were randomly selected to build a RAPToR reference with monotonous genes (Supplementary Table 1), on which all samples were then staged. Samples were also staged on the reference built from the other tissue (Supplementary Fig. 8).

Staging H. Sapiens single-cell embryo development . We filtered Petroupoulos et al. 39 single-cell counts to remove genes with a median expression of log(TPM + 1) &lt; 0.5 across all samples, leaving 8,482 genes. Twenty percent of cells were randomly sampled from each of the five time points (305 cells) to build an interpolated reference (Supplementary Table 1), on which all 1,529 cells were then staged (only non-reference cells are shown in Fig. 2l). Cell pseudo-times are from the original study 39 , provided as sample metadata in the ArrayExpress entry.

Probing robustness of reference interpolation. Robustness of reference interpolation to the choice of dimensionality reduction method and number of components was evaluated using either the C. elegans time series by Kim et al. 20 (as above), or the one by Meeuse et al. 21 as references.

Robustness was evaluated computing sum squared (SSQ) of gene-expression prediction error by reference models using PCA or ICA and 2-16 components with the Kim et al. time series, and 2-20 with the Meeuse et al. time series. The model formula was fixed to the one defined in Supplementary Table 1. The SSQ prediction error is defined as SSQerror = ( X n × m -X pred ) 2 n × m , with n samples, m genes.

For six conditions-ICA/PCA, each at three different numbers of components-we staged the reference samples as well as an independent

C. elegans time series 26 on the interpolated reference (only samples within reference boundaries were staged on the Kim et al. reference). We evaluated models built from 4, 9, and 14 PCA or ICA components for the Kim et al. reference and models built from 10, 20 and 25 PCA or ICA components for the Meeuse et al. reference.

We then reported the R ² value of a linear fit of RAPToR estimates by the chronological age of the samples in each condition (Supplementary Table 2), as well as correlation scores between samples and the interpolated reference at the estimate (Supplementary Fig. 5).

Estimating the impact of gene-set size on staging. The impact of gene-set size on staging was evaluated by staging the C. elegans larval time series by Hendriks et al. 26 on the reference built from the Kim et al. 20 samples, as above.

We staged the samples using 50 random gene sets of sizes 16,000, 12,000, 8,000, 4,000, 2,000, and 1,000. The resulting estimates were used to compute confidence intervals for varying bootstrap set sizes. We reported the median absolute deviation of estimates to the full gene set estimate plus interpolation resolution (that is, the size of half the confidence interval).

The same approach was repeated for smaller gene-set sizes of 2,000, 1,000, 500, and 250, this time staging the samples with and without priors (defined as 1.5 times the chronological age of the samples to account for the developmental speed difference with the reference; prior standard deviation was set to 10).

Tissue-specific staging and quantification of soma-germline heterochrony. Two hundred and eight RILs from a cross between N2 (Bristol) and CB4856 (Hawaii) strains of C. elegans were genotyped at 1,455 single-nucleotide polymorphism markers, and collected as young-adult hermaphrodites (originally intended as one time point) for profiling by microarray with one sample per RIL 11 .

Microarray intensities were first normalized within arrays with LOESS using the 'normalizeWithinArrays' function of the 'limma' library. Arrays corresponding to pooled mixed stage controls were then discarded. Samples were discarded when the 99 th percentile of the distribution of their Spearman correlation coefficients with other samples fell below 0.95, leaving 193 samples for analysis.

We staged the samples with the 'Cel\_larv\_YA ' reference 21 of the wormRef package (Supplementary Table 1), using 1,000 interpolated time points.

Samples were first staged using the entire available gene set to obtain the global estimates, then with soma and germline-specific gene sets to obtain the corresponding tissue-specific estimates: the soma gene set corresponds to the oscillatory genes denoted 'osc' in Hendriks et al. 26 . The 'germline' gene set corresponds to the union of 'germline\_intrinsic' , 'spermatogenesis\_enriched' , and 'oogenesis\_enriched' gene sets defined in Reinke et al. 26 . Estimating soma age required the use of the global estimate as prior (owing to oscillations in gene expression generating multiple correlation peaks), with the prior standard deviation set to ten for all samples. Germline age estimates required no prior.

To compare expression dynamics between reference and RILs, we kept the overlapping genes between the non-interpolated reference and the samples, quantile normalized both datasets together, and performed an ICA ('ica' function of 'icafast') extracting 46 components, explaining 95% of the variance in the joined data. For components 2-8, capturing developmental signal (IC1 captured batch effect, IC9 captured genetic variation), we defined contributing genes as those with an absolute loading above 1.96. We then tested for enrichment of soma, oogenesis, and spermatogenesis categories in these genes with a two-sided hypergeometric test, the P values of which were adjusted across all tests with the BenjaminiHochberg method.

To test heterochrony between RILs and the reference, we fit splines on the reference samples in IC2-IC8 (using the same model as the RAPToR reference; Supplementary Table 1) and computed root mean square error between the fit and RILs for each component. This was done for global, soma and germline age estimates of RILs (Fig. 3g), and for shifted values of global age estimates ( -5 h to + 5 h; Extended Data Fig. 4c).

To test the existence of heterochrony among the RILs, we fit identical models on the RIL expression data using the 'lmFit' function in limma with global, soma, or germline age values as predictors. We used natural cubic splines ('ns' function in the 'splines' library) on the age with four, six, or eight degrees of freedom. Choice between models (at equal spline degrees of freedom) was done per gene on the basis of highest R ² value.

QTL analysis on soma-germline heterochrony. The multivariate QTL analysis on soma-germline heterochrony among RILs defined as (soma age) -(germline age) was performed by random forest regression 55 with or without batch as a covariate. Each RIL was genotyped at 1,455 single-nucleotide polymorphism markers 11 . Redundant markers were filtered out from the selected 193 RILs, missing values for the remaining 1,105 markers are imputed with the 'rfImpute' function and random forest regression was fit with 5,000 trees using the 'randomForest' function; both functions are from the 'randomForest' package (v.4.6.14). The random forest selection frequency was used as importance measure, adjusted for selection bias 55 , which was estimated by fitting 500 forests of 10 trees to Gaussian noise.

We estimated the null probability distribution of random forest selection frequency through 100 trait permutations, calculated empirical P values and adjusted them for FDR.

## ARTICLES

Cross-species staging. Staging non-model Drosophila on D. melanogaster. We used the interpolated reference we built from D. melanogaster embryo development data 25 (Dme\_embryo reference of the drosoRef package; Supplementary Table 1) to stage developmental time series of six Drosophila species 41 : D. melanogaster , D. simulans , D. ananassae , D. pseudoobscura , D. permisilis and D. virilis profiled by microarrays. We used orthologs provided by the authors 41 . The number of interpolated time points in the reference was set to 500.

Developmental speed difference from D. melanogaster was determined with a linear model without intercept predicting RAPToR estimates with the chronological age of samples, with species as covariate and including interaction. A comparison with the original scaling factors 41 is shown in Supplementary Table 3.

To determine whether RAPToR estimates or the linearly scaled age from the study 41 is the better development indicator, we fit identical linear models on gene expression (lmFit function of limma) with either value as the predictor, and species as covariate. Age is modeled using a natural cubic spline with two to eight degrees of freedom (ns function of splines). For each gene, we use R ² to compare the goodness of fit of either model (Supplementary Fig. 10). No interaction between age and species coefficients was considered as temporal scaling of development between species is already applied.

We evaluated the effect of species distance on staging through the maximal correlation score between the samples and the reference (that is, at their age estimate).

Staging C. elegans on Drosophila. We staged a C. elegans embryo time series 27 on the interpolated reference we built from the D. melanogaster embryo development time series 25 ('Dme\_embryo' reference, drosoRef package). First, poor-quality C. elegans samples were discarded when the 99 th percentile of the distribution of their Spearman correlation coefficients with other samples fell below 0.67. Additionally, a sample (GSM1487346, or 'sample\_0029') was also excluded as it clearly appeared as an outlier on multiple ICA components (Extended Data Fig. 8). Four samples (GSM1487318, GSM1487319, GSM1487320, and GSM1487321, or 'sample\_0001' to '\_0004') were further removed owing to erroneous chronological age (Extended Data Fig. 8), leaving 127 samples.

We then performed the staging using a restricted worm-fly ortholog set 45 . We also did staging on a second reference interpolated as above but using the first two instead of eight components. For both references, the number of interpolated time points was set to 500.

Further analysis is restricted to the overlapping set of orthologs between worm and fly datasets (3,194 genes). We ranked genes by Spearman correlation between the C. elegans embryo time series and their matching timepoints in the second D. melanogaster reference. We then selected the 10% of genes with highest correlation (319 genes) and staged the C. elegans samples once more on the second D. melanogaster reference, evaluating staging performance with Spearman correlation and the R ² of a linear model between chronological age and estimated age.

Hierarchical clustering of the top 10% genes in the original D. melanogaster reference data 25 ('hclust' function on the euclidean distance matrix of gene-centered log(TPM + 1)), resulted in 3 clusters with over 20 genes. We then evaluated Gene Ontology enrichment in each cluster with gProfiler 56 using the 3,194 overlapping set of worm-fly orthologs as background (Supplementary Table 4-6).

Staging M. musculus single cells on H. sapiens single cells . We filtered Deng et al. 42 mouse early embryo single-cell counts to remove genes with a median expression of log(TPM + 1) &lt; 0.5 across all samples, leaving 6,506 genes. We then staged cells using all available mouse-human orthologs from ensembl (v.99), on the interpolated reference we built from H. sapiens embryo development single cells 39 (see above; Supplementary Table 1). We compared age estimates with chronological age, as well as with pseudo-time rankings similar to Petropoulos et al. 39 (Extended Data Fig. 7). We fit a principal curve (principal\_curve function of princurve library) on the first three components of a PCA on the 1,000 most variable genes.

As done above for worm-fly staging, we then restricted further analysis to the overlapping set of orthologs between mouse and human datasets (6,057 genes), ranked genes to select the 10% of genes with highest correlation between both species (509 genes), and staged the M. musculus single cells once more on the human reference, evaluating staging performance with Spearman correlation and the R ² of a linear model between chronological and estimated age.

Hierarchical clustering of the top 10% of genes in the original H. sapiens reference data (expression matrix aggregated per sampled time point, gene-centered log(TPM + 1)), resulted in 5 clusters with over 20 genes. We then evaluated Gene Ontology enrichment in each cluster with gProfiler 56 using the 6,057 overlapping set of mouse-human orthologs as background (Supplementary Tables 7-11).

Staging H. sapiens on B. taurus. We filtered an RNA-seq cow early embryo time series 44 to keep genes with median log(TPM + 1) expression &gt; 0, and built a reference with half of the samples. We similarly kept genes of a microarray human embryo time series 43 with log( X + 1) expression &gt; 2 and built a reference with half the samples (Supplementary Table 1). As only morphological stages (for example, two-cell, four-cell) and not timings were given in both datasets, we used timings from the literature 57 . The number of interpolated time points for both references was set to 100. We then staged all cow samples on the human reference and vice versa, using all available human-cow orthologs from ensembl (v.99) (Supplementary Fig. 11).

Exploiting RAPToR age estimates. Drug dose response on developmental delay in C. elegans. Expression profiles of young C. elegans adults exposed to drugs 46 were staged on the 'Cel\_larv\_YA ' reference 21  from the wormRef package (Supplementary Table 1), with 500 interpolated time points in the reference. We estimated global, soma-, and germline-specific ages (see Tissue-specific staging ). For each age type, we then subtracted the age of the control sample within each replicate of each drug assay to compute the developmental difference by treatment group. We fit a linear model with drug, dose, and interaction on the age differences to assess the significance of the effects.

Increasing statistical power in differential expression analyses. WT and pash-1 ts C. elegans samples 47 were staged on the 'Cel\_YA\_2' reference 22  from the wormRef package (Supplementary Table 1), with 500 interpolated time points in the reference. The second replicate of the first WT time point (wt\_h0.2) was omitted from further analysis owing to its extreme developmental displacement and lack of comparable mutant sample.

We fit identical linear models with the 'lmFit' function in the 'limma' library to test for differential expression, including either chronological or estimated age modeled with a natural cubic spline ('ns' function in 'splines' , degree of freedom = 2), strain, and their interaction.

Effect of strain and development was then assessed by considering the significance of appropriate model coefficients (interaction and strain coefficients for strain effect, spline and interaction coefficients for development effect), with the 'topTable' function in the 'limma' library. Differential expression was considered significant at 0.05 Benjamini-Hochberg FDR.

To test the effect of similar random age differences from chronological age, we generated 100 'random age' sets by sampling age differences from the distribution of (chronological age) -(estimated age) values, estimated with the 'density' function in R. Sampled age differences were then added to the chronological age, and the same model and analysis as above was applied. The goodness-of-fit per gene is assessed using R ².

Quantifying developmentally driven changes in gene expression. Given any two groups of expression profiling samples ' A ' and 'B' , we first stage them, then fit a linear model per gene on log2(TPM + 1) (or log 2 (intensity + 1) for microarray data) to compute the observed log 2 (fold changes) (logFCs) of ' A ' versus 'B' samples. Then we fit the same model on reference profiles at matching time points to compute logFCs expected from development only (Extended Data Fig. 9). We use squared Pearson correlation between observed and expected logFCs to quantify the variance explained by development in the observed logFC.

Control and post-dauer C. elegans samples 48 were germline-staged (see above) on the 'Cel\_larv\_YA ' reference 21 , and on the 'Cel\_YA\_2' reference 22 of the wormRef package for confirmation, as they landed near the edges of the first reference. The number of interpolated time points in the Cel\_larv\_YA and Cel\_YA\_2 references were set to 1,000 and 500 respectively. Using the method described above, we quantified the differential expression explained only by difference in developmental stages between the control and post-dauer samples.

We could not compare our results to the original results as we were unable to exactly reproduce the distribution of differential expression and P values of the original t -test analysis. We therefore recalculated differential gene expression using linear models (function 'lmFit' in 'limma' library in R).

Recovering direct perturbation effects using reference data. WT and xrn-2 time series of C. elegans late-larval development 49 were staged on the 'Cel\_larv\_YA ' reference 21 from the wormRef package (Supplementary Table 1), with 500 interpolated time points. We restricted further analysis to the genes with both at least five raw counts for at least one sample, and overlapping with the reference gene set (17,656 genes).

Defining the differential expression gold standard. To establish the gold standard of differentially expressed genes, we selected time points 8-10 of xrn-2 and WT, as they had the best (estimated) developmental match. We then calculate differential expression fitting a generalized linear model (GLM) on raw counts using the glmFit function of egdeR (v.3.36.0), including only the strain variable (model\_1), and considered genes differentially expressed with Benjamini-Hochberg adjusted P values &lt; 0.05 of a likelihood ratio test (glmLRT function of 'edgeR') on the strain coefficient.

Evaluating gold-standard gene detection decrease with age gap. To test how increasing mismatch in developmental time between xrn-2 and WT impacts differential expression analysis we apply the same GLM used for the gold standard (model\_1) to calculate differential expression between the mutant and WT samples shifted by -1, -2, -3, -5, and -7 time points and we estimated expression changes explained by development as detailed above. We then evaluated how well model\_1 P values detect gold-standard differentially expressed genes at increasing age gaps by precision-recall (PR) curves and area under PR curves (AUPRC) using the 'prediction' function of the 'ROCR' package (v.1.0.11).

## NATURE METHODS

Correcting expression changes from development. To accurately account for developmental changes we combine the samples of interest with the interpolated reference.

For each set of samples (including WT and mutant samples), we define the reference window to include as the range of age estimates widened by a 1-h margin on either side. For example, age estimates of the 'WT-1' set range 51.7-58.3 h. Thus, we include a 50.7-59.3 h window of interpolated reference.

We transform the interpolated reference data to artificial counts assuming a fixed library size of 25 × 10 6 counts per sample and a fixed number of reads 'per gene length' defined by the median of available gene lengths:

ArtificialCounts = (interpolatedTPM/(10 6 )) × ((25 × 10 6 )/ median(geneLengths)) × geneLengths

The artificial count matrix is then joined to the sample count matrix and a GLM is fit (' glmFit' in ' edgeR') , including batch (between reference and sample data), the variable of interest (strain) where reference data is grouped together with the control, and developmental time modeled with splines ('ns' function in 'splines'). To select the optimal spline degree of freedom for each window, we minimized the residual SSQ of a linear model fit on the reference window only (Extended Data Fig. 10h). Only model coefficients of the variable of interest (strain logFCs) are considered.

We first evaluated how well strain logFCs detects differentially expressed genes from the gold standard using PR curves and AUPRC ('prediction' function in 'ROCR'). We then defined an age-corrected classifier (ACC) as the weighted mean of the model\_1 P value and strain logFC of the model including the reference:

ACC = w × strainLogFC + (1 -w ) × ( -log 10 (model\_1Pval))

with w , the weight ratio of either classifier. We defined the optimal w as the value for which the AUPRC is maximal, and estimated it for each set of WT shifts. At optimal w , we then reported the AUPRC of our ACC and compared it to the standard model.

As the optimal w cannot usually be estimated in this way, we explored the relationship between optimal w and observed/expected logFC correlation (as defined above) calculated for a larger amount of WT three-sample sets (Supplementary Table 13).

Reporting summary. Further information on research design is available in the Nature Research Reporting Summary linked to this article.

## data availability

Source data for all figures is provided. Source data are provided with this paper.

## Code availability

The code to download and (pre)process the data, perform the analyses and generate the figures of this paper can be found at https://gitbio.ens-lyon.fr/LBMC/ qrg/raptor-analysis

## ARTICLES

。

。

Extended Data Fig. 1 | RAPToR estimates fit gene expression data better than chronological age. a , RAPtoR estimates of D. melanogaster single-embryo samples 27  staged on a reference built from bulk data 25 plotted against established BLIND ranks 27 . b , Percentage of genes better fitted by either RAPtoR estimates or chronological age modeled with splines using 2-8 degrees of freedom in otherwise identical models. c , R2 of models from ( b ) gene count in each half of the plot is indicated in the corners. d,e , Principal components plotted along chronological age ( d ), and RAPtoR estimates ( e ) (as in Fig. 2d-f).

。

## NAtURE MEtHODS

Extended Data Fig. 2 | Reference interpolation allows RAPToR estimates at high resolution. a , RAPtoR estimates of a zebrafish embryonic time-series from 9 spawns 28  staged on a reference built from Domazet et al. data 23 plotted against original developmental ranks 28 . b , First 2 principal components of the zebrafish time-series plotted against RAPtoR age estimates. Spawns are color-coded. c,d , RAPtoR estimates of the zebrafish time-series on the non-interpolated reference ( i.e the sampling time of the reference sample with the highest correlation) vs. original developmental ranks ( c ) and vs. standard RAPtoR estimates (as in a ) ( d ). In a,c,d , original reference time points within the plot area are shown on the right, in blue.

Extended Data Fig. 3 | Tissue-specific staging yields soma and germline ages. a , RAPtoR estimates of C. elegans Recombinant Inbred Lines (RILs) 11 staged on the larval to young-adult reference built from Meeuse et al. 21 vs. Francesconi &amp; Lehner 12 estimates. b-d , Comparison of RAPtoR estimates of global age vs. germline age ( b ), global age vs. soma age ( c ), and soma age vs. germline age ( d ). e , Distribution of soma-germline heterochrony.

Extended Data Fig. 4 | See next page for caption.

## ARTICLES

## NAtURE MEtHODS

Extended Data Fig. 4 | A delayed germline and an advanced soma. a , Independent Components from ICA on C. elegans Recombinant Inbred Lines (RILs) 11 joined to the (non-interpolated) reference data 21  plotted along chronological age and RAPtoR global estimates for the reference (orange) and RILs (black) respectively. b , Gene loadings on ICA components for all genes (n/thinspace = /thinspace14132), germline genes (oogen. n/thinspace = /thinspace582, sperm. n/thinspace = /thinspace596) and soma (n/thinspace = /thinspace2005) categories. Each box within violins spans the interquartile range (IQR), the central white dot denotes the median, and whiskers extend to 1.5 × IQR in either direction. Category enrichment p-values derive from a two-sided hypergeometric test on genes with absolute loadings above 1.96. From left to right, p-values are IC2: p/thinspace &gt; /thinspace0.99, p/thinspace &lt; /thinspace1e- 10, and p/thinspace &gt; /thinspace0.99; IC3: p/thinspace &lt; /thinspace1e- 10, p/thinspace = /thinspace2.66e-06, and p/thinspace = /thinspace0.022; IC4: p/thinspace &gt; /thinspace0.99, p/thinspace &gt; /thinspace0.99, and p/thinspace &lt; /thinspace1e- 10; IC5: p/thinspace &gt; /thinspace0.99, p/thinspace &gt; /thinspace0.99, and p/thinspace &lt; /thinspace1e- 10; IC6: p/thinspace &lt; /thinspace1e- 10, p/thinspace &gt; /thinspace0.99, and p/thinspace = /thinspace6.54e-04; IC7: p/thinspace &gt; /thinspace0.99, p/thinspace &gt; /thinspace0.99, and p/thinspace &lt; /thinspace1e- 10; IC8: p/thinspace &gt; /thinspace0.99, p/thinspace &gt; /thinspace0.99, and p/thinspace &lt; /thinspace1e- 10. c,d , Summed ( c ) and per-component ( d ) Root Mean Square Error (RMSE) between RILs and reference fit on IC2-IC8 when shifting RIL (global) age estimates. RMSE per-component shows heterochrony, with soma dynamics of RILs matching younger reference time and the reverse for germline dynamics. *: p/thinspace &lt; /thinspace0.05, **: p/thinspace &lt; /thinspace0.01, ***: p/thinspace &lt; /thinspace0.001.

Extended Data Fig. 5 | soma-germline heterochrony among C. elegans recombinant lines. Recombinant Inbred Lines (RILs) 11 are staged on the larval to young-adult reference built from Meeuse et al. samples 21 . a , Percentage of genes better fitted by either RAPtoR global, soma, or germline age estimates, modeled with splines with 4, 6, or 8 degrees of freedom in otherwise identical models. Genes are classified into spermatogenesis, oogenesis, somatic, or other (see methods). b , R2 per gene of models with global, soma, or germline age estimates as predictors for 4, 6, and 8 spline degrees of freedom.

Extended Data Fig. 6 | RAPToR age estimates synchronize expression dynamics across species. a-c , Principal components of Drosophila embryogenesis in 6 species 41 plotted along chronological age ( a ), linearly scaled chronological age 41 ( b ), and RAPtoR age estimates on a D. melanogaster reference 25 ( c ).

M.musculuschronologicalage(hpost-fertilization)

Extended Data Fig. 7 | staging M. musculus single cells on h. sapiens reference. Single cells from M. musculus embryos 42  were staged on a H. sapiens single-cell embryogenesis reference 39  using orthologs. a , First 2 principal components of a PCA done on the 1000 most variable genes. A principal curve is fit on the first 3 components. Cells are colored by RAPtoR age estimate on the H. sapiens reference. b , RAPtoR age estimates of M. musculus single cells on H. sapiens reference vs. cell ranks along principal curve (a). c , Chronological age of M. musculus single cells vs. RAPtoR age estimates on H. sapiens reference using top 10% most correlated genes between mouse and human for staging (see methods). d , H. sapiens (red) and M. musculus (black) clustered gene expression profiles (aggregated per time point) of highest-correlated genes between both species (see methods).

Extended Data Fig. 8 | staging C. elegans embryogenesis with d. melanogaster. a , C. elegans embryo samples from Levin et al. 27 staged on the D. melanogaster reference built from Graveley et al. 25  samples. Gaps appear in the estimates, likely at points where fly expression dynamics are incompatible with those of worms. b , As in ( a ), staging on the adjusted fly reference and using top 10% most correlated genes between fly and worm embryogenesis (see methods). c , D. melanogaster (red) and C. elegans (black) clustered gene expression profiles of highest-correlated genes between both species (see methods). d , ICA components of the C. elegans embryo time course plotted along sampling time. Both the red highlighted outlier and 4 samples with erroneous chronological age (circled in IC1) are omitted from analysis (see methods).

Extended Data Fig. 9 | estimating the impact of development by integrating reference data. a-c , Cartoon detailing how the log-fold-changes (logFCs) of a differential expression analysis between two sample groups ( a ) and the logFCs of their matching time points in the RAPtoR interpolated reference ( b ) can be compared to quantify the impact of development ( c ).

Extended Data Fig. 10 | Correcting the effect of development by integrating reference data. Samples from C. elegans time-course experiments of wildt-type (Wt) and xrn-2 mutants, profiled by Miki et al. 49 , and staged on the larval to young-adult reference built from Meeuse et al. samples 21 , are used to validate developmental correction approach (see also Fig. 5f-i). a , Cartoon of a model integrating a window of reference data, with Strain and Batch coefficients shown in blue. b , Number of DE genes found by a standard differential expression model (FDR/thinspace &lt; /thinspace0.05) increases with the age gaps between compared groups, with a quasi-constant fraction of truly DE genes. c , Area under PR curves (AUPRC) in detecting gold-standard DE genes for standard differential expression model p-value, age-corrected logFCs, or the age-corrected classifier for each shifted Wt subset. d , w parameter optimization for shifted Wt sets, by maximizing area under the PR curves. e , PR curves of gold-standard gene detection by the age-corrected classifier for each shifted Wt subset. f , Correlation of expected development logFCs and observed logFCs between the xrn-2 subset and combinations of 3-sample Wt sets (note these are not the 'Wt -n' subsets, see Supplementary table 13). g , Relationship between optimal w and sample-reference logFC correlation, as in ( f ). h , Optimal spline degree-of-freedom (df) selection for the different Wt shifted sets by reaching a residual Sum of Square (SSQ) plateau. the selected df increases with the shift, which is expected since the reference window to include gets larger and may thus contain more complex dynamics. DE, Differentially Expressed. logFC, log2 fold-change. FDR, false discovery rate, PR: Precision-Recall.

All code for data analysis and figure plotting is available at https://gitbio.ens-lyon.fr/LBMC/qrg/raptor-analysis

We share our tool as an open-source R package on github: https://github.com/LBMC/RAPToR

For manuscripts utilizing custom algorithms or software that are central to the research but not yet described in published literature, software must be made available to editors and reviewers. We strongly encourage code deposition in a community repository (e.g. GitHub). See the Nature Research guidelines for submitting code &amp; software for further information.

## Data

Policy information about availability of data

All manuscripts must include a data availability statement. This statement should provide the following information, where applicable:

- Accession codes, unique identifiers, or web links for publicly available datasets
- A list of figures that have associated raw data
- A description of any restrictions on data availability

The full list of datasets and accession numbers is given in Supplementary Table 12.

All the data used in this study were previously published and are available on GEO, GSE49043 (C. elegans larval development time course, 20C), GSE52861 (C. elegans late larval development time course, 25C), GSE60471 (D. melanogaster embryonic development time course), GSE60755 (C. elegans embryonic development time course), GSE60619 (D. rerio embryonic development time course), GSE83395 (D. rerio early embryonic samples), GSE39897 (M. musculus embryonic development time course), GSE76316 (M. musculus first molar embryonic development, replicate 1), GSE23857 (C. elegans late larval/young adult recombinant inbred lines), GSE696 (C. elegans young adult to adult development time course), GSE130811 (C. elegans larval to adult time course), GSE12298 (C. elegans young adult toxicity assay of 3 drugs), GSE97775 (C. elegans larval development time course of WT and xrn-2 mutant), GSE45719 (M. musculus single-cell early-embryo development), GSE93826 (C. elegans aging time course), GSE77110 (C. elegans aging time course), GSE12290 (C. elegans single-worm aging time course), GSE71620 (H. sapiens profiling of brain tissues), GSE178436 (B. taurus early-embryo development time course), and  GSE29397 (H. sapiens early-embryo development time course), available on ArrayExpress, E-MTAB-404 (Embryonic development time course of 6 Drosophila species, E-MTAB-1333 (C. elegans young adult to adult mutant vs. wild type experiment, and E-MTAB-3929 (H. sapiens single-cell early-embryo development), available on fruitfly.org (D. melanogaster embryonic development time course), in the supplementary data of the original publication (M. musculus embryonic development time course), on request to the authors and now in our code repository (C. elegans post-dauer vs. control, and D. melanogaster aging time course). The data from Sémon et al. (M. musculus first molar embryonic development, replicate 2) is, at the time of writing, in the process of acquiring a GEO accession ID. Source data for all plots is provided.

## Field-specific reporting

Please select the one below that is the best fit for your research. If you are not sure, read the appropriate sections before making your selection.

- [x] Life sciences

- [ ] Behavioural &amp; social sciences

- [ ] Ecological, evolutionary &amp; environmental sciences

For a reference copy of the document with all sections, see nature.com/documents/nr-reporting-summary-flat.pdf

## Life sciences study design

All studies must disclose on these points even when the disclosure is negative.

Sample size

Sample sizes were determined by the previously published data used for this study. Datasets were chosen to showcase the various functionalities of RAPToR.

Data exclusions

As described in methods, poor quality samples were filtered out per dataset by spearman correlation with other samples. Please refer to methods for full details.

Further reasons for omitting samples are also detailed in the methods section of relevant analyses (e.g. keeping samples within appropriate developmental range for staging, or removing a clear outlier from ICA components).

Replication

All attempts at replication of data analysis were successful. No replication of experimental data was performed since we did not collect any experimental data.

Randomization

Randomization of samples is not applicable to our study as we do not collect any experimental data.

Blinding

Blinding was not relevant to our study as we report a new software as main finding

## Reporting for specific materials, systems and methods

We require information from authors about some types of materials, experimental systems and methods used in many studies. Here, indicate whether each material, system or method listed is relevant to your study. If you are not sure if a list item applies to your research, read the appropriate section before selecting a response.

## Materials &amp; experimental systems

- [x] n/a Involved in the study

- [ ] Antibodies

- [ ] Eukaryotic cell lines

- [ ] Palaeontology and archaeology

- [ ] Animals and other organisms

- [ ] Human research participants

- [ ] Clinical data

- [ ] Dual use research of concern

## Methods

n/a

Involved in the study

- [ ] ChIP-seq

- [ ] Flow cytometry

- [ ] MRI-based neuroimaging
