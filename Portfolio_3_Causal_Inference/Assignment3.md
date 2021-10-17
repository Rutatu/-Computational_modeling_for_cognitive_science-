Assignment 3 - Causal inference
================
RF
2/5/2020

## Assignment 3 - Exploring causal inference issues

In this assignment we explore some issues related to multiple
regressions (regressions with more than one predictor), and inferred
(causal) relations between variables. N.B. the data is simulated (to
make sure I know the actual mechanism generating it), but it’s based on
a real study. So bear with a longish introduction to get into the
details of what we are doing and why it is important.

### Altercentric intrusion in schizophrenia

People with schizophrenia often report altered control and distinction
of self-other representations: intrusive thoughts, hearing of voices,
delusions of mind reading, paranoia, etc (a substantial portion of the
psychotic symptoms experienced in schizophrenia). These have been
variously attributed to hypermentalizing (over attribution of mental
states to others), social impairment (over preoccupation with own
thought processes), hyper socialization (inability to inhibit
information from others), etc.

The current study investigates 1) whether schizophrenia is indeed
related to altered control and distinction of self-other
representations, in particular altercentric intrusions (inability to
inhibit social information), and 2) whether these are related to the
relevant psychotic symptoms. N.B. the actual study also investigates
egocentric intrusion, do check the papers below if interested.

The task is a slightly modified version of this:
<https://www.ncbi.nlm.nih.gov/pubmed/20731512> You look at a picture
with some dots visible to you, as well as with a different person with a
different set of dots visible to them. The number of dots you see and
that the other sees can be the same (congruent condition) or not
(incongruent condition). You are tasked to indicate whether a given
number (e.g. 3) matches the number of dots you see (and the dots visible
to the other person are irrelevant to the task).

The tasks investigates altercentric intrusion: will your reaction time
change according to whether the other person is seeing the same amount
of dots as you, or not? The idea is that if you correctly inhibit social
information, your reaction time should not change, as the information
about the other person is not relevant. On the contrary, if you
nevertheless use task irrelevant social information, you’ll be slower at
indicating whether 3 is the right number of dots when the other person
sees a different amount of dots than you (conflicting information). The
bigger the difference between RTs in the congruent and incongruent
condition the bigger the altercentric intrusion effect.

For each participant you have 6 variables: 1) ID, 2)
AltercentricIntrusion (continuous score), 3) Diagnosis (schizophrenia
vs. control), 4) VoiceHearing (severity of voice hearing symptoms,
continuous score of the severity of the symptom as measured by a
clinician), 5) MindReading (severity of delusions of mind reading,
continuous score of the severity of the symptom as measured by a
clinician); 6) Apathy (severity of lack of motivation in taking care of
oneself, from washing to showing up at work, continuous score of the
severity of the symptom as measured by a clinician).

The research questions you have to answer are the following:

## First part

Q1.1) Does schizophrenia involved altercentric intrusion? Define model
and priors. Test the implications of your priors (prior predictive
checks) and if needed adjust them. Run the model. Test the quality of
the fitted model (posterior predictive checks). Assess the evidence in
favor of an increased altercentric intrusion in schizophrenia. Report
the model and the results, including plots.

``` r
pacman::p_load(tidyverse, brms, PerformanceAnalytics, parallel, vcov, GridExtra)
```

    ## Installing package into 'D:/Users/thram_000/Documents/R/win-library/3.6'
    ## (as 'lib' is unspecified)

    ## Warning: package 'GridExtra' is not available (for R version 3.6.1)

    ## Warning: Perhaps you meant 'gridExtra' ?

    ## Warning: unable to access index for repository http://www.stats.ox.ac.uk/pub/RWin/bin/windows/contrib/3.6:
    ##   kan ikke åbne adresse 'http://www.stats.ox.ac.uk/pub/RWin/bin/windows/contrib/3.6/PACKAGES'

    ## Warning in p_install(package, character.only = TRUE, ...):

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'GridExtra'

    ## Warning in pacman::p_load(tidyverse, brms, PerformanceAnalytics, parallel, : Failed to install/load:
    ## GridExtra

``` r
set.seed(123)

cores<- detectCores()#for parallel processing
# Prepare the data

d <- read_csv("Ass3.csv")
```

    ## Parsed with column specification:
    ## cols(
    ##   VoiceHearing = col_double(),
    ##   MindReading = col_double(),
    ##   Apathy = col_double(),
    ##   AltercentricIntrusion = col_double(),
    ##   ID = col_double(),
    ##   Diagnosis = col_double()
    ## )

``` r
summary(d)
```

    ##   VoiceHearing      MindReading          Apathy        AltercentricIntrusion
    ##  Min.   :-1.8754   Min.   :-1.4875   Min.   :-1.4747   Min.   :1.494        
    ##  1st Qu.: 0.5295   1st Qu.: 0.4483   1st Qu.: 0.5856   1st Qu.:3.322        
    ##  Median : 1.1920   Median : 1.1987   Median : 1.1123   Median :4.046        
    ##  Mean   : 1.1780   Mean   : 1.1320   Mean   : 1.1911   Mean   :3.952        
    ##  3rd Qu.: 1.8583   3rd Qu.: 1.7438   3rd Qu.: 1.8822   3rd Qu.:4.611        
    ##  Max.   : 3.6905   Max.   : 3.7404   Max.   : 3.5015   Max.   :6.312        
    ##        ID           Diagnosis   
    ##  Min.   :  1.00   Min.   :0.00  
    ##  1st Qu.: 75.75   1st Qu.:0.00  
    ##  Median :150.50   Median :0.00  
    ##  Mean   :150.50   Mean   :0.25  
    ##  3rd Qu.:225.25   3rd Qu.:0.25  
    ##  Max.   :300.00   Max.   :1.00

``` r
d$Diagnosis <- plyr::revalue(as.character(d$Diagnosis), 
                             c("0"="Controls", "1"="Schizophrenia"))

d <- d %>%
  mutate(
    ID = as.factor(ID),
    Diagnosis = as.factor(Diagnosis)
  )


summary(d)
```

    ##   VoiceHearing      MindReading          Apathy        AltercentricIntrusion
    ##  Min.   :-1.8754   Min.   :-1.4875   Min.   :-1.4747   Min.   :1.494        
    ##  1st Qu.: 0.5295   1st Qu.: 0.4483   1st Qu.: 0.5856   1st Qu.:3.322        
    ##  Median : 1.1920   Median : 1.1987   Median : 1.1123   Median :4.046        
    ##  Mean   : 1.1780   Mean   : 1.1320   Mean   : 1.1911   Mean   :3.952        
    ##  3rd Qu.: 1.8583   3rd Qu.: 1.7438   3rd Qu.: 1.8822   3rd Qu.:4.611        
    ##  Max.   : 3.6905   Max.   : 3.7404   Max.   : 3.5015   Max.   :6.312        
    ##                                                                             
    ##        ID              Diagnosis  
    ##  1      :  1   Controls     :225  
    ##  2      :  1   Schizophrenia: 75  
    ##  3      :  1                      
    ##  4      :  1                      
    ##  5      :  1                      
    ##  6      :  1                      
    ##  (Other):294

``` r
# Define the formula
# Define the formula

AltercentricDiagnosis_f0 <- bf(
  AltercentricIntrusion ~ 1 + Diagnosis
)

AltercentricDiagnosis_f <- bf(
  AltercentricIntrusion ~ 0 + Diagnosis
)


# Design the priors

get_prior(AltercentricDiagnosis_f0, family = gaussian, d)
```

    ##                 prior     class                   coef group resp dpar nlpar
    ## 1                             b                                             
    ## 2                             b DiagnosisSchizophrenia                      
    ## 3 student_t(3, 4, 10) Intercept                                             
    ## 4 student_t(3, 0, 10)     sigma                                             
    ##   bound
    ## 1      
    ## 2      
    ## 3      
    ## 4

``` r
get_prior(AltercentricDiagnosis_f, family = gaussian, d)
```

    ##                 prior class                   coef group resp dpar nlpar bound
    ## 1                         b                                                   
    ## 2                         b      DiagnosisControls                            
    ## 3                         b DiagnosisSchizophrenia                            
    ## 4 student_t(3, 0, 10) sigma

``` r
priorDiagnosis <- c(
  prior(normal(4, 1), class = b), # mean and then 2SD - inspect the data
  prior(normal(1, 2), class = sigma)
) 

# Test the priors

AltercentricDiagnosis_PriorCheck_m <- brm(
  formula = AltercentricDiagnosis_f,
  data = d,
  family = gaussian,
  prior = priorDiagnosis,
  sample_prior = "only", #meaning we want to sample the prior
  cores = cores,
  file = "AltercentricDiagnosis_PriorCheck_m"
)
#
pp_check(AltercentricDiagnosis_PriorCheck_m, nsamples = 100)+ggtitle("AI~Diagnosis PriorCheck")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
## Fitting the model
AltercentricDiagnosis_m <- brm(
  formula = AltercentricDiagnosis_f,
  data = d,
  family = gaussian,
  prior = priorDiagnosis,
  sample_prior = T,
  cores = cores,
  file = "AltercentricDiagnosis_m"
)

#plot model
plot(AltercentricDiagnosis_m)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

``` r
# Posterior predictive check
pp_check(AltercentricDiagnosis_m, nsamples = 100)+ggtitle("AI~Diagnosis Posterior predictive check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-1-3.png)<!-- -->

``` r
## Check the model for warnings
AltercentricDiagnosis_m
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 0 + Diagnosis 
    ##    Data: d (Number of observations: 300) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##                        Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## DiagnosisControls          3.86      0.06     3.74     3.98 1.00     4163
    ## DiagnosisSchizophrenia     4.23      0.10     4.02     4.42 1.00     4498
    ##                        Tail_ESS
    ## DiagnosisControls          2945
    ## DiagnosisSchizophrenia     3124
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.92      0.04     0.85     0.99 1.00     4817     3104
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
# Hypothesis testing + updating check
plot(hypothesis(AltercentricDiagnosis_m,
           "DiagnosisSchizophrenia > DiagnosisControls"))[[1]]+ggtitle("Is AI higher in Schizophrenia patients?")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-1-4.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-1-5.png)<!-- -->

``` r
hypothesis(AltercentricDiagnosis_m,
           "DiagnosisSchizophrenia > DiagnosisControls")
```

    ## Hypothesis Tests for class b:
    ##                 Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio
    ## 1 (DiagnosisSchizop... > 0     0.36      0.12     0.17     0.56        999
    ##   Post.Prob Star
    ## 1         1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(conditional_effects(AltercentricDiagnosis_m), points=T,plot=F)[[1]]+ggtitle("Estimates of AI ~ Diagnosis")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-1-6.png)<!-- -->

``` r
#evidence ratio - the posterior is a bunch of samples and we want to count how many of them are in/compatimble with our hypothesis and we count them and gives us the ratio (21 min)
```

The model indicates a credible difference in altercentric intrusion in
the two groups supporting our hypothesis (b = 0.36, CIs = 0.16, 0.57, ER
= 1332). Controls showed on average an altercentric intrusion effect of
3.86 (CIs 3.74, 3.98), and schizophrenia of 4.22 (CIs = 4.01, 4.43).
\[Add plot of the effects\]

# SI

The model had no divergences, a Rhat of 1, and Effective Sample Sizes
above 2000 for both Bulk and Tail. \[Add prior and posterior checks
plots; add updating check plot\]

Q1.2) Is altercentric intrusion related to specific symptoms *in the
patients*? Identify which of the symptoms could be relevant. Should you
include more than one symptom? Build models, priors, predictive checks.
Assess the evidence and report models and results, including plots.
Discuss whether the results make sense.

``` r
## Isolate disorder group
d <-  d %>% 
  mutate(
    AltercentricIntrusion = scale(AltercentricIntrusion),
    VoiceHearing = scale(VoiceHearing),
    MindReading = scale(MindReading),
    Apathy = scale(Apathy)
  )

schizoData <- d %>% 
  filter(Diagnosis == "Schizophrenia")


## define formula and priors ##

(schizoSummary<- psych::describe(schizoData))
```

    ##                       vars  n   mean    sd median trimmed    mad   min    max
    ## VoiceHearing             1 75   0.99  0.70   0.96    1.01   0.75 -0.46   2.43
    ## MindReading              2 75   0.69  0.91   0.63    0.71   1.00 -1.58   2.69
    ## Apathy                   3 75   0.94  0.73   0.99    0.96   0.73 -0.66   2.40
    ## AltercentricIntrusion    4 75   0.30  0.88   0.21    0.28   0.84 -1.69   2.56
    ## ID*                      5 75 154.97 93.26 160.00  155.67 131.95  3.00 298.00
    ## Diagnosis*               6 75   2.00  0.00   2.00    2.00   0.00  2.00   2.00
    ##                        range  skew kurtosis    se
    ## VoiceHearing            2.90 -0.19    -0.76  0.08
    ## MindReading             4.26 -0.17    -0.54  0.10
    ## Apathy                  3.06 -0.30    -0.71  0.08
    ## AltercentricIntrusion   4.25  0.25    -0.20  0.10
    ## ID*                   295.00 -0.05    -1.45 10.77
    ## Diagnosis*              0.00   NaN      NaN  0.00

``` r
# voice hearing
AI_VoiceHearing_f1 <- bf(
  AltercentricIntrusion ~ 1 + VoiceHearing)


# MindReading
AI_MindReading_f1 <- bf(
  AltercentricIntrusion ~ 1 + MindReading)


# Apathy
AI_Apathy_f1 <- bf(
  AltercentricIntrusion ~ 1 + Apathy)

# VH + MR
AI_VH_MR_f1 <- bf(
  AltercentricIntrusion ~ 1 + VoiceHearing + MindReading)

# ALL
AI_VH_MR_A_f1 <- bf(
  AltercentricIntrusion ~ 1 + Apathy + MindReading + VoiceHearing)


## Design the priors ##

# getting quick summary plus SD


# define prior

priorR2 <- c(
  prior(normal(0, 1), class = Intercept), # mean and then 2SD - inspect the data
  prior(normal(1, 2), class = sigma),
  prior(normal(0, .3), class = b) # sigma = the average error we expect
)

## a function for automatizaying prior check

check_prior <- function(formula, data, prior, filename){
  model <- brm(
  formula = formula,
  data = data,
  family = gaussian,
  prior = prior,
  sample_prior = "only",
  control = list( adapt_delta = 0.95),
  core = cores,
  file = filename)
  
  p_check <- pp_check(model, nsamples = 100)
  return(p_check)
}

## prior checking

checkprior_VH <- check_prior(AI_VoiceHearing_f1, schizoData, priorR2, "checkprior_VH")
checkprior_VH+ggtitle("AI~VoiceHearing Prior check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
checkprior_MR <- check_prior(AI_MindReading_f1, schizoData, priorR2, "checkprior_MR")
checkprior_MR+ggtitle("AI~MindReading Prior check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
checkprior_Apathy <- check_prior(AI_Apathy_f1, schizoData, priorR2, "checkprior_Apathy")
checkprior_Apathy+ggtitle("AI~Apathy Prior check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
checkprior_VH_MR <- check_prior(AI_VH_MR_f1, schizoData, priorR2, "checkprior_VH_MR")
checkprior_VH_MR+ggtitle("AI~VoiceHearing+MindReading Prior check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
checkprior_VC_MR_A <- check_prior(AI_VH_MR_A_f1, schizoData, priorR2, "checkprior_VC_MR_A")
checkprior_VC_MR_A+ggtitle("AI~VoiceHearing+MindReading+Apathy Prior check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

## fitting

``` r
### voice hearing
VoiceHearing1 <- brm(
    formula = AI_VoiceHearing_f1,
    data = schizoData,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "VoiceHearing1"
    
  )

#plot model

plot(VoiceHearing1)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# a plot for a posterior check
VH_posch <- pp_check(VoiceHearing1, nsamples = 100)
VH_posch+ggtitle("Ai~VoiceHearing Posterior Check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
# test and plot your hypothesis
VH_hyp <- hypothesis(VoiceHearing1, "VoiceHearing > 0")
VH_hyp
```

    ## Hypothesis Tests for class b:
    ##           Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (VoiceHearing) > 0     0.06      0.13    -0.15     0.29       2.18      0.69
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(conditional_effects(VoiceHearing1), points = T)[[1]]+ggtitle("Ai~VoiceHearing")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

``` r
## mind reading ##
MindReading1 <- brm(
    formula = AI_MindReading_f1,
    data = schizoData,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "MindReading1"
    )

MindReading1
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + MindReading 
    ##    Data: schizoData (Number of observations: 75) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept       0.24      0.13    -0.01     0.49 1.00     3936     2879
    ## MindReading     0.07      0.11    -0.15     0.29 1.00     3843     2708
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.90      0.07     0.77     1.06 1.00     3516     2897
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
#plot model

plot(MindReading1)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-5.png)<!-- -->

``` r
# a plot for a posterior check
MR_posch <- pp_check(MindReading1, nsamples = 100)
MR_posch+ggtitle("Ai~MindReading Posterior Check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-6.png)<!-- -->

``` r
# test and plot your hypothesis
MR_hyp <- hypothesis(MindReading1, "MindReading > 0")
MR_hyp
```

    ## Hypothesis Tests for class b:
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (MindReading) > 0     0.07      0.11    -0.11     0.26          3      0.75
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(MR_hyp)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-7.png)<!-- -->

``` r
plot(conditional_effects(MindReading1), points = T)[[1]]+ggtitle("Ai~MindReading")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-8.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-9.png)<!-- -->

``` r
## apathy ##
Apathy1 <- brm(
    formula = AI_Apathy_f1,
    data = schizoData,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "Apathy1"
    )

Apathy1
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + Apathy 
    ##    Data: schizoData (Number of observations: 75) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept     0.48      0.16     0.17     0.79 1.00     3904     2577
    ## Apathy       -0.19      0.12    -0.43     0.05 1.00     4086     2754
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.89      0.08     0.75     1.06 1.00     4047     3131
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
#plot model

plot(Apathy1)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-10.png)<!-- -->

``` r
# a plot for a posterior check
A_posch <- pp_check(Apathy1, nsamples = 100)
A_posch+ggtitle("Ai~Apathy Posterior Check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-11.png)<!-- -->

``` r
# test and plot your hypothesis
A_hyp <- hypothesis(Apathy1, "Apathy < 0")
A_hyp
```

    ## Hypothesis Tests for class b:
    ##     Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
    ## 1 (Apathy) < 0    -0.19      0.12    -0.39     0.01      15.13      0.94     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(A_hyp)  
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-12.png)<!-- -->

``` r
plot(conditional_effects(Apathy1), points = T)[[1]]+ggtitle("Ai~Apathy")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-13.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-3-14.png)<!-- -->

``` r
## voice hearing, mind reading

AI_VH_MR_m1 <- brm(
  formula = AI_VH_MR_f1,
  data = schizoData,
  family = gaussian,
  prior = priorR2,
  sample_prior = T,
  cores = cores,
  file = "AI_VH_MR_m1"
)
AI_VH_MR_m1
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + VoiceHearing + MindReading 
    ##    Data: schizoData (Number of observations: 75) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept        0.12      0.21    -0.29     0.55 1.00     3170     2911
    ## VoiceHearing     0.11      0.15    -0.19     0.39 1.00     3374     2888
    ## MindReading      0.10      0.11    -0.13     0.32 1.00     3387     2943
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.90      0.08     0.77     1.07 1.00     3874     3139
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
#plot model

plot(AI_VH_MR_m1)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
# Posterior predictive check
pp_check(AI_VH_MR_m1, nsamples = 100)+ggtitle("AI~VoiceHearing+MindReading Posterior Predictive Check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
# Hypothesis testing + updating check

hypothes_VHMR_VH <- hypothesis(AI_VH_MR_m1,
           "VoiceHearing > 0")

hypothes_VHMR_MR <- hypothesis(AI_VH_MR_m1,
           "MindReading > 0")

hypothes_VHMR_VH
```

    ## Hypothesis Tests for class b:
    ##           Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (VoiceHearing) > 0     0.11      0.15    -0.14     0.35       3.42      0.77
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
hypothes_VHMR_MR
```

    ## Hypothesis Tests for class b:
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (MindReading) > 0      0.1      0.11    -0.09     0.29       4.46      0.82
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(hypothes_VHMR_VH,plot=F)[[1]]+ggtitle("AI~VoiceHearing+MindReading")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

``` r
plot(hypothes_VHMR_MR,plot=F)[[1]]+ggtitle("AI~VoiceHearing+MindReading")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-4.png)<!-- -->

``` r
## conditional effects
plot(conditional_effects(AI_VH_MR_m1), points = T)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-5.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-6.png)<!-- -->

``` r
## Voice hearing, mind reading, apathy
AI_VH_MR_A_m1 <- brm(
  formula = AI_VH_MR_A_f1,
  data = schizoData,
  family = gaussian,
  prior = priorR2,
  sample_prior = T,
  cores = cores,
  file = "AI_VH_MR_A_m1"
)


AI_VH_MR_A_m1
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + Apathy + MindReading + VoiceHearing 
    ##    Data: schizoData (Number of observations: 75) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept        0.37      0.30    -0.21     0.96 1.00     3439     2955
    ## Apathy          -0.17      0.14    -0.46     0.11 1.00     3725     2838
    ## MindReading      0.04      0.13    -0.20     0.29 1.00     3553     3445
    ## VoiceHearing     0.05      0.15    -0.24     0.35 1.00     3980     3123
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.90      0.08     0.76     1.06 1.00     4185     2824
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
#plot model

plot(AI_VH_MR_A_m1)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-7.png)<!-- -->

``` r
# Posterior predictive check
pp_check(AI_VH_MR_A_m1, nsamples = 100)+ggtitle("AI~VoiceHearing+MindReading+Apathy Posterior Predictive Check")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-8.png)<!-- -->

``` r
# Hypothesis testing + updating check

hypothes3_VH <- hypothesis(AI_VH_MR_A_m1,
           "VoiceHearing > 0")

hypothes3_MR <- hypothesis(AI_VH_MR_A_m1,
           "MindReading > 0")

hypothes3_Apathy <- hypothesis(AI_VH_MR_A_m1,
           "Apathy < 0")


hypothes3_VH
```

    ## Hypothesis Tests for class b:
    ##           Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (VoiceHearing) > 0     0.05      0.15    -0.19      0.3       1.79      0.64
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
hypothes3_MR
```

    ## Hypothesis Tests for class b:
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (MindReading) > 0     0.04      0.13    -0.17     0.24        1.6      0.62
    ##   Star
    ## 1     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
hypothes3_Apathy
```

    ## Hypothesis Tests for class b:
    ##     Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
    ## 1 (Apathy) < 0    -0.17      0.14     -0.4     0.06       7.99      0.89     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(hypothes3_VH,plot=F)[[1]]+ggtitle("AI~VoiceHearing+MindReading+Apathy")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-9.png)<!-- -->

``` r
plot(hypothes3_MR,plot=F)[[1]]+ggtitle("AI~VoiceHearing+MindReading+Apathy")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-10.png)<!-- -->

``` r
plot(hypothes3_Apathy,plot=F)[[1]]+ggtitle("AI~VoiceHearing+MindReading+Apathy")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-11.png)<!-- -->

``` r
## conditional effects
plot(conditional_effects(AI_VH_MR_A_m1), points = T)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-12.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-13.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-4-14.png)<!-- -->

## Model comparison

``` r
## adding criterion loo
VoiceHearing1 <- add_criterion(VoiceHearing1, criterion = "loo")
```

    ## Automatically saving the model object in 'VoiceHearing1.rds'

``` r
MindReading1 <- add_criterion(MindReading1, criterion = "loo")
```

    ## Automatically saving the model object in 'MindReading1.rds'

``` r
Apathy1 <- add_criterion(Apathy1, criterion = "loo")
```

    ## Automatically saving the model object in 'Apathy1.rds'

``` r
AI_VH_MR_m1 <- add_criterion(AI_VH_MR_m1, criterion = "loo")
```

    ## Automatically saving the model object in 'AI_VH_MR_m1.rds'

``` r
AI_VH_MR_A_m1 <- add_criterion(AI_VH_MR_A_m1, criterion = "loo")
```

    ## Automatically saving the model object in 'AI_VH_MR_A_m1.rds'

``` r
loo_compare(VoiceHearing1,
           MindReading1,
           Apathy1,
           AI_VH_MR_m1,
           AI_VH_MR_A_m1)
```

    ##               elpd_diff se_diff
    ## Apathy1        0.0       0.0   
    ## MindReading1  -1.1       1.3   
    ## VoiceHearing1 -1.3       1.4   
    ## AI_VH_MR_A_m1 -1.6       0.4   
    ## AI_VH_MR_m1   -1.6       1.2

``` r
## adding wights
loo_model_weights(VoiceHearing1,
           MindReading1,
           Apathy1,
           AI_VH_MR_m1,
           AI_VH_MR_A_m1)
```

    ## Method: stacking
    ## ------
    ##               weight
    ## VoiceHearing1 0.000 
    ## MindReading1  0.000 
    ## Apathy1       1.000 
    ## AI_VH_MR_m1   0.000 
    ## AI_VH_MR_A_m1 0.000

## Second part

Q2.1) However, we know that the diagnosis is based on symptom
assessment: if the overall sum of symptoms is severe enough, the
participant gets a diagnosis. In other words, by selecting the patients,
and including the symptoms in the model we might have inadvertently
introduced an issue in our inference. Do try to draw a causal graph
(Directed Acyclical Graph) of the variables and compare it with the
types of causal graphs presented in the slides. Discuss which biases you
might have introduced.

Q2.2.) Redesign your analysis following the graph and report how the
results change

``` r
# formula
champion_A_f <- bf(
  AltercentricIntrusion ~ 1 + Apathy)

# prior check
checkprior_champion <- check_prior(champion_A_f, d, priorR2,"AI~Apathy_allData")
checkprior_champion+ggtitle("AI~Apathy prior check using all data")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
# model
champion_A_m <- brm(
    formula = champion_A_f,
    data = d,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "Champion_Apathy_m"
    )

champion_A_m
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + Apathy 
    ##    Data: d (Number of observations: 300) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept     0.00      0.06    -0.12     0.11 1.00     3849     3250
    ## Apathy        0.08      0.06    -0.03     0.20 1.00     4019     2656
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     1.00      0.04     0.92     1.09 1.00     3590     2518
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
#plot model

plot(champion_A_m)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
# a plot for a posterior check
champion_posch <- pp_check(champion_A_m, nsamples = 100)
champion_posch+ggtitle("AI~Apathy Posterior Pedictive Check using all data")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
# test and plot your hypothesis
champion_hyp <- hypothesis(champion_A_m, "Apathy < 0")
champion_hyp
```

    ## Hypothesis Tests for class b:
    ##     Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
    ## 1 (Apathy) < 0     0.08      0.06    -0.01     0.18       0.08      0.07     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(champion_hyp,plot = F)[[1]]+ggtitle("Using all data") 
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
plot(conditional_effects(champion_A_m), points = T,plot=F)[[1]]+ggtitle("Estimates of AI ~ Apathy using all data")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->

``` r
# formula
TheoryFormula <- bf(
  AltercentricIntrusion ~ 1 + MindReading+VoiceHearing)

# prior check
checkprior_TheoryFormula <- check_prior(TheoryFormula, d, priorR2,"AI~VH_MR_AllData")
checkprior_TheoryFormula
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-1.png)<!-- -->

``` r
# model
TheoryFormulaModel<- brm(
    formula = TheoryFormula,
    data = d,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "TheoryModel_VH_MR"
    )

TheoryFormulaModel
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + MindReading + VoiceHearing 
    ##    Data: d (Number of observations: 300) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept        0.00      0.06    -0.11     0.11 1.00     4395     3042
    ## MindReading      0.16      0.05     0.06     0.27 1.00     4194     3243
    ## VoiceHearing     0.16      0.06     0.06     0.27 1.00     4524     3325
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.97      0.04     0.90     1.06 1.00     4808     3011
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
#plot model

plot(TheoryFormulaModel)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-2.png)<!-- -->

``` r
# a plot for a posterior check
theory_posch <- pp_check(TheoryFormulaModel, nsamples = 100)
theory_posch +ggtitle("AI~MindReading+VoiceHearing Prior Check using all data")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-3.png)<!-- -->

``` r
# test and plot your hypothesis
theory_hyp_MR <- hypothesis(TheoryFormulaModel, "MindReading > 0")
theory_hyp_MR
```

    ## Hypothesis Tests for class b:
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (MindReading) > 0     0.16      0.05     0.08     0.25        799         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
theory_hyp_VH <- hypothesis(TheoryFormulaModel, "VoiceHearing > 0")
theory_hyp_VH
```

    ## Hypothesis Tests for class b:
    ##           Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (VoiceHearing) > 0     0.16      0.06     0.07     0.26     570.43         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(theory_hyp_MR,plot=F)[[1]]+ggtitle("AI~MindReading+VoiceHearing Hypothesis Test using all data") 
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-4.png)<!-- -->

``` r
plot(theory_hyp_VH,plot=F)[[1]]+ggtitle("AI~MindReading+VoiceHearing Hypothesis Test using all data")
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-5.png)<!-- -->

``` r
plot(conditional_effects(TheoryFormulaModel), points = T)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-6.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Theory%20correction-7.png)<!-- -->

``` r
## adding criterion loo
TheoryFormulaModel <-
  add_criterion(TheoryFormulaModel, criterion = "loo")
```

    ## Automatically saving the model object in 'TheoryModel_VH_MR.rds'

``` r
champion_A_m  <- add_criterion(champion_A_m, criterion = "loo")
```

    ## Automatically saving the model object in 'Champion_Apathy_m.rds'

``` r
loo_compare(TheoryFormulaModel,
            champion_A_m)
```

    ##                    elpd_diff se_diff
    ## TheoryFormulaModel  0.0       0.0   
    ## champion_A_m       -8.2       4.3

``` r
## adding wights
loo_model_weights(TheoryFormulaModel,
                  champion_A_m)
```

    ## Method: stacking
    ## ------
    ##                    weight
    ## TheoryFormulaModel 0.953 
    ## champion_A_m       0.047

``` r
### voice hearing
VoiceHearing2 <- brm(
    formula = AI_VoiceHearing_f1,
    data = d,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "VoiceHearing2"
    
  )
  VoiceHearing2
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + VoiceHearing 
    ##    Data: d (Number of observations: 300) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept        0.00      0.06    -0.11     0.10 1.00     4180     3032
    ## VoiceHearing     0.19      0.06     0.07     0.30 1.00     4154     2539
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.99      0.04     0.91     1.07 1.00     4545     3014
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
# a plot for a posterior check
VH_posch2 <- pp_check(VoiceHearing2, nsamples = 100)
VH_posch2
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-1.png)<!-- -->

``` r
# test and plot your hypothesis
VH_hyp2 <- hypothesis(VoiceHearing2, "VoiceHearing > 0")
VH_hyp2
```

    ## Hypothesis Tests for class b:
    ##           Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (VoiceHearing) > 0     0.19      0.06     0.09     0.28        999         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(VH_hyp2)  
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-2.png)<!-- -->

``` r
plot(conditional_effects(VoiceHearing2), points = T)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-3.png)<!-- -->

``` r
## mind reading

MindReading2 <- brm(
    formula = AI_MindReading_f1,
    data = d,
    family = gaussian,
    prior = priorR2,
    sample_prior = T,
    refresh = 0,
    cores = cores,
    file = "MindReading2"
    )

MindReading2
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + MindReading 
    ##    Data: d (Number of observations: 300) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept      -0.00      0.06    -0.11     0.11 1.00     4536     3000
    ## MindReading     0.19      0.05     0.08     0.29 1.00     4611     3050
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.99      0.04     0.91     1.07 1.00     4277     3158
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
# a plot for a posterior check
MR_posch2 <- pp_check(MindReading2, nsamples = 100)
MR_posch2
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-4.png)<!-- -->

``` r
# test and plot your hypothesis
MR_hyp2 <- hypothesis(MindReading2, "MindReading > 0")
MR_hyp2
```

    ## Hypothesis Tests for class b:
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (MindReading) > 0     0.19      0.05      0.1     0.28       3999         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(MR_hyp2)  
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-5.png)<!-- -->

``` r
plot(conditional_effects(MindReading2), points = T)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-6.png)<!-- -->

``` r
# VH + MR + A
AI_VH_MR_A_m2 <- brm(
  formula = AI_VH_MR_A_f1,
  data = d,
  family = gaussian,
  prior = priorR2,
  sample_prior = T,
  cores = cores,
  file = "AI_VH_MR_A_m2"
)

AI_VH_MR_A_m2
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: AltercentricIntrusion ~ 1 + Apathy + MindReading + VoiceHearing 
    ##    Data: d (Number of observations: 300) 
    ## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup samples = 4000
    ## 
    ## Population-Level Effects: 
    ##              Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept        0.00      0.05    -0.11     0.11 1.00     4903     2917
    ## Apathy           0.02      0.06    -0.10     0.13 1.00     4171     3087
    ## MindReading      0.16      0.06     0.06     0.27 1.00     5466     3458
    ## VoiceHearing     0.16      0.06     0.05     0.27 1.00     4826     3277
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     0.97      0.04     0.90     1.06 1.00     4739     3021
    ## 
    ## Samples were drawn using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
# Posterior predictive check
pp_check(AI_VH_MR_A_m2, nsamples = 100)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-7.png)<!-- -->

``` r
# Hypothesis testing + updating check

hypothes3_VH2 <- hypothesis(AI_VH_MR_A_m2,
           "VoiceHearing > 0")

hypothes3_MR2 <- hypothesis(AI_VH_MR_A_m2,
           "MindReading > 0")

hypothes3_Apathy2 <- hypothesis(AI_VH_MR_A_m2,
           "Apathy > 0")


hypothes3_VH2
```

    ## Hypothesis Tests for class b:
    ##           Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (VoiceHearing) > 0     0.16      0.06     0.07     0.26     332.33         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
hypothes3_MR2
```

    ## Hypothesis Tests for class b:
    ##          Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob
    ## 1 (MindReading) > 0     0.16      0.06     0.07     0.25        799         1
    ##   Star
    ## 1    *
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
hypothes3_Apathy2
```

    ## Hypothesis Tests for class b:
    ##     Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
    ## 1 (Apathy) > 0     0.02      0.06    -0.08     0.11       1.58      0.61     
    ## ---
    ## 'CI': 90%-CI for one-sided and 95%-CI for two-sided hypotheses.
    ## '*': For one-sided hypotheses, the posterior probability exceeds 95%;
    ## for two-sided hypotheses, the value tested against lies outside the 95%-CI.
    ## Posterior probabilities of point hypotheses assume equal prior probabilities.

``` r
plot(hypothes3_VH2)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-8.png)<!-- -->

``` r
plot(hypothes3_MR2)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-9.png)<!-- -->

``` r
plot(hypothes3_Apathy2)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-10.png)<!-- -->

``` r
## conditional effects
plot(conditional_effects(AI_VH_MR_A_m2), points = T)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-11.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-12.png)<!-- -->![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/Building%20other%20models-13.png)<!-- -->

``` r
## adding criterion loo
VoiceHearing2 <- add_criterion(VoiceHearing2, criterion = "loo")
```

    ## Automatically saving the model object in 'VoiceHearing2.rds'

``` r
MindReading2 <- add_criterion(MindReading2, criterion = "loo")
```

    ## Automatically saving the model object in 'MindReading2.rds'

``` r
champion_A_m  <- add_criterion(champion_A_m, criterion = "loo")
```

    ## Automatically saving the model object in 'Champion_Apathy_m.rds'

``` r
TheoryFormulaModel <- add_criterion(TheoryFormulaModel, criterion = "loo")
```

    ## Automatically saving the model object in 'TheoryModel_VH_MR.rds'

``` r
AI_VH_MR_A_m2 <- add_criterion(AI_VH_MR_A_m2, criterion = "loo")
```

    ## Automatically saving the model object in 'AI_VH_MR_A_m2.rds'

``` r
loo_compare(VoiceHearing2,
            MindReading2,
            champion_A_m,
            TheoryFormulaModel,
            AI_VH_MR_A_m2)
```

    ##                    elpd_diff se_diff
    ## TheoryFormulaModel  0.0       0.0   
    ## AI_VH_MR_A_m2      -0.9       0.3   
    ## MindReading2       -3.4       3.0   
    ## VoiceHearing2      -3.6       2.9   
    ## champion_A_m       -8.2       4.3

``` r
## adding wights
loo_model_weights(VoiceHearing2,
            MindReading2,
            champion_A_m,
            TheoryFormulaModel,
            AI_VH_MR_A_m2)
```

    ## Method: stacking
    ## ------
    ##                    weight
    ## VoiceHearing2      0.066 
    ## MindReading2       0.127 
    ## champion_A_m       0.000 
    ## TheoryFormulaModel 0.806 
    ## AI_VH_MR_A_m2      0.000

## Third part

These issues are very difficult to think through, and not knowing the
causal mechanisms generating the data in advance makes our inferences
even more unreliable. To explore these issues, I recommend using
simulations. In other words, defining a “true” model, generating data
from it and assessing what different analyses would lead you to infer
(and therefore which biases they might introduce). You can find the code
I used to simulate your data below.

Q3.1) Look through the code and identify whether the results you have
match the underlying truth. Discuss what you have learned.

Q3.2) OPTIONAL: is this a general pattern? Try varying the parameters
(e.g. correlation values) and assess whether the new dataset(s) leads to
the same biases in your analysis.

``` r
pacman::p_load(MASS, tidyverse, psych)

seed <- 1981 # Defining a seed so the results are always the same
n <- 300 # Defining the amount of participants

SymptomCorr <- .2 # Defining the correlation of symptoms (as they tend to co-occur)
EffectCorrRel <- .2 # Defining the correlation between relevant symptoms and effect (Some symptoms are positively correlated with the effect)
EffectCorrIrrel <- 0 # Defining the correlation between irrelevant symptoms and effect (none)

# Creating the variance-covariance matrix for the variables we want to generate (3 symptoms, 1 effect)
Sigma <- matrix(data=c(1,SymptomCorr,SymptomCorr,EffectCorrRel,
                       SymptomCorr,1,SymptomCorr,EffectCorrRel,
                       SymptomCorr,SymptomCorr,1,EffectCorrIrrel,
                       EffectCorrRel,EffectCorrRel,EffectCorrIrrel,1),
                       nrow=4,ncol=4)

## Generate data from a multivariate (mvr) normal (n) distribution
data <- mvrnorm(n = n, # number of participant
        mu = c(1.2, 1.2, 1.2, 4), # mean of each variable
        Sigma) # variance co-variance matrix

# Giving meaningful names to variables and add ID
data <- data.frame(
  VoiceHearing = data[,1], 
  MindReading =  data[,2],
  Apathy =  data[,3], 
  AltercentricIntrusion = data[,4],
  ID = seq(nrow(data)))

# Assessing whether the participant has schizophrenia (high enough sum of symptoms)
# Here we choose participants scoring above 75% percentile (the most severe ones)
data$Diagnosis <- 0
data$Diagnosis[(data$VoiceHearing + data$MindReading + data$Apathy) > 
              quantile(data$VoiceHearing + data$MindReading + data$Apathy, .75)] <-1

## Plotting the relation between variables in schizophrenia
data1 <- data %>% subset(Diagnosis==1) %>% dplyr::select(-Diagnosis, -ID)
pairs.panels(data1)
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
## Plotting the relation between variables all participants
pairs.panels(dplyr::select(data,-Diagnosis, -ID))
```

![](Assignment3_Peter-_Bella_Jakub_Bianka_Ruta_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

``` r
write_csv(data, "datAss3.csv")
```
