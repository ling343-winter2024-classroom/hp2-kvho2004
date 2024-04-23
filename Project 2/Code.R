# Clear the workspace
rm(list=ls())

# Load all data points of the four critical regions from 60 participants
rtAll<-read.csv(file.choose(),header=TRUE)
View(rtAll)

# Summarize the dataset
summary(rtAll)
str(rtAll)

# Sum of the RTs data points
length(rtAll$RT)

# Define the range of qualified data points
up<-rtAll$RT==2000 
low<-rtAll$RT==100

# Remove outliers of RT > 2000 and RT < 100
rest<-subset(rtAll,(rtAll$RT>100)&(rtAll$RT<2000))

# Save the rest of the data points in a separate file
write.csv(rest,"rest_bet100and2000.csv")

# Length of the rest data points
length(rest$RT)

# Calculate the proportion of the removed data points
len_rem=length(rtAll$RT)-length(rest$RT)
len_rem
prop_rem=len_rem/length(rtAll$RT)
prop_rem

# install package "plyr"
install.packages("plyr")

# Load "plyr"
library(plyr)

# Get the means for each verb type
means <- tapply( rest[,"RT"], rest[,c("Position","VerbType")], mean )
means



### Formula for LMEM-based non-significance intervals, copied from Politzer-Ahles (2017).
LMEMinterval <- function(
  formula,
  data,
  boot.type="percentile",
  conf=.95,
  nsim=NULL){
  
  # This convenience function calculates LMEM-based "confidence" intervals for
  #	a given design and dataset. 
  # Parameters:
  #	formula: a usual model formula, with one DV and one or more IVs. Currently
  #		this function is able to handle functions in DVs (e.g., using log(RT)
  #		rather than RT as a DV), but not in IVs. And I haven't done much testing
  #		of DVs with functions so there may be bugs; I prefer just creating a 
  #		new column in the data frame (e.g., creating a logRT column).
  #		Also note that this is currently only implemented for single-factor
  #		designs. If you have a factorial (e.g. 2x2) design, this function will
  #		collapse it into a single-factor (e.g. 1x4) design.
  #	data: a data frame with repeated measures data
  #	conf: The confidence level (between 0 and 1) for the CI. Defaults .95.
  #	boot.type: which type of bootstrap to use. Defaults to "percentile". If set 
  #		to anything else, it will instead use normal bootstrap. 
  #		Percentile bootstrap is more accurate but slower, as it requires more  
  #     iterations to get accurate.
  #	nsim: Number of bootstrap replicates to use. By default this will be 2000 if
  #		boot.type=="percentile" and 200 otherwise, but you can set `nsim` to 
  #		override that.
  
  # Load the lme4 and boot packages
  require( lme4 )
  require( boot )
  
  # Figure out the DV and the IVs.
  # This doesn't use all.var() because that strips away functions,
  #	whereas sometimes your formula might be like log(DV) ~ rather than just DV.
  vars <- rownames(attr(terms(formula),"factors"))
  
  # Figure out the DV
  DV <- vars[1]
  
  # Figure out what the random effects are. The first line finds which
  #	IVs look like random effects terms (which ones have pipes), and 
  #	the next line grabs the stuff after the pipe
  ranef_idx <- which( unlist( lapply( vars, function(x){ length( grep("|", x, fixed=T ) ) } ) )>0 )
  grouping.vars <- unlist( lapply( vars[ranef_idx],
                                   function(x){ strsplit( x, " | ", fixed=T )[[1]][2] } ) )
  
  # Figure out the fixed IVs
  IVs <- vars[-c(1,ranef_idx)]
  
  # handles cases where the DV has a function around it (e.g. when the DV
  #	is `log(RT)` rather than just `RT`
  realDV <- all.vars(formula)[1] 
  if( DV != realDV ){
    func <- gsub( paste0("(",realDV,")"), "", DV, fixed=T )
    DV <- realDV
    data[,DV] <- unlist(lapply( data[,DV], func ) )
  }
  
  ### A function to do the scaling. It first fits an intercept-only model to the 
  ### data, then subtracts the residuals and adds the intercept (the grand mean)
  LMEscale <- function( formula, data ){
    model <- lmer( formula, data )
    data$LMEscale <- as.numeric( resid(model) + fixef(model)["(Intercept)"] )
    return(data)
  }
  
  # Scale the data, using a model with only a fixed intercept
  # and random intercepts
  lmerformula <- paste( DV, " ~ 1 + ", paste( "(1|", grouping.vars, ")", collapse=" + " ) ) 
  data <- LMEscale( lmerformula, data )
  
  ### The rest of the code handles making bootstrap CIs of the scaled data. The 
  ###	general procedure is as follows: to get the bootstrap CIs we have to fit 
  ###	an lmer	model to the scaled data. To make the models more likely to 
  ###	converge, we want to fit the models without random correlation 
  ###	parameters; to do this use a hack from https://rpubs.com/Reinhold/22193,
  ###	which requires first calculating a temporary model [which may not 
  ###	converge] and then extracting dummy coefficients directly from
  ###	its model matrix to construct the good model.
  ###	Finally, we bootstrap the good model a lot of times and use the bootstrap
  ###	coefficients to get "confidence" intervals.
  
  # Collapse design into one factor (just treating each condition as its 
  #	own condition, without considering main effects, interactions, etc.)
  data$Condition <- factor( do.call( paste0, lapply( IVs, function(IV){ data[,IV] } ) ) )
  
  # Create the temporary model, which may not converge, it doesn't matter
  lmerformula <- paste( "LMEscale ~ 0 + Condition + ", paste( "(1|", grouping.vars, ")", collapse=" + " ) )
  junkmodel <- lmer( lmerformula, data )
  
  # Pull out dummy variables from model matrix https://rpubs.com/Reinhold/22193
  mydummies <- list()
  for ( i in 1:dim( model.matrix(junkmodel) )[2] ) {
    data[,paste0("c",i)] <- model.matrix(junkmodel)[,i]	
  }
  
  # Make random effect terms using the dummy variables rather than the big 
  #	'Condition' variable. Per https://rpubs.com/Reinhold/22193, this ensures 
  #	that random correlations between the random effects will not 
  #  be used.
  #	We also specify no random intercepts; because the data are scaled, every
  #	subject/item/whatever should already have a mean of 0, so the random
  #	intercept is meaningless and in fact often cannot be fit anyway.
  ran <- paste( "0 +", paste( "c", 1:dim( model.matrix(junkmodel) )[2],
                              collapse="+", sep="" ) )
  
  # Now fit the good model. Because there is no fixed-effect intercept, it will 
  #	estimate a coefficient for each condition, rather than estimating 
  #	comparisons between conditions.
  lmerformula <- paste( "LMEscale ~ 0 + Condition + ", paste( "(", ran, "||", grouping.vars, ")", collapse=" + " ) )
  model <- lmer( lmerformula, data )
  
  # A function that gets the fixed-effect estimates from a model; this will 
  #	be used for bootstrapping
  getcoef <- function(.){ getME(., "beta" ) }
  
  # Print a message so we know the function is going
  message( "Bootstrapping LME-scaled values, may be very slow..." )
  
  # Figures out the number of bootstrap replicates to use (unless the user
  #	already specified how many)
  if( is.null(nsim) ) {
    nsim <- ifelse( boot.type=="percentile", 2000, 200 )
  }
  
  # Sets some variables that will be needed depending on whether we do a 
  #	percentile or normal bootstrap
  if (boot.type=="percentile" ) {
    ci_idx <- 4:5
    ci.type <- "percent"
  } else {
    ci_idx <- 2:3
    ci.type <- "normal"
  }
  
  # Bootstrap the model. This line is what takes time
  bigboo <- bootMer( model, getcoef, nsim=nsim )
  
  # Extracts the requested CIs from the bootstrap samples
  CIs <- do.call(rbind, lapply( 1:length(fixef(model)), function(x){ boot.ci( bigboo, index=x, type=substr(ci.type, 1, 4), conf=conf )[[ci.type]][ci_idx] } ) )
  
  # Gives human-friendly row names and column names
  rownames(CIs) <- substr( names( fixef(model) ), nchar("Condition")+1, nchar(names(fixef(model))) )
  colnames(CIs) <- c("Lower", "Upper")
  
  # Returns the CIs
  return(CIs)
}

# Calculate the LMEM-based non-significance intervals
system.time( CIs <- LMEMinterval( RT ~ VerbType*Position + (1|Subject) + (1|Item), rest, boot.type="percentile" ) )

# Difference-adjust them (widen each margin of error by sqrt(2), and update each interval to be the middle of the interval +- the difference-adjusted margin of error)
daME <- (( CIs[,2] - CIs[,1] )/2) * sqrt(2)
intervals <- cbind( rowMeans(CIs)-daME, rowMeans(CIs)+daME )

# Get the intervals 
intervals

# The codes for Figure 1 (and Figure 2 mentioned in the manuscript) are presented at the end of this file. 



### Statistical analyses were conducted by performing separate linear mixed-effects models (lmer) 
###   with the lme4 package (Bates et al., 2015) on the data from each of the four critical regions:
###   Verb, NP, NP+1, NP+2. As the dependent variable, the reading times were log-transformed to 
###   yield a model with approximately normal residuals. The Verb Type was dummy coded. 
###   We started with the full structure of random effects supported by the design, which included 
###   crossed random intercepts for both participants and items as well as random slope parameters
###   for the main effects of Verb Type and NP Predictability. We then simplified the random effects 
###   structure via model comparisons to obtain the maximal fitting model for the data, 
###   using α = .2 (Matuschek et al., 2017). The main effects of Verb Type * NP Predictability, Verb Type, 
###   and NP Predictability were evaluated by means of likelihood ratio test, by comparing a model with 
###   one factor to the same model but without the factor. As for the main effect of Verb Type, 
###   pairwise comparisons were conducted by summarizing the maximal model (with the summary() function),
###   and then releveling the baseline condition. 

library(lme4)
library(lmerTest)

# The categorical variable Verb Type was dummy coded. 
contr.treatment(3)
contrasts(rest$VerbType) = contr.treatment(3)
str(rest$VerbType)

# Subset the rest of the data by positions: Veb, NP, NP+1, NP+2 
Verb<-subset(rest,Position=="Verb")
View(Verb)
str(Verb)
write.csv(Verb,"Verb.csv")
Verb$RT=log(Verb$RT)

NP<-subset(rest,Position=="NP")
View(NP)
str(NP)
write.csv(NP,"NP.csv")
NP$RT=log(NP$RT)

NP1<-subset(rest,Position=="NP_1")
View(NP1)
str(NP1)
write.csv(NP1,"NP1.csv")
NP1$RT=log(NP1$RT)

NP2<-subset(rest,Position=="NP_2")
View(NP2)
str(NP2)
write.csv(NP2,"NP2.csv")
NP2$RT=log(NP2$RT)



## Find the maximal fitting model (using α = .2 (Matuschek et al., 2017))
# Model comparisons with the data at the Verb regionn.
VerbFullMdl<- lmer(RT ~ VerbType * NPPredictability +
                     (1 + VerbType * NPPredictability |Subject) + 
                     (1 + VerbType * NPPredictability |Item), 
                   data = Verb,REML = FALSE, 
                   control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                       optCtrl=list(maxfun=1e5))) # Fail to converge
VerbMdl_inter1<- lmer(RT ~ VerbType * NPPredictability +
                        (1 + VerbType * NPPredictability |Subject) + 
                        (1 + VerbType + NPPredictability|Item), 
                      data = Verb,REML = FALSE, 
                      control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                          optCtrl=list(maxfun=1e5))) 
VerbMdl_inter2<- lmer(RT ~ VerbType * NPPredictability +
                        (1 + VerbType + NPPredictability |Subject) + 
                        (1 + VerbType * NPPredictability |Item), 
                      data = Verb,REML = FALSE, 
                      control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                          optCtrl=list(maxfun=1e5))) # Fail to converge
VerbMdl<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType + NPPredictability |Subject) + 
                 (1 + VerbType + NPPredictability |Item), 
               data = Verb,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) 
anova(VerbFullMdl,VerbMdl_inter1) 
anova(VerbFullMdl,VerbMdl_inter2) 
anova(VerbMdl_inter1,VerbMdl) 
anova(VerbMdl_inter2,VerbMdl)  

VerbMdl1<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + VerbType + NPPredictability|Subject) + 
                  (1 + VerbType |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5))) 
VerbMdl2<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + VerbType |Subject) + 
                  (1 + VerbType + NPPredictability |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbMdl3<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + VerbType |Subject) + 
                  (1 + VerbType |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbMdl4<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + VerbType + NPPredictability|Subject) + 
                  (1 + NPPredictability |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbMdl5<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + NPPredictability|Subject) + 
                  (1 + VerbType + NPPredictability |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbMdl6<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + NPPredictability |Subject) + 
                  (1 + NPPredictability |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5))) # Fail to converge
VerbMdl7<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + VerbType |Subject) + 
                  (1 + NPPredictability |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbMdl8<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + NPPredictability |Subject) + 
                  (1 + VerbType |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbMdl9<- lmer(RT ~ VerbType * NPPredictability +
                  (1 + VerbType |Subject) + 
                  (1 |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
anova(VerbMdl,VerbMdl1)  
anova(VerbMdl,VerbMdl2)  
anova(VerbMdl1,VerbMdl3)     
anova(VerbMdl2,VerbMdl3)  
anova(VerbMdl,VerbMdl4) 
anova(VerbMdl,VerbMdl5) 
anova(VerbMdl4,VerbMdl6)
anova(VerbMdl5,VerbMdl6)  
anova(VerbMdl4,VerbMdl7) 
anova(VerbMdl1,VerbMdl8) 
anova(VerbMdl3,VerbMdl9) 
anova(VerbMdl7,VerbMdl9) 

# Model comparisons with the data of NP.
NPFullMdl<- lmer(RT ~ VerbType * NPPredictability +
                   (1 + VerbType * NPPredictability|Subject) + 
                   (1 + VerbType * NPPredictability|Item), 
                 data = NP,REML = FALSE, 
                 control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                     optCtrl=list(maxfun=1e5))) # Fail to converge
NPMdl_inter1<- lmer(RT ~ VerbType * NPPredictability +
                      (1 + VerbType * NPPredictability|Subject) + 
                      (1 + VerbType + NPPredictability |Item), 
                    data = NP,REML = FALSE, 
                    control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                        optCtrl=list(maxfun=1e5)))
NPMdl_inter2<- lmer(RT ~ VerbType * NPPredictability +
                      (1 + VerbType + NPPredictability |Subject) + 
                      (1 + VerbType * NPPredictability|Item), 
                    data = NP,REML = FALSE, 
                    control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                        optCtrl=list(maxfun=1e5))) # Fail to converge
NPMdl<- lmer(RT ~ VerbType * NPPredictability +
               (1 + VerbType + NPPredictability|Subject) + 
               (1 + VerbType + NPPredictability |Item), 
             data = NP,REML = FALSE, 
             control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                 optCtrl=list(maxfun=1e5)))
anova(NPFullMdl,NPMdl_inter1)
anova(NPFullMdl,NPMdl_inter2)
anova(NPMdl_inter1,NPMdl)
anova(NPMdl_inter2,NPMdl)

NPMdl1<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType + NPPredictability|Subject) + 
                (1 + VerbType |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl2<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + VerbType + NPPredictability |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl3<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + VerbType |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl4<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType + NPPredictability|Subject) + 
                (1 + NPPredictability |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl5<- lmer(RT ~ VerbType * NPPredictability +
                (1 + NPPredictability|Subject) + 
                (1 + VerbType + NPPredictability |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl6<- lmer(RT ~ VerbType * NPPredictability +
                (1 + NPPredictability|Subject) + 
                (1 + NPPredictability |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl7<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + NPPredictability |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl8<- lmer(RT ~ VerbType * NPPredictability +
                (1 + NPPredictability |Subject) + 
                (1 + VerbType |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPMdl9<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType |Subject) + 
                (1 |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
anova(NPMdl,NPMdl1)
anova(NPMdl,NPMdl2)
anova(NPMdl1,NPMdl3)
anova(NPMdl2,NPMdl3)
anova(NPMdl,NPMdl4)
anova(NPMdl,NPMdl5)
anova(NPMdl4,NPMdl6)  
anova(NPMdl5,NPMdl6) 
anova(NPMdl4,NPMdl7) 
anova(NPMdl1,NPMdl8) 
anova(NPMdl3,NPMdl9) 
anova(NPMdl7,NPMdl9) 

# Model comparisons with the data of NP+1.
NP1FullMdl<- lmer(RT ~ VerbType * NPPredictability +
                    (1 + VerbType * NPPredictability|Subject) + 
                    (1 + VerbType * NPPredictability|Item), 
                  data = NP1,REML = FALSE, 
                  control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                      optCtrl=list(maxfun=1e5))) # Fail to converge
NP1Mdl_inter1<- lmer(RT ~ VerbType * NPPredictability +
                       (1 + VerbType * NPPredictability|Subject) + 
                       (1 + VerbType + NPPredictability |Item), 
                     data = NP1,REML = FALSE, 
                     control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                         optCtrl=list(maxfun=1e5)))
NP1Mdl_inter2<- lmer(RT ~ VerbType * NPPredictability +
                       (1 + VerbType + NPPredictability |Subject) + 
                       (1 + VerbType * NPPredictability|Item), 
                     data = NP1,REML = FALSE, 
                     control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                         optCtrl=list(maxfun=1e5)))
NP1Mdl<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType + NPPredictability |Subject) + 
                (1 + VerbType + NPPredictability |Item), 
              data = NP1,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
anova(NP1FullMdl,NP1Mdl_inter1)
anova(NP1FullMdl,NP1Mdl_inter2)
anova(NP1Mdl_inter1,NP1Mdl)
anova(NP1Mdl_inter2,NP1Mdl)

NP1Mdl1<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType + NPPredictability|Subject) + 
                 (1 + VerbType |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP1Mdl2<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + VerbType + NPPredictability |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) # Fail to converge
NP1Mdl3<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + VerbType |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP1Mdl4<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType + NPPredictability|Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP1Mdl5<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + NPPredictability|Subject) + 
                 (1 + VerbType + NPPredictability |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) # Fail to converge
NP1Mdl6<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + NPPredictability|Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP1Mdl7<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP1Mdl8<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + NPPredictability|Subject) + 
                 (1 + VerbType |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) # Fail to converge
NP1Mdl9<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
anova(NP1Mdl,NP1Mdl1)
anova(NP1Mdl,NP1Mdl2)
anova(NP1Mdl1,NP1Mdl3)
anova(NP1Mdl2,NP1Mdl3)
anova(NP1Mdl,NP1Mdl4)
anova(NP1Mdl,NP1Mdl5)
anova(NP1Mdl4,NP1Mdl6)  
anova(NP1Mdl5,NP1Mdl6)  
anova(NP1Mdl4,NP1Mdl7) 
anova(NP1Mdl1,NP1Mdl8) 
anova(NP1Mdl3,NP1Mdl9) 
anova(NP1Mdl7,NP1Mdl9) 

# Model comparisons with the data of NP+2.
NP2FullMdl<- lmer(RT ~ VerbType * NPPredictability +
                    (1 + VerbType * NPPredictability|Subject) + 
                    (1 + VerbType * NPPredictability |Item), 
                  data = NP2,REML = FALSE, 
                  control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                      optCtrl=list(maxfun=1e5))) # Fail to converge
NP2Mdl_inter1<- lmer(RT ~ VerbType * NPPredictability + 
                       (1 + VerbType * NPPredictability|Subject) + 
                       (1 + VerbType + NPPredictability |Item), 
                     data = NP2,REML = FALSE, 
                     control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                         optCtrl=list(maxfun=1e5))) # Fail to converge
NP2Mdl_inter2<- lmer(RT ~ VerbType * NPPredictability +
                       (1 + VerbType + NPPredictability |Subject) + 
                       (1 + VerbType * NPPredictability|Item), 
                     data = NP2,REML = FALSE, 
                     control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                         optCtrl=list(maxfun=1e5)))
NP2Mdl<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType + NPPredictability |Subject) + 
                (1 + VerbType + NPPredictability |Item), 
              data = NP2,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
anova(NP2FullMdl,NP2Mdl_inter1)
anova(NP2FullMdl,NP2Mdl_inter2)
anova(NP2Mdl_inter1,NP2Mdl)
anova(NP2Mdl_inter2,NP2Mdl)

NP2Mdl1<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType + NPPredictability|Subject) + 
                 (1 + VerbType |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP2Mdl2<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + VerbType + NPPredictability |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) # Fail to converge
NP2Mdl3<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + VerbType |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP2Mdl4<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType + NPPredictability|Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP2Mdl5<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + NPPredictability |Subject) + 
                 (1 + VerbType + NPPredictability |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) 
NP2Mdl6<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + NPPredictability |Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) 
NP2Mdl7<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP2Mdl8<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + NPPredictability |Subject) + 
                 (1 + VerbType |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) 
NP2Mdl9<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5))) 
anova(NP2Mdl,NP2Mdl1)
anova(NP2Mdl,NP2Mdl2)
anova(NP2Mdl1,NP2Mdl3)
anova(NP2Mdl2,NP2Mdl3)
anova(NP2Mdl,NP2Mdl4)
anova(NP2Mdl,NP2Mdl5)
anova(NP2Mdl4,NP2Mdl6)
anova(NP2Mdl5,NP2Mdl6) 
anova(NP2Mdl4,NP2Mdl7)  
anova(NP2Mdl1,NP2Mdl8) 
anova(NP2Mdl3,NP2Mdl9) 
anova(NP2Mdl7,NP2Mdl9) # α = .1789 < .2

## The maximal random structure was determined using α = .2 (Matuschek et al., 2017).
##   Based on the above model comparisons, the maximal model fitting the data is
##   RT ~ VT * NPP + (1 + VT|Subject) + (1 + NPP|Item).
##   Following analyses are to evaluate the main effects of Verb Type * NP Predictability,
##   Verb Type, and NP Predictability by means of likelihood ratio test,
##   i.e., comparing a model with one factor to the same model but without the factor.

# Evaluate the main effects with the data for Verb
VerbMax<- lmer(RT ~ VerbType * NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + NPPredictability |Item), 
               data = Verb,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
VerbVtNp<- lmer(RT ~ VerbType + NPPredictability +
                  (1 + VerbType |Subject) + 
                  (1 + NPPredictability |Item), 
                data = Verb,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
VerbNp<- lmer(RT ~ NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + NPPredictability |Item), 
              data = Verb,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
VerbVt<- lmer(RT ~ VerbType +
                (1 + VerbType |Subject) + 
                (1 + NPPredictability |Item), 
              data = Verb,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
anova(VerbMax,VerbVtNp)
anova(VerbVtNp,VerbNp)
anova(VerbVtNp,VerbVt)

# Evaluate the main effects with the data for NP
NPMax<- lmer(RT ~ VerbType * NPPredictability +
               (1 + VerbType |Subject) + 
               (1 + NPPredictability |Item), 
             data = NP,REML = FALSE, 
             control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                 optCtrl=list(maxfun=1e5)))
NPVtNp<- lmer(RT ~ VerbType + NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + NPPredictability |Item), 
              data = NP,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NPNp<- lmer(RT ~ NPPredictability +
              (1 + VerbType |Subject) + 
              (1 + NPPredictability |Item), 
            data = NP,REML = FALSE, 
            control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                optCtrl=list(maxfun=1e5)))
NPVt<- lmer(RT ~ VerbType +
              (1 + VerbType |Subject) + 
              (1 + NPPredictability |Item), 
            data = NP,REML = FALSE, 
            control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                optCtrl=list(maxfun=1e5)))
anova(NPMax,NPVtNp)
anova(NPVtNp,NPNp)
anova(NPVtNp,NPVt)

# Evaluate the main effects with the data for NP1
NP1Max<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + NPPredictability |Item), 
              data = NP1,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NP1Max_VT<- lmer(RT ~ NPPredictability +
                   (1 + VerbType |Subject) + 
                   (1 + NPPredictability |Item), 
                 data = NP1,REML = FALSE, 
                 control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                     optCtrl=list(maxfun=1e5)))
NP1VtNp<- lmer(RT ~ VerbType + NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP1,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP1Np<- lmer(RT ~ NPPredictability +
               (1 + VerbType |Subject) + 
               (1 + NPPredictability |Item), 
             data = NP1,REML = FALSE, 
             control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                 optCtrl=list(maxfun=1e5)))
NP1Vt<- lmer(RT ~ VerbType +
               (1 + VerbType |Subject) + 
               (1 + NPPredictability |Item), 
             data = NP1,REML = FALSE, 
             control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                 optCtrl=list(maxfun=1e5)))
anova(NP1Max,NP1VtNp)
anova(NP1VtNp,NP1Np)
anova(NP1VtNp,NP1Vt)

# Evaluate the main effects with the data for NP2
NP2Max<- lmer(RT ~ VerbType * NPPredictability +
                (1 + VerbType |Subject) + 
                (1 + NPPredictability |Item), 
              data = NP2,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
NP2VtNp<- lmer(RT ~ VerbType + NPPredictability +
                 (1 + VerbType |Subject) + 
                 (1 + NPPredictability |Item), 
               data = NP2,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
NP2Np<- lmer(RT ~ NPPredictability +
               (1 + VerbType |Subject) + 
               (1 + NPPredictability |Item), 
             data = NP2,REML = FALSE, 
             control=lmerControl(optimizer = "nloptwrap",
                                 optCtrl = list(algorithm="NLOPT_LN_COBYLA",
                                                xtol_abs=1e-5,ftol_abs=1e-5)))
NP2Vt<- lmer(RT ~ VerbType +
               (1 + VerbType |Subject) + 
               (1 + NPPredictability |Item), 
             data = NP2,REML = FALSE, 
             control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                 optCtrl=list(maxfun=1e5)))
anova(NP2Max,NP2VtNp)
anova(NP2VtNp,NP2Np)
anova(NP2VtNp,NP2Vt)



## Examine the main effect of Verb Type
# Pairwise comparisons: Verb
summary(VerbMax)
VerbMax_rel<- lmer(RT ~ relevel(VerbType,ref="PrefV") * NPPredictability +
                     (1 + VerbType |Subject) + 
                     (1 + NPPredictability |Item), 
                   data = Verb,REML = FALSE, 
                   control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                       optCtrl=list(maxfun=1e5)))
summary(VerbMax_rel)

# Pairwise comparisons: NP
summary(NPMax)
NPMax_rel<- lmer(RT ~ relevel(VerbType,ref="PrefV") * NPPredictability +
                   (1 + VerbType |Subject) + 
                   (1 + NPPredictability |Item), 
                 data = NP,REML = FALSE, 
                 control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                     optCtrl=list(maxfun=1e5)))
summary(NPMax_rel)

# Pairwise comparisons: NP+1
summary(NP1Max)
NP1Max_rel<- lmer(RT ~ relevel(VerbType,ref="PrefV") * NPPredictability +
                    (1 + VerbType |Subject) + 
                    (1 + NPPredictability |Item), 
                  data = NP1,REML = FALSE, 
                  control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                      optCtrl=list(maxfun=1e5)))
summary(NP1Max_rel)

# Pairwise comparisons: NP+2
summary(NP2Max)
NP2Max_rel<- lmer(RT ~ relevel(VerbType,ref="PrefV") * NPPredictability +
                    (1 + VerbType |Subject) + 
                    (1 + NPPredictability |Item), 
                  data = NP2,REML = FALSE, 
                  control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                      optCtrl=list(maxfun=1e5)))
summary(NP2Max_rel)



### Since the interaction between Verb Type and NP Predictability (see Table 4) has also been found,
###   we conducted Additional statistical analyses to examine how the NP Predictability 
###   interacted with Verb Type to affect the reading times. We examined in which verb type
###   the NP predictability may affect the reading times. All the data were subset into three
###   groups in terms of Verb Type. The separate linear mixed-effects analyses were performed 
###   with these groups, with the log-transformed reading times as the dependent variable, 
###   NP Predictability as a fixed factor, and by-participants and by-items random slope for 
###   the main effects of NP Predictability. The statistical analyses were conducted on the data 
###   from each of the two critical regions (i.e., NP+1, NP+2) where the interaction effect of 
###   Verb Type * NP Predictability was found. 

# Subset the data by Verb Type
AspData <- subset(rest,VerbType=="AspV")
NprefData <- subset(rest,VerbType=="NonPrefV")
PrefData <- subset(rest,VerbType=="PrefV")

write.csv(AspData,"AspVData.csv")
write.csv(NprefData,"NonprefData.csv")
write.csv(PrefData,"PrefData.csv")

# Subset the AspData by positions: NP+1, NP+2
NP1AspData<-subset(AspData ,Position=="NP_1")
NP2AspData<-subset(AspData,Position=="NP_2")

# Subset the NprefData by positions: NP+1, NP+2
NP1NprefData<-subset(NprefData ,Position=="NP_1")
NP2NprefData<-subset(NprefData,Position=="NP_2")

# Subset the PrefData by positions: NP+1, NP+2
NP1PrefData<-subset(PrefData ,Position=="NP_1")
NP2PrefData<-subset(PrefData,Position=="NP_2")



## Examine the main effect of NP Predictability on reading times with the AspData 
#  NP+1 
NP1Asp<- lmer(RT ~ NPPredictability  +
                (1 + NPPredictability|Subject) + 
                (1 + NPPredictability|Item), 
              data = NP1AspData,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
summary(NP1Asp)

# NP+2
NP2Asp<- lmer(RT ~ NPPredictability  +
                (1 + NPPredictability|Subject) + 
                (1 + NPPredictability|Item), 
              data = NP2AspData,REML = FALSE, 
              control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                  optCtrl=list(maxfun=1e5)))
summary(NP2Asp)



## Examine the main effect of NP Predictability on reading times with the PrefData 
# NP+1
NP1Pref<- lmer(RT ~ NPPredictability  +
                 (1 + NPPredictability|Subject) + 
                 (1 + NPPredictability|Item), 
               data = NP1PrefData,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
summary(NP1Pref)

# NP+2
NP2Pref<- lmer(RT ~ NPPredictability  +
                 (1 + NPPredictability|Subject) + 
                 (1 + NPPredictability|Item), 
               data = NP2PrefData,REML = FALSE, 
               control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                   optCtrl=list(maxfun=1e5)))
summary(NP2Pref)



## Examine the main effect of NP Predictability on reading times with the NprefData 
# NP+1
NP1Npref<- lmer(RT ~ NPPredictability  +
                  (1 + NPPredictability|Subject) + 
                  (1 + NPPredictability|Item), 
                data = NP1NprefData,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
summary(NP1Npref)

# NP+2
NP2Npref<- lmer(RT ~ NPPredictability  +
                  (1 + NPPredictability|Subject) + 
                  (1 + NPPredictability|Item), 
                data = NP2NprefData,REML = FALSE, 
                control=lmerControl(optimizer="bobyqa",calc.derivs = FALSE,
                                    optCtrl=list(maxfun=1e5)))
summary(NP2Npref)



###----Figure 1: Reading Times by Verb Type and Position----
# Load the data in which outliers (>2000ms and <100ms) have been removed.
rest<-read.csv(file.choose(),header=TRUE)
View(rest)

### Formula for LMEM-based non-significance intervals, copied from Politzer-Ahles (2017).
LMEMinterval <- function(
  formula,
  data,
  boot.type="percentile",
  conf=.95,
  nsim=NULL){
  
  # This convenience function calculates LMEM-based "confidence" intervals for
  #	a given design and dataset. 
  # Parameters:
  #	formula: a usual model formula, with one DV and one or more IVs. Currently
  #		this function is able to handle functions in DVs (e.g., using log(RT)
  #		rather than RT as a DV), but not in IVs. And I haven't done much testing
  #		of DVs with functions so there may be bugs; I prefer just creating a 
  #		new column in the data frame (e.g., creating a logRT column).
  #		Also note that this is currently only implemented for single-factor
  #		designs. If you have a factorial (e.g. 2x2) design, this function will
  #		collapse it into a single-factor (e.g. 1x4) design.
  #	data: a data frame with repeated measures data
  #	conf: The confidence level (between 0 and 1) for the CI. Defaults .95.
  #	boot.type: which type of bootstrap to use. Defaults to "percentile". If set 
  #		to anything else, it will instead use normal bootstrap. 
  #		Percentile bootstrap is more accurate but slower, as it requires more  
  #     iterations to get accurate.
  #	nsim: Number of bootstrap replicates to use. By default this will be 2000 if
  #		boot.type=="percentile" and 200 otherwise, but you can set `nsim` to 
  #		override that.
  
  # Load the lme4 and boot packages
  require( lme4 )
  require( boot )
  
  # Figure out the DV and the IVs.
  # This doesn't use all.var() because that strips away functions,
  #	whereas sometimes your formula might be like log(DV) ~ rather than just DV.
  vars <- rownames(attr(terms(formula),"factors"))
  
  # Figure out the DV
  DV <- vars[1]
  
  # Figure out what the random effects are. The first line finds which
  #	IVs look like random effects terms (which ones have pipes), and 
  #	the next line grabs the stuff after the pipe
  ranef_idx <- which( unlist( lapply( vars, function(x){ length( grep("|", x, fixed=T ) ) } ) )>0 )
  grouping.vars <- unlist( lapply( vars[ranef_idx],
                                   function(x){ strsplit( x, " | ", fixed=T )[[1]][2] } ) )
  
  # Figure out the fixed IVs
  IVs <- vars[-c(1,ranef_idx)]
  
  # handles cases where the DV has a function around it (e.g. when the DV
  #	is `log(RT)` rather than just `RT`
  realDV <- all.vars(formula)[1] 
  if( DV != realDV ){
    func <- gsub( paste0("(",realDV,")"), "", DV, fixed=T )
    DV <- realDV
    data[,DV] <- unlist(lapply( data[,DV], func ) )
  }
  
  ### A function to do the scaling. It first fits an intercept-only model to the 
  ### data, then subtracts the residuals and adds the intercept (the grand mean)
  LMEscale <- function( formula, data ){
    model <- lmer( formula, data )
    data$LMEscale <- as.numeric( resid(model) + fixef(model)["(Intercept)"] )
    return(data)
  }
  
  # Scale the data, using a model with only a fixed intercept
  # and random intercepts
  lmerformula <- paste( DV, " ~ 1 + ", paste( "(1|", grouping.vars, ")", collapse=" + " ) ) 
  data <- LMEscale( lmerformula, data )
  
  ### The rest of the code handles making bootstrap CIs of the scaled data. The 
  ###	general procedure is as follows: to get the bootstrap CIs we have to fit 
  ###	an lmer	model to the scaled data. To make the models more likely to 
  ###	converge, we want to fit the models without random correlation 
  ###	parameters; to do this use a hack from https://rpubs.com/Reinhold/22193,
  ###	which requires first calculating a temporary model [which may not 
  ###	converge] and then extracting dummy coefficients directly from
  ###	its model matrix to construct the good model.
  ###	Finally, we bootstrap the good model a lot of times and use the bootstrap
  ###	coefficients to get "confidence" intervals.
  
  # Collapse design into one factor (just treating each condition as its 
  #	own condition, without considering main effects, interactions, etc.)
  data$Condition <- factor( do.call( paste0, lapply( IVs, function(IV){ data[,IV] } ) ) )
  
  # Create the temporary model, which may not converge, it doesn't matter
  lmerformula <- paste( "LMEscale ~ 0 + Condition + ", paste( "(1|", grouping.vars, ")", collapse=" + " ) )
  junkmodel <- lmer( lmerformula, data )
  
  # Pull out dummy variables from model matrix https://rpubs.com/Reinhold/22193
  mydummies <- list()
  for ( i in 1:dim( model.matrix(junkmodel) )[2] ) {
    data[,paste0("c",i)] <- model.matrix(junkmodel)[,i]	
  }
  
  # Make random effect terms using the dummy variables rather than the big 
  #	'Condition' variable. Per https://rpubs.com/Reinhold/22193, this ensures 
  #	that random correlations between the random effects will not 
  #  be used.
  #	We also specify no random intercepts; because the data are scaled, every
  #	subject/item/whatever should already have a mean of 0, so the random
  #	intercept is meaningless and in fact often cannot be fit anyway.
  ran <- paste( "0 +", paste( "c", 1:dim( model.matrix(junkmodel) )[2],
                              collapse="+", sep="" ) )
  
  # Now fit the good model. Because there is no fixed-effect intercept, it will 
  #	estimate a coefficient for each condition, rather than estimating 
  #	comparisons between conditions.
  lmerformula <- paste( "LMEscale ~ 0 + Condition + ", paste( "(", ran, "||", grouping.vars, ")", collapse=" + " ) )
  model <- lmer( lmerformula, data )
  
  # A function that gets the fixed-effect estimates from a model; this will 
  #	be used for bootstrapping
  getcoef <- function(.){ getME(., "beta" ) }
  
  # Print a message so we know the function is going
  message( "Bootstrapping LME-scaled values, may be very slow..." )
  
  # Figures out the number of bootstrap replicates to use (unless the user
  #	already specified how many)
  if( is.null(nsim) ) {
    nsim <- ifelse( boot.type=="percentile", 2000, 200 )
  }
  
  # Sets some variables that will be needed depending on whether we do a 
  #	percentile or normal bootstrap
  if (boot.type=="percentile" ) {
    ci_idx <- 4:5
    ci.type <- "percent"
  } else {
    ci_idx <- 2:3
    ci.type <- "normal"
  }
  
  # Bootstrap the model. This line is what takes time
  bigboo <- bootMer( model, getcoef, nsim=nsim )
  
  # Extracts the requested CIs from the bootstrap samples
  CIs <- do.call(rbind, lapply( 1:length(fixef(model)), function(x){ boot.ci( bigboo, index=x, type=substr(ci.type, 1, 4), conf=conf )[[ci.type]][ci_idx] } ) )
  
  # Gives human-friendly row names and column names
  rownames(CIs) <- substr( names( fixef(model) ), nchar("Condition")+1, nchar(names(fixef(model))) )
  colnames(CIs) <- c("Lower", "Upper")
  
  # Returns the CIs
  return(CIs)
}

# Calculate the LMEM-based non-significance intervals
system.time( CIs <- LMEMinterval( RT ~ VerbType*Position + (1|Subject) + (1|Item), rest, boot.type="percentile" ) )

# Difference-adjust them (widen each margin of error by sqrt(2), and update each interval to be the middle of the interval +- the difference-adjusted margin of error)
daME <- (( CIs[,2] - CIs[,1] )/2) * sqrt(2)
intervals <- cbind( rowMeans(CIs)-daME, rowMeans(CIs)+daME )

# Get the means for each verb type
means <- tapply( rest[,"RT"], rest[,c("Position","VerbType")], mean )

# Adjust the margins
par(mar=c(5.1,5.1,2.1,2.1), mgp=c(3,1,0))

# Set up the plot area
plot(1,1,type="n", xlim=c(-1.7,1.7), ylim=range(c(380,500)), ylab="Reading time (ms)", xlab="Region", xaxt="n", main="Mean reading times by verb type and region", cex.main=1.5, cex.lab=1.3, cex.axis=1.4)

# x-axis
axis(1, at=c(-1.5, -.5, .5, 1.5), labels=c("Verb","NP","NP+1","NP+2"), cex.axis=1.2)

# For each verb type, plot mean reading time at each region
lines( c(-1.55,-.55,.45,1.45), means[c(4,1,2,3),1],type="o", lwd=2, pch=0, col="#CC0000")
lines( c(-1.5,-.5,.5,1.5), means[c(4,1,2,3),2],type="o", lwd=2, pch=1, col="black")
lines( c(-1.45,-.45,.55,1.55), means[c(4,1,2,3),3],type="o", lwd=2, pch=2, col="#3300FF")

# Plot difference-adjusted non-significant intervals of each verb type at the verb region
arrows( -1.55, intervals[4,1], -1.55, intervals[4,2], lwd=2, angle=90, length=.15, code=3, col = "#CC0000")
arrows( -1.5, intervals[8,1], -1.5, intervals[8,2], lwd=2, angle=90, length=.15, code=3, col = "black")
arrows( -1.45, intervals[12,1], -1.45, intervals[12,2], lwd=2, angle=90, length=.15, code=3, col = "#3300FF")

# Plot difference-adjusted non-significant intervals of each verb type at the NP region
arrows( -.55, intervals[1,1], -.55, intervals[1,2], lwd=2, angle=90, length=.15, code=3, col = "#CC0000")
arrows( -.5, intervals[5,1], -.5, intervals[5,2], lwd=2, angle=90, length=.15, code=3, "black")
arrows( -.45, intervals[9,1], -.45, intervals[9,2], lwd=2, angle=90, length=.15, code=3, col = "#3300FF")

# Plot difference-adjusted non-significant intervals of each verb type at the NP+1 region
arrows( .45, intervals[2,1], .45, intervals[2,2], lwd=2, angle=90, length=.15, code=3, col = "#CC0000")
arrows( .5, intervals[6,1], .5, intervals[6,2], lwd=2, angle=90, length=.15, code=3, col = "black")
arrows( .55, intervals[10,1], .55, intervals[10,2], lwd=2, angle=90, length=.15, code=3, col = "#3300FF" )

# Plot difference-adjusted non-significant intervals of each verb type at the NP+2 region
arrows( 1.45, intervals[3,1], 1.45, intervals[3,2], lwd=2, angle=90, length=.15, code=3, col = "#CC0000")
arrows( 1.5, intervals[7,1], 1.5, intervals[7,2], lwd=2, angle=90, length=.15, code=3, col = "black")
arrows( 1.55, intervals[11,1], 1.55, intervals[11,2], lwd=2, angle=90, length=.15, code=3, col = "#3300FF" )

# Add a legend
legend("topright", inset=.01, cex=1, horiz=F, legend=c("AspV","NonPrefV","PrefV"), pch=c(0,1,2),  col=c("#CC0000","black","#3300FF" ), lty=1)



###----Figure 2: Relation between RT and NP Predictability by Verb Type----
# Install the package "ggplot2"
library(ggplot2)

# Load the subsetted data at the NP+1 region
NP1<-read.csv(file.choose(),header=TRUE)


## Summarize data. (Copied from Cookbook for R, Retrieved from http://www.cookbook-r.com/Manipulating_data/Summarizing_data/, March 10,2020 )
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

# Summarize the data at NP+1 region
ds_NP1<- summarySE(NP1,  measurevar="RT",  groupvars=c("NPPredictability","VerbType", "Item"), na.rm=TRUE)  

# Plot the relation between reading time and NP predictability
f2 <- ggplot(ds_NP1,aes(x = NPPredictability, y = RT )) +
  geom_point(aes(shape = VerbType, color = VerbType), size = 3) +
  geom_smooth(aes(color = VerbType),method = "lm", se = FALSE, formula = y ~ x) +
  labs(x = "NP predictability", y = "Reading time (ms)") +
  scale_x_continuous(limits=c(-0.01,1)) +
  scale_shape_manual(values=c(0,1,2)) +
  scale_color_manual(values=c("#CC0000", "black", "#3300FF")) +
  ggtitle("Relation between reading time and NP predictability") +
  theme(panel.background = element_rect(fill="white", color="black"),
        panel.grid = element_blank(), 
        plot.title = element_text(lineheight=.8, face="bold",hjust = .5),
        legend.position=c(.90,.85),
        legend.background = element_rect(colour = "black"),
        text = element_text(size=15))
f2






