
# eXERCISE

# 1. Randomized complete block design: Cotton

# This exercise demonstrates how to apply statistical methods that you already know, like the
# ANOVA and multiple comparisons, to standard experimental designs.

library("lme4"); library("emmeans");
library("multcomp"); library("multcompView")

# We need to modifz some options
options ( contrasts = c("contr.sum","contr.sum") )

cotton <- read.table ("v14-u41-01.csv", header=T, sep=";", dec=",",
                      stringsAsFactors = T)
str(cotton)

# There is one treatment factor (fertilizer) with five factor levels. The design is balanced: All
#treatments, i.e. all five fertilizers, occur equally often. The design is complete: Each treatment
#occurs in each block.

#Analysis of variance

boxplot( Strength ~ Fert , data=cotton )


# analysis of a RCBD looks pretty similar to a two-factor ANOVA. However, there are no
#interactions between the blocks and the main factor. In the classical analysis, both the main
#factor and the block are treated as fixed.

m1 <- aov ( Strength ~ Blk + Fert, data=cotton )
summary ( m1 )

#The different fertilizers seem to result in very different yields and indeed, the fertilizer effect is
#significant in the ANOVA whereas the block effect is not.

# means and effects

model.tables ( m1, "means" )
model.tables ( m1, "effects" )

# Least significant differences

va <- anova ( m1 ) # Analysis of variance
b <- 3; t <- 5 # Number of blocks and treatments
sed <- sqrt ( 2 * va[3,3] / b) # Standard error of a difference of means
dfe <- va[3,1] # Degrees of freedom residuals
lsd <- qt(0.975,dfe) * sed # Least significant difference
hsd <- qtukey(0.95,t,dfe) * sed / sqrt(2) # Honestly significant diff.
sed; lsd; hsd


# Pairwise comparisons, unadjusted
e.F <- emmeans ( m1, ~Fert )
cld ( e.F, adjust="none" )
test( contrast (e.F , "pairwise"), adjust="none")

# Pairwise comparisons, adjusted

cld ( e.F, adjust="tukey" )
test( contrast (e.F , "pairwise"), adjust="tukey")

# It is also possible to choose a different significance level for the CLD:
cld ( e.F, adjust="tukey", alpha=0.1 )

##############################
### 2 Latin Square: Sugar beets exercise
##############################################

#If there is structured environmental variation in two directions, we also do blocking in two
#directions. In a Latin Square, each treatment occurs in each row and each column. The row and
#column blocks are therefore complete. The design is balanced.


sugar <- read.table ( "v14-u41-02.csv", header=T, sep=";", dec=",",
                      stringsAsFactors = T)
str(sugar)
sugar

# There is one treatment factor (variety) with six factor levels.

# visualise the data
boxplot ( Yield ~ Variety, data=sugar )

# ANOVA
m2 <- aov ( Yield ~ Row + Col + Variety, data=sugar )
summary ( m2 )

# The variety e???ect is signi???cant, row and column e???ects are not.

# LSD and HSD
r <- 6
va <- anova ( m2 )
sed <- sqrt ( 2* va[4,3] / r)
dfe <- va[4,1]
lsd <- qt(0.975,dfe) * sed
hsd <- qtukey(0.95,r,dfe) * sed / sqrt(2)
sed; lsd; hsd

# r is the number of treatments (and also the number of rows and columns - why?).

# Means and pairwise comparisons
e.V <- emmeans(m2, ~Variety)
cld ( e.V, adjust="none" )
test( contrast ( e.V, "pairwise"), adjust="none" )

###############################################################
#### 3 Simple lattice: Soy beans
###############################################################

# Library

library("lme4"); library("emmeans");
library("multcomp"); library("multcompView")
options ( contrasts = c("contr.sum","contr.sum") )

#The simple lattice (Zweisatzgitter) is a partially balanced, incomplete design. It contains the
#first two replications of a balanced lattice, meaning that only some treatment pairs occur together
#in the same block.

# data 

soy <- read.table ( "v14-u41-03.csv", header=T, sep=";", dec=",",
                    stringsAsFactors = T)
str(soy)

# ANOVA

#The first line creates a nowcolumn in the data frame which helps to identify each block correctly:
#Block 1 in rep 1 is not the same as block 1 in rep 2.

soy$Block.in.Rep <- soy$Block:soy$Rep # Sequence of model terms!!
m3 <- aov ( Yield ~ Rep + Block.in.Rep + Variety , data=soy)
summary(m3)

#Apparently, the replication and the block have a significant influence on the yield.

# LSD
  
r <- 2; k <- 5
va <- anova ( m3 )
msb <- va[2,3]; 
mse <- va[4,3]; dfe <- va[4,1]
w <- (msb-mse) / (k*(r-1)*msb) # w is a special equation in the lecture.

# Different values for first and second associates:

sed.l1 <- sqrt( 2*mse/r*(1+(r-1)*w)); lsd.l1 <- qt(0.975,dfe) * sed.l1
sed.l0 <- sqrt( 2*mse/r*(1+r*w)) ; lsd.l0 <- qt(0.975,dfe) * sed.l0

sed.l1; sed.l0
lsd.l1; lsd.l0

# For the estimation of the SED, we need the value w

##We estimate two SEDs because we have so-called first and second associates. First associates
#occurred together in one block in the trial ( D 1), second associates did not ( D 0).
#Consequently, we can estimate the differences between first associates with greater precision
#than between second associates. This in turn means that the entry means of second associates
#have to show a greater difference before we consider them to be significantly different.
#Look at the field plan to find out which value you need for the comparisons.

##########
#Calculation of adjusted entry means using a mixed model##
#####

# Adjusted entry means

#"Adjusted entry means" (also: "adjusted treatment means") mean that the observed means are
#adjusted by the effect for the replication and the block. We can estimate these means by using
#the equation from the lecture slides or we can use the function lmer() that fits a linear mixed
#model on the data. The output is then handed over to emmeans() which will return the adjusted
#entry means.

m3a <- lmer(Yield ~ (1|Rep) + (1|Block.in.Rep) + Variety , data=soy)
emmeans(m3a, ~Variety)
 
# pairwise comparisons

e.v <- emmeans(m3a, ~Variety)
cld(e.v, adjust= "none" )
test (contrast (e.v, "pairwise"), adjust= "none")


################################################################
## 4 Randomized complete block design with two factors: Cowpeas
#############################################################################

# A randomized complete block design can also have two factors. In this case, each treatment
#combination occurs in all blocks.

peas <- read.table ( "v14-u41-04.csv", header=T, sep=";", dec=",",
                     stringsAsFactors = T)
str(peas)


# There are 3 varieties (variety: treatment factor A), 3 spacings (spacing: treatment factor B) and
#4 blocks.

# Analysis of variance
# The interaction between spacing and variety is included here.

boxplot ( Yield ~ Spacing:Variety , data=peas )
m4 <- aov ( Yield ~ Block + Spacing + Variety + Spacing:Variety ,
            data=peas)
summary(m4)

# All factors, including the interactions, are significant in the analysis.

# Plot interaction

with ( peas, interaction.plot ( Spacing, Variety, Yield, type="b"))

# The interaction plot shows crossing lines, pointing to interactions between variety and spacing.

# Least significant differences

r <- 4; a <- 3; b <- 3
va <- anova ( m4 )
mse <- va[5,3]; mse
dfe <- va[5,1]; dfe
sed.a <- sqrt( 2*mse/(r*b))
sed.b <- sqrt( 2*mse/(r*a))
sed.ab <- sqrt( 2*mse/(r))
lsd.a <- qt(0.975,dfe) * sed.a
lsd.b <- qt(0.975,dfe) * sed.b
lsd.ab <- qt(0.975,dfe) * sed.ab

# Like in the two-factor ANOVA from the previous chapter, we have to estimate 
#different LSDs # for each factor and for the interactions. r is the number of
#blocks. They are complete so it is also the number of replications. a is the
#number of varieties and b is the number of spacings.

# Means and pairwise comparisons
#Means for each combination of spacing and variety:
e.V1 <- emmeans(m4, ~Spacing:Variety)

#The results are averaged over the levels of block.

#Means of the spacings within the varieties

e.V2 <- emmeans(m4, ~Spacing|Variety)
cld ( e.V2, adjust="none", sort=F )
test( contrast ( e.V2, "pairwise"), adjust="none" )

#Comparisons are made between the different spacings within the levels of variety.
#Means of the varieties within the spacings:

e.V3 <- emmeans(m4, ~Variety|Spacing)
cld ( e.V3, adjust="none",sort=F )
test( contrast ( e.V3, "pairwise"), adjust="none" )

#Comparisons are made between the different varieties within the levels of spacing.


###############################################
# 5 Split plot: Alfalfa
####################################

#Split plot designs (Spaltanlagen) can be used if one of the two factors under investigation is
#restricted to application to larger plots (e.g. tillage) whereas the other factor can be applied to
#smaller plots and/or if there is more precision needed for one factor than for the other.

alf <- read.table ( "v14-u41-05.csv", header=T, sep=";", dec=",",
                    stringsAsFactors = T)
str(alf)

# There are 3 varieties, 4 dates for final cutting and 6 blocks. The main plots (Großparzellen)
#are the varieties and the sub-plots (Kleinparzellen) are the cutting dates. Main plots as well as
#sub-plots occur within the blocks!

#Research question: Do the different varieties have different yields when harvested at different
#final cutting dates?

# Analysis of variance
#Variety and cutting date are included in the analysis as fixed factors. Their interaction is fixed,
#too.

m5 <- aov(Yield ~ Variety + Date + Date:Variety
          + Error(Block/Variety), data=alf)
summary (m5)

# we calcluate error for block and variety hier.
# The split plot is a balanced design. For this situation, Error() within a formula allows the
# specification of random effects. The block is random here.

# Interaction plot

par (mfrow= c(1, 2))
with ( alf, interaction.plot ( Variety, Date, Yield, type="b" ))
with (alf, interaction.plot ( Date, Variety, Yield, type="b" ))

# Least significant differences

#The degrees of freedom are needed for the LSD. For the interactions, they cannot be taken
#from the summary()-Output but have to be approximated. One possibility is the Satterthwaite
#Approximation.

r <- 6 ; a <- 3 ; b <- 4
ms.me <- summary(m5)[[2]][[1]][2,3]
ms.se <- summary(m5)[[3]][[1]][3,3]
df.me <- summary(m5)[[2]][[1]][2,1]
df.se <- summary(m5)[[3]][[1]][3,1]
df.ia <- ( ms.me + (b-1)*ms.se ) ^2 / # Satterthwaite Approximation
  ( ms.me ^2 / df.me + ((b-1)*ms.se)^2 / df.se )
sed.1 <- sqrt(2*ms.me/r/b); sed.2 <- sqrt(2*ms.se/r/ a)
sed.3 <- sqrt(2*ms.se/r ); sed.4 <- sqrt(2*(ms.me+(b-1)*ms.se)/r/b)
lsd.1 <- qt(0.975,df.me) * sed.1 ; lsd.2 <- qt(0.975,df.se) * sed.2
lsd.3 <- qt(0.975,df.se) * sed.3 ; lsd.4 <- qt(0.975,df.ia) * sed.4 

sed.1; sed.2; sed.3; sed.4
lsd.1; lsd.2; lsd.3; lsd.4

# Pairwise comparisons
#SED estimates from above can be found in the contrast() outputs. The standard error in the
#cld() output is the SEM.

e.V <- emmeans(m5, ~Variety)
cld ( e.V, adjust="none",sort=F)
test( contrast ( e.V, "pairwise"), adjust="none" )

e.D <- emmeans(m5, ~Date)
cld ( e.D, adjust="none",sort=F)
test( contrast ( e.D, "pairwise"), adjust="none" )


e.DV <- emmeans(m5, ~Date|Variety)
cld ( e.DV, adjust="none",sort=F )
test( contrast ( e.DV, "pairwise"), adjust="none" )

e.VD <- emmeans(m5, ~Variety|Date)
cld ( e.VD, adjust="none",sort=F )
test( contrast ( e.VD, "pairwise"), adjust="none" )