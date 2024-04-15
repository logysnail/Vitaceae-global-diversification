##参考http://bamm-project.org/

setwd("E:/VitaceaeDiversification_YYC/analysis/diversification/bamm_v16")

#树文件准备=============================================================================================================================================
install.packages('ape')
install.packages('phytools')
install.packages('BAMMtools')
install.packages('plot')


library(ape)
library(phytools)
library(BAMMtools)
tre<-read.tree("dated_Vitaceae_505_1-e12_median.tre")
#Drop outgroups
outgroup<-c("Leea_rubra_MAC05716_A577","Leea_guineensis_Wen2010091","Leea_monticola_Wen9569","Leea_cuspidifera_Wen9621","Leea_aequata_Wen8382","Leea_indica_Wen15011","Leea_macrophylla_Ren55105_R24","Leea_gonioptera_Wen10711","Leea_aculeata_Wen7058","Leea_spinea_Wen9575")
tre<-drop.tip(tre,outgroup)
write.tree(tre,"dated_Vitaceae_505_1-e12_dropout.tre")
#Your phylogenetic tree must be ultrametric, it must be fully bifurcating (no polytomies), and all branch lengths must be greater than 0
is.ultrametric(tre)
#if FALSE
#强行ultramatric
library(phytools)
tre<-force.ultrametric(tre,method="extend")
write.tree(tre,"dated_Vitaceae_508_1e-13_dropout_ultra.tre")
# 排列好
tre<-ladderize(tre)
write.tree(tre,"dated_Vitaceae_508_1e-13_dropout_ultra_ladd.tre")
# read the tree
tre<-read.tree("dated_Vitaceae_508_1e-13_dropout_ultra_ladd.tre")
is.binary(tre)
min(tre$edge.length)

#Prior设置=============================================================================================================================================
packageDescription('BAMMtools') #check version of BAMMtools
library(BAMMtools)
setBAMMpriors(tre)
##用myPriors.txt中参数替换ControlFile.txt中部分参数
#可选额外设置测试==========================================================================================================================================
#(Alternative setting in control file)
#The frequency in which a time-flip proposal occurs, relative to other proposals, is given by updateRateLambdaTimeMod
#under development?
updateRateLambdaTimeMode = 0/1
# Relative frequency of MCMC moves that flip the time mode
# (time-constant <=> time-variable)
# METROPOLIS COUPLED MCMC
#The \Delta T value should be set such that the probability of accepting a chain swap proposal is between 20% and 60%.
#For small to medium sized trees (< 1,000 taxa), we have found that \Delta T = 0.1 works well.
#For large trees, \Delta T = 0.05 or \Delta T = 0.01 works better.
deltaT = 0.01/0.05/0.1
#Determine the best settings
#In the tools directory of the BAMM GitHub repository, we have provided an R script, chainSwapPercent.R
#print out the percent acceptance of the chain swap proposals by testing all combinations of the given values of swap period, \Delta T, and number of chains
#assume that the scripts are located in ~/bamm/tools and that the BAMM executable can be run as bamm
#run in R
#source("~/bamm/tools/chainSwapPercent.R")
source("E:/TEST/deltaT/chainSwapPercent.R")
chainSwapPercent(bammPath = 'bamm', controlfile = 'divcontrol_Vitaceae.txt',
    nChains = c(2, 4, 6), deltaT = c(0.01, 0.05,0.1),
    swapPeriod = c(100, 1000), nGenerations = 20000,
    burnin = 0.2, deleteTempFiles = TRUE)
#测试结果percent在20%~60%为最优？然后改control file中的设置
#每次测试结果不同————取numberOfChains = 4，\Delta T = 0.1，swapPeriod = 1000

#运行bamm==============================================================================================================================================
#command-line to run bamm
#controlfile,和树和取样比例文件和bamm程序放在同一个文件夹
bamm -c divcontrol_Vitaceae.txt

#生成文件：
#The /run_info.txt/ file, containing a summary of your parameters/settings
#An /mcmc_out.txt/ or equivalent file, containing raw MCMC information useful in diagnosing convergence
#An /event_data.txt/ file or equivalent, containing all of evolutionary rate parameters and their topological mappings
#A /prior.txt/ file or equivalent, giving the prior expectation on the number of shift events (this is optional and can be turned off).
#A /chain_swap.txt/ file, containing data about each chain swap proposal (when a proposal occured, which chains might be swapped, and whether the swap was accepted).

##运行结束后，copy “event_data.txt" and "mcmc_out.txt"至R working directory

#========================================================================================================================================================
###1. Assessing MCMC convergence
mcmcout<-read.csv("mcmc_out.txt",header=T) #读取mcmcout
plot(mcmcout$logLik ~ mcmcout$generation) #plot the log-likelihood trace of your MCMC output file.
  ##PDF: This can give you a ballpark idea of whether your run has converged

burnstart <- floor(0.1 * nrow(mcmcout)) #discard some as burnin (10%)
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#ESS检测
library(coda) #Use coda library to check the effective sample sizes of the log-likelihood and the number of shift events present in each sample
effectiveSize(postburn$N_shifts) #at least 200
effectiveSize(postburn$logLik) #at least 200
#additional test for convergence by analyzing multiple independent BAMM runs: whether the runs are converging on similar distributions by analyzing the branch-specific marginal rate shift probabilities (see marginalShiftProbsTree).

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###2. The bammdata object
library(BAMMtools)
tre <- read.tree("dated_Vitaceae_505_1-e12_dropout.tre")
tre <- ladderize(tre, right = FALSE)
edata<-getEventData(tre,eventdata="event_data.txt",burnin=0.1) #读取edata
summary(edata)

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###3. Mean phylorate plot平均分化速率作图
#visualizing mean, model-averaged diversification rates at any point along every branch of a phylogenetic tree
#The mean phylorate plot summarizes rate dynamics across the entire posterior
#colors for any branch segment reflect the mean of the marginal density of evolutionary rates at any particular point in time
#不同颜色展示模式
q <- plot(edata, breaksmethod='linear', lwd=2)
addBAMMshifts(edata, cex=2)
addBAMMlegend(q, location='topleft')
#plot函数不再支持plot.bammdata.base包3.6开始就不支持plot.*函数了，统一使用plot()就好了。 forecast包从8.1开始也不支持forecast.*了，改为统一使用forecast
q <- plot(edata,breaksmethod='linear', logcolor = TRUE, lwd=2)
addBAMMshifts(edata, cex=2)
addBAMMlegend(q, location='topleft')

q <- plot(edata, breaksmethod='jenks', lwd=2)
addBAMMshifts(edata, cex=2)
addBAMMlegend(q, location='topleft')


pdffn = "Mean_phylorate_jenks_tiplabel.pdf"
pdf(pdffn, width=50, height=150)
q <- plot(edata, labels = TRUE, breaksmethod='jenks', lwd=1) #相应参数可改变
addBAMMshifts(edata, cex=6)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

pdffn = "Mean_phylorate_jenks.pdf"
pdf(pdffn, width=14, height=20)
q <- plot(edata, breaksmethod='jenks', lwd=1) #相应参数可改变
addBAMMshifts(edata, cex=3)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###4. A phylorate plot of any sample任一样本的分化速率作图
e20 <- subsetEventData(edata, index = 20) #the 20th sample from the posterior
plot(e20, lwd=2)
addBAMMshifts(e20, cex=2)

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###5. How many rate shifts速率转变次数的概率分布
post_probs <- table(postburn$N_shifts) / nrow(postburn) #a vector of model posterior probabilities
names(post_probs) #To see which models are part of the set that were sampled；没有“0”说明有速率转变的明确证据
post_probs[1] / post_probs[2] #to compute the posterior odds ratio for (say) two models ‘X’ and ‘Y’ (X and Y must be integers)
shift_probs <- summary(edata) #a dataframe giving the posterior probabilities of each rate shift count observed during simulation of the posterior
shift_probs #查看number of rate shifts。每一行代表一个模型，第一列为模型速度转变的次数，第二列为相应模型的relative probability。该数据可在excel中做条形图

#是否发生转变看number of shifts的后验概率/k shifts模型的Bayes factor评估------------------------------------
#The evidence for rate heterogeneity comes from considering the posterior probabilities on the number of shifts
#or - even better - the Bayes factor evidence in favor of model with k shifts (M_k) relative to a model with 0 shifts (M_0)

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###6. Bayes factors for model comparison速率转变次数的模型选择
#Bayes factors – unlike model posterior probabilities – are theoretically invariant to the prior on the number of shifts
#We suggest that (usually) the overall best model from a BAMM analysis is the model with the highest Bayes factor relative to the null model, M_0
#compute Bayes factor evidence in favor of one model relative to another
bfmat <- computeBayesFactors(mcmcout, expectedNumberOfShifts=1, burnin=0.1) #a pairwise matrix of Bayes factors
#In general, the first column of this output matrix is the comparison of all the models relative to the model with the lowest number of supported shifts (often zero)
#>12 to be consistent with at least some effect; > 20 generally imply strong evidence for one model over another; > 50 are very strong evidence

#先验与后验分布的差异visualizing the prior and posterior simultaneously------------------------------------
#to see what models are not being sampled in the posterior, and to evaluate how far from the prior the posterior has moved
plotPrior(mcmcout, expectedNumberOfShifts=1)

#树上各分支发生转变的先验概率prior shift probabilities-----------------------------------------------------
branch_priors <- getBranchShiftPriors(tre, expectedNumberOfShifts = 1)
plot.phylo(branch_priors)
#The object branch_priors is now a copy of our phylogenetic tree, but where each branch length is equal to the prior probability of a rate shift
#The prior probability of a shift is proportional to the branch length

#树上各分支发生转变的概率（边际——后验？）marginal shift probabilities--------------------------------------
marg_probs <- marginalShiftProbsTree(edata)
plot.phylo(marg_probs)
#The object marg_probs is a copy of the original phylogenetic tree, but where the branch lengths have been replaced by the branch-specific marginal shift probabilities
#the length of a given branch is equal to the percentage of samples from the posterior that contain a rate shift on that particular branch
#massive evidence for rate heterogeneity can be in your dataset, but marginal shift probabilities will be a function of the frequency distribution of distinct alternative shift configurations

#树上各分支发生转变的后验可能性posterior odds of a rate shift normalized by the prior expectation----------
mo <- marginalOddsRatioBranches(edata, branch_priors)
plot.phylo(mo)#找不到对象mo
#The object mo is a copy of our phylogenetic tree where the branch lengths have been scaled to equal the corresponding marginal odds ratio
#The marginal odds ratios tell you the weight of evidence supporting a shift along a particular branch after normalizing by the branch length itself
#What the marginal odds do provide is an estimate of the “density” of shifts on a particular branch, independent of the length of the branch
#When there is a conflict, we suggest that the marginal shift probabilities are your overall best measure of the location of a shift
#But the marginal odds ratio is critical for distinguishing core and non-core rate shifts

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###7. Bayesian credible sets of shift configurations筛选可信集
#marginal odds ratios less than 5 are fairly weak, excluding them to focus on branches with marginal probabilities that are substantially elevated relative to the prior
#summarize the credible set of macroevolutionary rate configurations using this marginal odds ratio criterion
#The 95% credible set is the set of distinct shift configurations that account for 95% of the probability of the data
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct #the number of distinct shift configurations in the data
summary(css)
plot(css)
plot.credibleshiftset(css) #这步不需要？plots for each of the N shift configurations with the highest posterior probabilities

pdffn = "0.95Credible_Shift_Set.pdf"
pdf(pdffn, width=14, height=20)
plot(css, breaksmethod='jenks', legend = TRUE)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

css$indices[[x]] #第x可能的转变模式
index <- css$indices[[x]][y] #第x可能的转变模式中的第y个样本
rsample <- subsetEventData(edata, index=index)
plot(rsample)
addBAMMshifts(rsample, cex=2)
#某一结果展示

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###8. The single best shift configuration展示单个最佳的速率转变模式
#8.1最优/最可能的转变模式Overall best shift configuration【建议展示MAP】-----------------------------------
#show the maximum a posteriori probability (MAP) shift configuration
#This is the distinct shift configuration with the highest posterior probability - e.g., the one that was sampled most often
#if you show just a single shift configuration estimated with BAMM for publication, we recommend showing the MAP configuration
#which should also be the most frequent shift configuration in your credible shift set
#should match the first plot from the panel of plots we obtained with plot.credibleshiftset
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
plot(best, lwd=1.25)
addBAMMshifts(best, cex=2)

pdffn = "MAP_best_config.pdf"
pdf(pdffn, width=14, height=20)
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
q <- plot(best, breaksmethod='jenks', lwd=1)
addBAMMshifts(best, cex=3)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

pdffn = "MAP_best_config_tip.pdf"
pdf(pdffn, width=60, height=70)
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
q <- plot(best, labels = TRUE, breaksmethod='jenks', lwd=1.25)
addBAMMshifts(best, cex=6)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#8.2最大可信转变模式Maximum shift credibility configuration【often same as MAP】----------------------------
#An alternative estimate of the most likely shift configuration is the maximum shift credibility configuration (MSC)
#to extract the shift configuration that maximizes the marginal probability of rate shifts along individual branches
#This concept is analogous to the maximum clade credibility tree in a Bayesian phylogenetic analysis
#We generally recommend using the MAP shift configuration (getBestShiftConfiguration) over the MSC configuration
#except for very large phylogenies (e.g. thousands of taxa, all shift configurations may have low probability)
#Often, however, the two approaches will estimate the same shift configuration
msc.set <- maximumShiftCredibility(edata, maximize='product') #A number of samples from the posterior potentially have identical credibility scores
msc.config <- subsetEventData(edata, index = msc.set$sampleindex) #pull out a single representative
plot(msc.config, lwd=2)
addBAMMshifts(msc.config, cex = 2)

pdffn = "MSC_best_config_alternative.pdf"
pdf(pdffn, width=14, height=20)
q <- plot(msc.config, breaksmethod='jenks', lwd=1)
addBAMMshifts(msc.config, cex=3)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

pdffn = "MSC_best_config_alternative_tip.pdf"
pdf(pdffn, width=50, height=70)
q <- plot(msc.config, labels = TRUE, breaksmethod='jenks', lwd=1)
addBAMMshifts(msc.config, cex=6)
addBAMMlegend(q, location='topleft')
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###9. Rate-through-time analysis速率时间曲线
pdffn = "RTTplot_Vitaceae.pdf"
pdf(pdffn, width=8, height=20)
par(mfrow=c(3,1))
plotRateThroughTime(edata, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,0.5), cex.axis=2)
plotRateThroughTime(edata, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,0.5), cex.axis=2)
plotRateThroughTime(edata, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,0.5), cex.axis=2)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#clade RTT (Tetrastigma; Cyphostemma; Cissus; Vitis-Ampelocissus)
pdffn = "tree_nodelabel.pdf" #check nodelabel
pdf(pdffn, width=50, height=70)
plot.phylo(tre)
nodelabels()
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

Temrca <- getMRCA(tre, tip = c("Tetrastigma_erubescens_var_monophyllum_APG27514_Dang522", "Tetrastigma_trifoliolatum_Wen10758_V004"))
Cymrca <- getMRCA(tre, tip = c("Cyphostemma_sp_RG6878_R088", "Cyphostemma_dehongense_jianfei_22"))
Cimrca <- getMRCA(tre, tip = c("Cissus_sp_JGR_2013_Lombardi8021", "CissusII_penninervis_Jackes_9850"))
VAmrca <- getMRCA(tre, tip = c("Vitis_sp_Wen12693_MG664802", "Ampelocissus_divaricata_ATB12_Dang605"))
Vimrca <- getMRCA(tre, tip = c("Vitis_sp_Wen12693_MG664802", "Ampelocissus_acapulcensis_Wen8696"))
Ammrca <- getMRCA(tre, tip = c("Ampelocissus_sp_nov_Wen8343_A428", "Ampelocissus_divaricata_ATB12_Dang605"))

pdffn = "RTTplots_Vitaceae&clades.pdf"
pdf(pdffn, width=14, height=20)
par(mfrow=c(7,3))
plotRateThroughTime(edata, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #all Vitaceae
text(x=60, y= 1.5, label="Vitaceae speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitaceae extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitaceae netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Temrca, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #Tetrastigma
text(x=60, y= 1.5, label="Tetrastigma speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Temrca, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Tetrastigma extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Temrca, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Tetrastigma netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Cymrca, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #Cyphostemma
text(x=60, y= 1.5, label="Cyphostemma speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Cymrca, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Cyphostemma extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Cymrca, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Cyphostemma netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Cimrca, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #Cissus
text(x=60, y= 1.5, label="Cissus speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Cimrca, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Cissus extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Cimrca, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Cissus netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=VAmrca, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #Vitis-Ampelocissus
text(x=60, y= 1.5, label="Vitis-Ampelocissus speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=VAmrca, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitis-Ampelocissus extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=VAmrca, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitis-Ampelocissus netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Vimrca, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #Vitis
text(x=60, y= 1.5, label="Vitis speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Vimrca, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitis extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Vimrca, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitis netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Ammrca, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2) #Ampelocissus
text(x=60, y= 1.5, label="Ampelocissus speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Ammrca, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Ampelocissus extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Ammrca, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1.5), xlim=c(80,0), cex.axis=2)
text(x=60, y= 1.5, label="Ampelocissus netdiv", font=4, cex=1.5, pos=4)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it
###################################
Tetrastigma <- getMRCA(tre, tip = c("Tetrastigma_sp_Wen8256", "Tetrastigma_sp_Wen10768_V006"))
Cyphostemma <- getMRCA(tre, tip = c("Cyphostemma_sp_RG6878_R088", "Cyphostemma_dehongense_jianfei_22"))
Vitis <- getMRCA(tre, tip = c("Vitis_vinifera_NC_007957", "Vitis_sp_Wen12676_MG664805"))
Ampelocissus <- getMRCA(tre, tip = c("Ampelocissus_sp_nov_Wen8343_A428", "Nothocissus_spicifera_Wen11675_N1128"))

pdffn = "RTTplots_acc_clade.pdf"
pdf(pdffn, width=14, height=20)
par(mfrow=c(7,3))
plotRateThroughTime(edata, node=Tetrastigma, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1), xlim=c(25,0), cex.axis=2) #Tetrastigma
text(x=60, y= 1.5, label="Tetrastigma speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Tetrastigma, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Tetrastigma extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Tetrastigma, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Tetrastigma netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Cyphostemma, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1), xlim=c(25,0), cex.axis=2) #Cyphostemma
text(x=60, y= 1.5, label="Cyphostemma speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Cyphostemma, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Cyphostemma extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Cyphostemma, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Cyphostemma netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Vitis, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1), xlim=c(25,0), cex.axis=2) #Cissus
text(x=60, y= 1.5, label="Vitis speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Vitis, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitis extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Vitis, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Vitis netdiv", font=4, cex=1.5, pos=4)

plotRateThroughTime(edata, node=Ampelocissus, ratetype="speciation", intervalCol="blue", avgCol="blue", ylim=c(0,1), xlim=c(25,0), cex.axis=2) #Cissus
text(x=60, y= 1.5, label="Ampelocissus speciation", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Ampelocissus, ratetype="extinction", intervalCol="red", avgCol="red", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Ampelocissus extinction", font=4, cex=1.5, pos=4)
plotRateThroughTime(edata, node=Ampelocissus, ratetype="netdiv", intervalCol="green", avgCol="green", ylim=c(0,1), xlim=c(25,0), cex.axis=2)
text(x=60, y= 1.5, label="Ampelocissus netdiv", font=4, cex=1.5, pos=4)


dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###10. Clade-specific evolutionary rates分支分化速率
allrates <- getCladeRates(edata) #the mean rate across all species for each sample in the posterior, averaging over any rate heterogeneity that occurs within your focal clade

#mean rate and the 95% highest posterior density (HPD) for all Vitaceae
mean(allrates$lambda)
quantile(allrates$lambda, c(0.05, 0.95))
#mean rate and the 95% highest posterior density (HPD) for a specific clade
Tetrastigmarates <- getCladeRates(edata, node=Temrca)
mean(Tetrastigmarates$lambda)
quantile(Tetrastigmarates$lambda, c(0.05, 0.95))
Cyphostemmarates <- getCladeRates(edata, node=Cymrca)
mean(Cyphostemmarates$lambda)
quantile(Cyphostemmarates$lambda, c(0.05, 0.95))
Cissusrates <- getCladeRates(edata, node=Cimrca)
mean(Cissusrates$lambda)
quantile(Cissusrates$lambda, c(0.05, 0.95))
VArates <- getCladeRates(edata, node=Cimrca)
mean(VArates$lambda)
quantile(VArates$lambda, c(0.05, 0.95))
mean(VArates$mu)
quantile(VArates$mu, c(0.05, 0.95))

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###11. Macroevolutionary cohort analysis群组分析
#provides a way of summarizing the extent to which species share correlated macroevolutionary dynamics
#visualize the pairwise probabilities that any two species share a common macroevolutionary rate regime
cmat <- getCohortMatrix(edata)
pdffn = "Macroevolutionary_cohort_matrix.pdf"
pdf(pdffn, width=14, height=20)
cohorts(cmat, edata, lwd=3, pal="temperature", use.plot.bammdata=TRUE)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#--------------------------------------------------------------------------------------------------------------------------------------------------------
###12. Cumulative shift probabilities累积转变概率？
#a potentially useful complement to cohort analysis
#shows the cumulative probability, on each branch, that a shift occurred somewhere between the focal branch and the root of the tree
#occurrence of such a shift implies that evolutionary dynamics on the focal branch are decoupled from the “background” diversification at the root of the tree
cst <- cumulativeShiftProbsTree(edata)
plot.phylo(cst)

pdffn = "0.95Cumulative_shift_probability_tree.pdf"
pdf(pdffn, width=50, height=70)
cst <- cumulativeShiftProbsTree(edata)
edgecols <- rep('black', length(tre$edge.length))
is_highprobshift <- cst$edge.length >= 0.95
edgecols[ is_highprobshift ] <- "red"
plot.phylo(tre, edge.color = edgecols)
dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#=========================================================================================================================================================
#rate heterogeneity in the form of a histogram（rate的频率分布与颜色呈现）
library(BAMMtools)
#load data
data(primates, events.primates)
ed <- getEventData(primates, events.primates, burnin=0.25, type = 'trait')

#create phylorate plot to generate output
q <- plot(edata, breaksmethod='jenks') #linear; quantile; logcolor; color.interval
ratesHistogram(q, plotBrks = TRUE, xlab = 'trait rates')
title(main='jenks', cex.main=1)

q <- plot(edata, breaksmethod='linear') #linear; quantile; logcolor; color.interval
