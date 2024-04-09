library(diversitree)
library(phytools)
setwd("###")
phy = read.tree(file = 'dated_Vitaceae_nooutgroup.tre')
phy = force.ultrametric(phy)
state_file = read.csv('###')
Life form = c(state_file$Life form)
names(Life form) = c(state_file$Taxon)
phy$tip.state = Life form
sampling.f <- c(0.5,0.5)
names(sampling.f) <- c(0,1)
p <- starting.point.bisse(phy)
lik_bi <- make.bisse(phy, phy$tip.state, sampling.f=sampling.f)
fit_bi <- find.mle(lik_bi, p, method="subplex")
fit_bi$lnLik
lik.1 <- constrain(lik_bi, lambda1 ~ lambda0)
fit.1 <- find.mle(lik.1, p[argnames(lik.1)])
fit.1$lnLik
lik.2 <- constrain(lik_bi, q10 ~ q01)
fit.2 <- find.mle(lik.2, p[argnames(lik.2)])
fit.2$lnLik
lik.3 <- constrain(lik_bi, mu1 ~ mu0)
fit.3 <- find.mle(lik.3, p[argnames(lik.3)])
fit.3$lnLik
lik.4 <- constrain(lik_bi, mu1 ~ mu0, q10 ~ q01)
fit.4 <- find.mle(lik.4, p[argnames(lik.4)])
fit.4$lnLik
model = round(rbind(full=coef(fit_bi), equal.1=coef(fit.1, TRUE), equal.2=coef(fit.2, TRUE), equal.3=coef(fit.3, TRUE), equal.4=coef(fit.4, TRUE)), 3)
write.csv(model, file = "Life form_BiSSE_model.csv")
anv = anova(fit_bi, equal.1=fit.1, equal.2=fit.2, equal.3=fit.3, equal.4=fit.4)
write.csv(anv, file = "Life form_BiSSE_anova.csv")
prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))
set.seed(1)
tmp <- mcmc(lik_bi, fit_bi$par, nsteps=100, prior=prior,
lower=0, w=rep(1, 6), print.every=0)
w <- diff(sapply(tmp[2:7], range))
samples <- mcmc(lik_bi, fit_bi$par, nsteps=10000, w=w, lower=0, prior=prior,
print.every=1000)
write.csv(samples, file = "Life form_BiSSE_samples.csv")
samples = read.csv(file = "Life form_BiSSE_samples.csv")
