library(SurvPT)

close_to <- function(x, y, tolerance=1e-10) {
    stopifnot(isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tolerance)))
}
summary_columns <- c('lower','upper','mean','median')
summarize <- function(draws) t(apply(draws,2,function(z) {
    c(quantile(z,c(0.025,0.975),names=FALSE),mean(z),median(z))
}))
native_predict <- getFromNamespace('SurvPTPredict','SurvPT')
native_sample <- getFromNamespace('SurvPTSample','SurvPT')
X <- rbind(c(0.2,0.2),c(0,0.7),c(0.5,1.3),c(1,Inf),c(2,2),c(2,Inf))

for (family in 1:3) {
    pars <- switch(family,NULL,0.6,c(0.6,1.4))
    fit <- survPT(X,prior=family,prior_params=pars,max_M=3,verbose=FALSE)
    before <- serialize(fit,NULL)
    tau <- fit$Tau
    x <- c(tau+1,0,tau,2*tau,Inf,0.5,tau+1)
    outside <- x > tau
    p <- fit$Tail_Prob
    rate <- -log(p)/tau
    point <- predict(fit,x,sample=FALSE)
    density <- predict(fit,x,type='density',sample=FALSE)
    native <- native_predict(fit$Model,x,matrix(numeric(),0L,0L))
    close_to(point$estimate[!outside],native$survival_probs[!outside])
    close_to(density$estimate[!outside],native$density_probs[!outside])
    expected <- p*pexp(x[outside]-tau,rate,lower.tail=FALSE)
    close_to(point$estimate[outside],expected)
    close_to(density$estimate[outside],rate*expected)
    stopifnot(!anyNA(point),!anyNA(density),identical(point$x,x))
    close_to(point$estimate[3],p)
    close_to(point$estimate[4],p^2)
    close_to(attr(point,'tail')$estimate,p)
    close_to(predict(fit,tau*(1+1e-10),sample=FALSE)$estimate,p,1e-8)

    # The density over the entire tail must integrate to the fitted tail mass.
    tail_integral <- integrate(function(t)
        predict(fit,t,type='density',sample=FALSE)$estimate,tau,Inf)$value
    close_to(tail_integral,p,1e-7)
    body_integral <- sum(vapply(setdiff(fit$Terminal_Nodes,fit$Tail_Node),function(j) {
        integrate(function(t) predict(fit,t,type='density',sample=FALSE)$estimate,
                  fit$Left[j],fit$Right[j])$value
    },numeric(1)))
    close_to(body_integral+tail_integral,1,1e-7)

    # Independently reconstruct every posterior curve from the native leaf
    # draws. Each tail gets its own rate, not a rate based on the MAP or mean.
    set.seed(710+family)
    pred <- predict(fit,x,n_samples=300,burn_in=50,thin=3)
    set.seed(710+family)
    post_density <- predict(fit,x,type='density',n_samples=300,burn_in=50,thin=3)
    set.seed(710+family)
    draws <- native_sample(fit$Model,300L,50L,3L,FALSE)$terminal_probs
    tail_draws <- draws[,match(fit$Tail_Node,fit$Terminal_Nodes)]
    rates <- -log(tail_draws)/tau
    curves <- native_predict(fit$Model,x,draws)
    expected_survival <- vapply(x[outside],function(t)
        tail_draws*pexp(t-tau,rates,lower.tail=FALSE),numeric(nrow(draws)))
    curves$survival_draws[,outside] <- expected_survival
    curves$density_draws[,outside] <- expected_survival*rates
    close_to(as.matrix(pred[,summary_columns]),summarize(curves$survival_draws))
    close_to(as.matrix(post_density[,summary_columns]),summarize(curves$density_draws))
    close_to(attr(pred,'tail')$mean,mean(tail_draws))
    close_to(pred$mean[4],mean(tail_draws^2))
    stopifnot(identical(attr(pred,'tail'),attr(post_density,'tail')),
              !anyNA(pred),!anyNA(post_density),
              all(unlist(pred[5,-1])==0),all(unlist(post_density[5,-1])==0),
              is.null(attr(pred,'samples')),is.null(attr(pred,'tail_samples')),
              identical(before,serialize(fit,NULL)))
    close_to(unlist(pred[1,-1]),unlist(pred[7,-1]))
    set.seed(710+family)
    reversed <- predict(fit,rev(x),n_samples=300,burn_in=50,thin=3)
    close_to(as.matrix(pred),as.matrix(reversed[rev(seq_along(x)),]))
    stopifnot(nrow(predict(fit,numeric(),n_samples=1,burn_in=0))==0)
}

# Zero or tiny mass and infinite times must not create Inf * 0 or 0 * Inf.
kernel <- getFromNamespace('pt_exponential_tail','SurvPT')
for (density in c(FALSE,TRUE)) {
    values <- kernel(c(0,1e-300,0.2),c(2,1e300,Inf),1,density)
    stopifnot(all(is.finite(values)),all(values[1,]==0),all(values[,3]==0))
}
stopifnot(inherits(tryCatch(kernel(1,2,1,FALSE),error=identity),'error'))

# Initialization-only and legacy overlap fits also permit tail evaluation.
for (args in list(list(max_EM_iter=0),list(allocated=0),list(flat=TRUE))) {
    fit <- do.call(survPT,c(list(X=X,max_M=2,verbose=FALSE),args))
    pred <- predict(fit,c(fit$Tau,2*fit$Tau,Inf),sample=FALSE)
    close_to(pred$estimate,c(fit$Tail_Prob,fit$Tail_Prob^2,0))
}
cat('Exponential tail prediction checks passed.\n')
