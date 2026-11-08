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
    point <- predict(fit,x,sample=FALSE)
    density <- predict(fit,x,type='density',sample=FALSE)
    native <- native_predict(fit$Model,x,matrix(numeric(),0L,0L))
    close_to(point$estimate[!outside],native$survival_probs[!outside])
    close_to(density$estimate[!outside],native$density_probs[!outside])
    stopifnot(all(is.na(point$estimate[outside])),
              all(is.na(density$estimate[outside])),
              !anyNA(point$estimate[!outside]),
              !anyNA(density$estimate[!outside]),identical(point$x,x))
    close_to(point$estimate[3],p)
    close_to(attr(point,'tail')$estimate,p)
    stopifnot(is.na(predict(fit,tau*(1+1e-10),sample=FALSE)$estimate))

    # Body density plus the separately estimated tail mass accounts for all mass.
    body_integral <- sum(vapply(setdiff(fit$Terminal_Nodes,fit$Tail_Node),function(j) {
        integrate(function(t) predict(fit,t,type='density',sample=FALSE)$estimate,
                  fit$Left[j],fit$Right[j])$value
    },numeric(1)))
    close_to(body_integral+p,1,1e-7)

    # Check posterior summaries within the body against native leaf draws.
    # Beyond tau, every curve summary must remain NA.
    set.seed(710+family)
    pred <- predict(fit,x,n_samples=300,burn_in=50,thin=3)
    set.seed(710+family)
    post_density <- predict(fit,x,type='density',n_samples=300,burn_in=50,thin=3)
    set.seed(710+family)
    draws <- native_sample(fit$Model,300L,50L,3L,FALSE)$terminal_probs
    tail_draws <- draws[,match(fit$Tail_Node,fit$Terminal_Nodes)]
    curves <- native_predict(fit$Model,x,draws)
    close_to(as.matrix(pred[!outside,summary_columns]),
             summarize(curves$survival_draws[,!outside,drop=FALSE]))
    close_to(as.matrix(post_density[!outside,summary_columns]),
             summarize(curves$density_draws[,!outside,drop=FALSE]))
    close_to(attr(pred,'tail')$mean,mean(tail_draws))
    close_to(unlist(pred[3,-1]),unlist(attr(pred,'tail')[1,-1]))
    stopifnot(identical(attr(pred,'tail'),attr(post_density,'tail')),
              all(is.na(pred[outside,-1])),all(is.na(post_density[outside,-1])),
              !anyNA(pred[!outside,]),!anyNA(post_density[!outside,]),
              is.null(attr(pred,'samples')),is.null(attr(pred,'tail_samples')),
              identical(before,serialize(fit,NULL)))
    close_to(unlist(pred[1,-1]),unlist(pred[7,-1]))
    set.seed(710+family)
    reversed <- predict(fit,rev(x),n_samples=300,burn_in=50,thin=3)
    close_to(as.matrix(pred),as.matrix(reversed[rev(seq_along(x)),]))
    # Empty and entirely out-of-range grids still provide tail summaries.
    for (type in c('survival','density')) {
        for (times in list(numeric(),Inf,c(tau+1,Inf))) {
            for (n in c(0L,1L,10L)) {
                result <- predict(fit,times,type=type,n_samples=n,burn_in=0)
                stopifnot(identical(result$x,times),all(is.na(result[,-1,drop=FALSE])),
                          !anyNA(attr(result,'tail')))
                close_to(attr(result,'tail')$estimate,p)
            }
        }
    }
}

# Initialization-only, legacy overlap and flat-prior fits use the same boundary.
for (args in list(list(max_EM_iter=0),list(allocated=0),list(flat=TRUE))) {
    fit <- do.call(survPT,c(list(X=X,max_M=2,verbose=FALSE),args))
    pred <- predict(fit,c(fit$Tau,2*fit$Tau,Inf),sample=FALSE)
    close_to(pred$estimate,c(fit$Tail_Prob,NA_real_,NA_real_))
}
cat('NA tail prediction checks passed.\n')
