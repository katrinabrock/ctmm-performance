orig_dexp2 <- function(x,Exp=exp(-x)) { ifelse(Exp<0.7071068,1-Exp^2,2*Exp*sinh(x)) }
if_else_dexp2 <- function(x,Exp=exp(-x)) { if(Exp<0.7071068) return(1-Exp^2) else return(2*Exp*sinh(x))}

byidx_dexp2 <- function(x,Exp=exp(-x)) {
    test <- Exp<0.7071068
    ans <- test
    ans[test] <- 1-Exp[test]^2
    ans[!test] <- 2*Exp[!test]*sinh(x[!test])
    return(ans)
}

Rcpp::cppFunction('
double cpp_dexp2(double x, NumericVector Exp_in = NumericVector::create()){
  double Exp = Exp_in.size() == 0 ? exp(-x) : Exp_in[0];
  if(Exp<0.7071068){
      return 1 - pow(Exp, 2.0);
  } else {
      return 2.0 * Exp * sinh(x);
  }
}')

run_over_inputs <- function(fn, precompute_Exp = FALSE, vector_input=FALSE){
    input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -10, -100)
    if(precompute_Exp & vector_input) return(fn(input, exp(-input)))
    if(vector_input) return(fn(input))
    if(precompute_Exp) return(sapply(input, function(x) return(fn(x, exp(-x)))))
    return(sapply(input, fn))
}


result <- microbenchmark::microbenchmark(
    orig00 = run_over_inputs(orig_dexp2),
    if_else00 = run_over_inputs(if_else_dexp2),
    cpp00 = run_over_inputs(cpp_dexp2),
    byidx00 = run_over_inputs(byidx_dexp2),
    orig10 = run_over_inputs(orig_dexp2, TRUE),
    if_else10 = run_over_inputs(if_else_dexp2, TRUE),
    cpp10 = run_over_inputs(cpp_dexp2, TRUE),
    byidx10 = run_over_inputs(byidx_dexp2, TRUE),
    orig01 = run_over_inputs(orig_dexp2, FALSE, TRUE),
    byidx01 = run_over_inputs(byidx_dexp2, FALSE, TRUE),
    orig11 = run_over_inputs(orig_dexp2, TRUE, TRUE),
    byidx11 = run_over_inputs(byidx_dexp2, TRUE, TRUE),
    check = 'identical'
)


library(ggplot2)
last_nth_char<-function(s, n) substr(s, nchar(s)-n+1, nchar(s)-n+1)

result$exp_precalc <- as.logical(as.numeric(last_nth_char(as.character(result$expr), 2)))
result$vector_input<- as.logical(as.numeric(last_nth_char(as.character(result$expr), 1)))
result$fn <- gsub('[0-9]', '', result$expr)

png('expbench.png')
ggplot(
        result,
        aes(x=fn, y=time)
        )+geom_boxplot() + facet_grid(vector_input ~ exp_precalc, labeller = 'label_both')
dev.off()