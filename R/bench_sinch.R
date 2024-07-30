orig_sinch <- Vectorize( function(x, SINH=sinh(x))
{
  if(x==0)
  { return(1) }
  else
  { return(SINH/x) }
} )

ifelse_sinch <- function(x, SINH=sinh(x)) {
  return(ifelse(x==0, 1, SINH/x))
}

byidx_sinch <- function(x, SINH=sinh(x)) {
  test <- x==0
  ans <- test
  ans[test] <- 1
  ans[!test] <- SINH[!test]/x[!test]
  return(ans)
}


input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -10, -100)

run_over_inputs <- function(fn, precompute = FALSE, vector_input=FALSE){
    input <- c(100, 10, 1, .1, .01, 0, -.01, -.1, -10, -100)
    if(precompute & vector_input) return(fn(input, sinh(input)))
    if(vector_input) return(fn(input))
    if(precompute) return(sapply(input, function(x) return(fn(x, sinh(x)))))
    return(sapply(input, fn))
}

result <- microbenchmark::microbenchmark(
    orig00 = run_over_inputs(orig_sinch),
    ifelse00 = run_over_inputs(ifelse_sinch),
    byidx00 = run_over_inputs(byidx_sinch),
    orig10 = run_over_inputs(orig_sinch, TRUE),
    ifelse10 = run_over_inputs(ifelse_sinch, TRUE),
    byidx10 = run_over_inputs(byidx_sinch, TRUE),
    orig01 = run_over_inputs(orig_sinch, FALSE, TRUE),
    ifelse01 = run_over_inputs(ifelse_sinch, FALSE, TRUE),
    byidx01 = run_over_inputs(byidx_sinch, FALSE, TRUE),
    orig11 = run_over_inputs(orig_sinch, TRUE, TRUE),
    ifelse11 = run_over_inputs(ifelse_sinch, TRUE, TRUE),
    byidx11 = run_over_inputs(byidx_sinch, TRUE, TRUE),
    check = 'identical'
)


library(ggplot2)
last_nth_char<-function(s, n) substr(s, nchar(s)-n+1, nchar(s)-n+1)

result$sinh_precalc <- as.logical(as.numeric(last_nth_char(as.character(result$expr), 2)))
result$vector_input<- as.logical(as.numeric(last_nth_char(as.character(result$expr), 1)))
result$fn <- gsub('[0-9]', '', result$expr)

png('sincbench.png')
ggplot(
        result,
        aes(x=fn, y=time)
        )+geom_boxplot() + facet_grid(vector_input ~ sinh_precalc, labeller = 'label_both')
dev.off()
