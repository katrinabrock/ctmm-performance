Rcpp::cppFunction('float cpp_dexp2(float x){
  double Exp = exp(-x);
  return 1 - pow(Exp, 2.0);
}')

Rcpp::cppFunction('
double cpp_dexp2(double x, NumericVector Exp_in = NumericVector::create()){
  double Exp = Exp_in.size() == 0 ? exp(-x) : Exp_in[0];
  if(Exp<0.7071068){
      return 1 - pow(Exp, 2.0);
  } else {
      return 2.0 * Exp * sinh(x);
  }
}')

print(cpp_dexp2(1))
print(cpp_dexp2(1, 1))

if_else_dexp2(1)
