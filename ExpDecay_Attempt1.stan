data {
  int LENGTH;
  vector[LENGTH] Y;
}
parameters {
  real lambda;
}
model {
  lambda ~ gamma(1, 1);
  for(i in 1:LENGTH)
    Y[i] ~ poisson(lambda[i]);
}
generated quantities {
  real pred;
  pred = poisson_rng(lambda);
}
