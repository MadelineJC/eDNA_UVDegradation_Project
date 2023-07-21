
data {
  int LENGTH;
  vector[LENGTH] Y;
}
parameters {
  real lambda;
}
model {
  real alpha;
  real beta;
  alpha = 1;
  beta = 1;
  lambda ~ gamma(alpha, beta);
  Y ~ exponential(lambda);
}
generated quantities {
  real pred;
  pred = exponential_rng(lambda);
}

