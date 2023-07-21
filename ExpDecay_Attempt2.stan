// https://link.springer.com/article/10.3758/s13428-016-0746-9#Sec1 

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
