%priors for translog

%default:
priors.ex = zeros(beta_length,1);
priors.cov = diag(ones(beta_length,1)).*100;