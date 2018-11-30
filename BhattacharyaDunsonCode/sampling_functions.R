## MCMC Functions for Latent factor models in Bhattacharya and Dunson (2011)

update_latent_factors <- function(mc, Y){
  # Lmsg = bsxfun(@times,Lambda,ps);
  # Veta1 = eye(k) + Lmsg'*Lambda;
  #       T = cholcov(Veta1); [Q,R] = qr(T);
  #       S = inv(R); Veta = S*S';                   % Veta = inv(Veta1)
  #       Meta = Y*Lmsg*Veta;                        % n x k 
  #       eta = Meta + normrnd(0,1,[n,k])*S';        % update eta in a block
  
  # Lambda_prec3 = with(mc, crossprod(Lambda, diag(1/sigma)))
  # Lambda_prec2 = t(with(mc, apply(Lambda, 2, function(l) l*(1/sigma))))
  k = mc$k

  Lambda_prec = t(mc$Lambda * 1/mc$sigma)
  Prec = with(mc, diag(k) +  Lambda_prec %*% Lambda)
  Prec_chol = chol(Prec)
  R = qr.R(qr(Prec_chol))
  S = solve(R)
  cov = tcrossprod(S)
  
  # Update all of the latent factors at once
  mc$eta = t(tcrossprod(cov %*% Lambda_prec, Y) + t(S) %*% matrix(rnorm(n*k), nrow = k, ncol = n))
  mc
}

# Update the factor loadings, Lambda
update_factor_loadings <- function(mc, Y){
#           eta2 = eta'*eta;
# for j = 1:p
# Qlam = diag(Plam(j,:)) + ps(j)*eta2; blam = ps(j)*(eta'*Y(:,j));
#                                                    Llam = chol(Qlam,'lower'); zlam = normrnd(0,1,k,1);
#                                                    vlam = Llam\blam; mlam = Llam'\vlam; ylam = Llam'\zlam;
#                                                    Lambda(j,:) = (ylam + mlam)';  
#                                                    end

  Y_cols = lapply(1:p, function(i) Y[,i])
  eta_t_eta = with(mc,crossprod(eta))
  phitau = t(with(mc, t(phi) * tauh))
  # phitau = with(mc, tcrossprod(phi, diag(tauh)))
  phitau_rows = lapply(1:p, function(i) phitau[i,])
  # sigma_j = as.vector(mc$sigma, "list")
  mc$Lambda = t(mapply(function(y, phitau, sigma) update_factor_loading_row(y, phitau, sigma, mc$eta, eta_t_eta, mc$k), Y_cols, phitau_rows, mc$sigma))
  mc
}

update_factor_loading_row <- function(y, phitau, sigma, eta, eta_t_eta, k){
  D = diag(phitau)
  prec = D + 1/sigma * eta_t_eta
  L = chol(prec)
  mean = solve(prec, 1/sigma * t(eta) %*% y)
  mean + solve(L, matrix(rnorm(k), nrow=k))
}

update_phi <- function(mc){
  # psijh = gamrnd(df/2 + 0.5,1./(df/2 + bsxfun(@times,Lambda.^2,tauh')));
  mc$phi = t(apply(mc$Lambda, 1, function(L) rgamma(mc$k, (df + 1)/2, (df + L^2*mc$tauh)/2)))
  mc
}

update_delta <- function(mc){
  # mat = bsxfun(@times,psijh,Lambda.^2);
  # ad = ad1 + 0.5*p*k; bd = bd1 + 0.5*(1/delta(1))*sum(tauh.*sum(mat)');
  #                                                     delta(1) = gamrnd(ad,1/bd);
  #                                                     tauh = cumprod(delta);
  #                                                     
  #                                                     for h = 2:k
  #                                                     ad = ad2 + 0.5*p*(k-h+1); bd = bd2 + 0.5*(1/delta(h))*sum(tauh(h:end).*sum(mat(:,h:end))');
  # delta(h) = gamrnd(ad,1/bd); tauh = cumprod(delta);
  # end
  k = mc$k

  Lambdasq_phi_sum = apply(mc$Lambda^2 * mc$phi, 2, sum)
  mc$delta[1] = rgamma(1, ad1 + p*k/2, bd1 + (1/2) * (1/mc$delta[1]) * sum(mc$tauh*Lambdasq_phi_sum))
  mc$tauh = cumprod(mc$delta)
  
  for(h in 2:k){
    # This is the version from the original paper
    # mc$delta[h] = rgamma(1, ad2 + p*(k - h + 1)/2, bd2 + (1/2)*(1/mc$delta[h]) * sum(mc$tauh[h:k] * Lambdasq_phi_sum[h:k]))
    
    # This is modified based on the Durante note
    mc$delta[h] = rgamma(1, ad2 + p*(k - h + 1)/2, bd2 * (1/mc$delta[h-1]) + (1/2)*(1/mc$delta[h]) * sum(mc$tauh[h:k] * Lambdasq_phi_sum[h:k]))
    mc$tauh = cumprod(mc$delta)
  }
  
  mc
}

update_sigma <- function(mc, Y){
  # Ytil = Y - eta*Lambda'; 
  #       ps=gamrnd(as + 0.5*n,1./(bs+0.5*sum(Ytil.^2)))';
  # Sigma=diag(1./ps);
  err = Y - tcrossprod(mc$eta, mc$Lambda)
  err_sum = apply(abs(err^2), 2, sum)
  mc$sigma = 1/rgamma(p, as + n/2, bs + err_sum/2)
  mc
}

update_k <- function(mc){
  # prob = 1/exp(b0 + b1*i);                % probability of adapting
  # uu = rand;
  # lind = sum(abs(Lambda) < epsilon)/p;    % proportion of elements in each column less than eps in magnitude
  # vec = lind >=prop;num = sum(vec);       % number of redundant columns
  # 
  # if uu < prob
  # if  i > 20 && num == 0 && all(lind < 0.995)

  # elseif num > 0

  # end
  # end
  # nofout(i+1) = k;
  # 
  k = mc$k
  
  prob = 1/exp(b0 + b1*i)
  u = runif(1)
  proportion_redundant = apply(abs(mc$Lambda) < epsilon, 2, sum) / p
  redundant = proportion_redundant >= proportion
  nredundant = sum(redundant)
  
  if( u < prob){
    if((i > 20) & (nredundant == 0) & all(proportion_redundant < 0.995)){
      k = k + 1
      mc = add_column(mc)
    }
    else if(nredundant > 0){
      k = max(1, k - nredundant)
      mc = drop_columns(mc, redundant)
    }
  }
  
  mc$k = k
  return(mc)
}

add_column <- function(mc){
  #     k = k + 1;
  #     Lambda(:,k) = zeros(p,1);
  #     eta(:,k) = normrnd(0,1,[n,1]);
  #     psijh(:,k) = gamrnd(df/2,2/df,[p,1]);
  #     delta(k) = gamrnd(ad2,1/bd2);
  #     tauh = cumprod(delta);
  #     Plam = bsxfun(@times,psijh,tauh');
  mc$Lambda = cbind(mc$Lambda, 0)
  mc$eta = cbind(mc$eta, rnorm(n))
  mc$phi = cbind(mc$phi, rgamma(p, df/2, df/2))
  mc$delta = c(mc$delta, rgamma(1, ad2, bd2))
  mc$tauh = cumprod(mc$delta)
  
  mc
}

drop_columns <- function(mc, redundant){
  #     nonred = setdiff(1:k,find(vec)); % non-redundant loadings columns
  #     k = max(k - num,1);
  #     Lambda = Lambda(:,nonred);
  #     psijh = psijh(:,nonred);
  #     eta = eta(:,nonred);
  #     delta = delta(nonred);
  #     tauh = cumprod(delta);
  #     Plam = bsxfun(@times,psijh,tauh');
  keep = !redundant
  mc$Lambda = mc$Lambda[,keep]
  mc$phi = mc$phi[,keep]
  mc$eta = mc$eta[,keep]
  mc$delta = mc$delta[keep]
  mc$tauh = cumprod(mc$delta)
  mc
}

save_samples <- function(stores, mc, i){
#   % -- save sampled values (after thinning) -- %
#     if mod(i,thin)==0 && i > burn
#   Omega = Lambda*Lambda' + Sigma;
#   Omega1 = Omega .* sqrt(VY'*VY);
# Omegaout = Omegaout + Omega(:)/sp;
# Omega1out = Omega1out + Omega1(:)/sp;
# nof1out((i-burn)/thin) = (nofout((i-burn)/thin) - num)*(num > 0);
# end
  
  stores$lambda_store[[i]] = mc$Lambda
  stores$sigma_store[[i]] = mc$sigma
  stores$eta_store[[i]] = mc$eta
  stores$phi_store[[i]] = mc$phi
  stores$delta_store[[i]] = mc$delta
  stores$tauh_store[[i]] = mc$tauh
  # stores$omega_store[[i]] = (mc$Lambda %*% t(mc$Lambda) + diag(mc$sigma)) * scaleMat
  stores$omega_store[[i]] = mc$Lambda %*% t(mc$Lambda) + diag(mc$sigma)
  stores$k_store[[i]] = mc$k
  
  stores
}


