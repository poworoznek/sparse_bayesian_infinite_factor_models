# build interactions matrix Psi from strung out vector csi
# see sampler bits

build_Psi = function(csi, k){
  diags = cumsum(rev(1:k)) - rev(1:k)+1
  Psi = diag(k)
  Psi[upper.tri(Psi)] = csi[-diags] / 2
  Psi = Psi + t(Psi)
  diag(Psi) = csi[diags]
  return(Psi)
}
