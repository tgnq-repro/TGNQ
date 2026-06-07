# =============================================================================
# additive_parameter.R
# -----------------------------------------------------------------------------
# Alternative parameter functions for the ADDITIVE TGNQ model
#
#   theta_{g h}(tau) = alpha_g(tau) + beta_h(tau).
#
# When the simulation flag `additive == TRUE`, the simulation driver sources
# this file *after* `generate_data.R`, so that the definitions below OVERRIDE
# alpha.func / beta.func / nu.func / gamma.func of the baseline DGP.
#
# Identification restrictions used here (see Remark 2 of the paper):
#   * alpha_1(tau) = 0 for all tau     -> fix the row-group baseline.
# =============================================================================

# -----------------------------------------------------------------------------
# alpha.func (row-group / susceptibility effect)
#   g = 1 : baseline group, fixed at 0 for identifiability
#   g = 2,3 : monotone functions of the latent quantile noise u
# -----------------------------------------------------------------------------
alpha.func<-function(u, g){
  if (g == 1)
    return(0)
  if (g == 2)
    return((pt(u+1,df=4)*0.1+0.05))
  if (g == 3 )
    return((pt(u-1,df=4)*0.1+0.2))
}

# -----------------------------------------------------------------------------
# beta.func (column-group / influence effect)
#   h = 2 chosen as a "zero" reference, so that h = 1 and h = 3 are
#   clearly distinguished from each other in the additive structure.
# -----------------------------------------------------------------------------
beta.func<-function(u, h){
  if (h == 1)
    return((punif(u,-12,12)/12+0.2)*2)
  if (h == 3)
    return(((punif(u,-12,12)/3)))
  if (h==2)
    return(u*0)
}

# -----------------------------------------------------------------------------
# nu.func (group-specific autoregressive effect)
# -----------------------------------------------------------------------------
nu.func<-function(u, g, G){
  if (g == 1)
    return(pnorm(u-0.5)*0.2)
  if (g == 2)
    return(pnorm(u)*0.2)
  if (g == 3)
    return(pnorm(u+0.5)*0.2)
}

# -----------------------------------------------------------------------------
# gamma.func (group-specific covariate effects for the p-th covariate)
# -----------------------------------------------------------------------------
gamma.func<-function(u, g, p){
  if (g == 1 & p==2)
    return(pnorm(u)*0.2)
  if (g == 2& p==2)
    return(pnorm(u+0.4)*0.2)
  if (g == 3&p==2)
    return(pnorm(u-0.2)*0.2)
  
  if (g == 1 & p==3)
    return(pt(u,df=2)*0.2)
  if (g == 2& p==3)
    return(pt(u+0.4,df=2)*0.2)
  if (g == 3&p==3)
    return(pt(u-0.2,df=2)*0.2)
}
