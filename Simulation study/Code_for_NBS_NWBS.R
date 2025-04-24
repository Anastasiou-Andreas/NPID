## We also include the method explained in the new paper by Wang, Yu, Rinaldo.

## Utils from that paper
sizeToCPS = function(mCi, dist) {
  # function adapted from Alex Reinhart's code
  # Convert a source size in mCi (miliCuries) to an approximate number of counts
  # at dist (meters), using calibrations with Cs-137 sources of known sizes.
  
  knownSize = 0.000844 # mCi
  knownDist = 0.05 # meters
  knownCounts = 630 # cps
  mu = 0.0100029 # attenuation coefficient, in units of m^{-1}, for 660 keV in air
  
  (mCi / knownSize) * knownCounts * (knownDist / dist)^2 * exp(-mu * (knownDist + dist))
}


inject_source = function(background_spectrum, background_cps, source_spectrum, source_size, distance) {
  # Creates a new spectrum based on a combination of the background and source spectra
  # background_spectrum and source_spectrum: spectral densities, i.e. normalized to 1
  # background_cps: how many counts per second (on average) does the background contribute to the detector?
  # source_size: strength of source in mCi
  # distance: distance to source in meters
  source_cps = sizeToCPS(source_size, distance)
  # take linear combination of background and source densities, weighted by CPS
  new_spectrum = background_cps*background_spectrum + source_cps*source_spectrum
  #plot(source_spectrum)
  return(	new_spectrum)
}


cumulative_frac = function(x) cumsum(x)/sum(x)

## Utils_functions from that paper
Delta_se_t = function(y,s,e,t,N)
{
  #T =   dim(y)[2]
  n =  dim(y)[2]
  
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te =sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux =  as.vector(y[s:t,])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y =  as.vector(y[s:e,])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st =  temp(vec_y)# temp(grid)
  
  aux = y[(t+1):e,]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te =  temp(vec_y)# temp(grid)
  
  temp =  sqrt( n_st*n_te / n_se   ) *max(abs(Fhat_te - Fhat_st  ))
  
  return(temp )
}

compute_F0 =   function(grid, p0 )
{
  cdf  = cumsum(p0)  
  loc  = (1:length(p0))/length(p0)
  F0 =  rep(0, length(grid))
  
  for(j in 1:length(grid))
  {
    ind =  which.min(abs(grid[j] - loc   )) 
    F0[j] =   cdf[ind[1]]
  }
  
  return(F0)
}

######################################

dist_change_points =  function(Shat,S0)
{
  
  if(length(Shat)==0)
  {return(Inf)}
  
  if(length(S0)==0)
  {return(-Inf)}
  
  temp =rep(0,length(S0))
  for(j in 1:length(S0))
  {
    temp[j] = min(abs(S0[j] - Shat))
  }
  return( max(temp) )
}

#############################################################3

#############################################################################
########################################################################################

new_BS = function(y, gam=1,s,e,flag=0,S = NULL,Dval=NULL,pos=1, N)
{
  if(e-s <   1+  3*gam  ||  flag ==1)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = c(pos,NULL)
    return(S);
  }
  else{
    e = e- gam
    s = s + gam 
    ###  calculate  statistics
    a =  rep(0,e-s+1  )
    
    for(t  in    (s+1):(e-1)  )
    {
      
      a[t-s ] =  Delta_se_t(y,s,e,t,N)
    }
    
    
    best_t  =  which.max( a)
    
    #  if(a[best_t]  <tau  )
    #  {
    #3   return(S)
    #  }
    #   print(s+best_t)
    best_t =   s+ best_t 
    pos1 =  pos
    pos2 = pos
    pos1[length(pos)] = pos[length(pos)]+1
    pos2[length(pos)] = pos[length(pos)]+1
    temp1 = new_BS(y,gam,s,best_t-1,flag,S,Dval,pos1,N)
    temp2 = new_BS(y, gam,best_t+1,e,flag, S,Dval,pos2,N)
    S1 = temp1$S
    Dval1 = temp1$Dval    
    pos1 = temp1$pos 
    S2 = temp2$S
    Dval2 = temp2$Dval 
    pos2 = temp2$pos 
    
    #pos = c(pos,pos[length(pos)]+1) 
    S =  c(S,best_t)
    Dval = c(Dval,a[best_t-s])
    
    
    S  =   c(S,S1,S2)
    Dval =  c(Dval,Dval1,Dval2)
    pos = c(pos,pos1,pos2)
    return(list(S=S,Dval = Dval,pos=pos))
  }
}

new_BS_threshold =  function(temp,tau,p)
{
  ind = which(temp$Dval >  tau)
  
  Shat =  c()
  
  if(length(ind)==0)
  {
    return(NULL)
  }
  
  for( j in 1:length(ind))
  {
    if(p[ind[j]]==0)
    {
      Shat = c(Shat,temp$S[ind[j]])
    }
    if(p[ind[j]] > 0  && min(abs(Shat - temp$S[p[ind[j]]] ))==0   )
    {
      Shat = c(Shat,temp$S[ind[j]])
    }
  }
  
  return(Shat)
}
parent =  function(temp)
{
  p= rep(1,length(temp$pos))
  p[1] = 0
  for(i  in 2:length(temp$pos))
  {
    ind = which(temp$pos[1:(i-1)] == (temp$pos[i]-1))
    p[i] =  ind[length(ind)]
  }
  return(p)
}


#######################################3

new_WBS =  function(y, gam=1,s,e,flag=0,S=NULL,Dval=NULL,pos=1,alpha,beta,N)
{
  
  
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  # 
  ind = which( beta_new- alpha_new > 1+  3*gam)
  #  alpha_new =  alpha_new[ind]h( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind] + gam
  beta_new=   beta_new[ind] - gam
  M =  length(alpha_new)
  
  xi = 1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  # print(S)
  
  #  beta_new=   beta_new[ind]
  # # 
  #  ind = whic
  if(M ==  0 ||  pos[length(pos)]>12)
  {
    S =  c(S,NULL)
    Dval =  c(Dval, NULL)
    pos  = NULL#c(pos,NULL)
    return(list(S=S,pos = pos,Dval = Dval))
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
  #  print()
  for( m in 1:M  )
  {
    temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
    for(t  in    (alpha_new[m]+1):(beta_new[m]-1 )  )
    {
      temp[t-(alpha_new[m]) ] =  Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
    }
    best_ind  =  which.max(temp)
    a[m] =  alpha_new[m] +  best_ind
    b[m] =   temp[best_ind]
  }
  best_ind =  which.max(b)
  
  
  #  print(b)
  #  if(b[best_ind]  <tau  )
  # {
  #  return(S)
  #}
  # print(S)
  #best_t =   s+ best_t 
  pos1 =  pos
  pos2 = pos
  pos1[length(pos)] = pos[length(pos)]+1
  pos2[length(pos)] = pos[length(pos)]+1
  
  #print(c(a[best_ind],b[best_ind]))
  temp2 = new_WBS(y,gam,a[best_ind]+1,e,flag, S,Dval,pos2,alpha, beta,N)
  temp1 = new_WBS(y,gam,s,a[best_ind]-1,flag, S,Dval,pos1,alpha, beta,N)
  S1 = temp1$S 
  Dval1 = temp1$Dval     
  pos1 = temp1$pos  
  S2 = temp2$S 
  Dval2 = temp2$Dval 
  pos2 = temp2$pos 
  
  S =  c(S,a[best_ind])
  Dval = c(Dval,b[best_ind])
  
  
  S  =   c(S,S1,S2)
  Dval =  c(Dval,Dval1,Dval2)
  pos = c(pos,pos1,pos2)
  
  return(list(S=S,Dval = Dval,pos=pos))
  # S  =   c(S,S1,S2)
  # 
  # 
  # return(S)
  
}

#######################



#############################################################3
################################################################


NBS_full =  function(y,z,gam,N)
{
  n = max(N)
  T = dim(y)[1]
  
  temp1 = new_BS(z, gam,1,T,0,NULL,NULL,1, N)  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:min(20,length(Dval))]-10^{-5})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    if(length(aux)==0)
    {
      break;
    }
    S[[j]] = sort(aux)
  }
  #}
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }  
  lamb =log(sum(N))/1.5#2.5#1.5#2#2.555#
  for(j in 1:length(S))#)
  {
    if(length(S[[j]])==0)
    {
      j = j+1;
      
      if(j>length(S))
        break;
    }
    
    B2  =  S[[j]]
    if(j==length(S))
    {
      B1 = NULL
    }
    if(j< length(S))
    {
      B1 = S[[j+1]]
    }
    temp = setdiff(B2,B1)
    
    st =  -10^15
    #Delta_se_t(z,eta1,eta2,eta,N)^2 
    for(l in 1:length(temp))
    {
      eta =  temp[l]
      
      if(j == length(S))
      {
        eta1 = 1
        eta2 = T
      }
      if(j < length(S))
      {
        for(k in 1:length(S[[j+1]]))
        {
          if(S[[j+1]][k]> eta  )
            break;
        }
        if(S[[j+1]][k]> eta )
        {
          eta2 = S[[j+1]][k]
          
          if(k ==1)
            eta1 = 1
          
          if(k > 1)
            eta1 = S[[j+1]][k-1]+1
        }
        if(S[[j+1]][k]< eta )
        {
          eta1 = S[[j+1]][k]+1
          eta2 = T
        }
      }
      st_aux = Delta_se_t(y,eta1,eta2,eta,N)^2
      # print(st_au)
      if(st_aux> st)
      {
        st = st_aux
      }
    }###  close for defining  eta1 and eta2
    
    
    # print(c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb)
    if(st >   lamb)
    {
      #B1 = B2
      return(B2)
    }
    # print(st)
  }
  #c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb
  
  return(B1)
  # 
  
}
#############################################3

arg_max_Delta_se_t = function(y,s,e,t,N)
{
  #T =   dim(y)[2]
  n =  dim(y)[2]
  
  n_st = sum(N[s:t])  #n*(t-s+1)
  n_se = sum(N[s:e])  #n*(e-s+1)
  n_te =sum(N[(t+1):e]) #n*(e-(t+1) +1)
  
  aux =  as.vector(y[s:t,])
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  vec_y =  as.vector(y[s:e,])
  vec_y = vec_y[which(is.na(vec_y)==FALSE)]
  Fhat_st =  temp(vec_y)# temp(grid)
  
  aux = y[(t+1):e,]
  aux = aux[which(is.na(aux)==FALSE)]
  temp = ecdf(aux)
  Fhat_te =  temp(vec_y)# temp(grid)
  
  #temp =  sqrt( n_st*n_te / n_se   ) *max(abs(Fhat_te - Fhat_st  ))
  ind =  which.max(abs(Fhat_te - Fhat_st  ))
  z_hat = min(vec_y[ ind])
  #print(min(abs(Fhat_te - Fhat_st  )))
  return(z_hat)
}
################################################3



NWBS_full =  function(y,z,gam,N,alpha,beta)
{
  
  n = max(N)
  T = dim(y)[1]
  #grid = 
  temp1 = new_WBS(z, gam,1,T,0,NULL,NULL,1,alpha,beta,N)
  
  Dval = temp1$Dval
  p1 =  parent(temp1)
  aux = sort(Dval,decreasing = TRUE)
  tau_grid = rev(aux[1:min(20,length(Dval))]-10^{-4})
  tau_grid = c(tau_grid,10)
  
  S =  c()
  for( j in 1:length(tau_grid))
  {
    aux = new_BS_threshold(temp1,tau_grid[j],p1)
    
    S[[j]] = sort(aux)
  }
  
  T= dim(y)[1]
  S = unique(S)
  if(length(S)==0)
  {
    return(NULL)
  }  
  
  lamb = log(sum(N))/1.5
  for(j in 1:(length(S)))
  {
    if(length(S[[j]]) == 0)
    {
      j = 2
    }
    B2  =  S[[j]]
    if(j < length(S))
    {
      B1 = S[[j+1]]
    }
    if(j  == length(S))
    {
      B1 = NULL
    }
    temp = setdiff(B2,B1)
    
    st =  -10^15
    #Delta_se_t(z,eta1,eta2,eta,N)^2 
    count = 0
    for(l in 1:length(temp))
    {
      eta =  temp[l]
      
      if( length(B1)==0)
      {
        eta1 = 1
        eta2 = T
      }
      if( length(B1)>0)
      {
        for(k in 1:length(B1))
        {
          if(B1[k]> eta  )
            break;
        }
        if(B1[k]> eta )
        {
          eta2 = B1[k]
          
          if(k ==1)
            eta1 = 1
          
          if(k > 1)
            eta1 = B1[k-1]+1
        }
        if(B1[k]< eta )
        {
          eta1 = B1[k]+1
          eta2 = T
        }
      }
      if( length(temp) > 1     )
      {
        st_aux = max(wbs_Delta_se_t(y,eta1,eta2,eta,N,alpha,beta)^2)  
        if(st_aux == 0)
        {
          st_aux = Delta_se_t(y,eta1,eta2,eta,N)^2
        }
      }
      if(length(temp)==1){st_aux = max(Delta_se_t(y,eta1,eta2,eta,N)^2)  }
      
      if(st_aux> st)
      {
        st = st_aux
        #    B1 = B2
      }
    }  
    
    if(st >  lamb)#  st > lamb
    {
      # B1 = B2
      return(B2)
    }
    #     print(st)
  }#
  #c1 - c2 - Delta_se_t(z,eta1+1,eta2,eta,N)^2 + lamb
  return(B1)
  # 
  
}

wbs_Delta_se_t  =  function(y,s,e,t,N,alpha, beta)
{
  alpha_new =  pmax(alpha,s)
  beta_new = pmin(beta,e)
  # 
  #  ind = which( beta_new- alpha_new > 1+  3*gam)
  #  alpha_new =  alpha_new[ind]
  #  beta_new=   beta_new[ind]
  # # 
  ind = which( beta_new- alpha_new >20)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  xi = 1/8
  alpha_new2 =alpha_new
  beta_new2  =  beta_new
  alpha_new =   ceiling((1-xi)*alpha_new2+  xi*beta_new2)
  beta_new =    ceiling((1-xi)*beta_new2 +  xi*alpha_new2)
  ind = which( beta_new- alpha_new >1)
  alpha_new =  alpha_new[ind]
  beta_new=   beta_new[ind]
  M =  length(alpha_new)
  
  # print(S)
  
  if(M ==  0)
  {
    return(0) 
  }
  
  b  =  rep(0,M)
  a =  rep(0,M)
  
  #print(beta_new)
  #print(alpha_new)
  #  print()
  
  for( m in 1:M  )
  {
    #  temp  =  rep(0,beta_new[m]-alpha_new[m]+1)
    #  for(t  in    (alpha_new[m]+1):(beta_new[m]-1 )  )
    #  {
    if(alpha_new[m]<t &&  t  < beta_new[m]  )
    {
      b[m] =  Delta_se_t(y,alpha_new[m],beta_new[m],t,N)
    }
    #  }
    # best_ind  =  which.max(temp)
    #  a[m] =  alpha_new[m] +  best_ind
    #  b[m] =   temp[best_ind]
  }
  #best_ind =  which.max(b)
  return(max(b))
  #return(b[best_ind])
}

