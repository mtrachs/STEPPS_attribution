########################################################################################################################
# 
########################################################################################################################
library(rstan)
library(VGAM)
library(gtools)


#######################################################################################################################
#
#######################################################################################################################
library(rstan)
library(RColorBrewer)
library(fields)
setwd('~/workflow_stepps_prediction/prediction/')
help.fun.loc <- 'utils/'
data.loc <- 'data/'
plot.loc <- 'plots/'
output.loc <-'/media/mathias/Seagate Expansion Drive/workflow_stepps_prediction/prediction/output/ensemble/' #'output/'

#z <- readRDS(paste0(data.loc,'z.RDS'))

#########################################################################################################################
source(paste(help.fun.loc,'pred_helper_funs.r',sep=''))
source(paste(help.fun.loc,'pred_plot_funs.r',sep=''))
########################################################################################################################

setwd('~/workflow_stepps_prediction/attribution/')
data.loc <- 'data/'
output.loc <- 'output/'
plot.loc <- 'plots/'

load(paste0(data.loc,'attribution.RData'))


model_names <- c('coord','coord_moisture','coord_climate','coord_climate_sq','coord_climate_autocor','coord_autocor',
                 'autocor','climate_autocor','climate','intercept')



for(mod_nam in model_names){




  fit <- read_stan_csv(paste(output.loc,'test_attribution_',mod_nam,'.csv',sep=''))
  #saveRDS(fit,paste0(data.loc,'fit_attribution_model_coord_mositure.RDS'))
  #fit <- readRDS(paste0(data.loc,'fit_attribution_model_coord_climate.RDS'))
  post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
  rm(fit)
  var_names<- colnames(post[,1,])
  par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
  nrun <- rowSums(post)!=0
  post <- array(dim=c(sum(nrun),1,dim(post)[3]),post[nrun,,])
  colnames(post[,1,]) <- var_names

  eta <- post[,,grep('eta',par_names)]
  rho <- post[,,grep('rho',par_names)]

  alpha_zero <- post[,,grep('alpha_zero',par_names)]
  alpha_one <- post[,,grep('alpha_one',par_names)]
  alpha_two <- post[,,grep('alpha_two',par_names)]
  alpha_three <- post[,,grep('alpha_three',par_names)]
  alpha_four <- post[,,grep('alpha_four',par_names)]
  alpha_five <- post[,,grep('alpha_five',par_names)]
  alpha_six <- post[,,grep('alpha_six',par_names)]
  mu_s_k <- post[,,grep('mu_s_k',par_names)]
  eps_s_t_k <- post[,,grep('eps_s_t_k',par_names)]


# 

# sapply(1:ncol(eta),function(x){
#   hist(eta[,x],main=taxa[x])
# })
# 
# sapply(1:ncol(rho),function(x){
#   hist(rho[,x],main=taxa[x])
# })

  eta_mean <- colMeans(eta)
  rho_mean <- colMeans(rho)
  alpha_zero_mean <- colMeans(alpha_zero)
  alpha_one_mean <- colMeans(alpha_one)
  alpha_two_mean <- colMeans(alpha_two)
  alpha_three_mean <- colMeans(alpha_three)
  alpha_four_mean <- colMeans(alpha_four)
  alpha_five_mean <- colMeans(alpha_five)
  alpha_six_mean <- colMeans(alpha_six)
  mu_s_k_mean <- colMeans(mu_s_k) ## K*S t1,spec1,t1,spec2,t1spe3,t1specK,
  eps_s_t_k_mean <- colMeans(eps_s_t_k) # S1,t1,k1, S2,t1,k1,S50,t1,k1;,S1,t2,k1,;S2 

#########################################################################################################################
  # if coords are used in this model then run this part
  #######################################################################################################################
  if(length(grep('autocor',mod_nam))>0){
    spat_autocor_resid <-  matrix(ncol=K,nrow=S,mu_s_k_mean,byrow=TRUE)
  }
  
  system.part.time <- lapply(1:S,function(s) matrix(ncol=K,nrow=time_slice))
  
  if(length(grep('coord',mod_nam))>0){
    
    system.part.space <- sapply(1:K,function(k){
      alpha_zero_mean[k] +
        alpha_one_mean[k]*coord_x +
        alpha_two_mean[k]*coord_y
    })
  

    


     for(s in 1:S){
      for (k in 1:K){
       
        if(length(alpha_three_mean)==0){
          system.part.time[[s]][,k] <- rep(system.part.space[s,k],time_slice)
        }
        
        if(length(grep('moisture',mod_nam))>0){
          system.part.time[[s]][,k] <- rep(system.part.space[s,k],time_slice) +
              alpha_three_mean[k]*M 
        }
        if((length(grep('climate',mod_nam))>0)&(length(alpha_five_mean)==0)){
          system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
           alpha_three_mean[k]*T + 
          alpha_four_mean[k]*M
        }
        if((length(grep('climate',mod_nam))>0)&(length(alpha_five_mean)>0)){
          system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
           alpha_three_mean[k]*T +
          alpha_four_mean[k]*M +
           alpha_five_mean[k]*T2+
          alpha_six_mean[k]*M2
        }
        if((length(grep('autocor',mod_nam))>0)&&(length(alpha_three_mean)>0)){
          system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
          alpha_three_mean[k]*T +
          alpha_four_mean[k]*M +
          rep(spat_autocor_resid[s,k],time_slice)
        }
        if((length(grep('autocor',mod_nam))>0)&(length(alpha_three_mean)==0)){
          system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
            rep(spat_autocor_resid[s,k],time_slice)
        }

        # +alpha_four_mean[k]*M

        #+ rep(spat_autocor_resid[s,k],time_slice)
      }
    }
  }
  
  
  
  
#####################################################################################################################
  # analyse models that do not include coordinates 'climate_autocor','climate'
  ###################################################################################################################
  if(length(grep('coord',mod_nam))==0){
   for(s in 1:S){
      for(k in 1:K){
        if(length(grep('autocor',mod_nam))==0){
          system.part.time[[s]][,k] <- alpha_zero_mean[k]+
                    alpha_one_mean[k]*T +
                    alpha_two_mean[k]*M
        }
        if((length(grep('autocor',mod_nam))==0)&&(length(alpha_three_mean)==0)){
          system.part.time[[s]][,k] <- rep(alpha_zero_mean[k],time_slice)
        }
        if((length(grep('autocor',mod_nam))>0)&(length(alpha_three_mean)>0)){
          system.part.time[[s]][,k] <- alpha_zero_mean[k]+
            alpha_three_mean[k]*T +
            alpha_four_mean[k]*M +
            rep(spat_autocor_resid[s,k],time_slice)
        }
        if((length(grep('autocor',mod_nam))>0)&(length(alpha_three_mean)==0)){
          system.part.time[[s]][,k] <- alpha_zero_mean[k]+ 
                                      rep(spat_autocor_resid[s,k],time_slice)
      }
    }
   }
  }  
 # saveRDS(system.part.time,paste0(data.loc,'alpha_coord_climate.RDS'))

#####################################################################################################################
# likelihood for each site and time
#####################################################################################################################

log_likelihood <- 
  sapply(1:S,function(x){
    sum(log(ddirichlet(t(r_mean_array[x,,2:(time_slice+1)]),exp(system.part.time[[x]]))))
  })

#saveRDS(log_likelihood_coord,paste0(data.loc,'log_likelihood_coord.RDS'))
saveRDS(log_likelihood,paste0(data.loc,'log_likelihood_',mod_nam,'.RDS'))

}

#hm...
# log_likelihood_ideal <- 
#   sapply(1:S,function(x){
#     sum(log(ddirichlet(t(r_mean_array[x,,2:(time_slice+1)]),t(r_mean_array[x,,2:(time_slice+1)]))))
#   })
#####################################################################################################################
#####################################################################################################################
log_likelihood_coord_climate_sq <- readRDS(paste0(data.loc,'log_likelihood_coord_climate_sq.RDS'))
log_likelihood_coord_climate <- readRDS(paste0(data.loc,'log_likelihood_coord_climate.RDS'))
log_likelihood_coord_moisture <- readRDS(paste0(data.loc,'log_likelihood_coord_moisture.RDS'))
log_likelihood_coord <- readRDS(paste0(data.loc,'log_likelihood_coord.RDS'))
log_likelihood_climate <- readRDS(paste0(data.loc,'log_likelihood_climate.RDS'))
log_likelihood_climate_autocor <- readRDS(paste0(data.loc,'log_likelihood_climate_autocor.RDS'))
log_likelihood_coord_climate_autocor <- readRDS(paste0(data.loc,'log_likelihood_coord_climate_autocor.RDS'))
log_likelihood_coord_autocor <- readRDS(paste0(data.loc,'log_likelihood_coord_autocor.RDS'))
log_likelihood_autocor <- readRDS(paste0(data.loc,'log_likelihood_autocor.RDS'))
log_likelihood_intercept <- readRDS(paste0(data.loc,'log_likelihood_intercept.RDS'))


 comp_likelihood <- cbind(log_likelihood_intercept,
                          log_likelihood_coord,
                          log_likelihood_coord_moisture,
                          log_likelihood_coord_climate,
                          log_likelihood_coord_climate_sq,
                          log_likelihood_coord_climate_autocor,
                          log_likelihood_coord_autocor,
                          log_likelihood_autocor,
                          log_likelihood_climate,
                          log_likelihood_climate_autocor)
 
 comp_likelihood <- as.data.frame(comp_likelihood)
 colMeans(comp_likelihood)
 
 diff_clim_autocor <- comp_likelihood$log_likelihood_climate_autocor -comp_likelihood$log_likelihood_autocor
 hist(diff_clim_autocor)
 
 
 
 
 # log_likelihood_coord_clim-log_likelihood_coord_clim_sq
write.table(round(comp_likelihood),paste(data.loc,'comp_likelihood.txt'),row.names = FALSE)
write.csv(round(comp_likelihood),paste(data.loc,'comp_likelihood.csv'),row.names = FALSE)



#########################################################################################################################
# predicted composition
#########################################################################################################################
predicted_comp <-   sapply(1:length(model_names), function(x) list())  
  
for(i in 1:length(model_names)){
    
    mod_nam <- model_names[i]
    
    
    fit <- read_stan_csv(paste(output.loc,'test_attribution_',mod_nam,'.csv',sep=''))
    #saveRDS(fit,paste0(data.loc,'fit_attribution_model_coord_mositure.RDS'))
    #fit <- readRDS(paste0(data.loc,'fit_attribution_model_coord_climate.RDS'))
    post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
    rm(fit)
    var_names<- colnames(post[,1,])
    par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
    nrun <- rowSums(post)!=0
    post <- array(dim=c(sum(nrun),1,dim(post)[3]),post[nrun,,])
    colnames(post[,1,]) <- var_names
    
    eta <- post[,,grep('eta',par_names)]
    rho <- post[,,grep('rho',par_names)]
    
    alpha_zero <- post[,,grep('alpha_zero',par_names)]
    alpha_one <- post[,,grep('alpha_one',par_names)]
    alpha_two <- post[,,grep('alpha_two',par_names)]
    alpha_three <- post[,,grep('alpha_three',par_names)]
    alpha_four <- post[,,grep('alpha_four',par_names)]
    alpha_five <- post[,,grep('alpha_five',par_names)]
    alpha_six <- post[,,grep('alpha_six',par_names)]
    mu_s_k <- post[,,grep('mu_s_k',par_names)]
    eps_s_t_k <- post[,,grep('eps_s_t_k',par_names)]
    
    
    # 
    
    # sapply(1:ncol(eta),function(x){
    #   hist(eta[,x],main=taxa[x])
    # })
    # 
    # sapply(1:ncol(rho),function(x){
    #   hist(rho[,x],main=taxa[x])
    # })
    
    eta_mean <- colMeans(eta)
    rho_mean <- colMeans(rho)
    alpha_zero_mean <- colMeans(alpha_zero)
    alpha_one_mean <- colMeans(alpha_one)
    alpha_two_mean <- colMeans(alpha_two)
    alpha_three_mean <- colMeans(alpha_three)
    alpha_four_mean <- colMeans(alpha_four)
    alpha_five_mean <- colMeans(alpha_five)
    alpha_six_mean <- colMeans(alpha_six)
    mu_s_k_mean <- colMeans(mu_s_k) ## K*S t1,spec1,t1,spec2,t1spe3,t1specK,
    eps_s_t_k_mean <- colMeans(eps_s_t_k) # S1,t1,k1, S2,t1,k1,S50,t1,k1;,S1,t2,k1,;S2 
    
    #########################################################################################################################
    # if coords are used in this model then run this part
    #######################################################################################################################
    if(length(grep('autocor',mod_nam))>0){
      spat_autocor_resid <-  matrix(ncol=K,nrow=S,mu_s_k_mean,byrow=TRUE)
    }
    
    system.part.time <- lapply(1:S,function(s) matrix(ncol=K,nrow=time_slice))
    
    if(length(grep('coord',mod_nam))>0){
      
      system.part.space <- sapply(1:K,function(k){
        alpha_zero_mean[k] +
          alpha_one_mean[k]*coord_x +
          alpha_two_mean[k]*coord_y
      })
      
      
      
      
      
      for(s in 1:S){
        for (k in 1:K){
          
          if(length(alpha_three_mean)==0){
            system.part.time[[s]][,k] <- rep(system.part.space[s,k],time_slice)
          }
          
          if(length(grep('moisture',mod_nam))>0){
            system.part.time[[s]][,k] <- rep(system.part.space[s,k],time_slice) +
            alpha_three_mean[k]*M 
          }
          if((length(grep('climate',mod_nam))>0)&(length(alpha_five_mean)==0)){
            system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
              alpha_three_mean[k]*T + 
              alpha_four_mean[k]*M
          }
          if((length(grep('climate',mod_nam))>0)&(length(alpha_five_mean)>0)){
            system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
            alpha_three_mean[k]*T +
            alpha_four_mean[k]*M +
            alpha_five_mean[k]*T2+
            alpha_six_mean[k]*M2
          }
          if(length(grep('autocor',mod_nam))>0){
            system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
            alpha_three_mean[k]*T +
            alpha_four_mean[k]*M +
            rep(spat_autocor_resid[s,k],time_slice)
          }
          if((length(grep('autocor',mod_nam))>0)&(length(alpha_three_mean)==0)){
            system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
              rep(spat_autocor_resid[s,k],time_slice)
          }
          
          # +alpha_four_mean[k]*M
          
          #+ rep(spat_autocor_resid[s,k],time_slice)
        }
        system.part.time[[s]] <- exp(system.part.time[[s]])/rowSums(exp(system.part.time[[s]]))
      }
    }
    
    
    
    
    #####################################################################################################################
    # analyse models that do not include coordinates 'climate_autocor','climate'
    ###################################################################################################################
    if(length(grep('coord',mod_nam))==0){
      for(s in 1:S){
        for(k in 1:K){
          if(length(grep('autocor',mod_nam))==0){
            system.part.time[[s]][,k] <- alpha_zero_mean[k]+
              alpha_one_mean[k]*T +
              alpha_two_mean[k]*M
          }
          if((length(grep('autocor',mod_nam))==0)&&(length(alpha_three_mean)==0)){
            system.part.time[[s]][,k] <- rep(alpha_zero_mean[k],time_slice)
          }
          if((length(grep('autocor',mod_nam))>0)&(length(alpha_three_mean)>0)){
            system.part.time[[s]][,k] <- alpha_zero_mean[k]+
              alpha_three_mean[k]*T +
              alpha_four_mean[k]*M +
              rep(spat_autocor_resid[s,k],time_slice)
          }
          if((length(grep('autocor',mod_nam))>0)&(length(alpha_three_mean)==0)){
            system.part.time[[s]][,k] <- alpha_zero_mean[k]+ 
              rep(spat_autocor_resid[s,k],time_slice)
          }
        }
        system.part.time[[s]] <- exp(system.part.time[[s]])/rowSums(exp(system.part.time[[s]]))
      }
    }  
    # saveRDS(system.part.time,paste0(data.loc,'alpha_coord_climate.RDS'))
    
    #####################################################################################################################
    # likelihood for each site and time
    #####################################################################################################################
    predicted_comp[[i]] <- system.part.time

  }   
   
   

 
 


# quadratic has lower likelihood: could potentially mean that mean of parameter estimates is 
# a really bad estimate for the overall effect of a variable as individual draws of the parameter are highly correlated
# could lokk at the individual draws to see if the likelihood gets better 

####################################################################################################################
# ensemble predictions
####################################################################################################################
fit <- read_stan_csv(paste(output.loc,'test_attribution_coord_climate_autocor.csv',sep=''))
post <- rstan::extract(fit, permuted=FALSE, inc_warmup=FALSE)
rm(fit)
var_names<- colnames(post[,1,])
par_names <- sapply(var_names, function(x) unlist(strsplit(x,'[[]'))[1])
nrun <- rowSums(post)!=0
post <- array(dim=c(sum(nrun),1,dim(post)[3]),post[nrun,,])
colnames(post[,1,]) <- var_names

eta <- post[,,grep('eta',par_names)]
rho <- post[,,grep('rho',par_names)]

alpha_zero <- post[,,grep('alpha_zero',par_names)]
alpha_one <- post[,,grep('alpha_one',par_names)]
alpha_two <- post[,,grep('alpha_two',par_names)]
alpha_three <- post[,,grep('alpha_three',par_names)]
alpha_four <- post[,,grep('alpha_four',par_names)]
mu_s_k <- post[,,grep('mu_s_k',par_names)]

# eta_mean <- colMeans(eta)
# rho_mean <- colMeans(rho)
# alpha_zero_mean <- colMeans(alpha_zero)
# alpha_one_mean <- colMeans(alpha_one)
# alpha_two_mean <- colMeans(alpha_two)
# alpha_three_mean <- colMeans(alpha_three)
# alpha_four_mean <- colMeans(alpha_four)
# alpha_five_mean <- colMeans(alpha_five)
# alpha_six_mean <- colMeans(alpha_six)
# mu_s_k_mean <- colMeans(mu_s_k) ## K*S t1,spec1,t1,spec2,t1spe3,t1specK,
# eps_s_t_k_mean <- colMeans(eps_s_t_k) # S1,t1,k1, S2,t1,k1,S50,t1,k1;,S1,t2,k1,;S2 

# system.part.time[[s]][,k] <-rep(system.part.space[s,k],time_slice) +
#   alpha_three_mean[k]*T +
#   alpha_four_mean[k]*M +
#   rep(spat_autocor_resid[s,k],time_slice)

spat_autocor_resid_ens  <-  lapply(1:nrow(mu_s_k),function(e) matrix(ncol=K,nrow=S,mu_s_k[e,],byrow=TRUE))

system.part.space.ensemble <- lapply(1:nrow(alpha_zero),function(e){ 
  sapply(1:K,function(k){
    alpha_zero[e,k] + alpha_one[e,k]*coord_x + alpha_two[e,k]*coord_y
  })
})


system.part.time.ensemble <- lapply(1:S,function(s){
  lapply(1:nrow(alpha_zero),function(e){
    sapply(1:K,function(k){
      rep(system.part.space.ensemble[[e]][s,k],time_slice) +
      alpha_three[e,k]*T +
      alpha_four[e,k]*M  +
      rep(spat_autocor_resid_ens[[e]][s,k],time_slice)
    
    })
  })
})




#saveRDS(system.part.time.ensemble,paste0(data.loc,'ensemble_coord_climate.RDS'))
saveRDS(system.part.time.ensemble,paste0(data.loc,'ensemble_coord_climate_autocor.RDS'))
#####################################################################################################################
#
#####################################################################################################################
ensemble_log_likelihood <- 
  sapply(1:S,function(x){
    log_likelihood_ensemble <- sapply(1:nrow(alpha_zero),function(e){
        sum(log(ddirichlet(t(r_mean_array[x,,2:(time_slice+1)]),exp(system.part.time.ensemble[[x]][[e]]))))
    })
    mean_lle <- mean(log_likelihood_ensemble)
    list(mean_lle =mean_lle,lle =log_likelihood_ensemble)
  })

#spatially autocorrelated residual gives skill the rest is nonsense


cbind(unlist(ensemble_log_likelihood[1,]), comp_likelihood$log_likelihood_coord_climate_autocor)




#######################################################################################################################
# predict compositions and make residual maps
#######################################################################################################################


residual_pred <- lapply(1:length(predicted_comp),function(z){ 
  lapply(1:length(predicted_comp[[1]]),function(x){
    t(r_mean_array[x,,2:32])-predicted_comp[[z]][[x]] 
  })
})


for(zz in 1:length(residual_pred)){
  residual_pred1 <- rep(NA,K) 

  for(i in 1:length(residual_pred[[1]])){
    residual_pred1 <- rbind(residual_pred1,residual_pred[[zz]][[i]])
  }

  residual_pred1 <- residual_pred1[-1,]


##########################################################################################################################
  limits <- list(xlims=range(coord.agg.final$east),ylims=range(coord.agg.final$north))

#saveRDS(r_mean,paste(data.loc,'r_mean.RDS',sep=''))

  breaks = c(-1,-0.5,-0.25,-0.15,-0.05, 0, 0.05,0.15, 0.25, 0.50, 1)
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                    function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  
  weird <- cut(residual_pred1,breaks,include.lowest=TRUE)
  wt <- table(weird)
  #breaks <- breaks[-c(which(wt==0))]
  breaklabels = apply(cbind(breaks[1:(length(breaks)-1)], breaks[2:length(breaks)]), 1, 
                      function(r) { sprintf("%0.2f - %0.2f", r[1], r[2]) })
  # 
  taxa =c('Ash','Beech','Birch','Chestnut','Hemlock','Hickory','Maple','Oak','Conifer','Hardwood','Pine','Spruce','Tamarack')

  plot_pred_maps(r_mean=residual_pred1,
               centers = coord.agg.final,
               taxa = taxa,#c('Dec','Con'),#major_taxa,
               t = seq(500,8000,250),
               N = S,
               K = K,
               T = time_slice,
               thresh = breaks,
               limits = limits,
               type = 'test',
               suff = model_names[zz],#'major_taxa',
               save_plots =TRUE,
               fpath = plot.loc,
               height = 15,
               width = 8)

}





