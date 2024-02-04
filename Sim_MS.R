
library(scales)
library(tidyverse)
library(inlabru)
library(INLA)
library(gridExtra)



### simulate species habitats in 1D
quantiles = seq(0,1,.0001)
abundance_target = -dnorm(quantiles,mean=0.25,sd=.15) + dnorm(quantiles,mean=0.75,sd=.15) 
abundance_bycatch =1-abundance_target
abundance_bycatch2 = 1-(dnorm(quantiles,mean=0.45,sd=.25)+dnorm(quantiles,mean=0.55,sd=.1))

plot(quantiles,abundance_target)
plot(quantiles,abundance_bycatch)
plot(quantiles,abundance_bycatch2)


cov_effect = data.frame(covariate_space = quantiles*100,
                        biomass_target=rescale(abundance_target,to=c(-1,1)) + 15 + rnorm(length(quantiles)), #
                        TS = rescale(abundance_target,to=c(0.1,1)))
mean_biomass_target = mean(cov_effect$biomass_target)
cov_effect = cov_effect %>% mutate(Random = .2,
                                   biomass_bycatch = rescale(abundance_bycatch,to=c(-1,1)) + 5 + rnorm(length(quantiles)), #
                                   CBS = rescale(abundance_bycatch,to=c(0.1,1)),
                                   UBS = rescale(abundance_bycatch2,to=c(0.1,1)),
                                   Proportional = rescale(TS,to=c(0.1,1)),
                                   Squared = TS^2,
                                   "Power 4" = TS^4,
                                   idx=1:n(),
                                   biomass_target_FI = rescale(abundance_target,to=c(-1,1)) + 15 + rnorm(length(quantiles),sd=.25),
                                   biomass_target_FD = rescale(TS,to=c(-1,1)) + 20 + rnorm(length(quantiles),sd=.5),
                                   biomass_bycatch_FI = rescale(abundance_bycatch,to=c(-1,1)) + 15 + rnorm(length(quantiles),sd=.25),
                                   biomass_bycatch_FD = rescale(CBS,to=c(-1,1)) + 20 + rnorm(length(quantiles),sd=.5),
                                   biomass_bycatch2_FI = rescale(abundance_bycatch2,to=c(-1,1)) + 15 + rnorm(length(quantiles),sd=.25),
                                   biomass_bycatch2_FD = rescale(UBS,to=c(-1,1)) + 20 + rnorm(length(quantiles),sd=.5))

##### get data for plotting
abundance_plotD = cov_effect %>% pivot_longer(c(3,6,7),names_to = "Species",values_to = "Species_Abundance")
prob_plotD = cov_effect[seq(1,10001,40),]  %>% 
  pivot_longer(c(4,8:10),names_to = "sampling",values_to = "Preferentiality") 


#### Get mean biomas of each species to assess model performances
mean_biomass_bycatch = mean(cov_effect$biomass_bycatch)
mean_biomass_target_FD = mean(cov_effect$biomass_target_FD)
mean_biomass_byc_FD = mean(cov_effect$biomass_bycatch_FD)
mean_biomass_byc2_FD = mean(cov_effect$biomass_bycatch2_FD)



###########################
### Figure 1 #########
#####################
n = 100 
random_sample = cov_effect[sample(size = n,cov_effect$idx,replace=T),] %>% 
  select(1,2,3,6,11:14) %>% mutate(sampling="Random", Preferentiality = 1.1)
proportional_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$Proportional,replace=T),] %>%  
  select(1,2,3,6,11:14) %>% mutate(sampling="Proportional", Preferentiality = 1.2)
squared_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$Squared,replace=T),] %>%  
  select(1,2,3,6,11:14) %>% mutate(sampling="Squared", Preferentiality = 1.3)
fourth_power_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$"Power 4",replace=T),] %>%
  select(1,2,3,6,11:14) %>% mutate(sampling="Power 4", Preferentiality = 1.4)


sim_data_all = rbind(random_sample,proportional_sample,squared_sample,fourth_power_sample)
sim_data_all = rbind(sim_data_all[,c("covariate_space","sampling","Preferentiality")],
                     prob_plotD[,c("covariate_space","sampling","Preferentiality")])
sim_data_all$sampling = factor(sim_data_all$sampling,levels = c("Random","Proportional","Squared","Power 4"))


png("samp_proba_both_samp_dist.png",height=400,width = 700)
ggplot()+
  geom_line(data=abundance_plotD,aes(x=covariate_space,y=Species_Abundance,color=Species),linewidth=2) + 
  #geom_point(data=prob_plotD,aes(x=covariate_space,y=Preferentiality,shape=sampling),linewidth=1) +
  #ggtitle("Sampling probability bycatch species") +
  theme_bw() + ylab("Standardized values") + 
  geom_point(data=sim_data_all,aes(x=covariate_space,y=Preferentiality ,shape=sampling),alpha=.6,size=1) +
  xlab("1D habitat space") +
  scale_fill_manual(values = c("orange", "purple","green","blue")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4","yellow"),name = "Species abundance")+
  theme(#axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=12),
    legend.text.align = 0)+
  scale_linetype_discrete(name = "Sampling probability")+
  scale_shape_discrete(name = "Preferentiality",
                       #values = c(1, 2, 3, 4),
                       labels = c(expression(p(x) %~% U(0,100)),
                                  expression(p(x) %prop% Y),
                                  expression(p(x) %prop% Y^2),
                                  expression(p(x) %prop% Y^4))
  )
dev.off()






#########################
####### Loop fit different models
#####################

n = 300
n_sims = 50
ptm <- Sys.time()
for(i in 1:n_sims){
  ## simulate sampling for iteration
  random_sample = cov_effect[sample(size = n,cov_effect$idx,replace=T),] %>% 
    select(1,2,3,5,6,12:17) %>% mutate(preferential="Random")
  proportional_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$Proportional,replace=T),] %>%  
    select(1,2,3,5,6,12:17) %>% mutate(preferential="Proportional")
  squared_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$Squared,replace=T),] %>%  
    select(1,2,3,5,6,12:17) %>% mutate(preferential="Squared")
  fourth_power_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$"Power 4" ,replace=T),] %>%
    select(1,2,3,5,6,12:17) %>% mutate(preferential="Power 4")
  
  
  #### 1D latent structure for inferring the habitat
  x = seq(0, 100, by = 2) ### habitat space
  mesh1D <- inla.mesh.1d(x, boundary = "free")
  matern <- inla.spde2.pcmatern(mesh1D,
                                prior.range = c(10, 0.75),
                                prior.sigma = c(0.1, 0.75),
                                constr = TRUE
  )
  
  matern_zero  <- inla.spde2.pcmatern(mesh1D,
                                      prior.range = c(10, 0.75),
                                      prior.sigma = c(0.1, 0.1),
                                      constr = TRUE
  )
  
  #################################
  ### fit conventional models ####
  ##############################
  
  form_target = biomass_target_FD ~ 
    covariate(covariate_space,model=matern) + 
    Intercept(1)
  
  form_bycatch = biomass_bycatch_FD ~ 
    covariate(covariate_space,model=matern) + 
    Intercept(1)
  
  form_bycatch2 = biomass_bycatch2_FD ~ 
    covariate(covariate_space,model=matern) + 
    Intercept(1)
  
  
  #### random sampling ####
  fit_random = bru(components = form_target,
                   data = random_sample,
                   family ="gaussian",
                   options = list(verbose=F))
  
  a = mean(predict(fit_random,cov_effect,~(Intercept+covariate))$mean)
  
  fit_random_byc = bru(components = form_bycatch,
                       data = random_sample,
                       family ="gaussian",
                       options = list(verbose=F))
  
  aa = mean(predict(fit_random_byc,cov_effect,~(Intercept+covariate))$mean)
  
  fit_random_byc2 = bru(components = form_bycatch2,
                        data = random_sample,
                        family ="gaussian",
                        options = list(verbose=F))
  
  aaa = mean(predict(fit_random_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  
  #### Preferential sampling, proportional to abundance intensity ####
  fit_proportional = bru(components = form_target,
                         data = proportional_sample,
                         family ="gaussian",
                         options = list(verbose=F))
  
  b = mean(predict(fit_proportional,cov_effect,~(Intercept+covariate))$mean)
  
  fit_proportional_byc = bru(components = form_bycatch,
                             data = proportional_sample,
                             family ="gaussian",
                             options = list(verbose=F))
  
  bb = mean(predict(fit_proportional_byc,cov_effect,~(Intercept+covariate))$mean)
  
  fit_proportional_byc2 = bru(components = form_bycatch2,
                              data = proportional_sample,
                              family ="gaussian",
                              options = list(verbose=F))
  
  bbb = mean(predict(fit_proportional_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  #### Preferential sampling, squared abundance intensity ####
  fit_squared = bru(components = form_target,
                    data = squared_sample,
                    family ="gaussian",
                    options = list(verbose=F))
  
  c = mean(predict(fit_squared,cov_effect,~(Intercept+covariate))$mean)
  
  fit_squared_byc = bru(components = form_bycatch,
                        data = squared_sample,
                        family ="gaussian",
                        options = list(verbose=F))
  
  cc = mean(predict(fit_squared_byc,cov_effect,~(Intercept+covariate))$mean)
  
  
  fit_squared_byc2 = bru(components = form_bycatch2,
                         data = squared_sample,
                         family ="gaussian",
                         options = list(verbose=F))
  
  ccc = mean(predict(fit_squared_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  #### Preferential sampling, 4th power abundance intensity ####
  fit_fourth_power = bru(components = form_target,
                         data = fourth_power_sample,
                         family ="gaussian",
                         options = list(verbose=F))
  
  d = mean(predict(fit_fourth_power,cov_effect,~(Intercept+covariate))$mean)
  
  fit_fourth_power_byc = bru(components = form_bycatch,
                             data = fourth_power_sample,
                             family ="gaussian",
                             options = list(verbose=F))
  
  dd = mean(predict(fit_fourth_power_byc,cov_effect,~(Intercept+covariate))$mean)
  
  fit_fourth_power_byc2 = bru(components = form_bycatch2,
                              data = fourth_power_sample,
                              family ="gaussian",
                              options = list(verbose=F))
  
  ddd = mean(predict(fit_fourth_power_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  
  if(i==1){
    Estimate_random =a
    Estimate_proportional=b
    Estimate_squared = c
    Estimate_4th_power = d
    
    Estimate_byc_random =aa
    Estimate_byc_proportional=bb
    Estimate_byc_squared = cc
    Estimate_byc_4th_power = dd
    
    Estimate_byc2_random =aaa
    Estimate_byc2_proportional=bbb
    Estimate_byc2_squared = ccc
    Estimate_byc2_4th_power = ddd
  }else{
    Estimate_random=c(Estimate_random,a)
    Estimate_proportional=c(Estimate_proportional,b)
    Estimate_squared =c(Estimate_squared,c)
    Estimate_4th_power = c(Estimate_4th_power,d)
    
    Estimate_byc_random=c(Estimate_byc_random,aa)
    Estimate_byc_proportional=c(Estimate_byc_proportional,bb)
    Estimate_byc_squared =c(Estimate_byc_squared,cc)
    Estimate_byc_4th_power = c(Estimate_byc_4th_power,dd)
    
    Estimate_byc2_random=c(Estimate_byc2_random,aaa)
    Estimate_byc2_proportional=c(Estimate_byc2_proportional,bbb)
    Estimate_byc2_squared =c(Estimate_byc2_squared,ccc)
    Estimate_byc2_4th_power = c(Estimate_byc2_4th_power,ddd)
  }
  
  ################################################
  ############## preferential models ########
  ########################################
  
  ## prepare data to mdoel sampling intensity
  df_random_pp = left_join(cov_effect,random_sample)%>%group_by(covariate_space) %>% 
    mutate(n_obs =length(which(!is.na(preferential))))
  df_proportional_pp = left_join(cov_effect,proportional_sample)%>%group_by(covariate_space) %>% 
    mutate(n_obs =length(which(!is.na(preferential))))
  df_squared_pp = left_join(cov_effect,squared_sample)%>%group_by(covariate_space) %>% 
    mutate(n_obs =length(which(!is.na(preferential))))
  df_4th_power_pp = left_join(cov_effect,fourth_power_sample)%>%group_by(covariate_space) %>% 
    mutate(n_obs =length(which(!is.na(preferential))))
  
  
  #### set a prior that pulls towards zero 
  prior_to_zero <- list(prior = 'gaussian', param = c(0, .051))
  
  
  ######################################################################
  ### XXX        --> Preferntial model with fishers error included #####
  ### XXX_simple --> Traditional preferntial models                #####
  ######################################################################
  
  ### model components
  cmp_target =  ~ covariate_pp(covariate_space,model=matern) +
    covariate_pp_copy(covariate_space, copy = "covariate_pp", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    covariate_biomass_target(covariate_space,model=matern_zero) +
    Intercept_biomass_target(1) + Intercept_pp(1) # +
  
  cmp_target_simple =  ~ covariate_pp(covariate_space,model=matern) +
    covariate_pp_copy(covariate_space, copy = "covariate_pp", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    Intercept_biomass_target(1) + Intercept_pp(1) 
  
  cmp_bycatch =  ~ covariate_pp(covariate_space,model=matern) +
    covariate_pp_copy(covariate_space, copy = "covariate_pp", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    Intercept_pp(1) +
    covariate_biomass_bycatch(covariate_space,model=matern_zero) +
    Intercept_biomass_bycatch(1)
  
  cmp_bycatch_simple =  ~ covariate_pp(covariate_space,model=matern) +
    covariate_pp_copy(covariate_space, copy = "covariate_pp", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    Intercept_pp(1) +
    Intercept_biomass_bycatch(1)
  
  cmp_bycatch2 =  ~     covariate_biomass_bycatch(covariate_space,model=matern_zero) +
    covariate_pp(covariate_space,model=matern) +
    covariate_pp_copy(covariate_space, copy = "covariate_pp", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    Intercept_pp(1) +
    Intercept_biomass_bycatch2(1)
  
  cmp_bycatch2_simple =  ~ covariate_pp(covariate_space,model=matern) +
    covariate_pp_copy(covariate_space, copy = "covariate_pp", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    Intercept_pp(1) +
    Intercept_biomass_bycatch2(1)
  
  
  
  #### formulas
  
  # sampling intensity
  form_pp = n_obs  ~   covariate_pp +   Intercept_pp
  
  # preferntial model with fishers error
  form_biomass_target = biomass_target_FD ~   covariate_pp_copy +  covariate_biomass_target +  Intercept_biomass_target
  form_biomass_bycatch = biomass_bycatch_FD ~   covariate_pp_copy +  covariate_biomass_bycatch +  Intercept_biomass_bycatch
  form_biomass_bycatch2 = biomass_bycatch2_FD ~   covariate_pp_copy +  covariate_biomass_bycatch +  Intercept_biomass_bycatch2
  
  # traditional preferntial model 
  form_biomass_target_simple = biomass_target_FD ~   covariate_pp_copy  +  Intercept_biomass_target
  form_biomass_bycatch_simple = biomass_bycatch_FD ~   covariate_pp_copy  +  Intercept_biomass_bycatch
  form_biomass_bycatch2_simple = biomass_bycatch2_FD ~   covariate_pp_copy  +  Intercept_biomass_bycatch2
  
  ####################
  ##### likeluhoods 
  
  ### sampling intensity
  lik_pp_random <- like("poisson",
                        formula = form_pp,
                        data = df_random_pp
  )
  
  
  lik_pp_proportional <- like("poisson",
                              formula = form_pp,
                              data = df_proportional_pp
  )
  
  lik_pp_squared <- like("poisson",
                         formula = form_pp,
                         data = df_squared_pp
  )
  
  lik_pp_4th_power <- like("poisson",
                           formula = form_pp,
                           data = df_4th_power_pp
  )
  
  
  #### preferential model with fishers error
  lik_biomass_target_random <- like("gaussian",
                                    formula = form_biomass_target,
                                    data = random_sample
  )
  
  lik_biomass_target_proportional <- like("gaussian",
                                          formula = form_biomass_target,
                                          data = proportional_sample
  )
  lik_biomass_target_squared <- like("gaussian",
                                     formula = form_biomass_target,
                                     data = squared_sample
  )
  lik_biomass_target_4th_power <- like("gaussian",
                                       formula = form_biomass_target,
                                       data = fourth_power_sample
  )
  
  lik_biomass_bycatch_random <- like("gaussian",
                                     formula = form_biomass_bycatch,
                                     data = random_sample
  )
  lik_biomass_bycatch_proportional <- like("gaussian",
                                           formula = form_biomass_bycatch,
                                           data = proportional_sample
  )
  lik_biomass_bycatch_squared <- like("gaussian",
                                      formula = form_biomass_bycatch,
                                      data = squared_sample
  )
  lik_biomass_bycatch_4th_power <- like("gaussian",
                                        formula = form_biomass_bycatch,
                                        data = fourth_power_sample
  )
  
  lik_biomass_bycatch2_random <- like("gaussian",
                                      formula = form_biomass_bycatch2,
                                      data = random_sample
  )
  lik_biomass_bycatch2_proportional <- like("gaussian",
                                            formula = form_biomass_bycatch2,
                                            data = proportional_sample
  )
  lik_biomass_bycatch2_squared <- like("gaussian",
                                       formula = form_biomass_bycatch2,
                                       data = squared_sample
  )
  lik_biomass_bycatch2_4th_power <- like("gaussian",
                                         formula = form_biomass_bycatch2,
                                         data = fourth_power_sample
  )
  
  ##### traditional preferential model
  lik_biomass_target_random_simple <- like("gaussian",
                                           formula = form_biomass_target_simple,
                                           data = random_sample
  )
  
  lik_biomass_target_proportional_simple <- like("gaussian",
                                                 formula = form_biomass_target_simple,
                                                 data = proportional_sample
  )
  lik_biomass_target_squared_simple <- like("gaussian",
                                            formula = form_biomass_target_simple,
                                            data = squared_sample
  )
  lik_biomass_target_4th_power_simple <- like("gaussian",
                                              formula = form_biomass_target_simple,
                                              data = fourth_power_sample
  )
  
  lik_biomass_bycatch_random_simple <- like("gaussian",
                                            formula = form_biomass_bycatch_simple,
                                            data = random_sample
  )
  lik_biomass_bycatch_proportional_simple <- like("gaussian",
                                                  formula = form_biomass_bycatch_simple,
                                                  data = proportional_sample
  )
  lik_biomass_bycatch_squared_simple <- like("gaussian",
                                             formula = form_biomass_bycatch_simple,
                                             data = squared_sample
  )
  lik_biomass_bycatch_4th_power_simple <- like("gaussian",
                                               formula = form_biomass_bycatch_simple,
                                               data = fourth_power_sample
  )
  
  lik_biomass_bycatch2_random_simple <- like("gaussian",
                                             formula = form_biomass_bycatch2_simple,
                                             data = random_sample
  )
  lik_biomass_bycatch2_proportional_simple <- like("gaussian",
                                                   formula = form_biomass_bycatch2_simple,
                                                   data = proportional_sample
  )
  lik_biomass_bycatch2_squared_simple <- like("gaussian",
                                              formula = form_biomass_bycatch2_simple,
                                              data = squared_sample
  )
  lik_biomass_bycatch2_4th_power_simple <- like("gaussian",
                                                formula = form_biomass_bycatch2_simple,
                                                data = fourth_power_sample
  )
  
  
  
  ###################
  ### Fit models ####
  ###################
  
  #### preferential models with fishers error
  pref_random_fit <- bru(cmp_target, lik_pp_random, lik_biomass_target_random)
  pref_proportional_fit <- bru(cmp_target, lik_pp_proportional, lik_biomass_target_proportional)
  pref_squared_fit <- bru(cmp_target, lik_pp_squared, lik_biomass_target_squared)
  pref_4th_power_fit <- bru(cmp_target, lik_pp_4th_power, lik_biomass_target_4th_power)
  
  pref_byc_random_fit <- bru(cmp_bycatch, lik_pp_random, lik_biomass_bycatch_random)
  pref_byc_proportional_fit <- bru(cmp_bycatch, lik_pp_proportional, lik_biomass_bycatch_proportional)
  pref_byc_squared_fit <- bru(cmp_bycatch, lik_pp_squared, lik_biomass_bycatch_squared)
  pref_byc_4th_power_fit <- bru(cmp_bycatch, lik_pp_4th_power, lik_biomass_bycatch_4th_power)
  
  pref_byc2_random_fit <- bru(cmp_bycatch2, lik_pp_random, lik_biomass_bycatch2_random)
  pref_byc2_proportional_fit <- bru(cmp_bycatch2, lik_pp_proportional, lik_biomass_bycatch2_proportional)
  pref_byc2_squared_fit <- bru(cmp_bycatch2, lik_pp_squared, lik_biomass_bycatch2_squared)
  pref_byc2_4th_power_fit <- bru(cmp_bycatch2, lik_pp_4th_power, lik_biomass_bycatch2_4th_power)
  

  
  #### traditional preferential models 
  simple_pref_random_fit <- bru(cmp_target_simple, lik_pp_random, lik_biomass_target_random_simple)
  simple_pref_proportional_fit <- bru(cmp_target_simple, lik_pp_proportional, lik_biomass_target_proportional_simple)
  simple_pref_squared_fit <- bru(cmp_target_simple, lik_pp_squared, lik_biomass_target_squared_simple)
  simple_pref_4th_power_fit <- bru(cmp_target_simple, lik_pp_4th_power, lik_biomass_target_4th_power_simple)
  
  simple_pref_byc_random_fit <- bru(cmp_bycatch_simple, lik_pp_random, lik_biomass_bycatch_random_simple)
  simple_pref_byc_proportional_fit <- bru(cmp_bycatch_simple, lik_pp_proportional, lik_biomass_bycatch_proportional_simple)
  simple_pref_byc_squared_fit <- bru(cmp_bycatch_simple, lik_pp_squared, lik_biomass_bycatch_squared_simple)
  simple_pref_byc_4th_power_fit <- bru(cmp_bycatch_simple, lik_pp_4th_power, lik_biomass_bycatch_4th_power_simple)
  
  simple_pref_byc2_random_fit <- bru(cmp_bycatch2_simple, lik_pp_random, lik_biomass_bycatch2_random_simple)
  simple_pref_byc2_proportional_fit <- bru(cmp_bycatch2_simple, lik_pp_proportional, lik_biomass_bycatch2_proportional_simple)
  simple_pref_byc2_squared_fit <- bru(cmp_bycatch2_simple, lik_pp_squared, lik_biomass_bycatch2_squared_simple)
  simple_pref_byc2_4th_power_fit <- bru(cmp_bycatch2_simple, lik_pp_4th_power, lik_biomass_bycatch2_4th_power_simple)
  
  
  
#### calculate quantiles@zero
  q0_beta_a =inla.pmarginal(c(0),pref_random_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_b =inla.pmarginal(c(0),pref_proportional_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_c =inla.pmarginal(c(0),pref_squared_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_d =inla.pmarginal(c(0),pref_4th_power_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  
  q0_beta_aa =inla.pmarginal(c(0),pref_byc_random_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_bb =inla.pmarginal(c(0),pref_byc_proportional_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_cc =inla.pmarginal(c(0),pref_byc_squared_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_dd =inla.pmarginal(c(0),pref_byc_4th_power_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  
  q0_beta_aaa =inla.pmarginal(c(0),pref_byc2_random_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_bbb =inla.pmarginal(c(0),pref_byc2_proportional_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_ccc =inla.pmarginal(c(0),pref_byc2_squared_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  q0_beta_ddd =inla.pmarginal(c(0),pref_byc2_4th_power_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
  
  
  ######################################
  ###### Estimated mean abundance #####
  ####################################
  
  ######## preferential model wirh fishers error
  p_a = mean(predict(pref_random_fit,cov_effect,
                     ~(Intercept_biomass_target+covariate_pp_copy + covariate_biomass_target),
                     include=c("Intercept_biomass_target","covariate_pp_copy","covariate_biomass_target"))$mean)
  
  p_b = mean(predict(pref_proportional_fit,cov_effect,
                     ~(Intercept_biomass_target+covariate_pp_copy + covariate_biomass_target),
                     include=c("Intercept_biomass_target","covariate_pp_copy","covariate_biomass_target"))$mean)
  
  p_c = mean(predict(pref_squared_fit,cov_effect,
                     ~(Intercept_biomass_target+covariate_pp_copy + covariate_biomass_target),
                     include=c("Intercept_biomass_target","covariate_pp_copy","covariate_biomass_target"))$mean)
  
  p_d = mean(predict(pref_4th_power_fit,cov_effect,
                     ~(Intercept_biomass_target+covariate_pp_copy + covariate_biomass_target),
                     include=c("Intercept_biomass_target","covariate_pp_copy","covariate_biomass_target"))$mean)
  
  p_aa = mean(predict(pref_byc_random_fit,cov_effect,
                      ~(Intercept_biomass_bycatch+covariate_pp_copy + covariate_biomass_bycatch),
                      include=c("Intercept_biomass_bycatch","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  p_bb = mean(predict(pref_byc_proportional_fit,cov_effect,
                      ~(Intercept_biomass_bycatch+covariate_pp_copy + covariate_biomass_bycatch),
                      include=c("Intercept_biomass_bycatch","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  p_cc = mean(predict(pref_byc_squared_fit,cov_effect,
                      ~(Intercept_biomass_bycatch+covariate_pp_copy + covariate_biomass_bycatch),
                      include=c("Intercept_biomass_bycatch","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  p_dd = mean(predict(pref_byc_4th_power_fit,cov_effect,
                      ~(Intercept_biomass_bycatch+covariate_pp_copy + covariate_biomass_bycatch),
                      include=c("Intercept_biomass_bycatch","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  
  p_aaa = mean(predict(pref_byc2_random_fit,cov_effect,
                       ~(Intercept_biomass_bycatch2+covariate_pp_copy + covariate_biomass_bycatch),
                       include=c("Intercept_biomass_bycatch2","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  p_bbb = mean(predict(pref_byc2_proportional_fit,cov_effect,
                       ~(Intercept_biomass_bycatch2+covariate_pp_copy + covariate_biomass_bycatch),
                       include=c("Intercept_biomass_bycatch2","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  p_ccc = mean(predict(pref_byc2_squared_fit,cov_effect,
                       ~(Intercept_biomass_bycatch2+covariate_pp_copy + covariate_biomass_bycatch),
                       include=c("Intercept_biomass_bycatch2","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  p_ddd = mean(predict(pref_byc2_4th_power_fit,cov_effect,
                       ~(Intercept_biomass_bycatch2+covariate_pp_copy + covariate_biomass_bycatch),
                       include=c("Intercept_biomass_bycatch2","covariate_pp_copy","covariate_biomass_bycatch"))$mean)
  
  
  
  
  ######### Tradional preferential models
  a_simple = mean(predict(simple_pref_random_fit,cov_effect,
                          ~(Intercept_biomass_target+covariate_pp_copy ), 
                          include=c("Intercept_biomass_target","covariate_pp_copy"))$mean)
  
  b_simple = mean(predict(simple_pref_proportional_fit,cov_effect,
                          ~(Intercept_biomass_target+covariate_pp_copy ), 
                          include=c("Intercept_biomass_target","covariate_pp_copy"))$mean)
  
  c_simple = mean(predict(simple_pref_squared_fit,cov_effect,
                          ~(Intercept_biomass_target+covariate_pp_copy ), 
                          include=c("Intercept_biomass_target","covariate_pp_copy"))$mean)
  
  d_simple = mean(predict(simple_pref_4th_power_fit,cov_effect,
                          ~(Intercept_biomass_target+covariate_pp_copy ), 
                          include=c("Intercept_biomass_target","covariate_pp_copy"))$mean)
  
  aa_simple = mean(predict(simple_pref_byc_random_fit,cov_effect,
                           ~(Intercept_biomass_bycatch+covariate_pp_copy ), 
                           include=c("Intercept_biomass_bycatch","covariate_pp_copy"))$mean)
  
  bb_simple = mean(predict(simple_pref_byc_proportional_fit,cov_effect,
                           ~(Intercept_biomass_bycatch+covariate_pp_copy ), 
                           include=c("Intercept_biomass_bycatch","covariate_pp_copy"))$mean)
  
  cc_simple = mean(predict(simple_pref_byc_squared_fit,cov_effect,
                           ~(Intercept_biomass_bycatch+covariate_pp_copy ), 
                           include=c("Intercept_biomass_bycatch","covariate_pp_copy"))$mean)
  
  dd_simple = mean(predict(simple_pref_byc_4th_power_fit,cov_effect,
                           ~(Intercept_biomass_bycatch+covariate_pp_copy ), 
                           include=c("Intercept_biomass_bycatch","covariate_pp_copy"))$mean)
  
  
  aaa_simple = mean(predict(simple_pref_byc2_random_fit,cov_effect,
                            ~(Intercept_biomass_bycatch2+covariate_pp_copy ), 
                            include=c("Intercept_biomass_bycatch2","covariate_pp_copy"))$mean)
  
  bbb_simple = mean(predict(simple_pref_byc2_proportional_fit,cov_effect,
                            ~(Intercept_biomass_bycatch2+covariate_pp_copy ), 
                            include=c("Intercept_biomass_bycatch2","covariate_pp_copy"))$mean)
  
  ccc_simple = mean(predict(simple_pref_byc2_squared_fit,cov_effect,
                            ~(Intercept_biomass_bycatch2+covariate_pp_copy ), 
                            include=c("Intercept_biomass_bycatch2","covariate_pp_copy"))$mean)
  
  ddd_simple = mean(predict(simple_pref_byc2_4th_power_fit,cov_effect,
                            ~(Intercept_biomass_bycatch2+covariate_pp_copy ), 
                            include=c("Intercept_biomass_bycatch2","covariate_pp_copy"))$mean)
  
  
  
  if(i==1){
    Pref_error_random = mean_biomass_target_FD - p_a
    Pref_error_proportional= mean_biomass_target_FD-p_b
    Pref_error_squared = mean_biomass_target_FD - p_c
    Pref_error_4th_power = mean_biomass_target_FD - p_d
    
    Pref_byc_error_random = mean_biomass_byc_FD - p_aa
    Pref_byc_error_proportional= mean_biomass_byc_FD - p_bb
    Pref_byc_error_squared = mean_biomass_byc_FD - p_cc
    Pref_byc_error_4th_power = mean_biomass_byc_FD - p_dd
    
    Pref_byc2_error_random = mean_biomass_byc2_FD - p_aaa
    Pref_byc2_error_proportional= mean_biomass_byc2_FD - p_bbb
    Pref_byc2_error_squared = mean_biomass_byc2_FD - p_ccc
    Pref_byc2_error_4th_power = mean_biomass_byc2_FD - p_ddd
    
    ##### pref simple
    simple_Pref_error_random = mean_biomass_target_FD - a_simple
    simple_Pref_error_proportional= mean_biomass_target_FD-b_simple
    simple_Pref_error_squared = mean_biomass_target_FD - c_simple
    simple_Pref_error_4th_power = mean_biomass_target_FD - d_simple
    
    simple_Pref_byc_error_random = mean_biomass_byc_FD - aa_simple
    simple_Pref_byc_error_proportional= mean_biomass_byc_FD - bb_simple
    simple_Pref_byc_error_squared = mean_biomass_byc_FD - cc_simple
    simple_Pref_byc_error_4th_power = mean_biomass_byc_FD - dd_simple
    
    simple_Pref_byc2_error_random = mean_biomass_byc2_FD - aaa_simple
    simple_Pref_byc2_error_proportional= mean_biomass_byc2_FD - bbb_simple
    simple_Pref_byc2_error_squared = mean_biomass_byc2_FD - ccc_simple
    simple_Pref_byc2_error_4th_power = mean_biomass_byc2_FD - ddd_simple
    
    # beta_random = beta_a
    # beta_proportional = beta_b
    # beta_squared = beta_c
    # beta_4th_power = beta_d
    # 
    # beta_byc_random = beta_aa
    # beta_byc_proportional = beta_bb
    # beta_byc_squared = beta_cc
    # beta_byc_4th_power = beta_dd
    # 
    # beta_byc2_random = beta_aaa
    # beta_byc2_proportional = beta_bbb
    # beta_byc2_squared = beta_ccc
    # beta_byc2_4th_power = beta_ddd
    
    
    
    ##### quantile at 0
    q0_beta_random = q0_beta_a
    q0_beta_proportional = q0_beta_b
    q0_beta_squared = q0_beta_c
    q0_beta_4th_power = q0_beta_d
    
    q0_beta_byc_random = q0_beta_aa
    q0_beta_byc_proportional = q0_beta_bb
    q0_beta_byc_squared = q0_beta_cc
    q0_beta_byc_4th_power = q0_beta_dd
    
    q0_beta_byc2_random = q0_beta_aaa
    q0_beta_byc2_proportional = q0_beta_bbb
    q0_beta_byc2_squared = q0_beta_ccc
    q0_beta_byc2_4th_power = q0_beta_ddd
    
  }else{
    Pref_error_random=c(Pref_error_random,mean_biomass_target_FD - p_a)
    Pref_error_proportional=c(Pref_error_proportional,mean_biomass_target_FD - p_b)
    Pref_error_squared =c(Pref_error_squared,mean_biomass_target_FD - p_c)
    Pref_error_4th_power = c(Pref_error_4th_power,mean_biomass_target_FD - p_d)
    
    Pref_byc_error_random=c(Pref_byc_error_random,mean_biomass_byc_FD -p_aa)
    Pref_byc_error_proportional=c(Pref_byc_error_proportional,mean_biomass_byc_FD -p_bb)
    Pref_byc_error_squared =c(Pref_byc_error_squared,mean_biomass_byc_FD -p_cc)
    Pref_byc_error_4th_power = c(Pref_byc_error_4th_power,mean_biomass_byc_FD -p_dd)
    
    Pref_byc2_error_random=c(Pref_byc2_error_random,mean_biomass_byc2_FD -p_aaa)
    Pref_byc2_error_proportional=c(Pref_byc2_error_proportional,mean_biomass_byc2_FD -p_bbb)
    Pref_byc2_error_squared =c(Pref_byc2_error_squared,mean_biomass_byc2_FD -p_ccc)
    Pref_byc2_error_4th_power = c(Pref_byc2_error_4th_power,mean_biomass_byc2_FD -p_ddd)
    
    #### simple preferential
    simple_Pref_error_random=c(simple_Pref_error_random,mean_biomass_target_FD - a_simple)
    simple_Pref_error_proportional=c(simple_Pref_error_proportional,mean_biomass_target_FD - b_simple)
    simple_Pref_error_squared =c(simple_Pref_error_squared,mean_biomass_target_FD - c_simple)
    simple_Pref_error_4th_power = c(simple_Pref_error_4th_power,mean_biomass_target_FD - d_simple)
    
    simple_Pref_byc_error_random=c(simple_Pref_byc_error_random,mean_biomass_byc_FD -aa_simple)
    simple_Pref_byc_error_proportional=c(simple_Pref_byc_error_proportional,mean_biomass_byc_FD -bb_simple)
    simple_Pref_byc_error_squared =c(simple_Pref_byc_error_squared,mean_biomass_byc_FD -cc_simple)
    simple_Pref_byc_error_4th_power = c(simple_Pref_byc_error_4th_power,mean_biomass_byc_FD -dd_simple)
    
    simple_Pref_byc2_error_random=c(simple_Pref_byc2_error_random,mean_biomass_byc2_FD -aaa_simple)
    simple_Pref_byc2_error_proportional=c(simple_Pref_byc2_error_proportional,mean_biomass_byc2_FD -bbb_simple)
    simple_Pref_byc2_error_squared =c(simple_Pref_byc2_error_squared,mean_biomass_byc2_FD -ccc_simple)
    simple_Pref_byc2_error_4th_power = c(simple_Pref_byc2_error_4th_power,mean_biomass_byc2_FD -ddd_simple)
    
    ##### quantile at 0
    q0_beta_random = c(q0_beta_random,q0_beta_a)
    q0_beta_proportional = c(q0_beta_proportional,q0_beta_b)
    q0_beta_squared = c(q0_beta_squared,q0_beta_c)
    q0_beta_4th_power = c(q0_beta_4th_power,q0_beta_d)
    
    q0_beta_byc_random = c(q0_beta_byc_random,q0_beta_aa)
    q0_beta_byc_proportional = c(q0_beta_byc_proportional,q0_beta_bb)
    q0_beta_byc_squared = c(q0_beta_byc_squared,q0_beta_cc)
    q0_beta_byc_4th_power = c(q0_beta_byc_4th_power,q0_beta_dd)
    
    q0_beta_byc2_random = c(q0_beta_byc2_random,q0_beta_aaa)
    q0_beta_byc2_proportional = c(q0_beta_byc2_proportional,q0_beta_bbb)
    q0_beta_byc2_squared = c(q0_beta_byc2_squared,q0_beta_ccc)
    q0_beta_byc2_4th_power = c(q0_beta_byc2_4th_power,q0_beta_ddd)
    
  }
  
  
  ######################################
  ######## Integrated models ###########
  ######################################
  
  ### get 1/3 of randomly collected data to mimmic a FI survey
  FI_data = cov_effect[c(1,sample(size = round(n/3)-2,cov_effect$idx,replace=T),nrow(cov_effect)),] %>%  #### allow for first and last observations for corners of distribution. Otherwise very influential
    select(1,2,3,5,6,12:17) %>% mutate(preferential="Random",Data_type = "FI")
  
  ### get remaining 2/3 of data using the different sampling schems 
  FD_data_random = random_sample[sample(1:n,round(n/3)*2), ] %>% mutate(Data_type = "FD")
  FD_data_proportional = proportional_sample[sample(1:n,round(n/3)*2), ] %>% mutate(Data_type = "FD")
  FD_data_squared = squared_sample[sample(1:n,round(n/3)*2), ] %>% mutate(Data_type = "FD")
  FD_data_fourth = fourth_power_sample[sample(1:n,round(n/3)*2), ] %>% mutate(Data_type = "FD")
  
  
  ###### cmps
  cmp_target =  ~ covariate_biomass_target(covariate_space,model=matern) +
    covariate_biomass_copy(covariate_space, copy = "covariate_biomass_target", fixed = T, hyper = list(beta = prior_to_zero)) +
    Intercept_FI_target(1) + Intercept_FD_target(1) # +
  
  cmp_bycatch =  ~ covariate_biomass_bycatch(covariate_space,model=matern) +
    covariate_biomass_copy(covariate_space, copy = "covariate_biomass_bycatch", fixed = T, hyper = list(beta = prior_to_zero)) +
    Intercept_FI_bycatch(1) + Intercept_FD_bycatch(1)
  
  cmp_bycatch2 =  ~ covariate_biomass_bycatch(covariate_space,model=matern) +
    covariate_biomass_copy(covariate_space, copy = "covariate_biomass_bycatch", fixed = T, hyper = list(beta = prior_to_zero)) +
    Intercept_FI_bycatch2(1) + Intercept_FD_bycatch2(1)
  
  
  ###### formulas
  form_FI_target = biomass_target_FI ~   covariate_biomass_target +  Intercept_FI_target
  form_FI_bycatch = biomass_bycatch_FI ~   covariate_biomass_bycatch +  Intercept_FI_bycatch
  form_FI_bycatch2 = biomass_bycatch2_FI ~   covariate_biomass_bycatch +  Intercept_FI_bycatch2
  
  form_FD_target = biomass_target_FD ~   covariate_biomass_copy  +  Intercept_FD_target
  form_FD_bycatch = biomass_bycatch_FD ~   covariate_biomass_copy  +  Intercept_FD_bycatch
  form_FD_bycatch2 = biomass_bycatch2_FD ~   covariate_biomass_copy  +  Intercept_FD_bycatch2
  
  
  ###### likelihoods
  
  ## FI
  lik_biomass_FI <- like("gaussian",
                         formula = form_FI_target,
                         data = FI_data
  )
  
  lik_biomass_FI_bycatch <- like("gaussian",
                                 formula = form_FI_bycatch,
                                 data = FI_data
  )
  
  lik_biomass_FI_bycatch2 <- like("gaussian",
                                  formula = form_FI_bycatch2,
                                  data = FI_data
  )
  
  ## FD
  lik_biomass_FD_random <- like("gaussian",
                                formula = form_FD_target,
                                data = FD_data_random
  )
  
  lik_biomass_FD_proportional <- like("gaussian",
                                      formula = form_FD_target,
                                      data = FD_data_proportional
  )
  
  lik_biomass_FD_squared <- like("gaussian",
                                 formula = form_FD_target,
                                 data = FD_data_squared
  )
  
  lik_biomass_FD_fourth <- like("gaussian",
                                formula = form_FD_target,
                                data = FD_data_fourth
  )
  
  lik_biomass_FD_random_bycatch <- like("gaussian",
                                        formula = form_FD_bycatch,
                                        data = FD_data_random
  )
  
  lik_biomass_FD_proportional_bycatch <- like("gaussian",
                                              formula = form_FD_bycatch,
                                              data = FD_data_proportional
  )
  
  lik_biomass_FD_squared_bycatch <- like("gaussian",
                                         formula = form_FD_bycatch,
                                         data = FD_data_squared
  )
  
  lik_biomass_FD_fourth_bycatch <- like("gaussian",
                                        formula = form_FD_bycatch,
                                        data = FD_data_fourth
  )
  
  lik_biomass_FD_random_bycatch2 <- like("gaussian",
                                         formula = form_FD_bycatch2,
                                         data = FD_data_random
  )
  
  lik_biomass_FD_proportional_bycatch2 <- like("gaussian",
                                               formula = form_FD_bycatch2,
                                               data = FD_data_proportional
  )
  
  lik_biomass_FD_squared_bycatch2 <- like("gaussian",
                                          formula = form_FD_bycatch2,
                                          data = FD_data_squared
  )
  
  lik_biomass_FD_fourth_bycatch2 <- like("gaussian",
                                         formula = form_FD_bycatch2,
                                         data = FD_data_fourth
  )
  
  #### target species
  comb_target_random_fit <- bru(cmp_target, lik_biomass_FI, lik_biomass_FD_random)
  comb_target_proportional_fit <- bru(cmp_target, lik_biomass_FI, lik_biomass_FD_proportional)
  comb_target_squared_fit <- bru(cmp_target, lik_biomass_FI, lik_biomass_FD_squared)
  comb_target_4th_power_fit <- bru(cmp_target, lik_biomass_FI, lik_biomass_FD_fourth)
  
  #### bycatch species
  comb_bycatch_random_fit <- bru(cmp_bycatch, lik_biomass_FI_bycatch, lik_biomass_FD_random_bycatch)
  comb_bycatch_proportional_fit <- bru(cmp_bycatch, lik_biomass_FI_bycatch, lik_biomass_FD_proportional_bycatch)
  comb_bycatch_squared_fit <- bru(cmp_bycatch, lik_biomass_FI_bycatch, lik_biomass_FD_squared_bycatch)
  comb_bycatch_4th_power_fit <- bru(cmp_bycatch, lik_biomass_FI_bycatch, lik_biomass_FD_fourth_bycatch)
  
  #### bycatch2 species
  comb_bycatch2_random_fit <- bru(cmp_bycatch2, lik_biomass_FI_bycatch2, lik_biomass_FD_random_bycatch2)
  comb_bycatch2_proportional_fit <- bru(cmp_bycatch2, lik_biomass_FI_bycatch2, lik_biomass_FD_proportional_bycatch2)
  comb_bycatch2_squared_fit <- bru(cmp_bycatch2, lik_biomass_FI_bycatch2, lik_biomass_FD_squared_bycatch2)
  comb_bycatch2_4th_power_fit <- bru(cmp_bycatch2, lik_biomass_FI_bycatch2, lik_biomass_FD_fourth_bycatch2)
  
  
  a_comb = mean(predict(comb_target_random_fit,cov_effect,
                        ~(Intercept_FD_target+ covariate_biomass_copy), 
                        include=c("Intercept_FD_target","covariate_biomass_copy"))$mean)
  
  b_comb = mean(predict(comb_target_proportional_fit,cov_effect,
                        ~(Intercept_FD_target+ covariate_biomass_copy), 
                        include=c("Intercept_FD_target","covariate_biomass_copy"))$mean)
  
  c_comb = mean(predict(comb_target_squared_fit,cov_effect,
                        ~(Intercept_FD_target+ covariate_biomass_copy), 
                        include=c("Intercept_FD_target","covariate_biomass_copy"))$mean)
  
  d_comb = mean(predict(comb_target_4th_power_fit,cov_effect,
                        ~(Intercept_FD_target+ covariate_biomass_copy), 
                        include=c("Intercept_FD_target","covariate_biomass_copy"))$mean)
  
  
  aa_comb = mean(predict(comb_bycatch_random_fit,cov_effect,
                         ~(Intercept_FD_bycatch+ covariate_biomass_copy), 
                         include=c("Intercept_FD_bycatch","covariate_biomass_copy"))$mean)
  
  bb_comb = mean(predict(comb_bycatch_proportional_fit,cov_effect,
                         ~(Intercept_FD_bycatch+ covariate_biomass_copy), 
                         include=c("Intercept_FD_bycatch","covariate_biomass_copy"))$mean)
  
  cc_comb = mean(predict(comb_bycatch_squared_fit,cov_effect,
                         ~(Intercept_FD_bycatch+ covariate_biomass_copy), 
                         include=c("Intercept_FD_bycatch","covariate_biomass_copy"))$mean)
  
  dd_comb = mean(predict(comb_bycatch_4th_power_fit,cov_effect,
                         ~(Intercept_FD_bycatch+ covariate_biomass_copy), 
                         include=c("Intercept_FD_bycatch","covariate_biomass_copy"))$mean)
  
  
  aaa_comb = mean(predict(comb_bycatch2_random_fit,cov_effect,
                          ~(Intercept_FD_bycatch2+ covariate_biomass_copy), 
                          include=c("Intercept_FD_bycatch2","covariate_biomass_copy"))$mean)
  
  bbb_comb = mean(predict(comb_bycatch2_proportional_fit,cov_effect,
                          ~(Intercept_FD_bycatch2+ covariate_biomass_copy), 
                          include=c("Intercept_FD_bycatch2","covariate_biomass_copy"))$mean)
  
  ccc_comb = mean(predict(comb_bycatch2_squared_fit,cov_effect,
                          ~(Intercept_FD_bycatch2+ covariate_biomass_copy), 
                          include=c("Intercept_FD_bycatch2","covariate_biomass_copy"))$mean)
  
  ddd_comb = mean(predict(comb_bycatch2_4th_power_fit,cov_effect,
                          ~(Intercept_FD_bycatch2+ covariate_biomass_copy), 
                          include=c("Intercept_FD_bycatch2","covariate_biomass_copy"))$mean)
  
  
  
  if(i==1){
    comb_error_random = mean_biomass_target_FD - a_comb
    comb_error_proportional= mean_biomass_target_FD-b_comb
    comb_error_squared = mean_biomass_target_FD - c_comb
    comb_error_4th_power = mean_biomass_target_FD - d_comb
    
    comb_byc_error_random = mean_biomass_byc_FD - aa_comb
    comb_byc_error_proportional= mean_biomass_byc_FD - bb_comb
    comb_byc_error_squared = mean_biomass_byc_FD - cc_comb
    comb_byc_error_4th_power = mean_biomass_byc_FD - dd_comb
    
    comb_byc2_error_random = mean_biomass_byc2_FD - aaa_comb
    comb_byc2_error_proportional= mean_biomass_byc2_FD - bbb_comb
    comb_byc2_error_squared = mean_biomass_byc2_FD - ccc_comb
    comb_byc2_error_4th_power = mean_biomass_byc2_FD - ddd_comb
    
    
  }else{
    comb_error_random=c(comb_error_random,mean_biomass_target_FD - a_comb)
    comb_error_proportional=c(comb_error_proportional,mean_biomass_target_FD - b_comb)
    comb_error_squared =c(comb_error_squared,mean_biomass_target_FD - c_comb)
    comb_error_4th_power = c(comb_error_4th_power,mean_biomass_target_FD - d_comb)
    
    comb_byc_error_random=c(comb_byc_error_random,mean_biomass_byc_FD -aa_comb)
    comb_byc_error_proportional=c(comb_byc_error_proportional,mean_biomass_byc_FD -bb_comb)
    comb_byc_error_squared =c(comb_byc_error_squared,mean_biomass_byc_FD -cc_comb)
    comb_byc_error_4th_power = c(comb_byc_error_4th_power,mean_biomass_byc_FD -dd_comb)
    
    comb_byc2_error_random=c(comb_byc2_error_random,mean_biomass_byc2_FD -aaa_comb)
    comb_byc2_error_proportional=c(comb_byc2_error_proportional,mean_biomass_byc2_FD -bbb_comb)
    comb_byc2_error_squared =c(comb_byc2_error_squared,mean_biomass_byc2_FD -ccc_comb)
    comb_byc2_error_4th_power = c(comb_byc2_error_4th_power,mean_biomass_byc2_FD -ddd_comb)
    
  }
  

  
  print(paste(i,i,i,i,Sys.time()-ptm))
  save.image("Simulation_results_Gaussian.RData")
}





############################################
######### Figure 2 #########################
##############################################


##############
### TS
sim_conventional_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                      Estimates = c(Estimate_random,Estimate_proportional,Estimate_squared,Estimate_4th_power))
sim_conventional_results$Error = sim_conventional_results$Estimates - mean_biomass_target_FD
sim_conventional_results$Sampling_scheme = factor(sim_conventional_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

##############
### CBS species
sim_conventional_byc_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                          Estimates = c(Estimate_byc_random,Estimate_byc_proportional,Estimate_byc_squared,Estimate_byc_4th_power))
sim_conventional_byc_results$Error = sim_conventional_byc_results$Estimates - mean_biomass_byc_FD
sim_conventional_byc_results$Sampling_scheme = factor(sim_conventional_byc_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

##############
### UBS species
sim_conventional_byc2_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                           Estimates = c(Estimate_byc2_random,Estimate_byc2_proportional,Estimate_byc2_squared,Estimate_byc2_4th_power))
sim_conventional_byc2_results$Error = sim_conventional_byc2_results$Estimates - mean_biomass_byc2_FD
sim_conventional_byc2_results$Sampling_scheme = factor(sim_conventional_byc2_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))



z = ggplot(sim_conventional_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank()) 

v = ggplot(sim_conventional_byc_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank()) 

y = ggplot(sim_conventional_byc2_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank())  


######################################
### PREFERENTIAL RESUKLTS

##############
### TS
sim_preferential_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                      Error = c(Pref_error_random,Pref_error_proportional,Pref_error_squared,Pref_error_4th_power))
sim_preferential_results$Sampling_scheme = factor(sim_preferential_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

##############
### CBS species
sim_byc_preferential_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                          Error = c(Pref_byc_error_random,Pref_byc_error_proportional,Pref_byc_error_squared,Pref_byc_error_4th_power))
sim_byc_preferential_results$Sampling_scheme = factor(sim_byc_preferential_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

##############
### UBS #####
sim_byc2_preferential_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                           Error = c(Pref_byc2_error_random,Pref_byc2_error_proportional,Pref_byc2_error_squared,Pref_byc2_error_4th_power))
sim_byc2_preferential_results$Sampling_scheme = factor(sim_byc2_preferential_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))


zz = ggplot(sim_preferential_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) +
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank()) 

vv = ggplot(sim_byc_preferential_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) +
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("")+
  theme(axis.text.x=element_blank()) 

yy = ggplot(sim_byc2_preferential_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank())  


##################################
###### traditional preferential 
sim_simple_preferential_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                             Error = c(simple_Pref_error_random,simple_Pref_error_proportional,simple_Pref_error_squared,simple_Pref_error_4th_power))
sim_simple_preferential_results$Sampling_scheme = factor(sim_simple_preferential_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))


sim_byc_simple_preferential_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                                 Error = c(simple_Pref_byc_error_random,simple_Pref_byc_error_proportional,simple_Pref_byc_error_squared,simple_Pref_byc_error_4th_power))
sim_byc_simple_preferential_results$Sampling_scheme = factor(sim_byc_simple_preferential_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))


sim_byc2_simple_preferential_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                                  Error = c(simple_Pref_byc2_error_random,simple_Pref_byc2_error_proportional,simple_Pref_byc2_error_squared,simple_Pref_byc2_error_4th_power))
sim_byc2_simple_preferential_results$Sampling_scheme = factor(sim_byc2_simple_preferential_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))




simple_zz = ggplot(sim_simple_preferential_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank()) 

simple_vv = ggplot(sim_byc_simple_preferential_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank())  

simple_yy = ggplot(sim_byc2_simple_preferential_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") +
  theme(axis.text.x=element_blank()) 





###############################################
####### combine models ######################
############################################

##############
### target species
combine_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                             Error = c(comb_error_random,comb_error_proportional,comb_error_squared,comb_error_4th_power))
combine_results$Sampling_scheme = factor(combine_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))


##############
### bycatch species
combine_byc_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                 Error = c(comb_byc_error_random,comb_byc_error_proportional,comb_byc_error_squared,comb_byc_error_4th_power))
combine_byc_results$Sampling_scheme = factor(combine_byc_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))


##############
### bycatch2 species
combine_byc2_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                  Error = c(comb_byc2_error_random,comb_byc2_error_proportional,comb_byc2_error_squared,comb_byc2_error_4th_power))
combine_byc2_results$Sampling_scheme = factor(combine_byc2_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))



comb_z = ggplot(combine_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") 

comb_v = ggplot(combine_byc_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") 

comb_y = ggplot(combine_byc2_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error)) + 
  geom_hline(yintercept = 0,color="red",linewidth=1) +
  geom_point(aes(x=Sampling_scheme,y=Error)) +
  ylab("")+ xlab("") 


####### Alltogether

col.titles = c("TS","CBS","UBS")
row.titles = c("Conventional SDM", 
               paste("Conventional", "preferential model", sep="\n"),
               "Extended \n preferential model",
               "Integrated SDM")

pl = list(z,v,y,
          simple_zz,simple_vv,simple_yy,
          zz,vv,yy,
          comb_z,comb_v,comb_y)

combine <- rbind(tableGrob(t(col.titles), theme = ttheme_minimal(), rows = ""), 
                 cbind(tableGrob(row.titles, theme = ttheme_minimal()), 
                       arrangeGrob(grobs = pl),  size = "last"), size = "last")

library(grid)

title <- textGrob(
  label = "Level of preferential sampling",
  x = .45, y = 0, hjust = 0, vjust = -2,
  gp = gpar(fontsize = 14, fontface = "bold", col = "black")
)
png("Result_boxplots.png",height=500,width = 900)
grid.arrange(combine, title, nrow = 2, heights = c(0.9, 0.1))
dev.off()





#########################################
############ Figure 3 #################
########################################

q0_pref_beta = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                          Scaling = c(q0_beta_random,q0_beta_proportional,q0_beta_squared,q0_beta_4th_power))
q0_pref_beta$Sampling_scheme = factor(q0_pref_beta$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

q0_pref_byc_beta = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                              Scaling = c(q0_beta_byc_random,q0_beta_byc_proportional,q0_beta_byc_squared,q0_beta_byc_4th_power))
q0_pref_byc_beta$Sampling_scheme = factor(q0_pref_byc_beta$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

q0_pref_byc2_beta = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                               Scaling = c(q0_beta_byc2_random,q0_beta_byc2_proportional,q0_beta_byc2_squared,q0_beta_byc2_4th_power))
q0_pref_byc2_beta$Sampling_scheme = factor(q0_pref_byc2_beta$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))

q0_beta_z = ggplot(q0_pref_beta) + geom_boxplot(aes(x=Sampling_scheme,y=Scaling)) + 
  geom_point(aes(x=Sampling_scheme,y=Scaling)) +
  ggtitle("TS")+ xlab("") + ylim(0,1)  + ylab("") +
  theme(axis.text=element_text(size=12),
        plot.title = element_text(size=18),
        axis.title=element_text(size=12,face="bold"))

q0_beta_v = ggplot(q0_pref_byc_beta) + geom_boxplot(aes(x=Sampling_scheme,y=Scaling)) + 
  geom_point(aes(x=Sampling_scheme,y=Scaling)) +
  ggtitle("CBS")+ xlab("") + ylab("") +
  ylim(0,1)+theme(axis.text=element_text(size=12),
                  plot.title = element_text(size=18),
                  axis.title=element_text(size=12,face="bold")) 

q0_beta_y = ggplot(q0_pref_byc2_beta) + geom_boxplot(aes(x=Sampling_scheme,y=Scaling)) + 
  geom_point(aes(x=Sampling_scheme,y=Scaling)) +
  ggtitle("UBS")+ xlab("") + ylab("") +
  ylim(0,1)+theme(axis.text=element_text(size=12),
                  plot.title = element_text(size=18),
                  axis.title=element_text(size=12,face="bold"))



png("q0_beta.png",height=350,width = 1000)
grid.arrange(q0_beta_z,q0_beta_v,q0_beta_y,
             bottom=textGrob("Level of preferential sampling", gp=gpar(fontsize=17, fontface = "bold")),
             left=textGrob("Quantile at zero", rot = 90, gp=gpar(fontsize=17, fontface = "bold")),
             #bottom="Level of preferential sampling",
             #left = "Quantile at zero",
             ncol=3)
dev.off()





#########################################
############ Figure 4 #################
########################################
model = pref_byc2_squared_fit
combine_covar = predict(model,cov_effect,
                        ~(Intercept_biomass_bycatch2+covariate_pp_copy + covariate_biomass_bycatch),
                        include=c("Intercept_biomass_bycatch2","covariate_pp_copy","covariate_biomass_bycatch"))

FX_covar = predict(model,cov_effect,
                   ~(Intercept_biomass_bycatch2+covariate_pp_copy ),
                   include=c("Intercept_biomass_bycatch2","covariate_pp_copy"))


LX_covar = predict(model,cov_effect,
                   ~(Intercept_biomass_bycatch2 + covariate_biomass_bycatch),
                   include=c("Intercept_biomass_bycatch2","covariate_biomass_bycatch"))

ea = rbind(combine_covar,FX_covar,LX_covar)
ea$effect = rep(c("alpha * f(x) + l(x)","alpha * f(x)", "l(x)"),each=nrow(cov_effect))


squared_sample$shape = "Modelled Data"
cov_effect$shape = "Simulated abundances"


png("fig_combine.png",height=400,width = 800)
ggplot() + 
  geom_point(data=cov_effect,aes(x=covariate_space,y=biomass_bycatch2_FD,shape=shape,color=shape),alpha=.3) +
  geom_point(data= squared_sample,aes(x=covariate_space,y=biomass_bycatch2_FD,shape=shape,color=shape),size=2) +
  geom_line(data=ea,aes(x=covariate_space,y=mean,linetype=effect),linewidth=1.5) +
  ylab("UBS   abundance") + xlab("Covariate space") +
  theme_bw() +
  scale_linetype_manual(values = c(1,2,3),name = "Fitted covariate fields",
                        labels = c(expression(alpha * f(x)),
                                   expression(alpha * f(x) + l(x)),
                                   expression(l(x)))) +
  scale_color_manual(values = c("red","grey"),name = "FD abundance") +
  guides(shape = "none") +
  theme(#axis.text=element_text(size=12),
    axis.title=element_text(size=15,face="bold"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=14),
    legend.text.align = 0)
dev.off()












