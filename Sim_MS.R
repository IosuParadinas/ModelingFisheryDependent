
library(scales)
library(tidyverse)
library(inlabru)
library(INLA)
library(gridExtra)
library(patchwork)  # For combining plots
library(purrr)


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
                                   biomass_target_FD = (rescale(abundance_target,to=c(-1,1)) + 15)*2.5 + rnorm(length(quantiles),sd=.5),
                                   biomass_bycatch_FI = rescale(abundance_bycatch,to=c(-1,1)) + 15 + rnorm(length(quantiles),sd=.25),
                                   biomass_bycatch_FD = (rescale(abundance_bycatch,to=c(-1,1)) + 15)*2.5 + rnorm(length(quantiles),sd=.5),
                                   biomass_bycatch2_FI = rescale(abundance_bycatch2,to=c(-1,1)) + 15 + rnorm(length(quantiles),sd=.25),
                                   biomass_bycatch2_FD = (rescale(abundance_bycatch2,to=c(-1,1)) + 15)*2.5 + rnorm(length(quantiles),sd=.5))


cor(cov_effect[,c(3,6,7,9,10)])

plot(cov_effect$UBS,cov_effect$Squared)
plot(cov_effect$UBS,cov_effect$"Power 4")

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
  select(1,2,3,6,11:14) %>% mutate(sampling="Random", Preferentiality = 1.4)
proportional_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$Proportional,replace=T),] %>%  
  select(1,2,3,6,11:14) %>% mutate(sampling="Proportional", Preferentiality = 1.3)
squared_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$Squared,replace=T),] %>%  
  select(1,2,3,6,11:14) %>% mutate(sampling="Squared", Preferentiality = 1.2)
fourth_power_sample = cov_effect[sample(size = n,cov_effect$idx,prob = cov_effect$"Power 4",replace=T),] %>%
  select(1,2,3,6,11:14) %>% mutate(sampling="Power 4", Preferentiality = 1.1)


sim_data_all = rbind(random_sample,proportional_sample,squared_sample,fourth_power_sample)
sim_data_all = rbind(sim_data_all[,c("covariate_space","sampling","Preferentiality")],
                     prob_plotD[,c("covariate_space","sampling","Preferentiality")])
sim_data_all$sampling = factor(sim_data_all$sampling,levels = c("Random","Proportional","Squared","Power 4"))

text_df1 = data.frame(idx=4:1,y=c(1.1,1.2,1.3,1.4));text_df1$x=105
text_df2 = data.frame(idx=1:4,y=c(.25,0.16,.1,-.05),x=c(25,28,36,45))

#png("C:/use/Pref_revision_paper/img/samp_proba_both_samp_dist.png",height=400,width = 700)
png("C:/Users/ip30/OneDrive - University of St Andrews/Desktop/samp_proba_both_samp_dist.png",height=400,width = 700)
ggplot()+
  geom_line(data=abundance_plotD,aes(x=covariate_space,y=Species_Abundance,color=Species),linewidth=4) + 
  #scale_discrete_manual(aesthetic = "linewidth", values = c(TS=40,CBS=2,UBS=2)) +
  scale_size_manual( values = c(TS=4,CBS=2,UBS=2)) +
  geom_text(data=text_df1,aes(x=x,y=y,label=idx),size=5)+
  geom_text(data=text_df2,aes(x=x,y=y,label=idx),size=5)+
  #geom_point(data=prob_plotD,aes(x=covariate_space,y=Preferentiality,shape=sampling),linewidth=1) +
  #ggtitle("Sampling probability bycatch species") +
  theme_bw() + ylab("Standardized values") + 
  geom_point(data=sim_data_all,aes(x=covariate_space,y=Preferentiality ,shape=sampling),alpha=.6,size=2) +
  xlab("1D habitat space") +
  scale_fill_manual(values = c("orange", "purple","green","blue")) +
  scale_color_manual(values = c("#F8766D", "#00BFC4","yellow"),name = "Species")+
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
# load("D:/Simulation_results_MS3.RData")
#  
 
n = 300
n_sims = 50
ptm <- Sys.time()
error=list()
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
  mesh1D <- fm_mesh_1d(x, boundary = "free")
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
  ### fit conventional SDM ####
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
                   options = list(control.inla = list(int.strategy = "eb")))
  
  a = mean(predict(fit_random,cov_effect,~(Intercept+covariate))$mean)
  
  fit_random_byc = bru(components = form_bycatch,
                       data = random_sample,
                       family ="gaussian",
                       options = list(control.inla = list(int.strategy = "eb")))
  
  aa = mean(predict(fit_random_byc,cov_effect,~(Intercept+covariate))$mean)
  
  fit_random_byc2 = bru(components = form_bycatch2,
                        data = random_sample,
                        family ="gaussian",
                        options = list(control.inla = list(int.strategy = "eb")))
  
  aaa = mean(predict(fit_random_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  
  #### Preferential sampling, proportional to abundance intensity ####
  fit_proportional = bru(components = form_target,
                         data = proportional_sample,
                         family ="gaussian",
                         options = list(control.inla = list(int.strategy = "eb")))
  
  b = mean(predict(fit_proportional,cov_effect,~(Intercept+covariate))$mean)
  
  fit_proportional_byc = bru(components = form_bycatch,
                             data = proportional_sample,
                             family ="gaussian",
                             options = list(control.inla = list(int.strategy = "eb")))
  
  bb = mean(predict(fit_proportional_byc,cov_effect,~(Intercept+covariate))$mean)
  
  fit_proportional_byc2 = bru(components = form_bycatch2,
                              data = proportional_sample,
                              family ="gaussian",
                              options = list(control.inla = list(int.strategy = "eb")))
  
  bbb = mean(predict(fit_proportional_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  #### Preferential sampling, squared abundance intensity ####
  fit_squared = bru(components = form_target,
                    data = squared_sample,
                    family ="gaussian",
                    options = list(control.inla = list(int.strategy = "eb")))
  
  c = mean(predict(fit_squared,cov_effect,~(Intercept+covariate))$mean)
  
  fit_squared_byc = bru(components = form_bycatch,
                        data = squared_sample,
                        family ="gaussian",
                        options = list(control.inla = list(int.strategy = "eb")))
  
  cc = mean(predict(fit_squared_byc,cov_effect,~(Intercept+covariate))$mean)
  
  
  fit_squared_byc2 = bru(components = form_bycatch2,
                         data = squared_sample,
                         family ="gaussian",
                         options = list(control.inla = list(int.strategy = "eb")))
  
  ccc = mean(predict(fit_squared_byc2,cov_effect,~(Intercept+covariate))$mean)
  
  #### Preferential sampling, 4th power abundance intensity ####
  fit_fourth_power = bru(components = form_target,
                         data = fourth_power_sample,
                         family ="gaussian",
                         options = list(control.inla = list(int.strategy = "eb")))
  
  d = mean(predict(fit_fourth_power,cov_effect,~(Intercept+covariate))$mean)
  
  fit_fourth_power_byc = bru(components = form_bycatch,
                             data = fourth_power_sample,
                             family ="gaussian",
                             options = list(control.inla = list(int.strategy = "eb")))
  
  dd = mean(predict(fit_fourth_power_byc,cov_effect,~(Intercept+covariate))$mean)
  
  fit_fourth_power_byc2 = bru(components = form_bycatch2,
                              data = fourth_power_sample,
                              family ="gaussian",
                              options = list(control.inla = list(int.strategy = "eb")))
  
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
    
    error[["fit_random"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_random$logfile) )
    error[["fit_squared"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_squared$logfile) )
    error[["fit_proportional"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_proportional$logfile) )
    error[["fit_fourth_power"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_fourth_power$logfile) )
    
    error[["fit_random_byc"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_random_byc$logfile) )
    error[["fit_squared_byc"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_squared_byc$logfile) )
    error[["fit_proportional_byc"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_proportional_byc$logfile) )
    error[["fit_fourth_power_byc"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_fourth_power_byc$logfile) )
    
    error[["fit_random_byc2"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_random_byc2$logfile) )
    error[["fit_squared_byc2"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_squared_byc2$logfile) )
    error[["fit_proportional_byc2"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_proportional_byc2$logfile) )
    error[["fit_fourth_power_byc2"]] = any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_fourth_power_byc2$logfile) )
    
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
    
    error[["fit_random"]] = c( error[["fit_random"]],any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_random$logfile) ))
    error[["fit_squared"]] = c(error[["fit_squared"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_squared$logfile) ))
    error[["fit_proportional"]] = c(error[["fit_proportional"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_proportional$logfile) ))
    error[["fit_fourth_power"]] = c(error[["fit_fourth_power"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_fourth_power$logfile) ))
    
    error[["fit_random_byc"]] = c(error[["fit_random_byc"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_random_byc$logfile) ))
    error[["fit_squared_byc"]] = c(error[["fit_squared_byc"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_squared_byc$logfile) ))
    error[["fit_proportional_byc"]] = c(error[["fit_proportional_byc"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_proportional_byc$logfile) ))
    error[["fit_fourth_power_byc"]] = c(error[["fit_fourth_power_byc"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_fourth_power_byc$logfile) ))
    
    error[["fit_random_byc2"]] = c(error[["fit_random_byc2"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_random_byc2$logfile) ))
    error[["fit_squared_byc2"]] = c(error[["fit_squared_byc2"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_squared_byc2$logfile) ))
    error[["fit_proportional_byc2"]] = c(error[["fit_proportional_byc2"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_proportional_byc2$logfile) ))
    error[["fit_fourth_power_byc2"]] = c(error[["fit_fourth_power_byc2"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = fit_fourth_power_byc2$logfile) ))
    
    
  }
  
  #########################################################################
  ############## Conventional/extended  preferential models ########
  ##############################################################
  
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
  cmp_ePM =  ~ covariate_biomass_err(covariate_space,model=matern) +
    covariate_biomass_copy(covariate_space, copy = "covariate_biomass", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    covariate_biomass(covariate_space,model=matern_zero) +
    Intercept_biomass(1) + Intercept_pp(1) # +
  
  cmp_cPM =  ~ covariate_biomass(covariate_space,model=matern) +
    covariate_biomass_copy(covariate_space, copy = "covariate_biomass", fixed = FALSE, hyper = list(beta = prior_to_zero)) +
    Intercept_biomass(1) + Intercept_pp(1) 
  
  
  #### formulas
  
  # sampling intensity
  form_pp_ePM = n_obs  ~ covariate_biomass_err +  covariate_biomass_copy +   Intercept_pp
  form_pp_cPM = n_obs  ~   covariate_biomass_copy +   Intercept_pp
  
  # preferntial model with fishers error
  form_target = biomass_target_FD ~     covariate_biomass +  Intercept_biomass
  form_bycatch = biomass_bycatch_FD ~    covariate_biomass +  Intercept_biomass
  form_bycatch2 = biomass_bycatch2_FD ~     covariate_biomass +  Intercept_biomass
  

  ####################
  ##### likeluhoods 
  
  ### sampling intensity
  lik_pp_cPM_random <- bru_obs("poisson",
                        formula = form_pp_cPM,
                        data = df_random_pp
  )
  
  
  lik_pp_cPM_proportional <- bru_obs("poisson",
                              formula = form_pp_cPM,
                              data = df_proportional_pp
  )
  
  lik_pp_cPM_squared <- bru_obs("poisson",
                         formula = form_pp_cPM,
                         data = df_squared_pp
  )
  
  lik_pp_cPM_4th_power <- bru_obs("poisson",
                           formula = form_pp_cPM,
                           data = df_4th_power_pp
  )
  
  lik_pp_ePM_random <- bru_obs("poisson",
                            formula = form_pp_ePM,
                            data = df_random_pp
  )
  
  
  lik_pp_ePM_proportional <- bru_obs("poisson",
                                  formula = form_pp_ePM,
                                  data = df_proportional_pp
  )
  
  lik_pp_ePM_squared <- bru_obs("poisson",
                             formula = form_pp_ePM,
                             data = df_squared_pp
  )
  
  lik_pp_ePM_4th_power <- bru_obs("poisson",
                               formula = form_pp_ePM,
                               data = df_4th_power_pp
  )
  
  
  #### preferential model with fishers error
  lik_target_random <- bru_obs("gaussian",
                                formula = form_target,
                                data = random_sample
  )
  
  lik_target_proportional <- bru_obs("gaussian",
                                      formula = form_target,
                                      data = proportional_sample
  )
  lik_target_squared <- bru_obs("gaussian",
                                 formula = form_target,
                                 data = squared_sample
  )
  lik_target_4th_power <- bru_obs("gaussian",
                                   formula = form_target,
                                   data = fourth_power_sample
  )
  
  
  lik_bycatch_random <- bru_obs("gaussian",
                                 formula = form_bycatch,
                                 data = random_sample
  )
  lik_bycatch_proportional <- bru_obs("gaussian",
                                       formula = form_bycatch,
                                       data = proportional_sample
  )
  lik_bycatch_squared <- bru_obs("gaussian",
                                  formula = form_bycatch,
                                  data = squared_sample
  )
  lik_bycatch_4th_power <- bru_obs("gaussian",
                                    formula = form_bycatch,
                                    data = fourth_power_sample
  )
  
  lik_bycatch2_random <- bru_obs("gaussian",
                                  formula = form_bycatch2,
                                  data = random_sample
  )
  lik_bycatch2_proportional <- bru_obs("gaussian",
                                        formula = form_bycatch2,
                                        data = proportional_sample
  )
  lik_bycatch2_squared <- bru_obs("gaussian",
                                   formula = form_bycatch2,
                                   data = squared_sample
  )
  lik_bycatch2_4th_power <- bru_obs("gaussian",
                                     formula = form_bycatch2,
                                     data = fourth_power_sample
  )
  
  
  
  ###################
  ### Fit models ####
  ###################
  
  #### preferential models with fishers error
  ePM_random_fit <- bru(cmp_ePM, lik_pp_ePM_random, lik_target_random,
                        options = list(control.inla = list(int.strategy = "eb")))
  ePM_proportional_fit <- bru(cmp_ePM, lik_pp_ePM_proportional, lik_target_proportional,
                              options = list(control.inla = list(int.strategy = "eb")))
  ePM_squared_fit <- bru(cmp_ePM, lik_pp_ePM_squared, lik_target_squared,
                         options = list(control.inla = list(int.strategy = "eb")))
  ePM_4th_power_fit <- bru(cmp_ePM, lik_pp_ePM_4th_power, lik_target_4th_power,
                           options = list(control.inla = list(int.strategy = "eb")))
  
  ePM_byc_random_fit <- bru(cmp_ePM, lik_pp_ePM_random, lik_bycatch_random,
                            options = list(control.inla = list(int.strategy = "eb")))
  ePM_byc_proportional_fit <- bru(cmp_ePM, lik_pp_ePM_proportional, lik_bycatch_proportional,
                                  options = list(control.inla = list(int.strategy = "eb")))
  ePM_byc_squared_fit <- bru(cmp_ePM, lik_pp_ePM_squared, lik_bycatch_squared,
                             options = list(control.inla = list(int.strategy = "eb")))
  ePM_byc_4th_power_fit <- bru(cmp_ePM, lik_pp_ePM_4th_power, lik_bycatch_4th_power,
                               options = list(control.inla = list(int.strategy = "eb")))
  
  ePM_byc2_random_fit <- bru(cmp_ePM, lik_pp_ePM_random, lik_bycatch2_random,
                             options = list(control.inla = list(int.strategy = "eb")))
  ePM_byc2_proportional_fit <- bru(cmp_ePM, lik_pp_ePM_proportional, lik_bycatch2_proportional,
                                   options = list(control.inla = list(int.strategy = "eb")))
  ePM_byc2_squared_fit <- bru(cmp_ePM, lik_pp_ePM_squared, lik_bycatch2_squared,
                              options = list(control.inla = list(int.strategy = "eb")))
  ePM_byc2_4th_power_fit <- bru(cmp_ePM, lik_pp_ePM_4th_power, lik_bycatch2_4th_power,
                                options = list(control.inla = list(int.strategy = "eb")))
  
  
  
  #### traditional preferential models 
  cPM_random_fit <- bru(cmp_cPM, lik_pp_cPM_random, lik_target_random,
                        options = list(control.inla = list(int.strategy = "eb")))
  cPM_proportional_fit <- bru(cmp_cPM, lik_pp_cPM_proportional, lik_target_proportional,
                              options = list(control.inla = list(int.strategy = "eb")))
  cPM_squared_fit <- bru(cmp_cPM, lik_pp_cPM_squared, lik_target_squared,
                         options = list(control.inla = list(int.strategy = "eb")))
  cPM_4th_power_fit <- bru(cmp_cPM, lik_pp_cPM_4th_power, lik_target_4th_power,
                           options = list(control.inla = list(int.strategy = "eb")))
  
  cPM_byc_random_fit <- bru(cmp_cPM, lik_pp_cPM_random, lik_bycatch_random,
                            options = list(control.inla = list(int.strategy = "eb")))
  cPM_byc_proportional_fit <- bru(cmp_cPM, lik_pp_cPM_proportional, lik_bycatch_proportional,
                                  options = list(control.inla = list(int.strategy = "eb")))
  cPM_byc_squared_fit <- bru(cmp_cPM, lik_pp_cPM_squared, lik_bycatch_squared,
                             options = list(control.inla = list(int.strategy = "eb")))
  cPM_byc_4th_power_fit <- bru(cmp_cPM, lik_pp_cPM_4th_power, lik_bycatch_4th_power,
                               options = list(control.inla = list(int.strategy = "eb")))
  
  cPM_byc2_random_fit <- bru(cmp_cPM, lik_pp_cPM_random, lik_bycatch2_random,
                             options = list(control.inla = list(int.strategy = "eb")))
  cPM_byc2_proportional_fit <- bru(cmp_cPM, lik_pp_cPM_proportional, lik_bycatch2_proportional,
                                   options = list(control.inla = list(int.strategy = "eb")))
  cPM_byc2_squared_fit <- bru(cmp_cPM, lik_pp_cPM_squared, lik_bycatch2_squared,
                              options = list(control.inla = list(int.strategy = "eb")))
  cPM_byc2_4th_power_fit <- bru(cmp_cPM, lik_pp_cPM_4th_power, lik_bycatch2_4th_power,
                                options = list(control.inla = list(int.strategy = "eb")))
  

#   
# #### calculate quantiles@zero
  
  
  
#   q0_beta_a =inla.pmarginal(c(0),ePM_random_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_b =inla.pmarginal(c(0),ePM_proportional_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_c =inla.pmarginal(c(0),ePM_squared_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_d =inla.pmarginal(c(0),ePM_4th_power_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   
#   q0_beta_aa =inla.pmarginal(c(0),ePM_byc_random_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_bb =inla.pmarginal(c(0),ePM_byc_proportional_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_cc =inla.pmarginal(c(0),ePM_byc_squared_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_dd =inla.pmarginal(c(0),ePM_byc_4th_power_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   
#   q0_beta_aaa =inla.pmarginal(c(0),ePM_byc2_random_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_bbb =inla.pmarginal(c(0),ePM_byc2_proportional_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_ccc =inla.pmarginal(c(0),ePM_byc2_squared_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   q0_beta_ddd =inla.pmarginal(c(0),ePM_byc2_4th_power_fit$marginals.hyperpar$`Beta for covariate_pp_copy`)
#   
  
  ######################################
  ###### Estimated mean abundance #####
  ####################################
  
  ######## preferential model wirh fishers error
  ePM_a = mean(predict(ePM_random_fit,cov_effect,
                     ~(Intercept_biomass+ covariate_biomass ),
                     include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_b = mean(predict(ePM_proportional_fit,cov_effect,
                     ~(Intercept_biomass + covariate_biomass ),
                     include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_c = mean(predict(ePM_squared_fit,cov_effect,
                     ~(Intercept_biomass + covariate_biomass ),
                     include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_d = mean(predict(ePM_4th_power_fit,cov_effect,
                     ~(Intercept_biomass + covariate_biomass ),
                     include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_aa = mean(predict(ePM_byc_random_fit,cov_effect,
                      ~(Intercept_biomass + covariate_biomass ),
                      include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_bb = mean(predict(ePM_byc_proportional_fit,cov_effect,
                      ~(Intercept_biomass + covariate_biomass ),
                      include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_cc = mean(predict(ePM_byc_squared_fit,cov_effect,
                      ~(Intercept_biomass + covariate_biomass ),
                      include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_dd = mean(predict(ePM_byc_4th_power_fit,cov_effect,
                      ~(Intercept_biomass + covariate_biomass ),
                      include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  
  ePM_aaa = mean(predict(ePM_byc2_random_fit,cov_effect,
                       ~(Intercept_biomass + covariate_biomass ),
                       include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_bbb = mean(predict(ePM_byc2_proportional_fit,cov_effect,
                       ~(Intercept_biomass + covariate_biomass ),
                       include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_ccc = mean(predict(ePM_byc2_squared_fit,cov_effect,
                       ~(Intercept_biomass + covariate_biomass ),
                       include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ePM_ddd = mean(predict(ePM_byc2_4th_power_fit,cov_effect,
                       ~(Intercept_biomass + covariate_biomass ),
                       include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  
  
  
  ######### Conventional preferential models
  a_cPM = mean(predict(cPM_random_fit,cov_effect,
                          ~(Intercept_biomass+covariate_biomass ), 
                          include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  b_cPM = mean(predict(cPM_proportional_fit,cov_effect,
                          ~(Intercept_biomass+covariate_biomass ), 
                          include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  c_cPM = mean(predict(cPM_squared_fit,cov_effect,
                          ~(Intercept_biomass+covariate_biomass ), 
                          include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  d_cPM = mean(predict(cPM_4th_power_fit,cov_effect,
                          ~(Intercept_biomass+covariate_biomass ), 
                          include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  aa_cPM = mean(predict(cPM_byc_random_fit,cov_effect,
                           ~(Intercept_biomass+covariate_biomass ), 
                           include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  bb_cPM = mean(predict(cPM_byc_proportional_fit,cov_effect,
                           ~(Intercept_biomass+covariate_biomass ), 
                           include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  cc_cPM = mean(predict(cPM_byc_squared_fit,cov_effect,
                           ~(Intercept_biomass+covariate_biomass ), 
                           include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  dd_cPM = mean(predict(cPM_byc_4th_power_fit,cov_effect,
                           ~(Intercept_biomass+covariate_biomass ), 
                           include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  
  aaa_cPM = mean(predict(cPM_byc2_random_fit,cov_effect,
                            ~(Intercept_biomass+covariate_biomass ), 
                            include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  bbb_cPM = mean(predict(cPM_byc2_proportional_fit,cov_effect,
                            ~(Intercept_biomass+covariate_biomass ), 
                            include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ccc_cPM = mean(predict(cPM_byc2_squared_fit,cov_effect,
                            ~(Intercept_biomass+covariate_biomass ), 
                            include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  ddd_cPM = mean(predict(cPM_byc2_4th_power_fit,cov_effect,
                            ~(Intercept_biomass+covariate_biomass ), 
                            include=c("Intercept_biomass","covariate_biomass"))$mean)
  
  
  
  if(i==1){
    ePM_error_random = mean_biomass_target_FD - ePM_a
    ePM_error_proportional= mean_biomass_target_FD-ePM_b
    ePM_error_squared = mean_biomass_target_FD - ePM_c
    ePM_error_4th_power = mean_biomass_target_FD - ePM_d
    
    ePM_byc_error_random = mean_biomass_byc_FD - ePM_aa
    ePM_byc_error_proportional= mean_biomass_byc_FD - ePM_bb
    ePM_byc_error_squared = mean_biomass_byc_FD - ePM_cc
    ePM_byc_error_4th_power = mean_biomass_byc_FD - ePM_dd
    
    ePM_byc2_error_random = mean_biomass_byc2_FD - ePM_aaa
    ePM_byc2_error_proportional= mean_biomass_byc2_FD - ePM_bbb
    ePM_byc2_error_squared = mean_biomass_byc2_FD - ePM_ccc
    ePM_byc2_error_4th_power = mean_biomass_byc2_FD - ePM_ddd
    
    ##### pref simple
    cPM_error_random = mean_biomass_target_FD - a_cPM
    cPM_error_proportional= mean_biomass_target_FD-b_cPM
    cPM_error_squared = mean_biomass_target_FD - c_cPM
    cPM_error_4th_power = mean_biomass_target_FD - d_cPM
    
    cPM_byc_error_random = mean_biomass_byc_FD - aa_cPM
    cPM_byc_error_proportional= mean_biomass_byc_FD - bb_cPM
    cPM_byc_error_squared = mean_biomass_byc_FD - cc_cPM
    cPM_byc_error_4th_power = mean_biomass_byc_FD - dd_cPM
    
    cPM_byc2_error_random = mean_biomass_byc2_FD - aaa_cPM
    cPM_byc2_error_proportional= mean_biomass_byc2_FD - bbb_cPM
    cPM_byc2_error_squared = mean_biomass_byc2_FD - ccc_cPM
    cPM_byc2_error_4th_power = mean_biomass_byc2_FD - ddd_cPM
    
    ### alpha coefficient
    #cPM
    beta_target_random_cPM = cPM_random_fit$summary.hyperpar[4,c(1:3,5)]
    beta_target_proportional_cPM = cPM_proportional_fit$summary.hyperpar[4,c(1:3,5)]
    beta_target_squared_cPM = cPM_squared_fit$summary.hyperpar[4,c(1:3,5)]
    beta_target_4th_power_cPM = cPM_4th_power_fit$summary.hyperpar[4,c(1:3,5)]

    beta_byc_random_cPM = cPM_byc_random_fit$summary.hyperpar[4,c(1:3,5)]
    beta_byc_proportional_cPM = cPM_byc_proportional_fit$summary.hyperpar[4,c(1:3,5)]
    beta_byc_squared_cPM = cPM_byc_squared_fit$summary.hyperpar[4,c(1:3,5)]
    beta_byc_4th_power_cPM = cPM_byc_4th_power_fit$summary.hyperpar[4,c(1:3,5)]

    beta_byc2_random_cPM = cPM_byc2_random_fit$summary.hyperpar[4,c(1:3,5)]
    beta_byc2_proportional_cPM = cPM_byc2_proportional_fit$summary.hyperpar[4,c(1:3,5)]
    beta_byc2_squared_cPM = cPM_byc2_squared_fit$summary.hyperpar[4,c(1:3,5)]
    beta_byc2_4th_power_cPM = cPM_byc2_4th_power_fit$summary.hyperpar[4,c(1:3,5)]
    
    #ePM
    beta_target_random_ePM = ePM_random_fit$summary.hyperpar[6,c(1:3,5)]
    beta_target_proportional_ePM = ePM_proportional_fit$summary.hyperpar[6,c(1:3,5)]
    beta_target_squared_ePM = ePM_squared_fit$summary.hyperpar[6,c(1:3,5)]
    beta_target_4th_power_ePM = ePM_4th_power_fit$summary.hyperpar[6,c(1:3,5)]
    
    beta_byc_random_ePM = ePM_byc_random_fit$summary.hyperpar[6,c(1:3,5)]
    beta_byc_proportional_ePM = ePM_byc_proportional_fit$summary.hyperpar[6,c(1:3,5)]
    beta_byc_squared_ePM = ePM_byc_squared_fit$summary.hyperpar[6,c(1:3,5)]
    beta_byc_4th_power_ePM = ePM_byc_4th_power_fit$summary.hyperpar[6,c(1:3,5)]
    
    beta_byc2_random_ePM = ePM_byc2_random_fit$summary.hyperpar[6,c(1:3,5)]
    beta_byc2_proportional_ePM = ePM_byc2_proportional_fit$summary.hyperpar[6,c(1:3,5)]
    beta_byc2_squared_ePM = ePM_byc2_squared_fit$summary.hyperpar[6,c(1:3,5)]
    beta_byc2_4th_power_ePM = ePM_byc2_4th_power_fit$summary.hyperpar[6,c(1:3,5)]
    
    error[["cPM_random_fit"]] = c( any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_random_fit$logfile) ))
    error[["cPM_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_proportional_fit$logfile) ))
    error[["cPM_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_squared_fit$logfile) ))
    error[["cPM_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_4th_power_fit$logfile) ))
    
    error[["cPM_byc_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_random_fit$logfile) ))
    error[["cPM_byc_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_proportional_fit$logfile) ))
    error[["cPM_byc_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_squared_fit$logfile) ))
    error[["cPM_byc_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_4th_power_fit$logfile) ))
    
    error[["cPM_byc2_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_random_fit$logfile) ))
    error[["cPM_byc2_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_proportional_fit$logfile) ))
    error[["cPM_byc2_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_squared_fit$logfile) ))
    error[["cPM_byc2_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_4th_power_fit$logfile) ))
    
    error[["ePM_random_fit"]] = c( any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_random_fit$logfile) ))
    error[["ePM_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_proportional_fit$logfile) ))
    error[["ePM_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_squared_fit$logfile) ))
    error[["ePM_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_4th_power_fit$logfile) ))
    
    error[["ePM_byc_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_random_fit$logfile) ))
    error[["ePM_byc_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_proportional_fit$logfile) ))
    error[["ePM_byc_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_squared_fit$logfile) ))
    error[["ePM_byc_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_4th_power_fit$logfile) ))
    
    error[["ePM_byc2_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_random_fit$logfile) ))
    error[["ePM_byc2_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_proportional_fit$logfile) ))
    error[["ePM_byc2_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_squared_fit$logfile) ))
    error[["ePM_byc2_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_4th_power_fit$logfile) ))
    
    ##### quantile at 0
    # q0_beta_random = q0_beta_a
    # q0_beta_proportional = q0_beta_b
    # q0_beta_squared = q0_beta_c
    # q0_beta_4th_power = q0_beta_d
    # 
    # q0_beta_byc_random = q0_beta_aa
    # q0_beta_byc_proportional = q0_beta_bb
    # q0_beta_byc_squared = q0_beta_cc
    # q0_beta_byc_4th_power = q0_beta_dd
    # 
    # q0_beta_byc2_random = q0_beta_aaa
    # q0_beta_byc2_proportional = q0_beta_bbb
    # q0_beta_byc2_squared = q0_beta_ccc
    # q0_beta_byc2_4th_power = q0_beta_ddd
    
  }else{
    ePM_error_random=c(ePM_error_random,mean_biomass_target_FD - ePM_a)
    ePM_error_proportional=c(ePM_error_proportional,mean_biomass_target_FD - ePM_b)
    ePM_error_squared =c(ePM_error_squared,mean_biomass_target_FD - ePM_c)
    ePM_error_4th_power = c(ePM_error_4th_power,mean_biomass_target_FD - ePM_d)
    
    ePM_byc_error_random=c(ePM_byc_error_random,mean_biomass_byc_FD -ePM_aa)
    ePM_byc_error_proportional=c(ePM_byc_error_proportional,mean_biomass_byc_FD -ePM_bb)
    ePM_byc_error_squared =c(ePM_byc_error_squared,mean_biomass_byc_FD -ePM_cc)
    ePM_byc_error_4th_power = c(ePM_byc_error_4th_power,mean_biomass_byc_FD -ePM_dd)
    
    ePM_byc2_error_random=c(ePM_byc2_error_random,mean_biomass_byc2_FD -ePM_aaa)
    ePM_byc2_error_proportional=c(ePM_byc2_error_proportional,mean_biomass_byc2_FD -ePM_bbb)
    ePM_byc2_error_squared =c(ePM_byc2_error_squared,mean_biomass_byc2_FD -ePM_ccc)
    ePM_byc2_error_4th_power = c(ePM_byc2_error_4th_power,mean_biomass_byc2_FD -ePM_ddd)
    
    #### simple preferential
    cPM_error_random=c(cPM_error_random,mean_biomass_target_FD - a_cPM)
    cPM_error_proportional=c(cPM_error_proportional,mean_biomass_target_FD - b_cPM)
    cPM_error_squared =c(cPM_error_squared,mean_biomass_target_FD - c_cPM)
    cPM_error_4th_power = c(cPM_error_4th_power,mean_biomass_target_FD - d_cPM)
    
    cPM_byc_error_random=c(cPM_byc_error_random,mean_biomass_byc_FD -aa_cPM)
    cPM_byc_error_proportional=c(cPM_byc_error_proportional,mean_biomass_byc_FD -bb_cPM)
    cPM_byc_error_squared =c(cPM_byc_error_squared,mean_biomass_byc_FD -cc_cPM)
    cPM_byc_error_4th_power = c(cPM_byc_error_4th_power,mean_biomass_byc_FD -dd_cPM)
    
    cPM_byc2_error_random=c(cPM_byc2_error_random,mean_biomass_byc2_FD -aaa_cPM)
    cPM_byc2_error_proportional=c(cPM_byc2_error_proportional,mean_biomass_byc2_FD -bbb_cPM)
    cPM_byc2_error_squared =c(cPM_byc2_error_squared,mean_biomass_byc2_FD -ccc_cPM)
    cPM_byc2_error_4th_power = c(cPM_byc2_error_4th_power,mean_biomass_byc2_FD -ddd_cPM)
    
    ##### alphas
    # cPM
    beta_target_random_cPM = rbind(beta_target_random_cPM,cPM_random_fit$summary.hyperpar[4,c(1:3,5)])
    beta_target_proportional_cPM = rbind(beta_target_proportional_cPM,cPM_proportional_fit$summary.hyperpar[4,c(1:3,5)])
    beta_target_squared_cPM = rbind(beta_target_squared_cPM,cPM_squared_fit$summary.hyperpar[4,c(1:3,5)])
    beta_target_4th_power_cPM = rbind(beta_target_4th_power_cPM,cPM_4th_power_fit$summary.hyperpar[4,c(1:3,5)])
    
    beta_byc_random_cPM = rbind(beta_byc_random_cPM,cPM_byc_random_fit$summary.hyperpar[4,c(1:3,5)])
    beta_byc_proportional_cPM = rbind(beta_byc_proportional_cPM,cPM_byc_proportional_fit$summary.hyperpar[4,c(1:3,5)])
    beta_byc_squared_cPM = rbind(beta_byc_squared_cPM,cPM_byc_squared_fit$summary.hyperpar[4,c(1:3,5)])
    beta_byc_4th_power_cPM = rbind(beta_byc_4th_power_cPM,cPM_byc_4th_power_fit$summary.hyperpar[4,c(1:3,5)])
    
    beta_byc2_random_cPM = rbind(beta_byc2_random_cPM,cPM_byc2_random_fit$summary.hyperpar[4,c(1:3,5)])
    beta_byc2_proportional_cPM = rbind(beta_byc2_proportional_cPM,cPM_byc2_proportional_fit$summary.hyperpar[4,c(1:3,5)])
    beta_byc2_squared_cPM = rbind(beta_byc2_squared_cPM,cPM_byc2_squared_fit$summary.hyperpar[4,c(1:3,5)])
    beta_byc2_4th_power_cPM = rbind(beta_byc2_4th_power_cPM,cPM_byc2_4th_power_fit$summary.hyperpar[4,c(1:3,5)])
    
    #ePM
    beta_target_random_ePM = rbind(beta_target_random_ePM,ePM_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_proportional_ePM = rbind(beta_target_proportional_ePM,ePM_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_squared_ePM = rbind(beta_target_squared_ePM,ePM_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_4th_power_ePM = rbind(beta_target_4th_power_ePM,ePM_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    beta_byc_random_ePM = rbind(beta_byc_random_ePM,ePM_byc_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_proportional_ePM = rbind(beta_byc_proportional_ePM,ePM_byc_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_squared_ePM = rbind(beta_byc_squared_ePM,ePM_byc_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_4th_power_ePM = rbind(beta_byc_4th_power_ePM,ePM_byc_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    beta_byc2_random_ePM = rbind(beta_byc2_random_ePM,ePM_byc2_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_proportional_ePM = rbind(beta_byc2_proportional_ePM,ePM_byc2_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_squared_ePM = rbind(beta_byc2_squared_ePM,ePM_byc2_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_4th_power_ePM = rbind(beta_byc2_4th_power_ePM,ePM_byc2_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    error[["cPM_random_fit"]] = c( error[["cPM_random_fit"]] , any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_random_fit$logfile) ))
    error[["cPM_proportional_fit"]] = c( error[["cPM_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_proportional_fit$logfile) ))
    error[["cPM_squared_fit"]] = c(error[["cPM_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_squared_fit$logfile) ))
    error[["cPM_4th_power_fit"]] = c(error[["cPM_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_4th_power_fit$logfile) ))
    
    error[["cPM_byc_random_fit"]] = c(error[["cPM_byc_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_random_fit$logfile) ))
    error[["cPM_byc_proportional_fit"]] = c(error[["cPM_byc_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_proportional_fit$logfile) ))
    error[["cPM_byc_squared_fit"]] = c(error[["cPM_byc_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_squared_fit$logfile) ))
    error[["cPM_byc_4th_power_fit"]] = c(error[["cPM_byc_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc_4th_power_fit$logfile) ))
    
    error[["cPM_byc2_random_fit"]] = c(error[["cPM_byc2_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_random_fit$logfile) ))
    error[["cPM_byc2_proportional_fit"]] = c(error[["cPM_byc2_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_proportional_fit$logfile) ))
    error[["cPM_byc2_squared_fit"]] = c(error[["cPM_byc2_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_squared_fit$logfile) ))
    error[["cPM_byc2_4th_power_fit"]] = c(error[["cPM_byc2_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = cPM_byc2_4th_power_fit$logfile) ))
    
    error[["ePM_random_fit"]] = c(error[["ePM_random_fit"]] , any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_random_fit$logfile) ))
    error[["ePM_proportional_fit"]] = c(error[["ePM_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_proportional_fit$logfile) ))
    error[["ePM_squared_fit"]] = c(error[["ePM_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_squared_fit$logfile) ))
    error[["ePM_4th_power_fit"]] = c(error[["ePM_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_4th_power_fit$logfile) ))
    
    error[["ePM_byc_random_fit"]] = c(error[["ePM_byc_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_random_fit$logfile) ))
    error[["ePM_byc_proportional_fit"]] = c(error[["ePM_byc_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_proportional_fit$logfile) ))
    error[["ePM_byc_squared_fit"]] = c(error[["ePM_byc_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_squared_fit$logfile) ))
    error[["ePM_byc_4th_power_fit"]] = c(error[["ePM_byc_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc_4th_power_fit$logfile) ))
    
    error[["ePM_byc2_random_fit"]] = c(error[["ePM_byc2_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_random_fit$logfile) ))
    error[["ePM_byc2_proportional_fit"]] = c(error[["ePM_byc2_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_proportional_fit$logfile) ))
    error[["ePM_byc2_squared_fit"]] = c(error[["ePM_byc2_squared_fit"]]  ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_squared_fit$logfile) ))
    error[["ePM_byc2_4th_power_fit"]] = c(error[["ePM_byc2_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ePM_byc2_4th_power_fit$logfile) ))
    
    
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
  
  
  ###### cmp
  cmp_ISDM =  ~ covariate_biomass(covariate_space,model=matern) +
    covariate_biomass_copy(covariate_space, copy = "covariate_biomass", fixed = F, hyper = list(beta = prior_to_zero)) +
    Intercept_FI(1) + Intercept_FD(1) # +

  ###### formulas
  form_FI_target = biomass_target_FI ~   covariate_biomass +  Intercept_FI
  form_FI_bycatch = biomass_bycatch_FI ~   covariate_biomass +  Intercept_FI
  form_FI_bycatch2 = biomass_bycatch2_FI ~   covariate_biomass +  Intercept_FI
  
  form_FD_target = biomass_target_FD ~   covariate_biomass_copy  +  Intercept_FD
  form_FD_bycatch = biomass_bycatch_FD ~   covariate_biomass_copy  +  Intercept_FD
  form_FD_bycatch2 = biomass_bycatch2_FD ~   covariate_biomass_copy  +  Intercept_FD
  
  
  ###### likelihoods
  
  ## FI
  lik_biomass_FI <- bru_obs("gaussian",formula = form_FI_target,data = FI_data )
  lik_biomass_FI_bycatch <- bru_obs("gaussian",formula = form_FI_bycatch,data = FI_data  )
  lik_biomass_FI_bycatch2 <- bru_obs("gaussian", formula = form_FI_bycatch2,data = FI_data  )
  
  ## FD
  lik_biomass_FD_random <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_random  )
  lik_biomass_FD_proportional <- bru_obs("gaussian",formula = form_FD_target,  data = FD_data_proportional  )
  lik_biomass_FD_squared <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_squared  )
  lik_biomass_FD_fourth <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_fourth  )
  
  lik_biomass_FD_random_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_random  )
  lik_biomass_FD_proportional_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_proportional  )
  lik_biomass_FD_squared_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_squared  )
  lik_biomass_FD_fourth_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_fourth  )
  
  lik_biomass_FD_random_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2, data = FD_data_random)
  lik_biomass_FD_proportional_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2,data = FD_data_proportional )
  lik_biomass_FD_squared_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2,data = FD_data_squared)
  lik_biomass_FD_fourth_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2,     data = FD_data_fourth)
  
  #### target species
  ISDM_target_random_fit <- bru(cmp_ISDM, lik_biomass_FI, lik_biomass_FD_random,
                                options = list(control.inla = list(int.strategy = "eb")))
  ISDM_target_proportional_fit <- bru(cmp_ISDM, lik_biomass_FI, lik_biomass_FD_proportional,
                                      options = list(control.inla = list(int.strategy = "eb")))
  ISDM_target_squared_fit <- bru(cmp_ISDM, lik_biomass_FI, lik_biomass_FD_squared,
                                 options = list(control.inla = list(int.strategy = "eb")))
  ISDM_target_4th_power_fit <- bru(cmp_ISDM, lik_biomass_FI, lik_biomass_FD_fourth,
                                   options = list(control.inla = list(int.strategy = "eb")))
  
  #### bycatch species
  ISDM_bycatch_random_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch, lik_biomass_FD_random_bycatch,
                                 options = list(control.inla = list(int.strategy = "eb")))
  ISDM_bycatch_proportional_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch, lik_biomass_FD_proportional_bycatch,
                                       options = list(control.inla = list(int.strategy = "eb")))
  ISDM_bycatch_squared_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch, lik_biomass_FD_squared_bycatch,
                                  options = list(control.inla = list(int.strategy = "eb")))
  ISDM_bycatch_4th_power_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch, lik_biomass_FD_fourth_bycatch,
                                    options = list(control.inla = list(int.strategy = "eb")))
  
  #### bycatch2 species
  ISDM_bycatch2_random_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch2, lik_biomass_FD_random_bycatch2,
                                  options = list(control.inla = list(int.strategy = "eb")))
  ISDM_bycatch2_proportional_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch2, lik_biomass_FD_proportional_bycatch2,
                                        options = list(control.inla = list(int.strategy = "eb")))
  ISDM_bycatch2_squared_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch2, lik_biomass_FD_squared_bycatch2,
                                   options = list(control.inla = list(int.strategy = "eb")))
  ISDM_bycatch2_4th_power_fit <- bru(cmp_ISDM, lik_biomass_FI_bycatch2, lik_biomass_FD_fourth_bycatch2,
                                     options = list(control.inla = list(int.strategy = "eb")))
  
  
  a_ISDM = mean(predict(ISDM_target_random_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass_copy), 
                        include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  b_ISDM = mean(predict(ISDM_target_proportional_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass_copy), 
                        include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  c_ISDM = mean(predict(ISDM_target_squared_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass_copy), 
                        include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  d_ISDM = mean(predict(ISDM_target_4th_power_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass_copy), 
                        include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  
  aa_ISDM = mean(predict(ISDM_bycatch_random_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass_copy), 
                         include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  bb_ISDM = mean(predict(ISDM_bycatch_proportional_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass_copy), 
                         include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  cc_ISDM = mean(predict(ISDM_bycatch_squared_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass_copy), 
                         include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  dd_ISDM = mean(predict(ISDM_bycatch_4th_power_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass_copy), 
                         include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  
  aaa_ISDM = mean(predict(ISDM_bycatch2_random_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass_copy), 
                          include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  bbb_ISDM = mean(predict(ISDM_bycatch2_proportional_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass_copy), 
                          include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  ccc_ISDM = mean(predict(ISDM_bycatch2_squared_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass_copy), 
                          include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  ddd_ISDM = mean(predict(ISDM_bycatch2_4th_power_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass_copy), 
                          include=c("Intercept_FD","covariate_biomass_copy"))$mean)
  
  
  if(i==1){
    ISDM_error_random = mean_biomass_target_FD - a_ISDM
    ISDM_error_proportional= mean_biomass_target_FD-b_ISDM
    ISDM_error_squared = mean_biomass_target_FD - c_ISDM
    ISDM_error_4th_power = mean_biomass_target_FD - d_ISDM
    
    ISDM_byc_error_random = mean_biomass_byc_FD - aa_ISDM
    ISDM_byc_error_proportional= mean_biomass_byc_FD - bb_ISDM
    ISDM_byc_error_squared = mean_biomass_byc_FD - cc_ISDM
    ISDM_byc_error_4th_power = mean_biomass_byc_FD - dd_ISDM
    
    ISDM_byc2_error_random = mean_biomass_byc2_FD - aaa_ISDM
    ISDM_byc2_error_proportional= mean_biomass_byc2_FD - bbb_ISDM
    ISDM_byc2_error_squared = mean_biomass_byc2_FD - ccc_ISDM
    ISDM_byc2_error_4th_power = mean_biomass_byc2_FD - ddd_ISDM
    
    
    error[["ISDM_target_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_random_fit$logfile) ))
    error[["ISDM_target_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_proportional_fit$logfile) ))
    error[["ISDM_target_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_squared_fit$logfile) ))
    error[["ISDM_target_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_4th_power_fit$logfile) ))
    
    error[["ISDM_bycatch_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_random_fit$logfile) ))
    error[["ISDM_bycatch_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_proportional_fit$logfile) ))
    error[["ISDM_bycatch_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_squared_fit$logfile) ))
    error[["ISDM_bycatch_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_4th_power_fit$logfile) ))
    
    error[["ISDM_bycatch2_random_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_random_fit$logfile) ))
    error[["ISDM_bycatch2_proportional_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_proportional_fit$logfile) ))
    error[["ISDM_bycatch2_squared_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_squared_fit$logfile) ))
    error[["ISDM_bycatch2_4th_power_fit"]] = c(any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_4th_power_fit$logfile) ))
    
    
  }else{
    ISDM_error_random=c(ISDM_error_random,mean_biomass_target_FD - a_ISDM)
    ISDM_error_proportional=c(ISDM_error_proportional,mean_biomass_target_FD - b_ISDM)
    ISDM_error_squared =c(ISDM_error_squared,mean_biomass_target_FD - c_ISDM)
    ISDM_error_4th_power = c(ISDM_error_4th_power,mean_biomass_target_FD - d_ISDM)
    
    ISDM_byc_error_random=c(ISDM_byc_error_random,mean_biomass_byc_FD -aa_ISDM)
    ISDM_byc_error_proportional=c(ISDM_byc_error_proportional,mean_biomass_byc_FD -bb_ISDM)
    ISDM_byc_error_squared =c(ISDM_byc_error_squared,mean_biomass_byc_FD -cc_ISDM)
    ISDM_byc_error_4th_power = c(ISDM_byc_error_4th_power,mean_biomass_byc_FD -dd_ISDM)
    
    ISDM_byc2_error_random=c(ISDM_byc2_error_random,mean_biomass_byc2_FD -aaa_ISDM)
    ISDM_byc2_error_proportional=c(ISDM_byc2_error_proportional,mean_biomass_byc2_FD -bbb_ISDM)
    ISDM_byc2_error_squared =c(ISDM_byc2_error_squared,mean_biomass_byc2_FD -ccc_ISDM)
    ISDM_byc2_error_4th_power = c(ISDM_byc2_error_4th_power,mean_biomass_byc2_FD -ddd_ISDM)
    
    
    error[["ISDM_target_random_fit"]] = c(error[["ISDM_target_random_fit"]] , any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_random_fit$logfile) ))
    error[["ISDM_target_proportional_fit"]] = c(error[["ISDM_target_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_proportional_fit$logfile) ))
    error[["ISDM_target_squared_fit"]] = c(error[["ISDM_target_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_squared_fit$logfile) ))
    error[["ISDM_target_4th_power_fit"]] = c(error[["ISDM_target_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_target_4th_power_fit$logfile) ))
    
    error[["ISDM_bycatch_random_fit"]] = c(error[["ISDM_bycatch_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_random_fit$logfile) ))
    error[["ISDM_bycatch_proportional_fit"]] = c(error[["ISDM_bycatch_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_proportional_fit$logfile) ))
    error[["ISDM_bycatch_squared_fit"]] = c(error[["ISDM_bycatch_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_squared_fit$logfile) ))
    error[["ISDM_bycatch_4th_power_fit"]] = c(error[["ISDM_bycatch_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch_4th_power_fit$logfile) ))
    
    error[["ISDM_bycatch2_random_fit"]] = c(error[["ISDM_bycatch2_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_random_fit$logfile) ))
    error[["ISDM_bycatch2_proportional_fit"]] = c(error[["ISDM_bycatch2_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_proportional_fit$logfile) ))
    error[["ISDM_bycatch2_squared_fit"]] = c(error[["ISDM_bycatch2_squared_fit"]]  ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_squared_fit$logfile) ))
    error[["ISDM_bycatch2_4th_power_fit"]] = c(error[["ISDM_bycatch2_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = ISDM_bycatch2_4th_power_fit$logfile) ))
    
    
  }
  
  
  
  
  
  
  
  
  
  
  ##############################################################
  ##### "mixture models"  HERE ########################
  ##################################################
  
  ### model components
  cmp_iPM =  ~ 
    covariate_biomass(covariate_space,model=matern) +
    covariate_biomass_copy_FI(covariate_space, copy = "covariate_biomass", fixed = F, hyper = list(beta = prior_to_zero)) +
    covariate_biomass_copy_pp(covariate_space, copy = "covariate_biomass", fixed = F, hyper = list(beta = prior_to_zero)) +
    Intercept_FD(1) + Intercept_pp(1) +  Intercept_FI(1) # +


  cmp_iEPM =  ~ 
    covariate_biomass(covariate_space,model=matern) +
    covariate_biomass_copy_FI(covariate_space, copy = "covariate_biomass", fixed = F, hyper = list(beta = prior_to_zero)) +
    covariate_biomass_copy_pp(covariate_space, copy = "covariate_biomass", fixed = F, hyper = list(beta = prior_to_zero)) +
    covariate_err(covariate_space,model=matern) +
    Intercept_FD(1) + Intercept_pp(1) +  Intercept_FI(1)# +
  
  ###### formulas
  form_pp_iPM = n_obs  ~   covariate_biomass_copy_pp +   Intercept_pp
  form_pp_iEPM = n_obs  ~   covariate_biomass_copy_pp +   Intercept_pp + covariate_err
  
  form_FI_target = biomass_target_FI ~   covariate_biomass_copy_FI +  Intercept_FI
  form_FI_bycatch = biomass_bycatch_FI ~   covariate_biomass_copy_FI +  Intercept_FI
  form_FI_bycatch2 = biomass_bycatch2_FI ~   covariate_biomass_copy_FI +  Intercept_FI
  
  form_FD_target = biomass_target_FD ~   covariate_biomass  +  Intercept_FD
  form_FD_bycatch = biomass_bycatch_FD ~   covariate_biomass  +  Intercept_FD
  form_FD_bycatch2 = biomass_bycatch2_FD ~   covariate_biomass  +  Intercept_FD
  
  ###############################
  ###### likelihoods ##########
  ###########################
  
  ### sampling intensity ##
  lik_pp_iPM_random <- bru_obs("poisson",formula = form_pp_iPM, data = df_random_pp )
  lik_pp_iPM_proportional <- bru_obs("poisson",formula = form_pp_iPM, data = df_proportional_pp  )
  lik_pp_iPM_squared <- bru_obs("poisson",  formula = form_pp_iPM, data = df_squared_pp  )
  lik_pp_iPM_4th_power <- bru_obs("poisson",  formula = form_pp_iPM, data = df_4th_power_pp  )
  
  lik_pp_iEPM_random <- bru_obs("poisson",formula = form_pp_iEPM, data = df_random_pp )
  lik_pp_iEPM_proportional <- bru_obs("poisson",formula = form_pp_iEPM, data = df_proportional_pp  )
  lik_pp_iEPM_squared <- bru_obs("poisson",  formula = form_pp_iEPM, data = df_squared_pp  )
  lik_pp_iEPM_4th_power <- bru_obs("poisson",  formula = form_pp_iEPM, data = df_4th_power_pp  )
  
  
  ## FI ##
  lik_biomass_target_FI <- bru_obs("gaussian",formula = form_FI_target,data = FI_data )
  lik_biomass_bycatch_FI <- bru_obs("gaussian",formula = form_FI_bycatch,data = FI_data )
  lik_biomass_bycatch2_FI <- bru_obs("gaussian",formula = form_FI_bycatch2,data = FI_data )
  
  ## FD ##
  ## target
  lik_biomass_FD_random <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_random )
  lik_biomass_FD_proportional <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_proportional)
  lik_biomass_FD_squared <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_squared  )
  lik_biomass_FD_fourth <- bru_obs("gaussian",formula = form_FD_target,data = FD_data_fourth  )
  
  ## bycatch
  lik_biomass_FD_random_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_random)
  lik_biomass_FD_proportional_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch, data = FD_data_proportional  )
  lik_biomass_FD_squared_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_squared  )
  lik_biomass_FD_fourth_bycatch <- bru_obs("gaussian",formula = form_FD_bycatch,data = FD_data_fourth  )
  

  #'bycatch2
  lik_biomass_FD_random_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2,data = FD_data_random)
  lik_biomass_FD_proportional_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2, data = FD_data_proportional  )
  lik_biomass_FD_squared_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2,data = FD_data_squared  )
  lik_biomass_FD_fourth_bycatch2 <- bru_obs("gaussian",formula = form_FD_bycatch2,data = FD_data_fourth  )
  

  #########################
  #### fitting ###########
  #######################
  
  #### target species
  iPM_target_random_fit <- bru(cmp_iPM, lik_biomass_target_FI, lik_biomass_FD_random,lik_pp_iPM_random,
                               options = list(control.inla = list(int.strategy = "eb")))
  iPM_target_proportional_fit <- bru(cmp_iPM, lik_biomass_target_FI, lik_biomass_FD_proportional,lik_pp_iPM_proportional,
                                     options = list(control.inla = list(int.strategy = "eb")))
  iPM_target_squared_fit <- bru(cmp_iPM, lik_biomass_target_FI, lik_biomass_FD_squared,lik_pp_iPM_squared,
                                options = list(control.inla = list(int.strategy = "eb")))
  iPM_target_4th_power_fit <- bru(cmp_iPM, lik_biomass_target_FI, lik_biomass_FD_fourth,lik_pp_iPM_4th_power,
                                  options = list(control.inla = list(int.strategy = "eb")))
  
  iEPM_target_random_fit <- bru(cmp_iEPM, lik_biomass_target_FI, lik_biomass_FD_random,lik_pp_iEPM_random,
                                options = list(control.inla = list(int.strategy = "eb")))
  iEPM_target_proportional_fit <- bru(cmp_iEPM, lik_biomass_target_FI, lik_biomass_FD_proportional,lik_pp_iEPM_proportional,
                                      options = list(control.inla = list(int.strategy = "eb")))
  iEPM_target_squared_fit <- bru(cmp_iEPM, lik_biomass_target_FI, lik_biomass_FD_squared,lik_pp_iEPM_squared,
                                 options = list(control.inla = list(int.strategy = "eb")))
  iEPM_target_4th_power_fit <- bru(cmp_iEPM, lik_biomass_target_FI, lik_biomass_FD_fourth,lik_pp_iEPM_4th_power,
                                   options = list(control.inla = list(int.strategy = "eb")))
  
  #### bycatch species
  iPM_bycatch_random_fit <- bru(cmp_iPM, lik_biomass_bycatch_FI, lik_biomass_FD_random_bycatch,lik_pp_iPM_random,
                                options = list(control.inla = list(int.strategy = "eb")))
  iPM_bycatch_proportional_fit <- bru(cmp_iPM, lik_biomass_bycatch_FI, lik_biomass_FD_proportional_bycatch,lik_pp_iPM_proportional,
                                      options = list(control.inla = list(int.strategy = "eb")))
  iPM_bycatch_squared_fit <- bru(cmp_iPM, lik_biomass_bycatch_FI, lik_biomass_FD_squared,lik_pp_iPM_squared,
                                 options = list(control.inla = list(int.strategy = "eb")))
  iPM_bycatch_4th_power_fit <- bru(cmp_iPM, lik_biomass_bycatch_FI, lik_biomass_FD_fourth,lik_pp_iPM_4th_power,
                                   options = list(control.inla = list(int.strategy = "eb")))
  
  iEPM_bycatch_random_fit <- bru(cmp_iEPM, lik_biomass_bycatch_FI, lik_biomass_FD_random_bycatch,lik_pp_iEPM_random,
                                 options = list(control.inla = list(int.strategy = "eb")))
  iEPM_bycatch_proportional_fit <- bru(cmp_iEPM, lik_biomass_bycatch_FI, lik_biomass_FD_proportional_bycatch,lik_pp_iEPM_proportional,
                                       options = list(control.inla = list(int.strategy = "eb")))
  iEPM_bycatch_squared_fit <- bru(cmp_iEPM, lik_biomass_bycatch_FI, lik_biomass_FD_squared_bycatch,lik_pp_iEPM_squared,
                                  options = list(control.inla = list(int.strategy = "eb")))
  iEPM_bycatch_4th_power_fit <- bru(cmp_iEPM, lik_biomass_bycatch_FI, lik_biomass_FD_fourth_bycatch,lik_pp_iEPM_4th_power,
                                    options = list(control.inla = list(int.strategy = "eb")))
  
  #### bycatch2 species
  iPM_bycatch2_random_fit <- bru(cmp_iPM, lik_biomass_bycatch2_FI, lik_biomass_FD_random_bycatch2,lik_pp_iPM_random,
                                 options = list(control.inla = list(int.strategy = "eb")))
  iPM_bycatch2_proportional_fit <- bru(cmp_iPM, lik_biomass_bycatch2_FI, lik_biomass_FD_proportional_bycatch2,lik_pp_iPM_proportional,
                                       options = list(control.inla = list(int.strategy = "eb")))
  iPM_bycatch2_squared_fit <- bru(cmp_iPM, lik_biomass_bycatch2_FI, lik_biomass_FD_squared_bycatch2,lik_pp_iPM_squared,
                                  options = list(control.inla = list(int.strategy = "eb")))
  iPM_bycatch2_4th_power_fit <- bru(cmp_iPM, lik_biomass_bycatch2_FI, lik_biomass_FD_fourth_bycatch2,lik_pp_iPM_4th_power,
                                    options = list(control.inla = list(int.strategy = "eb")))
  
  iEPM_bycatch2_random_fit <- bru(cmp_iEPM, lik_biomass_bycatch2_FI, lik_biomass_FD_random_bycatch2,lik_pp_iEPM_random,
                                  options = list(control.inla = list(int.strategy = "eb")))
  iEPM_bycatch2_proportional_fit <- bru(cmp_iEPM, lik_biomass_bycatch2_FI, lik_biomass_FD_proportional_bycatch2,lik_pp_iEPM_proportional,
                                        options = list(control.inla = list(int.strategy = "eb")))
  iEPM_bycatch2_squared_fit <- bru(cmp_iEPM, lik_biomass_bycatch2_FI, lik_biomass_FD_squared_bycatch2,lik_pp_iEPM_squared,
                                   options = list(control.inla = list(int.strategy = "eb")))
  iEPM_bycatch2_4th_power_fit <- bru(cmp_iEPM, lik_biomass_bycatch2_FI, lik_biomass_FD_fourth_bycatch2,lik_pp_iEPM_4th_power,
                                     options = list(control.inla = list(int.strategy = "eb")))
  
  
  #### here
  
  
  a_iPM = mean(predict(iPM_target_random_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean)
  b_iPM = mean(predict(iPM_target_proportional_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean)
  c_iPM = mean(predict(iPM_target_squared_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean)
  d_iPM = mean(predict(iPM_target_4th_power_fit,cov_effect,
                       ~(Intercept_FD+ covariate_biomass), 
                       include=c("Intercept_FD","covariate_biomass"))$mean)
  
  
  aa_iPM = mean(predict(iPM_bycatch_random_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean)
  bb_iPM = mean(predict(iPM_bycatch_proportional_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean)
  cc_iPM = mean(predict(iPM_bycatch_squared_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean)
  dd_iPM = mean(predict(iPM_bycatch_4th_power_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean)
  
  
  aaa_iPM = mean(predict(iPM_bycatch2_random_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass), 
                          include=c("Intercept_FD","covariate_biomass"))$mean)
  bbb_iPM = mean(predict(iPM_bycatch2_proportional_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass), 
                          include=c("Intercept_FD","covariate_biomass"))$mean)
  ccc_iPM = mean(predict(iPM_bycatch2_squared_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass), 
                          include=c("Intercept_FD","covariate_biomass"))$mean)
  ddd_iPM = mean(predict(iPM_bycatch2_4th_power_fit,cov_effect,
                          ~(Intercept_FD+ covariate_biomass), 
                          include=c("Intercept_FD","covariate_biomass"))$mean)
  
  
  a_iEPM = tryCatch(mean(predict(iEPM_target_random_fit,cov_effect,
                                 ~(Intercept_FD+ covariate_biomass), 
                                 include=c("Intercept_FD","covariate_biomass"))$mean),
                    error=function(e)NA)
  b_iEPM = tryCatch(mean(predict(iEPM_target_proportional_fit,cov_effect,
                                 ~(Intercept_FD+ covariate_biomass), 
                                 include=c("Intercept_FD","covariate_biomass"))$mean),
                    error=function(e)NA)
  c_iEPM = tryCatch(mean(predict(iEPM_target_squared_fit,cov_effect,
                                 ~(Intercept_FD+ covariate_biomass), 
                                 include=c("Intercept_FD","covariate_biomass"))$mean),
                    error=function(e)NA)
  d_iEPM = tryCatch(mean(predict(iEPM_target_4th_power_fit,cov_effect,
                                 ~(Intercept_FD+ covariate_biomass), 
                                 include=c("Intercept_FD","covariate_biomass"))$mean),
                    error=function(e)NA)
  
  
  aa_iEPM = tryCatch(mean(predict(iEPM_bycatch_random_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean),
                     error=function(e)NA)
  bb_iEPM = tryCatch(mean(predict(iEPM_bycatch_proportional_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean),
                     error=function(e)NA)
  cc_iEPM = tryCatch(mean(predict(iEPM_bycatch_squared_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean),
                     error=function(e)NA)
  dd_iEPM = tryCatch(mean(predict(iEPM_bycatch_4th_power_fit,cov_effect,
                        ~(Intercept_FD+ covariate_biomass), 
                        include=c("Intercept_FD","covariate_biomass"))$mean),
                     error=function(e)NA)
  
  
  aaa_iEPM = tryCatch(mean(predict(iEPM_bycatch2_random_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean),
                      error=function(e)NA)
  bbb_iEPM = tryCatch(mean(predict(iEPM_bycatch2_proportional_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean),
                      error=function(e)NA)
  ccc_iEPM = tryCatch(mean(predict(iEPM_bycatch2_squared_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean),
                      error=function(e)NA)
  ddd_iEPM = tryCatch(mean(predict(iEPM_bycatch2_4th_power_fit,cov_effect,
                         ~(Intercept_FD+ covariate_biomass), 
                         include=c("Intercept_FD","covariate_biomass"))$mean),
                      error=function(e)NA)
  
  

  
  if(i==1){
    iEPM_error_random = mean_biomass_target_FD - a_iEPM
    iEPM_error_proportional= mean_biomass_target_FD-b_iEPM
    iEPM_error_squared = mean_biomass_target_FD - c_iEPM
    iEPM_error_4th_power = mean_biomass_target_FD - d_iEPM
    
    iEPM_byc_error_random = mean_biomass_byc_FD - bb_iEPM
    iEPM_byc_error_proportional= mean_biomass_byc_FD - bb_iEPM
    iEPM_byc_error_squared = mean_biomass_byc_FD - cc_iEPM
    iEPM_byc_error_4th_power = mean_biomass_byc_FD - dd_iEPM
    
    iEPM_byc2_error_random = mean_biomass_byc2_FD - aaa_iEPM
    iEPM_byc2_error_proportional= mean_biomass_byc2_FD - bbb_iEPM
    iEPM_byc2_error_squared = mean_biomass_byc2_FD - ccc_iEPM
    iEPM_byc2_error_4th_power = mean_biomass_byc2_FD - ddd_iEPM
    
    ##### pref simple
    iPM_error_random = mean_biomass_target_FD - a_iPM
    iPM_error_proportional= mean_biomass_target_FD-b_iPM
    iPM_error_squared = mean_biomass_target_FD - c_iPM
    iPM_error_4th_power = mean_biomass_target_FD - d_iPM
    
    iPM_byc_error_random = mean_biomass_byc_FD - aa_iPM
    iPM_byc_error_proportional= mean_biomass_byc_FD - bb_iPM
    iPM_byc_error_squared = mean_biomass_byc_FD - cc_iPM
    iPM_byc_error_4th_power = mean_biomass_byc_FD - dd_iPM
    
    iPM_byc2_error_random = mean_biomass_byc2_FD - aaa_iPM
    iPM_byc2_error_proportional= mean_biomass_byc2_FD - bbb_iPM
    iPM_byc2_error_squared = mean_biomass_byc2_FD - ccc_iPM
    iPM_byc2_error_4th_power = mean_biomass_byc2_FD - ddd_iPM
    
    
    ##### quantile at 0
    # q0_beta_random = q0_beta_a
    # q0_beta_proportional = q0_beta_b
    # q0_beta_squared = q0_beta_c
    # q0_beta_4th_power = q0_beta_d
    # 
    # q0_beta_byc_random = q0_beta_aa
    # q0_beta_byc_proportional = q0_beta_bb
    # q0_beta_byc_squared = q0_beta_cc
    # q0_beta_byc_4th_power = q0_beta_dd
    # 
    # q0_beta_byc2_random = q0_beta_aaa
    # q0_beta_byc2_proportional = q0_beta_bbb
    # q0_beta_byc2_squared = q0_beta_ccc
    # q0_beta_byc2_4th_power = q0_beta_ddd
    
    
    ##### alphas iPM
    beta_target_random_iPM = rbind(iPM_target_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_proportional_iPM = rbind(iPM_target_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_squared_iPM = rbind(iPM_target_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_4th_power_iPM = rbind(iPM_target_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    beta_byc_random_iPM = rbind(iPM_bycatch_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_proportional_iPM = rbind(iPM_bycatch_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_squared_iPM = rbind(iPM_bycatch_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_4th_power_iPM = rbind(iPM_bycatch_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    beta_byc2_random_iPM = rbind(iPM_bycatch2_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_proportional_iPM = rbind(iPM_bycatch2_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_squared_iPM = rbind(iPM_bycatch2_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_4th_power_iPM = rbind(iPM_bycatch2_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    ##### alphas iEPM
    beta_target_random_iEPM = rbind(iEPM_target_random_fit$summary.hyperpar[8,c(1:3,5)])
    beta_target_proportional_iEPM = rbind(iEPM_target_proportional_fit$summary.hyperpar[8,c(1:3,5)])
    beta_target_squared_iEPM = rbind(iEPM_target_squared_fit$summary.hyperpar[8,c(1:3,5)])
    beta_target_4th_power_iEPM = rbind(iEPM_target_4th_power_fit$summary.hyperpar[8,c(1:3,5)])
    
    beta_byc_random_iEPM = rbind(iEPM_bycatch_random_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc_proportional_iEPM = rbind(iEPM_bycatch_proportional_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc_squared_iEPM = rbind(iEPM_bycatch_squared_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc_4th_power_iEPM = rbind(iEPM_bycatch_4th_power_fit$summary.hyperpar[8,c(1:3,5)])
    
    beta_byc2_random_iEPM = rbind(iEPM_bycatch2_random_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc2_proportional_iEPM = rbind(iEPM_bycatch2_proportional_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc2_squared_iEPM = rbind(iEPM_bycatch2_squared_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc2_4th_power_iEPM = rbind(iEPM_bycatch2_4th_power_fit$summary.hyperpar[8,c(1:3,5)])
    
    
  }else{
    iEPM_error_random=c(iEPM_error_random,mean_biomass_target_FD - a_iEPM)
    iEPM_error_proportional=c(iEPM_error_proportional,mean_biomass_target_FD - b_iEPM)
    iEPM_error_squared =c(iEPM_error_squared,mean_biomass_target_FD - c_iEPM)
    iEPM_error_4th_power = c(iEPM_error_4th_power,mean_biomass_target_FD - d_iEPM)
    
    iEPM_byc_error_random=c(iEPM_byc_error_random,mean_biomass_byc_FD -aa_iEPM)
    iEPM_byc_error_proportional=c(iEPM_byc_error_proportional,mean_biomass_byc_FD -bb_iEPM)
    iEPM_byc_error_squared =c(iEPM_byc_error_squared,mean_biomass_byc_FD -cc_iEPM)
    iEPM_byc_error_4th_power = c(iEPM_byc_error_4th_power,mean_biomass_byc_FD -dd_iEPM)
    
    iEPM_byc2_error_random=c(iEPM_byc2_error_random,mean_biomass_byc2_FD -aaa_iEPM)
    iEPM_byc2_error_proportional=c(iEPM_byc2_error_proportional,mean_biomass_byc2_FD -bbb_iEPM)
    iEPM_byc2_error_squared =c(iEPM_byc2_error_squared,mean_biomass_byc2_FD -ccc_iEPM)
    iEPM_byc2_error_4th_power = c(iEPM_byc2_error_4th_power,mean_biomass_byc2_FD -ddd_iEPM)
    
    #### simple preferential
    iPM_error_random=c(iPM_error_random,mean_biomass_target_FD - a_iPM)
    iPM_error_proportional=c(iPM_error_proportional,mean_biomass_target_FD - b_iPM)
    iPM_error_squared =c(iPM_error_squared,mean_biomass_target_FD - c_iPM)
    iPM_error_4th_power = c(iPM_error_4th_power,mean_biomass_target_FD - d_iPM)
    
    iPM_byc_error_random=c(iPM_byc_error_random,mean_biomass_byc_FD -aa_iPM)
    iPM_byc_error_proportional=c(iPM_byc_error_proportional,mean_biomass_byc_FD -bb_iPM)
    iPM_byc_error_squared =c(iPM_byc_error_squared,mean_biomass_byc_FD -cc_iPM)
    iPM_byc_error_4th_power = c(iPM_byc_error_4th_power,mean_biomass_byc_FD -dd_iPM)
    
    iPM_byc2_error_random=c(iPM_byc2_error_random,mean_biomass_byc2_FD -aaa_iPM)
    iPM_byc2_error_proportional=c(iPM_byc2_error_proportional,mean_biomass_byc2_FD -bbb_iPM)
    iPM_byc2_error_squared =c(iPM_byc2_error_squared,mean_biomass_byc2_FD -ccc_iPM)
    iPM_byc2_error_4th_power = c(iPM_byc2_error_4th_power,mean_biomass_byc2_FD -ddd_iPM)
    
    
    error[["iPM_target_random_fit"]] = c( error[["iPM_target_random_fit"]] , any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_target_random_fit$logfile) ))
    error[["iPM_target_proportional_fit"]] = c( error[["cPM_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_target_proportional_fit$logfile) ))
    error[["iPM_target_squared_fit"]] = c(error[["iPM_target_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_target_squared_fit$logfile) ))
    error[["iPM_target_4th_power_fit"]] = c(error[["iPM_target_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_target_4th_power_fit$logfile) ))
    
    error[["iPM_bycatch_random_fit"]] = c(error[["iPM_bycatch_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch_random_fit$logfile) ))
    error[["iPM_bycatch_proportional_fit"]] = c(error[["iPM_bycatch_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch_proportional_fit$logfile) ))
    error[["iPM_bycatch_squared_fit"]] = c(error[["iPM_bycatch_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch_squared_fit$logfile) ))
    error[["iPM_bycatch_4th_power_fit"]] = c(error[["iPM_bycatch_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch_4th_power_fit$logfile) ))
    
    error[["iPM_bycatch2_random_fit"]] = c(error[["iPM_bycatch2_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch2_random_fit$logfile) ))
    error[["iPM_bycatch2_proportional_fit"]] = c(error[["iPM_bycatch2_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch2_proportional_fit$logfile) ))
    error[["iPM_bycatch2_squared_fit"]] = c(error[["iPM_bycatch2_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch2_squared_fit$logfile) ))
    error[["iPM_bycatch2_4th_power_fit"]] = c(error[["iPM_bycatch2_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iPM_bycatch2_4th_power_fit$logfile) ))
    
    error[["iEPM_target_random_fit"]] = c(error[["iEPM_target_random_fit"]] , any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_target_random_fit$logfile) ))
    error[["iEPM_target_proportional_fit"]] = c(error[["iEPM_target_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_target_proportional_fit$logfile) ))
    error[["iEPM_target_squared_fit"]] = c(error[["iEPM_target_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_target_squared_fit$logfile) ))
    error[["iEPM_target_4th_power_fit"]] = c(error[["iEPM_target_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_target_4th_power_fit$logfile) ))
    
    error[["iEPM_bycatch_random_fit"]] = c(error[["iEPM_bycatch_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch_random_fit$logfile) ))
    error[["iEPM_bycatch_proportional_fit"]] = c(error[["iEPM_bycatch_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch_proportional_fit$logfile) ))
    error[["iEPM_bycatch_squared_fit"]] = c(error[["iEPM_bycatch_squared_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch_squared_fit$logfile) ))
    error[["iEPM_bycatch_4th_power_fit"]] = c(error[["iEPM_bycatch_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch_4th_power_fit$logfile) ))
    
    error[["iEPM_bycatch2_random_fit"]] = c(error[["iEPM_bycatch2_random_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch2_random_fit$logfile) ))
    error[["iEPM_bycatch2_proportional_fit"]] = c(error[["iEPM_bycatch2_proportional_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch2_proportional_fit$logfile) ))
    error[["iEPM_bycatch2_squared_fit"]] = c(error[["iEPM_bycatch2_squared_fit"]]  ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch2_squared_fit$logfile) ))
    error[["iEPM_bycatch2_4th_power_fit"]] = c(error[["iEPM_bycatch2_4th_power_fit"]] ,any(grepl(pattern = "iterative process seems to diverge, 'vb.correction' is aborted", x = iEPM_bycatch2_4th_power_fit$logfile) ))
    
    
    
    ##### quantile at 0
    # q0_beta_random = c(q0_beta_random,q0_beta_a)
    # q0_beta_proportional = c(q0_beta_proportional,q0_beta_b)
    # q0_beta_squared = c(q0_beta_squared,q0_beta_c)
    # q0_beta_4th_power = c(q0_beta_4th_power,q0_beta_d)
    # 
    # q0_beta_byc_random = c(q0_beta_byc_random,q0_beta_aa)
    # q0_beta_byc_proportional = c(q0_beta_byc_proportional,q0_beta_bb)
    # q0_beta_byc_squared = c(q0_beta_byc_squared,q0_beta_cc)
    # q0_beta_byc_4th_power = c(q0_beta_byc_4th_power,q0_beta_dd)
    # 
    # q0_beta_byc2_random = c(q0_beta_byc2_random,q0_beta_aaa)
    # q0_beta_byc2_proportional = c(q0_beta_byc2_proportional,q0_beta_bbb)
    # q0_beta_byc2_squared = c(q0_beta_byc2_squared,q0_beta_ccc)
    # q0_beta_byc2_4th_power = c(q0_beta_byc2_4th_power,q0_beta_ddd)
    
    ##### alphas iPM
    beta_target_random_iPM = rbind(beta_target_random_iPM,iPM_target_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_proportional_iPM = rbind(beta_target_proportional_iPM,iPM_target_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_squared_iPM = rbind(beta_target_squared_iPM,iPM_target_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_target_4th_power_iPM = rbind(beta_target_4th_power_iPM,iPM_target_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    beta_byc_random_iPM = rbind(beta_byc_random_iPM,iPM_bycatch_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_proportional_iPM = rbind(beta_byc_proportional_iPM,iPM_bycatch_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_squared_iPM = rbind(beta_byc_squared_iPM,iPM_bycatch_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc_4th_power_iPM = rbind(beta_byc_4th_power_iPM,iPM_bycatch_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    beta_byc2_random_iPM = rbind(beta_byc2_random_iPM,iPM_bycatch2_random_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_proportional_iPM = rbind(beta_byc2_proportional_iPM,iPM_bycatch2_proportional_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_squared_iPM = rbind(beta_byc2_squared_iPM,iPM_bycatch2_squared_fit$summary.hyperpar[6,c(1:3,5)])
    beta_byc2_4th_power_iPM = rbind(beta_byc2_4th_power_iPM,iPM_bycatch2_4th_power_fit$summary.hyperpar[6,c(1:3,5)])
    
    ##### alphas iEPM
    beta_target_random_iEPM = rbind(beta_target_random_iEPM,iEPM_target_random_fit$summary.hyperpar[8,c(1:3,5)])
    beta_target_proportional_iEPM = rbind(beta_target_proportional_iEPM,iEPM_target_proportional_fit$summary.hyperpar[8,c(1:3,5)])
    beta_target_squared_iEPM = rbind(beta_target_squared_iEPM,iEPM_target_squared_fit$summary.hyperpar[8,c(1:3,5)])
    beta_target_4th_power_iEPM = rbind(beta_target_4th_power_iEPM,iEPM_target_4th_power_fit$summary.hyperpar[8,c(1:3,5)])
    
    beta_byc_random_iEPM = rbind(beta_byc_random_iEPM,iEPM_bycatch_random_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc_proportional_iEPM = rbind(beta_byc_proportional_iEPM,iEPM_bycatch_proportional_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc_squared_iEPM = rbind(beta_byc_squared_iEPM,iEPM_bycatch_squared_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc_4th_power_iEPM = rbind(beta_byc_4th_power_iEPM,iEPM_bycatch_4th_power_fit$summary.hyperpar[8,c(1:3,5)])
    
    beta_byc2_random_iEPM = rbind(beta_byc2_random_iEPM,iEPM_bycatch2_random_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc2_proportional_iEPM = rbind(beta_byc2_proportional_iEPM,iEPM_bycatch2_proportional_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc2_squared_iEPM = rbind(beta_byc2_squared_iEPM,iEPM_bycatch2_squared_fit$summary.hyperpar[8,c(1:3,5)])
    beta_byc2_4th_power_iEPM = rbind(beta_byc2_4th_power_iEPM,iEPM_bycatch2_4th_power_fit$summary.hyperpar[8,c(1:3,5)])
    
  }
  
  
  ### prpe
  
  
  
  
  print(paste(i,i,i,i,Sys.time()-ptm))
  #save.image("C:/Users/ip30/OneDrive - University of St Andrews/Desktop/Preferential_OK/results/Simulation_results_Gamma_final.RData")
}

save.image("D:/Simulation_results_MS3.3.RData")

load("D:/Simulation_results_MS3.3.RData")

##############################################
##### Count convergence problems ##############
#############################################
df <- stack(error)

df %>% group_by(ind) %>% summarise(perc_conver_error= sum(values)/n_sims) %>% 
  filter(perc_conver_error>0)




############################################
######### Figure 2 #########################
##############################################


##############
### TS
sim_conventional_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                      Estimates = c(Estimate_random,Estimate_proportional,Estimate_squared,Estimate_4th_power))
sim_conventional_results$Error = sim_conventional_results$Estimates - mean_biomass_target_FD
sim_conventional_results$Sampling_scheme = factor(sim_conventional_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_conventional_results$Sp ="TS"
sim_conventional_results$Model ="SDM"

sim_conventional_results=sim_conventional_results[,-2]

##############
### CBS species
sim_conventional_byc_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                          Estimates = c(Estimate_byc_random,Estimate_byc_proportional,Estimate_byc_squared,Estimate_byc_4th_power))
sim_conventional_byc_results$Error = sim_conventional_byc_results$Estimates - mean_biomass_byc_FD
sim_conventional_byc_results$Sampling_scheme = factor(sim_conventional_byc_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_conventional_byc_results$Sp ="CBS"
sim_conventional_byc_results$Model ="SDM"

sim_conventional_byc_results=sim_conventional_byc_results[,-2]

##############
### UBS species
sim_conventional_byc2_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                           Estimates = c(Estimate_byc2_random,Estimate_byc2_proportional,Estimate_byc2_squared,Estimate_byc2_4th_power))
sim_conventional_byc2_results$Error = sim_conventional_byc2_results$Estimates - mean_biomass_byc2_FD
sim_conventional_byc2_results$Sampling_scheme = factor(sim_conventional_byc2_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_conventional_byc2_results$Sp ="UBS"
sim_conventional_byc2_results$Model ="SDM"

sim_conventional_byc2_results=sim_conventional_byc2_results[,-2]




######################################
### PREFERENTIAL RESUKLTS

##############
### TS
sim_ePM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                      Error = c(ePM_error_random,ePM_error_proportional,ePM_error_squared,ePM_error_4th_power))
sim_ePM_results$Sampling_scheme = factor(sim_ePM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_ePM_results$Sp ="TS"
sim_ePM_results$Model ="ePM"

##############
### CBS species
sim_byc_ePM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                          Error = c(ePM_byc_error_random,ePM_byc_error_proportional,ePM_byc_error_squared,ePM_byc_error_4th_power))
sim_byc_ePM_results$Sampling_scheme = factor(sim_byc_ePM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_byc_ePM_results$Sp ="CBS"
sim_byc_ePM_results$Model ="ePM"

##############
### UBS #####
sim_byc2_ePM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                           Error = c(ePM_byc2_error_random,ePM_byc2_error_proportional,ePM_byc2_error_squared,ePM_byc2_error_4th_power))
sim_byc2_ePM_results$Sampling_scheme = factor(sim_byc2_ePM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_byc2_ePM_results$Sp ="UBS"
sim_byc2_ePM_results$Model ="ePM"



##################################
###### traditional preferential 
sim_cPM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                             Error = c(cPM_error_random,cPM_error_proportional,cPM_error_squared,cPM_error_4th_power))
sim_cPM_results$Sampling_scheme = factor(sim_cPM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_cPM_results$Sp ="TS"
sim_cPM_results$Model ="cPM"


sim_byc_cPM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                                 Error = c(cPM_byc_error_random,cPM_byc_error_proportional,cPM_byc_error_squared,cPM_byc_error_4th_power))
sim_byc_cPM_results$Sampling_scheme = factor(sim_byc_cPM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_byc_cPM_results$Sp ="CBS"
sim_byc_cPM_results$Model ="cPM"

sim_byc2_cPM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                                  Error = c(cPM_byc2_error_random,cPM_byc2_error_proportional,cPM_byc2_error_squared,cPM_byc2_error_4th_power))
sim_byc2_cPM_results$Sampling_scheme = factor(sim_byc2_cPM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
sim_byc2_cPM_results$Sp ="UBS"
sim_byc2_cPM_results$Model ="cPM"







###############################################
####### ISDM models ######################
############################################

##############
### target species
combine_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                             Error = c(ISDM_error_random,ISDM_error_proportional,ISDM_error_squared,ISDM_error_4th_power))
combine_results$Sampling_scheme = factor(combine_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
combine_results$Sp ="TS"
combine_results$Model ="iSDM"

##############
### bycatch species
combine_byc_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                 Error = c(ISDM_byc_error_random,ISDM_byc_error_proportional,ISDM_byc_error_squared,ISDM_byc_error_4th_power))
combine_byc_results$Sampling_scheme = factor(combine_byc_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
combine_byc_results$Sp ="CBS"
combine_byc_results$Model ="iSDM"

##############
### bycatch2 species
combine_byc2_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                  Error = c(ISDM_byc2_error_random,ISDM_byc2_error_proportional,ISDM_byc2_error_squared,ISDM_byc2_error_4th_power))
combine_byc2_results$Sampling_scheme = factor(combine_byc2_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
combine_byc2_results$Sp ="UBS"
combine_byc2_results$Model ="iSDM"






###############################################
####### iPM models ######################
############################################

##############
### target species
iPM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                             Error = c(iPM_error_random,iPM_error_proportional,iPM_error_squared,iPM_error_4th_power))
iPM_results$Sampling_scheme = factor(iPM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
iPM_results$Sp ="TS"
iPM_results$Model ="iPM"

##############
### bycatch species
iPM_byc_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                 Error = c(iPM_byc_error_random,iPM_byc_error_proportional,iPM_byc_error_squared,iPM_byc_error_4th_power))
iPM_byc_results$Sampling_scheme = factor(iPM_byc_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
iPM_byc_results$Sp ="CBS"
iPM_byc_results$Model ="iPM"

##############
### bycatch2 species
iPM_byc2_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                  Error = c(iPM_byc2_error_random,iPM_byc2_error_proportional,iPM_byc2_error_squared,iPM_byc2_error_4th_power))
iPM_byc2_results$Sampling_scheme = factor(iPM_byc2_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
iPM_byc2_results$Sp ="UBS"
iPM_byc2_results$Model ="iPM"


###############################################
####### iEPM models ######################
############################################

##############
### target species
iEPM_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                             Error = c(iEPM_error_random,iEPM_error_proportional,iEPM_error_squared,iEPM_error_4th_power))
iEPM_results$Sampling_scheme = factor(iEPM_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
iEPM_results$Sp ="TS"
iEPM_results$Model ="iEPM"

##############
### bycatch species
iEPM_byc_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                 Error = c(iEPM_byc_error_random,iEPM_byc_error_proportional,iEPM_byc_error_squared,iEPM_byc_error_4th_power))
iEPM_byc_results$Sampling_scheme = factor(iEPM_byc_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
iEPM_byc_results$Sp ="CBS"
iEPM_byc_results$Model ="iEPM"

##############
### bycatch2 species
iEPM_byc2_results = data.frame(Sampling_scheme =rep(c("Random","Proportional","Squared","Power 4"),each=n_sims),
                                  Error = c(iEPM_byc2_error_random,iEPM_byc2_error_proportional,iEPM_byc2_error_squared,iEPM_byc2_error_4th_power))
iEPM_byc2_results$Sampling_scheme = factor(iEPM_byc2_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
iEPM_byc2_results$Sp ="UBS"
iEPM_byc2_results$Model ="iEPM"




######################################################
############ NEW Figure combined #################
############################################

all_results = rbind(sim_conventional_results,sim_conventional_byc_results,sim_conventional_byc2_results,
                    sim_ePM_results,sim_byc_ePM_results,sim_byc2_ePM_results,
                    sim_cPM_results,sim_byc_cPM_results,sim_byc2_cPM_results,
                    combine_results,combine_byc_results,combine_byc2_results,
                    iPM_results,iPM_byc_results,iPM_byc2_results,
                    iEPM_results,iEPM_byc_results,iEPM_byc2_results)


all_results$Sampling_scheme = factor(all_results$Sampling_scheme,levels = c("Random","Proportional","Squared","Power 4"))
all_results$Sp = factor(all_results$Sp,levels = c("TS","CBS","UBS"))
all_results$Model = rep(c("SDM" , "ePM" , "cPM" , "iSDM" ,"PiSDM" , "ePiSDM"),each=600)
all_results$Model = factor(all_results$Model,levels = c("SDM","cPM","ePM","iSDM","PiSDM","ePiSDM"))


png("C:/Users/ip30/OneDrive - University of St Andrews/Desktop/results_pref_Ok.png",height=500,width = 900)
ggplot(all_results) + geom_boxplot(aes(x=Sampling_scheme,y=Error,color=Model)) + 
  geom_hline(yintercept = 0,color="red",linewidth=.5,alpha=.4) +
  facet_wrap(~Sp) + theme_bw() +
  xlab("Level of preferential sampling") +
  ylab("Bias") + ylim(-1,2.5)+
  theme(#axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=12),
    legend.text.align = 0,strip.text.x = element_text(size = 14))
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



png("C:/use/Pref_revision_paper/img/q0_beta.png",height=350,width = 1000)
grid.arrange(q0_beta_z,q0_beta_v,q0_beta_y,
             bottom=textGrob("Level of preferential sampling", gp=gpar(fontsize=17, fontface = "bold")),
             left=textGrob("Quantile at zero", rot = 90, gp=gpar(fontsize=17, fontface = "bold")),
             #bottom="Level of preferential sampling",
             #left = "Quantile at zero",
             ncol=3)
dev.off()









####################################################
###### scaling parameters of joint effects ######
##############################################
beta_objects <- ls()[grepl("^beta_", ls())]
beta_objects
beta_objects[grepl("iEPM", beta_objects)]



library(tidyverse)
library(patchwork)  # For combining plots
library(grid)  # For textGrob

# Define the structure
models <- c("cPM", "ePM", "iPM", "iEPM")
species <- c("target", "byc", "byc2")
pref_sampling <- c("random", "proportional", "squared", "4th_power")

# Define nice labels for models
model_labels <- c(
  "cPM" = "cPM",
  "ePM" = "ePM",
  "iPM" = "PiSDM",
  "iEPM" = "ePiSDM"
)

# Define nice labels for species (rows)
species_labels <- c(
  "target" = "Target species",
  "byc" = "Correlated species",
  "byc2" = "Uncorrelated species"
)

# Define nice labels for preference sampling (columns)
pref_labels <- c(
  "random" = "Pref level: random",
  "proportional" = "Pref level: proportional",
  "squared" = "Pref level: squared",
  "4th_power" = "Pref level: 4th_power"
)

# Function to create a single panel plot
create_panel <- function(object_name, data, is_first_row, is_first_col, 
                         species_label, pref_label) {
  p <- data %>%
    mutate(sim = row_number()) %>%
    ggplot(aes(x = sim, y = mean)) +
    geom_point(size = 1.5, color = "blue", alpha = 0.7) +
    geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`), 
                  width = 0.2, color = "blue", alpha = 0.7) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 7),
      axis.title.x = element_blank(),  # Remove x-axis title from panels
      axis.title.y = element_blank(),  # Remove y-axis title from panels
      plot.margin = margin(2, 2, 2, 2)
    )
  
  # Add column title (only for first row)
  if (is_first_row) {
    p <- p + ggtitle(pref_label) +
      theme(plot.title = element_text(size = 9, hjust = 0.5))
  }
  
  # Add row label (only for first column) as a text annotation
  if (is_first_col) {
    p <- p + 
      annotate("text", x = -Inf, y = Inf, label = species_label, 
               hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold")
  }
  
  return(p)
}

# Function to create figure for one model
create_model_figure <- function(model_name) {
  
  # Create all 12 panels (3 species × 4 pref_sampling)
  plot_list <- list()
  
  for (i in seq_along(species)) {
    for (j in seq_along(pref_sampling)) {
      # Construct object name
      obj_name <- paste0("beta_", species[i], "_", pref_sampling[j], "_", model_name)
      
      # Get the data object
      plot_data <- get(obj_name)
      
      # Create panel
      p <- create_panel(
        obj_name, 
        plot_data, 
        is_first_row = (i == 1),
        is_first_col = (j == 1),
        species_label = species_labels[species[i]],
        pref_label = pref_labels[pref_sampling[j]]
      )
      
      # Store plot
      plot_list[[(i - 1) * length(pref_sampling) + j]] <- p
    }
  }
  
  # Combine all panels into one figure
  combined_plot <- wrap_plots(plot_list, ncol = 4, nrow = 3) +
    plot_annotation(
      title = model_labels[model_name],
      caption = "Simulation",  # Overall x-axis label at bottom
      theme = theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.caption = element_text(size = 11, hjust = 0.5)
      )
    )
  
  # Wrap with overall y-axis label
  final_plot <- wrap_elements(combined_plot) +
    labs(tag = expression(paste("Scaling parameter (", alpha, ")"))) +
    theme(
      plot.tag = element_text(size = 11, angle = 90),
      plot.tag.position = "left"
    )
  
  return(final_plot)
}

# Generate figures for all 4 models
figures <- map(models, create_model_figure)
names(figures) <- models

# Display or save figures
# To view:
figures$cPM
figures$ePM
figures$iPM
figures$iEPM

# To save:
walk2(figures, models, ~ggsave(
  filename = paste0("beta_estimates_", .y, ".png"),
  plot = .x,
  width = 12,
  height = 9,
  dpi = 300
))
