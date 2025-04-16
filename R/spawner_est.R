library(tidyverse)
library(tidybayes)
library(rstan)
library(loo)
library(patchwork)
library(ggrepel)
library(janitor)
library(scales)

theme_set(theme_bw())

#load data####
spawners <- read.csv("../../../okanagan_data/2025-04-07 Draft/Spawn Timing Data/spawn_timing.csv") %>% 
  mutate(date = as.Date(date), year = year(date), yday = yday(date))

spawners_old <- readxl::read_excel("../../../data/Okanagan/!CNAT 2020/CNAT_Escapement Units_2020.xlsx", sheet = "Escapement Units") %>% 
  janitor::clean_names()
spawners_old <- select(spawners_old, year, spawners = au_cindex) %>% 
  filter(!is.na(year)) %>% 
  arrange(year) %>% 
  mutate(spawners =ifelse(year == 2011, NA, spawners))

index_sk <- spawners %>% 
  group_by(location, year, yday, date) %>% 
  filter(location == "Index") %>% 
  filter(species == "sk") %>% 
  summarise(
    live = if (all(is.na(live))) NA_real_ else sum(live, na.rm = TRUE),
    dead = if (all(is.na(dead))) NA_real_ else sum(dead, na.rm = TRUE)) %>% 
  arrange(date) %>% 
  ungroup() %>% 
  mutate(day_ind = as.numeric(factor(yday)), 
         year_ind = as.numeric(factor(year))) %>% 
  filter(!is.na(live)) 
#filter(!(location == "Above McIntyre Dam" & year == 2010)) %>% 
#filter(live>0)

#initial plotting of data####
index_sk %>% 
  ggplot(aes(x = yday, y = live, color = location))+
  geom_point()+
  facet_wrap(~year)+
  scale_y_log10()

index_sk %>% 
  ggplot(aes(x = yday, y = live, color = location))+
  geom_point()+
  facet_wrap(~year)


#set up data for model####
priors <- data.frame(prior = c("log_runs_mu", "timings_mu", "timings_sigma", "spread","spread_sigma","residence_mu", "count_dispersion"),
                     v1 = c(10,280-240, 0,log(8-1), 0,log(11-1), 1), 
                     v2 = c(1,5,1,0.2, 10, 0.3, 0.2))


sp_dat = list(n_priors = nrow(priors), 
              priors = data.matrix(priors[,-1]),
              n_years = max(index_sk$year_ind),
              year = as.numeric(factor(index_sk$year)),
              day = index_sk$yday-240,
              n_obs = nrow(index_sk),
              live_counts = index_sk$live)

#stan model####

m <- stan_model(file = "./stan/spawners.stan")
fit <- sampling(m, data = sp_dat, chains = 4, cores = 4)

#model checks####
worst_Rhat <- summary(fit)$summary %>% 
  as.data.frame() %>% 
  mutate(Rhat = round(Rhat, 3)) %>% 
  arrange(desc(Rhat))

worst_Rhat %>% 
  filter(n_eff>3) %>% 
  ggplot(aes(x = n_eff, y = Rhat))+
  geom_point()+
  geom_hline(yintercept = 1.01, lty = 2)+
  geom_vline(xintercept = 400, lty = 2)

traceplot(fit, pars = rownames(worst_Rhat)[1:20])
traceplot(fit, pars = c("timing_sigma", "timing_mu", "spread_arrive_mu", "spread_arrive_sigma", "spread_death_mu", "spread_death_sigma"))

#posterior plots####
post <- extract(fit)

gather_draws(fit, residence) %>% 
  ggplot(aes(x = .variable, y = .value))+
  geom_hline(yintercept = 11, lty = 2)+
  stat_pointinterval()


spread_draws(fit, log_run[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = exp(log_run)))+
  stat_pointinterval(aes(color = "modelled"))+
  geom_point(data = spawners_old %>% filter(year > 1999), aes(x = year, y = spawners, color = "AUC"))+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), name = "Spawner abundance (thousands)", expand = expansion(mult = c(0, 0.02)))+
  scale_color_manual(values = c("red", 1), "")+
  xlab("")
ggsave("./figures/spawner_abund.pdf", height = 4, width = 7)
  
spread_draws(fit, timing[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = timing))+
  stat_pointinterval()+
  
  spread_draws(fit, spread_arrive[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = spread_arrive))+
  stat_pointinterval()+
  
  spread_draws(fit, spread_death[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = spread_death))+
  stat_pointinterval()+
  plot_layout(ncol = 1)

#plot estimated spawner curves####
spawn_curves.df <- data.frame()
for(j in 1:ncol(post$log_run)){
  for(i in 1:sp_dat$n_locations){
    if(sp_dat$is_valid[j,i] == 1){
      live <- sapply(1:max(sp_dat$day), FUN = function(x) {
        entered <- pnorm(x, post$timing[,j,i], post$spread_arrive[,j,i])
        exited <- pnorm(x, post$timing[,j,i] + post$residence[,i], post$spread_death[,j,i])
        
        exp(post$log_run[,j] + log(entered - (entered * exited))) * post$run_prop[,j,i]
      })
      
      spawn_curves.df <- bind_rows(spawn_curves.df, data.frame(year = j + min(spawners$year)-1, 
                                                               day = 1:max(sp_dat$day)+240,
                                                               location = location_levels$location[i],
                                                               fish = apply(live, 2, median),
                                                               l89 = apply(live, 2, quantile, probs = 0.065),
                                                               u89 = apply(live, 2, quantile, probs = 0.945)))
    }
  }
}

ggplot(spawn_curves.df, aes(x = day, y = fish, color = location, fill = location))+
  geom_line()+
  geom_ribbon(aes(ymin = l89, ymax = u89), color = NA, alpha = 0.2)+
  facet_wrap(~year, scales = "free_y")+
  geom_point(data = index_sk, aes(x = yday, y = live), size = 0.8)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")

ggplot(spawn_curves.df, aes(x = day, y = fish))+
  geom_line()+
  geom_ribbon(aes(ymin = l89, ymax = u89), color = NA, alpha = 0.2)+
  facet_grid(location~year, scales = "free_y")+
  geom_point(data = index_sk, aes(x = yday, y = live), size = 0.8)+
  scale_x_continuous(breaks = seq(240, 325, by = 35))

ggplot(spawn_curves.df, aes(x = day, y = fish, color = location, fill = location))+
  geom_line()+
  geom_ribbon(aes(ymin = l89, ymax = u89), color = NA, alpha = 0.2)+
  facet_wrap(~year)+
  geom_point(data = index_sk, aes(x = yday, y = live), size = 0.8)+
  scale_y_log10()+
  coord_cartesian(ylim = c(1,1e5))+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")

