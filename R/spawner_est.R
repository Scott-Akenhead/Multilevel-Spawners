library(tidyverse)
library(tidybayes)
library(rstan)
library(loo)
library(patchwork)
library(ggrepel)
library(janitor)
library(scales)
library(geomtextpath)

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
  ggplot(aes(x = year, y = timing + 240))+
  stat_pointinterval()+
  scale_y_continuous(breaks = seq(280, 292, by = 3))+
  ylab("Mean arrival day")+
  xlab("")+
  
  spread_draws(fit, spread_arrive[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = spread_arrive))+
  stat_pointinterval()+
  ylab("Arrival timing sd (days)")+
  xlab("")+
  
  spread_draws(fit, spread_death[year]) %>% 
  mutate(year = year + min(spawners$year)-1) %>% 
  ggplot(aes(x = year, y = spread_death))+
  stat_pointinterval()+
  plot_layout(ncol = 1)+
  ylab("Death timing sd (days)")+
  xlab("")+
  
  plot_annotation(tag_levels = "A")
ggsave("./figures/timing_spread.pdf", height = 8, width = 6)

#plot estimated spawner curves####
spawn_curves.df <- data.frame()
for(j in 1:ncol(post$log_run)){
  live <- sapply(1:max(sp_dat$day), FUN = function(x) {
    entered <- pnorm(x, post$timing[,j], post$spread_arrive[,j])
    exited <- pnorm(x, post$timing[,j] + post$residence, post$spread_death[,j])
    
    exp(post$log_run[,j] + log(entered - (entered * exited)))
  })
  
  spawn_curves.df <- bind_rows(spawn_curves.df, data.frame(year = j + min(spawners$year)-1, 
                                                           day = 1:max(sp_dat$day)+240,
                                                           fish = apply(live, 2, median),
                                                           l89 = apply(live, 2, quantile, probs = 0.065),
                                                           u89 = apply(live, 2, quantile, probs = 0.945)))
}


ggplot(spawn_curves.df, aes(x = day, y = fish))+
  geom_line()+
  geom_ribbon(aes(ymin = l89, ymax = u89), color = NA, alpha = 0.2)+
  facet_wrap(~year)+
  geom_point(data = index_sk, aes(x = yday, y = live), size = 0.8)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), name = "Spawner count (thousands)", expand = expansion(mult = c(0, 0.02)))+
  theme(strip.background = element_rect(color="NA", fill="NA"))+
  scale_x_continuous(labels = ~ format(as.Date(.x, origin = "2023-12-31"), "%b"), 
                     breaks = as.numeric(as.Date(c("2024-09-01", "2024-10-01", "2024-11-01")) - as.Date("2023-12-31")), name = "")
ggsave("./figures/spawner_curves.pdf",width = 8, height = 8)

j <- 25 # used 2024 as example
arrived <- sapply(15:max(sp_dat$day), FUN = function(x) {
  pnorm(x, post$timing[,j], post$spread_arrive[,j]) * exp(post$log_run[,j])
}) %>% 
  as_tibble() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "day", values_to = "spawners") %>%
  mutate(day = as.integer(gsub("V", "", day)) + 14) %>% 
  mutate(type = "arrived")

dead <- sapply(15:max(sp_dat$day), FUN = function(x) {
  pnorm(x, post$timing[,j] + post$residence, post$spread_death[,j]) * exp(post$log_run[,j])
}) %>% 
  as_tibble() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "day", values_to = "spawners") %>%
  mutate(day = as.integer(gsub("V", "", day))+14) %>% 
  mutate(type = "dead")

counts <- sapply(15:max(sp_dat$day), FUN = function(x) {
  entered <- pnorm(x, post$timing[,j], post$spread_arrive[,j])
  exited <- pnorm(x, post$timing[,j] + post$residence, post$spread_death[,j])
  
  exp(post$log_run[,j] + log(entered - (entered * exited)))
}) %>% 
  as_tibble() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(-draw, names_to = "day", values_to = "spawners") %>%
  mutate(day = as.integer(gsub("V", "", day)) + 14) %>% 
  mutate(type = "present")

example_data <- bind_rows(arrived, dead, counts) %>% 
  group_by(type, day) %>% 
  summarise(l89 = quantile(spawners, probs = 0.065), 
            u89 = quantile(spawners, probs = 0.945),
            spawners = median(spawners))

mean_arrival_day <- median(post$timing[,j])
median_residence <- median(post$residence)
y_value <- median(pnorm(mean_arrival_day, post$timing[,j], post$spread_arrive[,j]) * exp(post$log_run[,j]))
x_start <- mean_arrival_day + 240
x_end <- mean_arrival_day + median_residence + 240

residence_df <- data.frame(
  x = x_start,
  xend = x_end,
  y = y_value,
  yend = y_value,
  label = "residence"
)

ggplot(example_data, aes(x = day + 240, y = spawners, group = type))+
  geom_textline(size = 4, linewidth = 1, aes(label = type), hjust = 0.65, family = "sans") +
  scale_x_continuous(labels = ~ format(as.Date(.x, origin = "2023-12-31"), "%b ") %>% 
                       paste0(as.integer(format(as.Date(.x, origin = "2023-12-31"), "%d"))),
                     breaks = as.numeric(as.Date(c("2024-09-15", "2024-10-01", "2024-10-15", "2024-11-01", "2024-11-15")) - as.Date("2023-12-31")), 
                     name = "", expand = c(0, 0))+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), name = "Spawners (thousands)", expand = expansion(mult = c(0, 0.075)))+
  geom_textsegment(data = residence_df, aes(x = x_start, xend = x_end, y = y_value, yend = y_value, label = "residence"), linewidth = 0.5, inherit.aes = FALSE, size = 2.8, family = "sans")+
  annotate("text",
           x = Inf,
           y = max(example_data$spawners),  
           label = "total abundance",
           hjust = 1.03,
           vjust = -0.3,
           size = 4,
           family = "sans")+
  theme(panel.grid = element_blank())
ggsave("./figures/method_example_no_CI.pdf", width = 6, height = 4)  

ggplot(example_data, aes(x = day + 240, y = spawners, group = type))+
  geom_ribbon(data = example_data %>% filter(type == "present"), aes(ymin = l89, ymax = u89), color = NA, alpha = 0.2)+
  geom_textline(size = 4, linewidth = 1, aes(label = type), hjust = 0.65, family = "sans") +
  scale_x_continuous(labels = ~ format(as.Date(.x, origin = "2023-12-31"), "%b ") %>% 
                       paste0(as.integer(format(as.Date(.x, origin = "2023-12-31"), "%d"))),
                     breaks = as.numeric(as.Date(c("2024-09-15", "2024-10-01", "2024-10-15", "2024-11-01", "2024-11-15")) - as.Date("2023-12-31")), 
                     name = "", expand = c(0, 0))+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), name = "Spawners (thousands)", expand = expansion(mult = c(0, 0.075)))+
  geom_point(data = index_sk %>% filter(year == 2024), aes(x = yday, y = live, group = NA))+
  geom_textsegment(data = residence_df, aes(x = x_start, xend = x_end, y = y_value, yend = y_value, label = "residence"), linewidth = 0.5, inherit.aes = FALSE, size = 2.8, family = "sans")+
  annotate("text",
           x = Inf,
           y = max(example_data$spawners),  
           label = "total abundance",
           hjust = 1.03,
           vjust = -0.3,
           size = 4,
           family = "sans")+
  theme(panel.grid = element_blank())
ggsave("./figures/method_example_obs_CI.pdf", width = 6, height = 4)  

ggplot(example_data, aes(x = day + 240, y = spawners, group = type))+
  geom_ribbon(aes(ymin = l89, ymax = u89), color = NA, alpha = 0.2)+
  geom_textline(size = 4, linewidth = 1, aes(label = type), hjust = 0.65, family = "sans") +
  scale_x_continuous(labels = ~ format(as.Date(.x, origin = "2023-12-31"), "%b ") %>% 
                       paste0(as.integer(format(as.Date(.x, origin = "2023-12-31"), "%d"))),
                     breaks = as.numeric(as.Date(c("2024-09-15", "2024-10-01", "2024-10-15", "2024-11-01", "2024-11-15")) - as.Date("2023-12-31")), 
                     name = "", expand = c(0, 0))+
  scale_y_continuous(labels = label_number(scale = 1e-3, suffix = "k"), name = "Spawners (thousands)", expand = expansion(mult = c(0, 0.05)))+
  geom_point(data = index_sk %>% filter(year == 2024), aes(x = yday, y = live, group = NA))+
  annotate("text",
           x = Inf,
           y = max(example_data$spawners),  
           label = "total abundance",
           hjust = 1.03,
           vjust = -0.3,
           size = 4,
           family = "sans")+
  theme(panel.grid = element_blank())
ggsave("./figures/method_example.pdf", width = 6, height = 4)
