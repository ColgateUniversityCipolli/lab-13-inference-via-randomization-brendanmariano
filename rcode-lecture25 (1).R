################################################################################
# LECTURE 25 R CODE
# CIPOLLI
# MATH 240 - SPRING 2024
################################################################################
library(tidyverse)
library(patchwork)
library(boot)
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0256568

################################################################################
# Load Dataset
################################################################################
dat <- read_csv("powerplays.csv")

dat <- dat |>
  ###################################################
# FOLLOWING METHODS SECTION OF THE PAPER:
# Remove Edmonton Oilers and Toronto Maple Leafs 
#        (They played at home unlike other teams)
# Remove 12 Round Robin Games       
#        (These are low-stakes games)
# https://www.espn.com/nhl/schedule/_/date/20200801
####################################################
slice(-c(7,9,     # Edmonton Oilers V CHI 8/1
         27,29,   # Edmonton Oilers V CHI 8/3
         51,54,   # Edmonton Oilers V CHI 8/5
         72,74,   # Edmonton Oilers V CHI 8/7
         13,18,   # Toronto Maple Leafs V CBJ 8/2
         35,41,   # Toronto Maple Leafs V CBJ 8/4
         58,63,   # Toronto Maple Leafs V CBJ 8/6
         70,76,   # Toronto Maple Leafs V CBJ 8/7
         84,88,   # Toronto Maple Leafs V CBJ 8/9
         11,12,   # PHI V BOS 8/2
         15,16,   # STL V COL 8/2
         25,26,   # WSH V TB  8/3  
         30,32,   # DAL V VGK 8/3 (not adjacent)
         47,50,   # TB V BOS 8/5 (not adjacent)
         53,55,   # COL V DAL 8/5 (not adjacent)
         57,59,   # WSH V PHI 8/6 (not adjacent)
         60,66,   # VGK V STL 8/6 (not adjacent)
         81,82,   # VGK V COL 8/8
         79,80,   # PHI V TB  8/8
         83,85,   # BOS V WSH 8/9 (not adjacent)
         86,87    # DAL V STL 8/9
)) |>
  filter(Home == 1)
# Interested in `PP Opp`

# Mean Home Penalties (2015-2019) = 2.90
# (244+283+250+252+243)/(89+91+87+84+87)
mu0 <- 2.90

################################################################################
# Summarize the data
################################################################################
ggplot(data=dat) +
  geom_bar(aes(x=`PP Opp`, y=after_stat(count/sum(count))))+
  geom_hline(yintercept = 0)+
  theme_bw()+
  scale_x_continuous(name = "Penalties Called on Home Team",
                     breaks = 0:8) +
  ylab("Proportion")

dat |>
  summarize(mean = mean(`PP Opp`),
            sd   = sd(`PP Opp`))

################################################################################
# T Test and Interval
################################################################################
(ttest.res <- t.test(x = dat$`PP Opp`,
                     mu = mu0,
                     alternative = "two.sided"))

library(effectsize)
hedges_g(x = dat$`PP Opp`,
         mu = mu0,
         alternative = "two.sided")

# Resampling for plot
R <- 10000
resamples <- tibble(xbars=rep(NA, R))
s <- sd(dat$`PP Opp`)
n <- nrow(dat)
for(i in 1:R){
  curr.resample <- sample(dat$`PP Opp`,
                          size = nrow(dat),
                          replace = T)
  
  resamples$xbars[i] <- (mean(curr.resample)-mu0)/(s/sqrt(n))
}


################################################################################
# Plot the test results
################################################################################
ggdat.t <- tibble(t = seq(-4, 4, length.out=1000)) |>
  mutate(pdf = dt(t, df=n-1))
ggdat.obs <- tibble(t=ttest.res$statistic,
                    y=0)

t.breaks <- as.numeric(c(-ttest.res$statistic, 
                         qt(0.025, df=n-1),
                         0,
                         qt(0.975, df=n-1),
                         ttest.res$statistic))

xbar.breaks <- t.breaks*sd(dat$`PP Opp`)/sqrt(n) + mu0

t.plot <- ggplot() +
  # Plot Sampling Distribution under H0
  geom_line(data=ggdat.t,
            aes(x=t, y=pdf))+
  # Highlight Rejection Regions
  geom_ribbon(data=subset(ggdat.t, t<=qt(0.025, df=n-1)),
              aes(x=t, ymin=0, ymax=pdf),
              fill="grey", alpha=0.5)+
  geom_ribbon(data=subset(ggdat.t, t>=qt(0.975, df=n-1)),
              aes(x=t, ymin=0, ymax=pdf),
              fill="grey", alpha=0.5)+
  # Plot Resampling Distribution
  stat_density(data=resamples, aes(x=xbars), geom="line", color="lightgrey")+
  # Plot Observation
  geom_point(data=ggdat.obs,
             aes(x=t, y=y,color="Observation"))+
  geom_point(data=ggdat.obs, shape=19,
             aes(x=-t, y=y,color="Mirror Observation"))+
  scale_color_manual("",
                     values=c("black", "red"))+
  # Clean up
  geom_hline(yintercept=0)+
  theme_bw()+
  scale_x_continuous("t",
                     limits = (c(2.0, 4.5)-mu0)/(s/sqrt(n)),
                     breaks = round(t.breaks,2),
                     sec.axis = sec_axis(~.,
                                         bquote(bar(x)),
                                         breaks = t.breaks,
                                         labels = round(xbar.breaks,2)))+
  ylab("Density")+
  ggtitle("T Hypothesis Test")

t.plot

################################################################################
# Bootstrap Test Confidence Interval
################################################################################
R <- 10000
resamples <- tibble(xbars = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(dat$`PP Opp`,
                          size = nrow(dat),
                          replace = T)
  
  resamples$xbars[i] <- mean(curr.resample)
}

# Confidence Interval
quantile(resamples$xbars, c(0.025, 0.975))

library(boot)
boot.mean <- function(d, i){
  mean(d[i])
}
boots <- boot(data = dat$`PP Opp`,
              statistic = boot.mean,
              R = R)
boot.ci(boots, type="bca")

################################################################################
# Plot the Bootstrap Hypothesis Test
################################################################################
# shift so H0 is true
mean(resamples$xbars)
(delta <- mean(resamples$xbars) - 2.9)


resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low <- 2.9 - delta
high <- 2.9 + delta

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()

library(boot.pval)
boot.pval(boots, theta_null = mu0) # this does something slightly different

################################################################################
# Plot the Bootstrap Test
################################################################################
ggdat.obs <- tibble(xbar=mean(dat$`PP Opp`),
                    y=0)|>
  mutate(mirror = mu0 - (xbar-mu0))

bootstrap.plot <- ggplot() +
  # Plot Resampling Distribution
  stat_density(data=resamples, aes(x=xbars), geom="line", color="lightgrey")+
  stat_density(data=resamples, aes(x=xbars.shifted), geom="line", color="black")+
  geom_hline(yintercept=0)+
  # Plot Observation
  geom_point(data=ggdat.obs,
             aes(x=xbar, y=y,color="Observation"))+
  geom_point(data=ggdat.obs, shape=19,
             aes(x=mirror, y=y,color="Mirror Observation"))+
  scale_color_manual("",
                     values=c("black", "red"))+
  # Plot Observation
  # Clean up
  theme_bw()+
  scale_x_continuous(bquote(bar(x)),
                     limits = c(2.0, 4.5),
                     breaks = round(xbar.breaks,2))+
  ylab("Density")+
  ggtitle("Bootstrap Hypothesis Test")


(t.plot/bootstrap.plot) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


################################################################################
# Pivot Plot
################################################################################
nresamps <- 1000
resamples <- tibble(xbars=rep(NA, nresamps))

npivots <- 1000
pivots <- tibble(xbars=rep(NA, npivots))
pivot.plot.data <- tibble(p05=rep(NA, nresamps),
                          p25=rep(NA, nresamps),
                          p50=rep(NA, nresamps),
                          p75=rep(NA, nresamps),
                          p95=rep(NA, nresamps))

for(i in 1:nresamps){
  ###################################################
  # Resampling 
  ###################################################
  curr.resample <- sample(x = dat$`PP Opp`,
                          size = nrow(dat),
                          replace = T)
  resamples$xbars[i] <- mean(curr.resample)
  
  ###################################################
  #  Pivots
  ###################################################
  for(j in 1:npivots){
    curr.pivot <- sample(x = curr.resample,
                         size = nrow(dat),
                         replace = T)
    pivots$xbars[j] <- mean(curr.pivot)
  }
  # save the 5th, 25th, 50th, 75th, 95th percentiles of how far the pivots
  # differ from the current resample.
  pivot.plot.data$p05[i] <- quantile(pivots$xbars - mean(curr.resample), 0.05)
  pivot.plot.data$p25[i] <- quantile(pivots$xbars - mean(curr.resample), 0.25) 
  pivot.plot.data$p50[i] <- quantile(pivots$xbars - mean(curr.resample), 0.50)
  pivot.plot.data$p75[i] <- quantile(pivots$xbars - mean(curr.resample), 0.75)
  pivot.plot.data$p95[i] <- quantile(pivots$xbars - mean(curr.resample), 0.95)
  ###################################################
}

pivot.plot.data <- pivot.plot.data|>
  mutate(xbarresamp = resamples$xbars)   # Add the resampled means to the tibble

pivot.plot <- ggplot(data=pivot.plot.data)+
  # Plot differences at 5th percentile
  geom_smooth(aes(x=xbarresamp, y=p05, color="5%"))+
  geom_hline(yintercept = mean(pivot.plot.data$p05) + c(-2,2)*sd(pivot.plot.data$p05), 
             linetype = "dotted")+
  # Plot differences at 25th percentile
  geom_smooth(aes(x=xbarresamp, y=p25, color="25%"))+
  geom_hline(yintercept = mean(pivot.plot.data$p25) + c(-2,2)*sd(pivot.plot.data$p25), 
             linetype = "dotted")+
  # Plot differences at 50th percentile
  geom_smooth(aes(x=xbarresamp, y=p50, color="50%"))+
  geom_hline(yintercept = mean(pivot.plot.data$p50) + c(-2,2)*sd(pivot.plot.data$p50), 
             linetype = "dotted")+
  # Plot differences at 75th percentile
  geom_smooth(aes(x=xbarresamp, y=p75, color="75%"))+
  geom_hline(yintercept = mean(pivot.plot.data$p75) + c(-2,2)*sd(pivot.plot.data$p75), 
             linetype = "dotted")+
  # Plot differences at 95th percentile
  geom_smooth(aes(x=xbarresamp, y=p95, color="95%"))+
  geom_hline(yintercept = mean(pivot.plot.data$p95) + c(-2,2)*sd(pivot.plot.data$p95), 
             linetype = "dotted")+
  theme_bw() +
  xlab("Resampled Means")+
  ylab("Difference Between Resamples and Pivots")+
  scale_color_discrete("", 
                       breaks=rev(c("5%", "25%", "50%", "75%", "95%")))+
  labs(caption = "Not perfectly flat, but roughly within +/- 2SD")

pivot.plot

########################################################################
# Brendan's Plot
########################################################################
nresamps <- 250
npivots <- 1000
resamples <- tibble(xbars=rep(NA, nresamps))
pivots <- tibble(xbars=rep(NA, npivots))

resampling.plots <- ggplot()


for(i in 1:nresamps){
  ###################################################
  # Resampling 
  ###################################################
  curr.resample <- sample(x = dat$`PP Opp`,
                          size = nrow(dat),
                          replace = T)
  resamples$xbars[i] <- mean(curr.resample)
  
  ###################################################
  #  Pivots
  ###################################################
  for(j in 1:npivots){
    curr.pivot <- sample(x = curr.resample,
                         size = nrow(dat),
                         replace = T)
    pivots$xbars[j] <- mean(curr.pivot)
  }
  pivots <- pivots |> mutate(xbars = xbars - mean(curr.resample))
  resampling.plots <- resampling.plots +
    stat_density(data=pivots, 
                 aes(x=xbars),
                 color="black",
                 alpha = 0.25,
                 linewidth = 0.1,
                 geom="line")
}
resampling.plots <- resampling.plots +
  geom_hline(yintercept=0)+
  theme_bw()+
  xlab("Difference Between the Resamples and the Pivots")+
  ylab("Density")

resampling.plots

pivot.plot+resampling.plots # Close to perfect!

################################################################################
# Randomization Test
################################################################################
R <- 10000
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- dat$`PP Opp` - mu0
# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift *
    sample(x = c(-1, 1),
           size = length(x.shift),
           replace = T)
  
  rand$xbars[i] <- mean(curr.rand)
}
# Thinking is hard
rand <- rand |>
  mutate(xbars = xbars + mu0) # shifting back

# p-value 
(delta <- abs(mean(dat$`PP Opp`) - mu0))
(low <- mu0 - delta) # mirror
(high<- mu0 + delta)   # xbar

mean(rand$xbars <= low) +
  mean(rand$xbars >= high)


################################################################################
# Plot the Randomization Test
################################################################################
ggdat.obs <- tibble(xbar=mean(dat$`PP Opp`),
                    y=0)|>
  mutate(mirror = mu0 - (xbar-mu0))

randomization.plot <- ggplot() +
  # Plot Randomization Distribution
  stat_density(data=rand, aes(x=xbars), geom="line", color="black")+
  geom_hline(yintercept=0)+
  # Plot Observation
  geom_point(data=ggdat.obs,
             aes(x=xbar, y=y,color="Observation"))+
  geom_point(data=ggdat.obs, shape=19,
             aes(x=mirror, y=y,color="Mirror Observation"))+
  scale_color_manual("",
                     values=c("black", "red"))+
  # Plot Observation
  # Clean up
  theme_bw()+
  scale_x_continuous(bquote(bar(x)),
                     limits = c(2.0, 4.5),
                     breaks = round(xbar.breaks,2))+
  ylab("Density")+
  ggtitle("Randomization Hypothesis Test")


(t.plot/bootstrap.plot/randomization.plot) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


################################################################################
# Randomization Test Confidence Interval
################################################################################
R <- 1000
mu0.iterate <- 0.01
starting.point <- mean(dat$`PP Opp`)
mu.lower <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- dat$`PP Opp` - mu.lower
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.lower) # shifting back
  
  # p-value 
  (delta <- abs(mean(dat$`PP Opp`) - mu.lower))
  (low <- mu.lower - delta) # mirror
  (high<- mu.lower + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))

  if(p.val < 0.05){
    break
  }else{
    mu.lower <- mu.lower - mu0.iterate
  }
}

mu.upper <- starting.point
repeat{
  rand <- tibble(xbars = rep(NA, R))
  
  # PREPROCESSING: shift the data to be mean 0 under H0
  x.shift <- dat$`PP Opp` - mu.upper
  # RANDOMIZE / SHUFFLE
  for(i in 1:R){
    curr.rand <- x.shift *
      sample(x = c(-1, 1),
             size = length(x.shift),
             replace = T)
    
    rand$xbars[i] <- mean(curr.rand)
  }
  # Thinking is hard
  rand <- rand |>
    mutate(xbars = xbars + mu.upper) # shifting back
  
  # p-value 
  (delta <- abs(mean(dat$`PP Opp`) - mu.upper))
  (low <- mu.upper - delta) # mirror
  (high<- mu.upper + delta)   # xbar
  (p.val <- mean(rand$xbars <= low) +
      mean(rand$xbars >= high))
  
  if(p.val < 0.05){
    break
  }else{
    mu.upper <- mu.upper + mu0.iterate
  }
}

c(mu.lower, mu.upper)
