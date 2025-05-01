library(tidyverse)
library(e1071)
###########
#Lab 13
###########

#1.a
dat1 = read_csv("zebrafinches.csv")
dat = dat1$further
t.val = as.numeric(t.test(dat, alternative = "less", mu = 0)$statistic)
cdf.val = pnorm(t.val)
pdf.val = dnorm(t.val)
skew = skewness(dat)
n = length(dat)
error = cdf.val + (skew/sqrt(n)*(2*(t.val)^2+1)/6*pdf.val)
# What does the t-value tell us? 
#Tells us that the value is likely well outside of the null?
#1.b
total.error = vector()
seq = seq(-10,10,.01)
for(t in seq){
  curr.cdf.val = pnorm(t)
  curr.pdf.val = dnorm(t)
  skew = skewness(dat)
  n = length(dat)
  curr.error = (skew/sqrt(n)*(2*(t)^2+1)/6*curr.pdf.val)
  total.error = append(total.error, curr.error)
}
total.error = tibble(x = total.error) |>
  mutate(seq)

ggplot() +
  geom_line(data = total.error, aes(x = seq, y = x))

#As value increases beyond null hypothesis, total error increases?

#1.c
t.val.c = qnorm(.05) #Should degrees of freedom be 29 or 30?
pdf.val.c = dnorm(t.val.c)
n.c = (skew/(6*(.1*.05))*(2*t.val.c^2+1)*pdf.val.c)^2

#2.a
###################
# RESAMPLING
###################

#Farther
t.dist.dat = vector()
s = sd(dat)
for(i in 1:10000){
  curr.val = (mean(sample(dat, replace = TRUE)) - 0)/(s/sqrt(n))
  t.dist.dat = append(t.dist.dat, curr.val)
}
ncp = -1*(mean(t.dist.dat))
resamples.null.farther = t.dist.dat + ncp

#Closer
t.close.dat = vector()
s = sd(dat1$closer)
for(i in 1:10000){
  curr.val = (mean(sample(dat1$closer, replace = TRUE)) - 0)/(s/sqrt(n))
  t.close.dat = append(t.close.dat, curr.val)
}
ncp2 = -1*(mean(t.close.dat))
resamples.null.closer = t.close.dat + ncp2
#Difference
t.diff.dat = vector()
s = sd(dat1$diff)
for(i in 1:10000){
  curr.val = (mean(sample(dat1$diff, replace = TRUE)) - 0)/(s/sqrt(n))
  t.diff.dat = append(t.diff.dat, curr.val)
}
ncp3 = -1*(mean(t.diff.dat))
resamples.null.diff = t.diff.dat + ncp3
#2.b
##########################
# CALCULATING P-VALUES
##########################
#Farther p-value
x.bar.far = mean(t.dist.dat)
p.val.far = mean(resamples.null.farther <= x.bar.far)
#Closer p-value
x.bar.close = mean(t.close.dat)
p.val.close = mean(resamples.null.closer >=x.bar.close) #Why greater than or equal to?
#Difference p-value
x.bar.diff = mean(t.diff.dat)
x.bar.diff.lower = 0-x.bar.diff
x.bar.diff.upper = 0+x.bar.diff

p.val.diff = mean(resamples.null.diff <= x.bar.diff.lower) + mean(resamples.null.diff >= x.bar.diff.upper)

#2.c
#Farther 
quantile(resamples.null.farther, prob = .025)
#Closer
quantile(resamples.null.closer, prob = .025)
#Difference
quantile(resamples.null.diff, prob = .025)
#############################################
# Calculating Bootstrap confidence intervals
#############################################
#Farther 
n = length(dat1$further)
sd.f = sd(dat1$further)
lower.f = quantile(t.dist.dat, prob = .025)
x.bar.lower.f = lower.f/sqrt(n)*sd.f
upper.f = quantile(t.dist.dat, prob = .975)
x.bar.upper.f = upper.f/sqrt(n)*sd.f
#Closer
sd.c = sd(dat1$closer)
lower.c = quantile(t.close.dat, prob = .025)
x.bar.lower.c = lower.c/sqrt(n)*sd.c
upper.c = quantile(t.close.dat, prob = .975)
x.bar.upper.c = upper.c/sqrt(n)*sd.c

#Difference
sd.d = sd(dat1$diff)
lower.d = quantile(t.diff.dat, prob = .025)
x.bar.lower.d = lower.d/sqrt(n)*sd.d
upper.d = quantile(t.diff.dat, prob = .975)
x.bar.upper.d = upper.d/sqrt(n)*sd.d


#3.a
##########################
# RANDOMIZATION
##########################

############
# Further Data
############
mu0 = 0
R = 10000
rand = tibble(xbars = rep(NA, 10000))
samp.mean = mean(dat1$further)
x.shift.f = dat1$further - mu0
#Doing the randomization process
for(i in 1:R){
  curr.samp = x.shift * sample(x = c(-1,1), 
                                   size = length(x.shift.f),
                                   replace = T)
  rand$xbars[i] = mean(curr.samp)
}
#Not necessary
rand = rand |>
  mutate(xbars = xbars + mu0)
  
#P-value
furth.p = mean(rand<=samp.mean)

#############
# Closer Data
#############

rand.cl = tibble(xbars = rep(NA, 10000))
samp.mean.cl = mean(dat1$closer)
x.shift.c = dat1$further - mu0
#Doing the randomization process
for(i in 1:R){
  curr.samp = x.shift * sample(x = c(-1,1), 
                               size = length(x.shift.c),
                               replace = T)
  rand.cl$xbars[i] = mean(curr.samp)
}
#Not necessary
rand.cl = rand.cl |>
  mutate(xbars = xbars + mu0)

#P-value
closer.p = mean(rand>=samp.mean.cl)

#############
# Closer Data
#############

rand.diff = tibble(xbars = rep(NA, 10000))
samp.mean.diff = mean(dat1$diff)
x.shift.diff = dat1$diff - mu0
#Doing the randomization process
for(i in 1:R){
  curr.samp = x.shift * sample(x = c(-1,1), 
                               size = length(x.shift.diff),
                               replace = T)
  rand.diff$xbars[i] = mean(curr.samp)
}
#Not necessary
rand.diff = rand.diff |>
  mutate(xbars = xbars + mu0)

#P-value
delta = samp.mean.diff
low = mu0 - delta 
high = mu0 + delta
diff.p = mean(rand.diff>=high) + mean(rand.diff<=low)
#Why make the mirror value?



######################
# Confidence Intervals
######################

#############
# Farther
#############

# Lower val
R = 1000
starting.point.further = mean(dat1$further)
lower.val.further = starting.point.further
samp.mean = mean(dat1$further)
repeat{
  rand = tibble(xbars = rep(NA, 1000))
  x.shift.f = dat1$further - lower.val.further
  #Doing the randomization process
  for(i in 1:R){
    curr.samp = x.shift.f * sample(x = c(-1,1), 
                                 size = length(x.shift.f),
                                 replace = T)
    rand$xbars[i] = mean(curr.samp)
  }
  #Shifting back
  rand = rand |>
    mutate(xbars = xbars + lower.val.further)
  
  #P-value
  delta = abs(samp.mean) - lower.val.further
  lower = lower.val.further - delta #Mirror?
  upper = lower.val.further + delta
  further.p = mean(rand$xbars <= lower) + mean(rand$xbars >= upper )
  print(further.p)
  if(further.p < .05) {#Should it be .025 or .5?
    break
  }
  else{
    lower.val.further = lower.val.further - .01
  }
}
lower.val.further

#Upper val
R = 1000
starting.point = mean(dat1$further)
upper.val.further = starting.point.further
samp.mean = mean(dat1$further)
repeat{
  rand = tibble(xbars = rep(NA, 1000))
  x.shift.f = dat1$further - upper.val.further
  #Doing the randomization process
  for(i in 1:R){
    curr.samp = x.shift.f * sample(x = c(-1,1), 
                                   size = length(x.shift.f),
                                   replace = T)
    rand$xbars[i] = mean(curr.samp)
  }
  #Shifting back
  rand = rand |>
    mutate(xbars = xbars + upper.val.further)
  
  #P-value
  delta = abs(samp.mean) - upper.val.further
  lower = upper.val.further - delta #Mirror?
  upper = upper.val.further + delta
  further.p = mean(rand$xbars <= lower) + mean(rand$xbars >= upper )
  if(further.p < .05) {#Should it be .025 or .5?
    break
  }
  else{
    upper.val.further = upper.val.further + .01
  }
}
upper.val.further



#Extra work
rand = tibble(xbars = rep(NA, 1000))
lower.val.further = -.5
x.shift.f = dat1$further - lower.val.further
#Doing the randomization process
for(i in 1:R){
  curr.samp = x.shift.f * sample(x = c(-1,1), 
                                 size = length(x.shift.f),
                                 replace = T)
  rand$xbars[i] = mean(curr.samp)
}

rand = rand |>
  mutate(xbars = xbars + lower.val.further)


ggplot() + 
  geom_density(data = rand, aes(x = xbars))


