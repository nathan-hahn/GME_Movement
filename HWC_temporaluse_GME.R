## Plots of csum ag usage by ag usage class and individual
## 11-19-19

# Load packages
library(ggplot2)
theme_set(theme_bw() + theme(panel.border = element_rect(), panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.title=element_text(size=14), axis.text = element_text(size=12))) 
library(lubridate)
library(dplyr)
library(xts)
library(zoo)

Sys.setenv(TZ="Africa/Nairobi") 

## Cumulative Ag Usage Plots
#TODO: Add "speed" coloring to csum graphs -- show entropy of raiding. Steeper slope is hotter/faster
#TODO: Analyse HEC data to justify window size for rolling averages

## Load data - used points plus cluster classification (2019)
df <- readRDS('GMEcollars_001_used_2019-12-18.rds')
colnames(df)[colnames(df)=="Name"] <- "name"
colnames(df)[colnames(df)=="date"] <- "Fixtime"

df$ag.used <- ifelse(df$lc.estes == 1, 1, 0)
df$year <- year(df$Fixtime)
df <- droplevels(df)

## hist of dist2ag is interesting by individual. A lot of individuals in class 1 and 2 have bimodal distributions
## Not plotted

# p <- lapply(split(df, df$name),
#             function(x)
#               plot(cumsum(x$ag.used), main = x$name[1]))


## Cumulative Sum Graphs
csum <- lapply(split(df, df$id),
               function(x) cumsum(x$ag.used))

df$csum <- unlist(csum, use.names = FALSE)

# filter to 2017-2018 and to select eles
df.2018 <- droplevels(subset(df, year <= 2018))
# sample 3 individuals from each group for data clarity. Code for plotting all eles is below
sample <- c("Bonchugu", "Ukungu", "Mokomre", "Ukwaju", "Lowana", "Kimbizwa")
df.sample <- droplevels(subset(df.2018, df.2018$name %in% sample))

#TODO: Automate re-ordering of name factor level based on ag.class.mean
df.sample$name <- as.factor(df.sample$name)
df.sample$ag.class.mean <- as.factor(df.sample$ag.class.mean)
df.sample$name <- factor(df.sample$name, levels = sample)

# create color palettes for each ag class
pal1 <- colorRampPalette(c("grey", "black"))(2) #n
pal2 <- colorRampPalette(c("blue", "green"))(2) #n
pal3 <- colorRampPalette(c("orange", "red"))(2) #n

# ggplot version
ggplot(df.sample, aes(x = Fixtime, y = csum, color = interaction(as.factor(ag.class.mean), as.factor(name)))) + geom_line() +
  scale_colour_manual(values = c(pal1,pal2,pal3)) +
  labs(x = "Time", y = "cummulative ag usage", title = "Temporal Agriculture Usage - 2018", color = "agClass.Individual")

## Plotting of cumulative ag usage for all elephants
# plot all eles cum ag usage across time. See some eles stand out, but different ones start and stop at different times.
ggplot(df, aes(x = Fixtime, y = csum)) + geom_line(aes(colour = name))

# filter to first set of collar deployments (2017)
t <- droplevels(subset(df, Fixtime < "2018-08-01 00:00:00"))
t <- unique(t$id)

df.2018 <- droplevels(subset(df, id %in% t))
csum <- lapply(split(df.2018, df.2018$id),
               function(x) cumsum(x$ag.used))
df.2018$csum <- unlist(csum, use.names = FALSE)

# plot all eles cum ag usage for 2017/2018 - ggplot
names <- df.2018 %>%
  group_by(name, ag.class.mean) %>%
  tally() %>%
  arrange(ag.class.mean) 
df.2018$name <- as.factor(df.2018$name)
df.2018$ag.class.mean <- as.factor(df.2018$ag.class.mean)
df.2018$name <- factor(df.2018$name, levels = names$name)

# create color palettes for each ag class
pal1 <- colorRampPalette(c("grey", "black"))(5) #n
pal2 <- colorRampPalette(c("blue", "green"))(5) #n
pal3 <- colorRampPalette(c("orange", "red"))(1) #n

p1 <- ggplot(df.2018, aes(x = Fixtime, y = csum, color = interaction(as.factor(ag.class.mean), as.factor(name)))) + geom_line() +
  scale_colour_manual(values = c(pal1,pal2,pal3)) +
  labs(x = "Time", y = "cummulative ag usage", title = "Cummulative ag usage Grumeti: first collar deployments", color = "agClass.Individual")

p1

## Temporal use -- occupancy by day
# by day

fixes.d <- df %>%
  mutate(day = date(Fixtime)) %>% # simplify to date for grouping
  group_by(name, day) %>%
  tally()

ag.day <- df %>%
  mutate(day = date(Fixtime)) %>% # simplify to date for grouping
  group_by(name, ag.class.mean, day) %>%
  tally(lc.estes == 1) %>% # tally fixes in ag
  droplevels() %>%
  ungroup() %>%
  rename(n.ag = n) %>%
  mutate(n = fixes.d$n) %>%
  arrange(name, day) %>%
  mutate(day.occupancy = n.ag/n) %>%
  ungroup()

(ag.day) # the month for each individual where ag occupancy is highest (across all)

df$day.occupancy <- ag.day

# plot daily occupancy over time
names = c("Bonchugu", "Imara" ,"Ukwaju", "Lowana")
t <- filter(ag.day, name %in% names)
t$name <- factor(t$name, levels = names)

p2 <- ggplot(t, aes(x = day, y = day.occupancy)) + geom_line() + facet_grid(name~.)
p2

## Rolling stats
# rollstat function applied with variable window size
# res.x is a new data.frame with a rolling mean and upper/lower CI column


# function to return mean, sd, 95% conf intervals
rollstat <- function(x, na.rm = TRUE) {
  # x     = numeric vector
  # na.rm = boolean, whether or not to remove NA's
  
  m  <- mean(x, na.rm = na.rm)
  s  <- sd(x, na.rm = na.rm)
  hi <- m + 2*s
  lo <- m - 2*s
  
  ret <- c(mean = m, stdev = s, hi.95 = hi, lo.95 = lo) 
  return(ret)
}

# Apply rollstate. Converts to time series object and applies rollstat function
window <- 90*24*2 # set number of days and convert to 30 minute fixes
split <- split(df, df$name) # careful of what you split by. May cuase problems splitting by name if there are large gaps between recollars. Let na.rm = FALSE in rollstat

system.time({
ts <- lapply(split, as.ts, start = split$Fixtime[1], frequency = 1) # frequency set with 30min fixes
tx <- lapply(ts, as.xts)
res <- lapply(seq_along(tx), function(i) rollapply(tx[[i]]$ag.used, width = window, align = "center", by.column = FALSE,
              FUN = rollstat, na.rm = FALSE))
merge <- lapply(seq_along(res), function(i) cbind(split[[i]], as.data.frame(res[[i]])))
})

plot(merge[[28]]$Fixtime, merge[[28]]$mean) #check

# unlist 
res.90 <- do.call(rbind, merge)

# plot
res.90$ag.used <- ifelse(res.90$ag.used == 1, 0.75, NA) # adjust ag.used values for plotting for visibility
p90 <- ggplot(res.90, aes(x = Fixtime, color = name)) + #, color = name
  geom_point(aes(y = ag.used), color = "grey40", size = .2, alpha = 0.5) +
  #geom_ribbon(aes(ymin = lo.95, ymax = hi.95), size = .5, alpha = 0.4) +
  geom_line(aes(y = mean), size = 1, alpha = 0.7) +
  facet_wrap(~ name) +
  labs(title = "Ag Usage: Volatility and Trend",
       subtitle = "60-Day Moving Average with 95% CI Bands") +
  theme(legend.position="none")

p90






