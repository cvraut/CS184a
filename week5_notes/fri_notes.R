#install.packages('RCurl')
require('RCurl')
mlbdfurl <- 'http://goo.gl/rih9v9'

mlb.df <- textConnection(getURL(mlbdfurl, followlocation  = TRUE))

mlb.df <- read.table(mlb.df, header = TRUE)

# get summary of data

dim(mlb.df)
View(mlb.df)

# with() --> treat a df as if all cols where vectors w/ same names as field names
# kinda like localized load of a df
with(
  data = mlb.df,
  expr = { print(summary(Height)); print(summary(Weight)); }
)

# cut() --> creates a factor w/ lvls determined by numerical bins
# in this example make 2 bins >6ft or <6ft height
mlb.df$HeightBin <- with(
  data = mlb.df,
  expr = cut(x = Height, breaks = c(min(Height), 72, max(Height)), include.lowest = TRUE)
)
head(mlb.df$HeightBin)

# split() --> splits data by a factor & returns result
# good for turning 1df --> 2+ df's
mlb.height.list <- split(x = mlb.df, f = mlb.df$HeightBin)
View(mlb.height.list$`[67,72]`)

# do.call() --> applies all members of the args param to the funct specified in the whatt param
mlb.df.2 <- do.call(what = rbind, args = mlb.height.list)
dim(mlb.df.2)
View(mlb.df.2)

# apply() is kinda like map for cols
apply(X = mlb.df[,c('Height', 'Weight', 'Age')], MARGIN = 2, FUN = mean)
apply(X = mlb.df[,c('Height', 'Weight', 'Age')], MARGIN = 2, FUN = summary)
apply(X = mlb.df[,c('Position', 'HeightBin')], MARGIN = 2, FUN = table)

# lapply() is like apply for lists
lapply(X = mlb.df[c('Height', 'Weight', 'Age')], MARGIN = 2, FUN = summary)

# sapply() is like lapply but it returns things differently
sapply(X = mlb.df[c('Height', 'Weight', 'Age')], MARGIN = 2, FUN = summary)

# tapply() applies a function to data grouped by a factor
HeightTable <- with(data = mlb.df,
                    tapply(X = Position, INDEX = HeightBin, FUN = table)
)
do.call(
  what = rbind,
  args = HeightTable
)

# aggregate() --> apply a funct to data grouped by 1 or more factors
mlb.HxP.df <- with(
  
  data = mlb.df,
  
  aggregate(x = mlb.df[,4:6], by = list(Position = Position, HeightBin = HeightBin), FUN = median)
  
)
View(mlb.HxP.df)

# by() --> split a df by 1/+ indices & eval the expr on each subset
mlb.lmWxH.list <- with(
  data = mlb.df,
  expr = by(data = mlb.df, INDICES = Position, function(x) lm(Weight ~ Height, data = x)$coef)
)
head(x = mlb.lmWxH.list, n = 2)


do.call(
  what = rbind,
  args = mlb.lmWxH.list
)
