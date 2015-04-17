library(caTools)
set.seed(123)
spl <- sample.split(archive.df, .90)

evalTrain <- archive.df[spl==T,]
evalTest <- archive.df[spl==F,]

