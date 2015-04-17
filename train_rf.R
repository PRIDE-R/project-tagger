library(randomForest)
library(ROCR)

set.seed(123)
rfModelSpecies <- randomForest(project.tags ~ species + submissionType, data = evalTrain)
# Calculate accuracy on the test set
rfPredSpecies <- predict(rfModelSpecies, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredSpecies)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)


set.seed(123)
rfModelSpeciesTissues <- randomForest(project.tags ~ species + tissues + submissionType, data = evalTrain)
# Calculate accuracy on the test set
rfPredSpeciesTissues <- predict(rfModelSpeciesTissues, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredSpeciesTissues)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)


set.seed(123)
rfModelSpeciesTissuesPtm <- randomForest(project.tags ~ species + tissues + ptm.names + submissionType, data = evalTrain)
# Calculate accuracy on the test set
rfPredSpeciesTissuesPtm <- predict(rfModelSpeciesTissuesPtm, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredSpeciesTissuesPtm)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)


set.seed(123)
rfModelMeta <- randomForest(project.tags ~ species + tissues + ptm.names + instrument.names + submissionType, data = evalTrain)
# Calculate accuracy on the test set
rfPredMeta <- predict(rfModelMeta, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredMeta)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)

set.seed(123)
rfModelSpeciesTissuesInstrument <- randomForest(project.tags ~ species + tissues + instrument.names + submissionType, data = evalTrain)
# Calculate accuracy on the test set
rfPredSpeciesTissuesInstrument <- predict(rfModelSpeciesTissuesInstrument, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredSpeciesTissuesInstrument)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)


set.seed(123)
rfModelMetaText <- randomForest(project.tags ~ . - accession - publication.date - project.title - project.description, data = evalTrain, ntree = 1000)
# Calculate accuracy on the test set
rfPredMetaText <- predict(rfModelMetaText, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredMetaText)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)

set.seed(123)
rfModelMetaTextNoPtm <- randomForest(project.tags ~ . - accession - publication.date - project.title - project.description - ptm.names, data = evalTrain)
# Calculate accuracy on the test set
rfPredMetaTextNoPtm <- predict(rfModelMetaTextNoPtm, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredMetaTextNoPtm)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)

set.seed(123)
rfModelSpeciesTissueText <- randomForest(project.tags ~ . - accession - publication.date - project.title - project.description - ptm.names - instrument.names, data = evalTrain)
# Calculate accuracy on the test set
rfPredSpeciesTissueText <- predict(rfModelSpeciesTissueText, newdata=evalTest)
t <- table(evalTest$project.tags, rfPredSpeciesTissueText)
(t[1,1] + t[2,2] + t[3,3] + t[4,4] + t[5,5] + t[6,6] + t[7,7] + t[8,8]) / nrow(evalTest)

# do cv
set.seed(123)
x <- archive.df
x$project.tag <- NULL
x$accession <- NULL
x$project.title <- NULL
x$project.description <- NULL
x$publication.date <- NULL
x$num.assays <- NULL
y <- archive.df$project.tag
rf.cv <- rfcv(x, y, cv.fold=10)
with(rf.cv, plot(n.var, error.cv))
rf.cv
