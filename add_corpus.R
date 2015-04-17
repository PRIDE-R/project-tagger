library(tm)
evalTrain$AllText <- do.call(paste, evalTrain[,c("project.title","project.description")])
evalTest$AllText <- do.call(paste, evalTest[,c("project.title","project.description")])

corpusAll <- Corpus(VectorSource(c(evalTrain$AllText, evalTest$AllText)))
corpusAll <- tm_map(corpusAll, tolower)
corpusAll <- tm_map(corpusAll, PlainTextDocument)
corpusAll <- tm_map(corpusAll, removePunctuation)
corpusAll <- tm_map(corpusAll, removeWords, stopwords("english"))
corpusAll <- tm_map(corpusAll, stripWhitespace)
corpusAll <- tm_map(corpusAll, stemDocument)

# Generate term matrix
dtmAll <- DocumentTermMatrix(corpusAll)
sparseAll <- removeSparseTerms(dtmAll, 0.99)
allWords <- data.frame(as.matrix(sparseAll))

colnames(allWords) <- make.names(colnames(allWords))

# Find most significative terms
allWordsTrain2 <- head(allWords, nrow(evalTrain))
allWordsTrain2$Popular <- evalTrain$project.tag
logModelAllWords <- glm(Popular~., data=allWordsTrain2, family=binomial)
all_three_star_terms <- names(which(summary(logModelAllWords)$coefficients[,4]<0.001))
all_two_star_terms <- names(which(summary(logModelAllWords)$coefficients[,4]<0.01))
all_one_star_terms <- names(which(summary(logModelAllWords)$coefficients[,4]<0.05))

# Leave just those terms that are different between popular and unpopular articles
allWords <- subset(allWords, 
                   select=names(allWords) %in% all_one_star_terms)

# Split again
allWordsTrain <- head(allWords, nrow(evalTrain))
allWordsTest <- tail(allWords, nrow(evalTest))

# Add to dataframes
evalTrain <- cbind(evalTrain, allWordsTrain)
evalTest <- cbind(evalTest, allWordsTest)

# Remove original text variables
evalTrain$AllText <- NULL
evalTest$AllText <- NULL

