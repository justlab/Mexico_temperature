#function to create Cross validated 90% and 10% splits

splitdf <- function(dataframe, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, trunc(length(index)/10))
  testset <- dataframe[trainindex, ] # the stuff you leave behind to predict to afterwards (ie the 10% left out)
  trainset <- dataframe[-trainindex, ] # the stuff you run the model on (as in "train" the model) (ie the 90% left out)
  list(trainset=trainset,testset=testset)
}