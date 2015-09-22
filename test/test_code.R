
library(maptpx)
test_data <- read.table('../data/test_data.txt')
out <- topics(test_data,K=4,tol=0.01);


## gives error message
