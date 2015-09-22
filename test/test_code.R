
library(maptpx)
test_data <- read.table('test_data.txt')
head(test_data)
out <- topics(test_data,K=4,tol=0.01);

## error message in the Quasi Newton acceleration



