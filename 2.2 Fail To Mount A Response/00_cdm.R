

library(glue)
test <- data.frame(test=c(1,2,4),var=c("a","b","c"))
test

variable <- "test"

new_variable <- paste("Test",variable,sep="_")

test %>% mutate(!!new_variable := 2*!!as.name(variable))

