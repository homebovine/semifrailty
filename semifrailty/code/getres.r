obb <- 0.5
mba1 <- mba2 <- matrix(NA, 1000, 10)
mmbb <- matrix(NA, 1000, 2)
for(itr in 783: 1000){
    load(file = paste(paste("../res/", itr, sep = ""),  obb, "res", sep = "_"))
    mba1[itr, ] <- mba[, 1]
    mba2[itr, ] <- mba[, 2]
    mmbb[itr, ] <- mbb
}
