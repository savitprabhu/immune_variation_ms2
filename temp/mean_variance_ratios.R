
rm(list = ls())

# simulated df

df = data.frame(rep(NA, 50))
for(n in 1:30){
  df = cbind(df, name = rnorm(n = 50, mean = 3*n))
  names(df)[n+1] = paste("subset", n)
}
df[1] = NULL

pdf("mean_variance_relationship.pdf", width = 5, height = 5)
plot(apply(df, 2, mean), apply(df, 2, function(x){sd(x)/mean(x)}),
     xlab = "Mean", ylab = "CV = SD/Mean", 
     main = "Simulated data, 30 subsets,\n 50 observations,\n normally distributed")
text(50, 0.2, "negative correlation because CV has \n 'mean' in the denominator",
     col = "red")
plot(apply(df, 2, mean), apply(df, 2, function(x){sd(x)^2}),
     xlab = "Mean", ylab = "Var", 
     main = "Simulated data, 30 subsets,\n 50 observations,\n normally distributed")
text(20, 1.4, "no correlation, \nas expected", col = "red")
dev.off()

