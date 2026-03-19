##Applying methods from a previous analysis session
#mydata<-read.csv("association.data") #importing this dataset
#head(mydata) #visualization of data
#str(mydata) #structure of the dataset

##########################################
#This extended analysis combines previously learned statistical inference methods
#The data extension condenses methodology from multiple analysis sessions:
##Initial interpretation from Day 1 and Day 2
##Full workflows from Day 3, Day 4, Day 5
##########################################
library(DescTools) #mode
install.packages("psych") #geometric and harmonic mean
library(psych)

###########################################################################
###This was used in the 03-02-model notebook
mydata1 <- read.csv("~/Downloads/Data/Data/association-data.csv")
View(association.data)
head(mydata) #visualization of data
str(mydata) #structure of the dataset
##########################################################################
##########################
#FOLLOWING THE NOTEBOOK 03-02
##########################
population_height <- rnorm(1e6, 80, 10) # population
plot(density(population_height))
random_sample_height <- sample(population_height, 100) # D={x1, ..., xn} is a random sample
hist(random_sample_height)
##########################
#but what is the "parameter" \theta?
mean(mydata1$height) #[1] 179.9516 this is the true population mean!


###########################
#Working on the first dataset (mydata1): the disease score is associated with colour blindeness!!
##########################
#FOLLOWING THE DATA ANALYSIS IN THE NOTEBOOK 03-03
##########################
#####DESCRIPTIVE STATISTICS#####
#Typical quantitative summaries of data fall into several classes:
#• location or central tendency measures (mean, median, geometric mean, harmonic mean),
#• mode,
#• scale (or spread of the data): (min,max,range,sample variance, standard deviation, interquartile range, quartiles),
#• skewness (or asymmetry): (mean-median).

#HEIGHT#
mean_height<-mean(mydata1$height)
median_height<-median(mydata1$height)
geometric.mean_height<-geometric.mean(mydata1$height)
harmonic.mean_height<-harmonic.mean(mydata1$height)
min_height<-min(mydata1$height)
max_height<-max(mydata1$height)
range_height<-range(mydata1$height)
var_height<-var(mydata1$height)
sd_height<-sd(mydata1$height)
IQR_height<-IQR(mydata1$height)
skewness_height<-mean(mydata1$height)-median(mydata1$height)

#MASS#
mean_mass<-mean(mydata1$mass)
median_mass<-median(mydata1$mass)
geometric.mean_mass<-geometric.mean(mydata1$mass)
harmonic.mean_mass<-harmonic.mean(mydata1$mass)
min_mass<-min(mydata1$mass)
max_mass<-max(mydata1$mass)
range_mass<-range(mydata1$mass)
var_mass<-var(mydata1$mass)
sd_mass<-sd(mydata1$mass)
IQR_mass<-IQR(mydata1$mass)
skewness_mass<-mean(mydata1$mass)-median(mydata1$mass)

#METABOLITE_A#
#na.rm = TRUE allows for the NA variables to be ignored
mean_metabolite_A<-mean(mydata1$metabolite_A, na.rm = TRUE)
median_metabolite_A<-median(mydata1$metabolite_A, na.rm = TRUE)
geometric.mean_metabolite_A<-geometric.mean(mydata1$metabolite_A, na.rm = TRUE)
harmonic.mean_metabolite_A<-harmonic.mean(mydata1$metabolite_A, na.rm = TRUE)
min_metabolite_A<-min(mydata1$metabolite_A, na.rm = TRUE)
max_metabolite_A<-max(mydata1$metabolite_A, na.rm = TRUE)
range_metabolite_A<-range(mydata1$metabolite_A, na.rm = TRUE)
var_metabolite_A<-var(mydata1$metabolite_A, na.rm = TRUE)
sd_metabolite_A<-sd(mydata1$metabolite_A, na.rm = TRUE)
IQR_metabolite_A<-IQR(mydata1$metabolite_A, na.rm = TRUE)
skewness_metabolite_A<-mean(mydata1$metabolite_A, na.rm = TRUE)-median(mydata1$metabolite_A, na.rm = TRUE)

#METABOLITE_B#
mean_metabolite_B<-mean(mydata1$metabolite_B)
median_metabolite_B<-median(mydata1$metabolite_B)
geometric.mean_metabolite_B<-geometric.mean(mydata1$metabolite_B)
harmonic.mean_metabolite_B<-harmonic.mean(mydata1$metabolite_B)
min_metabolite_B<-min(mydata1$metabolite_B)
max_metabolite_B<-max(mydata1$metabolite_B)
range_metabolite_B<-range(mydata1$metabolite_B)
var_metabolite_B<-var(mydata1$metabolite_B)
sd_metabolite_B<-sd(mydata1$metabolite_B)
IQR_metabolite_B<-IQR(mydata1$metabolite_B)
skewness_metabolite_B<-mean(mydata1$metabolite_B)-median(mydata1$metabolite_B)

#METABOLITE_C#
mean_metabolite_C<-mean(mydata1$metabolite_C)
median_metabolite_C<-median(mydata1$metabolite_C)
geometric.mean_metabolite_C<-geometric.mean(mydata1$metabolite_C)
harmonic.mean_metabolite_C<-harmonic.mean(mydata1$metabolite_C)
min_metabolite_C<-min(mydata1$metabolite_C)
max_metabolite_C<-max(mydata1$metabolite_C)
range_metabolite_C<-range(mydata1$metabolite_C)
var_metabolite_C<-var(mydata1$metabolite_C)
sd_metabolite_C<-sd(mydata1$metabolite_C)
IQR_metabolite_C<-IQR(mydata1$metabolite_C)
skewness_metabolite_C<-mean(mydata1$metabolite_C)-median(mydata1$metabolite_C)

#METABOLITE_D#
mean_metabolite_D<-mean(mydata1$metabolite_D)
median_metabolite_D<-median(mydata1$metabolite_D)
geometric.mean_metabolite_D<-geometric.mean(mydata1$metabolite_D)
harmonic.mean_metabolite_D<-harmonic.mean(mydata1$metabolite_D)
min_metabolite_D<-min(mydata1$metabolite_D)
max_metabolite_D<-max(mydata1$metabolite_D)
range_metabolite_D<-range(mydata1$metabolite_D)
var_metabolite_D<-var(mydata1$metabolite_D)
sd_metabolite_D<-sd(mydata1$metabolite_D)
IQR_metabolite_D<-IQR(mydata1$metabolite_D)
skewness_metabolite_D<-mean(mydata1$metabolite_D)-median(mydata1$metabolite_D)

#METABOLITE_E#
mean_metabolite_E<-mean(mydata1$metabolite_E)
median_metabolite_E<-median(mydata1$metabolite_E)
geometric.mean_metabolite_E<-geometric.mean(mydata1$metabolite_E)
harmonic.mean_metabolite_E<-harmonic.mean(mydata1$metabolite_E)
min_metabolite_E<-min(mydata1$metabolite_E)
max_metabolite_E<-max(mydata1$metabolite_E)
range_metabolite_E<-range(mydata1$metabolite_E)
var_metabolite_E<-var(mydata1$metabolite_E)
sd_metabolite_E<-sd(mydata1$metabolite_E)
IQR_metabolite_E<-IQR(mydata1$metabolite_E)
skewness_metabolite_E<-mean(mydata1$metabolite_E)-median(mydata1$metabolite_E)

#DISEASE SCORE#
mean_disease_score<-mean(mydata1$disease_score)
median_disease_score<-median(mydata1$disease_score)
#geometric.mean_disease_score<-geometric.mean(mydata1$disease_score) ###negative values present so result is NaNs produced
#harmonic.mean_disease_score<-harmonic.mean(mydata1$disease_score) ###negative values present so result is NaNs produced
min_disease_score<-min(mydata1$disease_score)
max_disease_score<-max(mydata1$disease_score)
range_disease_score<-range(mydata1$disease_score)
var_disease_score<-var(mydata1$disease_score)
sd_disease_score<-sd(mydata1$disease_score)
IQR_disease_score<-IQR(mydata1$disease_score)
skewness_disease_score<-mean(mydata1$disease_score)-median(mydata1$disease_score)
########
#shortened
for (i in continuous_variables){
  variable_vectors<-mydata1[[i]]
  if (i=="metabolite_A"){
    variable_vectors<-na.omit(variable_vectors) #OpenAI for removal of NA variables in metabolite_A
  }
  results_con_variables<-c(
    mean(variable_vectors),
    median(variable_vectors),
    min(variable_vectors),
    max(variable_vectors),
    var(variable_vectors),
    sd(variable_vectors),
    IQR(variable_vectors),
    skewness<-mean(variable_vectors)-median(variable_vectors)
  )
  print(results_con_variables)
}


#visualization of the data on graphs
#some of the variables are continuous, while others are binary

#####CONTINUOUS NUMERICAL VARIABLES
par(mfrow=c(4,2))
hist(mydata1$height, breaks=50, col = "brown")
hist(mydata1$mass, breaks=50, col = "grey")
hist(mydata1$metabolite_A, breaks=50, col = "lightblue")
hist(mydata1$metabolite_B, breaks=50, col = "orange")
hist(mydata1$metabolite_C, breaks=50, col= "lightpink")
hist(mydata1$metabolite_D, breaks=50, col = "darkgreen")
hist(mydata1$metabolite_E, breaks=50, col = "violet")
hist(mydata1$disease_score, breaks=50, col = "red")


############################
#FOLLOWING THE NOTEBOOK 03-04 FOR ESTIMATORS
############################
#estimates 1 and 2 are calculated for all continuous variables as shown in Session 3
calculationforestimate2<-function(a,b) abs(a-b)/2+min(c(a,b)) #this was done by calculating the harmonic mean
#HEIGHT#
estimate1_height<-mean(mydata1$height) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_height) #179.9516
estimate2_height<-calculationforestimate2(max(mydata1$height),min(mydata1$height))
print(estimate2_height) #179.3979

#MASS#
estimate1_mass<-mean(mydata1$mass) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_mass) #74.03284
estimate2_mass<-calculationforestimate2(max(mydata1$mass),min(mydata1$mass))
print(estimate2_mass) #73.34385

#METABOLITE A#
estimate1_metabolite_A<-mean(mydata1$metabolite_A, na.rm = TRUE) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_metabolite_A) #0.8295026
estimate2_metabolite_A<-calculationforestimate2(max(mydata1$metabolite_A, na.rm = TRUE),min(mydata1$metabolite_A, na.rm = TRUE))
print(estimate2_metabolite_A) # 0.6467399

#METABOLITE B#
estimate1_metabolite_B<-mean(mydata1$metabolite_B) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_metabolite_B) #1.66165
estimate2_metabolite_B<-calculationforestimate2(max(mydata1$metabolite_B),min(mydata1$metabolite_B))
print(estimate2_metabolite_B) # 1.317804

#METABOLITE C#
estimate1_metabolite_C<-mean(mydata1$metabolite_C) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_metabolite_C) #0.2818546
estimate2_metabolite_C<-calculationforestimate2(max(mydata1$metabolite_C),min(mydata1$metabolite_C))
print(estimate2_metabolite_C) # 0.3950708

#METABOLITE D#
estimate1_metabolite_D<-mean(mydata1$metabolite_D) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_metabolite_D) #4.957974
estimate2_metabolite_D<-calculationforestimate2(max(mydata1$metabolite_D),min(mydata1$metabolite_D))
print(estimate2_metabolite_D) # 4.874373

#METABOLITE E#
estimate1_metabolite_E<-mean(mydata1$metabolite_E) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_metabolite_E) #10.19689
estimate2_metabolite_E<-calculationforestimate2(max(mydata1$metabolite_E),min(mydata1$metabolite_E))
print(estimate2_metabolite_E) # 10.08065

#DISEASE SCORES#
estimate1_disease_score<-mean(mydata1$disease_score) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_disease_score) #-0.2680843
estimate2_disease_score<-calculationforestimate2(max(mydata1$disease_score),min(mydata1$disease_score))
print(estimate2_disease_score) #2.968784
########
#shortened
continuous_variables<-c("height","mass","metabolite_A","metabolite_B","metabolite_C","metabolite_D","metabolite_E","disease_score")
for (i in continuous_variables){
  x<-mydata1[[i]]
  estimate1<-mean(x, na.rm=TRUE)
  estimate2<-calculationforestimate2(max(x, na.rm=TRUE), min(x, na.rm=TRUE))
print(estimate1) 
print(estimate2)
}




##############################
#Visualisation of the estimates
#############################
#HEIGHT#
tD1_height<-tD2_height<-c()
x1<-mydata1$height
n1<-50
theta_height <- mean(x1)
for (i in 1:1000) {
  samples1 <- sample(x1, size=n1, replace=TRUE)  # repeated sampling from data
  tD1_height<- c(tD1_height, mean(samples1))
  tD2_height<- c(tD2_height, (min(samples1) + max(samples1)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_height, col="blue", main = "Estimate1 height (mean): ")
abline(v = theta_height, col="red", lwd=2)
hist(tD2_height, col="grey", main = "Estimate2 height (mid-range)")
abline(v = theta_height, col="red", lwd=2)

#MASS#
tD1_mass<-tD2_mass<-c()
x2<-mydata1$mass
n2<-50 
theta_mass <- mean(x2)
for (i in 1:1000) {
  samples2 <- sample(x2, size=n2, replace=TRUE)  # repeated sampling from data
  tD1_mass<- c(tD1_mass, mean(samples2))
  tD2_mass <- c(tD2_mass, (min(samples2) + max(samples2)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_mass, col="green", main = "Estimate1 mass (mean)")
abline(v = theta_mass, col="red", lwd=2)
hist(tD2_mass, col="pink", main = "Estimate2 mass (mid-range)")
abline(v = theta_mass, col="red", lwd=2)


#METABOLITE A#
tD1_metabolite_A<-tD2_metabolite_A<-c()
x3<-mydata1$metabolite_A
n3<-100
theta_metabolite_A <- mean(x3, na.rm = TRUE)
for (i in 1:1000) {
  samples3 <- sample(x3, size=n3, replace=TRUE)  # repeated sampling from data
  tD1_metabolite_A <- c(tD1_metabolite_A, mean(samples3, na.rm = TRUE))
  tD2_metabolite_A <- c(tD2_metabolite_A, (min(samples3, na.rm = TRUE) + max(samples3, na.rm = TRUE)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_metabolite_A, col="brown", main = "Estimate1 metabolite A (mean)")
abline(v = theta_metabolite_A, col="red", lwd=2)
hist(tD2_metabolite_A, col="violet", main = "Estimate2 metabolite A (mid-range)")
abline(v = theta_metabolite_A, col="red", lwd=2)

#METABOLITE B#
tD1_metabolite_B<-tD2_metabolite_B<-c()
x4<-mydata1$metabolite_B
n4<-100
theta_metabolite_B <- mean(x4)
for (i in 1:1000) {
  samples4 <- sample(x4, size=n4, replace=TRUE)  # repeated sampling from data
  tD1_metabolite_B <- c(tD1_metabolite_B, mean(samples4))
  tD2_metabolite_B <- c(tD2_metabolite_B, (min(samples4) + max(samples4)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_metabolite_B, col="beige", main = "Estimate1 metabolite B (mean)")
abline(v = theta_metabolite_B, col="red", lwd=2)
hist(tD2_metabolite_B, col="lightblue", main = "Estimate2 metabolite B (mid-range)")
abline(v = theta_metabolite_B, col="red", lwd=2)

#METABOLITE C#
tD1_metabolite_C<-tD2_metabolite_C<-c()
x5<-mydata1$metabolite_C
n5<-100
theta_metabolite_C <- mean(x5)
for (i in 1:1000) {
  samples5 <- sample(x5, size=n5, replace=TRUE)  # repeated sampling from data
  tD1_metabolite_C <- c(tD1_metabolite_C, mean(samples5))
  tD2_metabolite_C <- c(tD2_metabolite_C, (min(samples5) + max(samples5)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_metabolite_C, col="orange", main = "Estimate1 metabolite C (mean)")
abline(v = theta_metabolite_C, col="red", lwd=2)
hist(tD2_metabolite_C, col="yellow", main = "Estimate2 metabolite C (mid-range)")
abline(v = theta_metabolite_C, col="red", lwd=2)

#METABOLITE D#
tD1_metabolite_D<-tD2_metabolite_D<-c()
x6<-mydata1$metabolite_D
n6<-100
theta_metabolite_D <- mean(x6)
for (i in 1:1000) {
  samples6 <- sample(x6, size=n6, replace=TRUE)  # repeated sampling from data
  tD1_metabolite_D <- c(tD1_metabolite_D, mean(samples6))
  tD2_metabolite_D <- c(tD2_metabolite_D, (min(samples6) + max(samples6)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_metabolite_D, col="darkmagenta", main = "Estimate1 metabolite D (mean)")
abline(v = theta_metabolite_D, col="red", lwd=2)
hist(tD2_metabolite_D, col="lightgreen", main = "Estimate2 metabolite D (mid-range)")
abline(v = theta_metabolite_D, col="red", lwd=2)

#METABOLITE E#
tD1_metabolite_E<-tD2_metabolite_E<-c()
x7<-mydata1$metabolite_E
n7<-100
theta_metabolite_E <- mean(x7)
for (i in 1:1000) {
  samples7 <- sample(x7, size=n7, replace=TRUE)  # repeated sampling from data
  tD1_metabolite_E <- c(tD1_metabolite_E, mean(samples7))
  tD2_metabolite_E <- c(tD2_metabolite_E, (min(samples7) + max(samples7)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_metabolite_E, col="coral", main = "Estimate1  metabolite E (mean)")
abline(v = theta_metabolite_E, col="red", lwd=2)
hist(tD2_metabolite_E, col="deeppink", main = "Estimate2 metabolite E(mid-range)")
abline(v = theta_metabolite_E, col="red", lwd=2)

#DISEASE SCORES#
tD1_disease_score<-tD2_disease_score<-c()
x8<-mydata1$disease_score
n8<-100
theta_disease_score <- mean(x8)
for (i in 1:1000) {
  samples8 <- sample(x8, size=n8, replace=TRUE)  # repeated sampling from data
  tD1_disease_score <- c(tD1_disease_score, mean(samples8))
  tD2_disease_score <- c(tD2_disease_score, (min(samples8) + max(samples8)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_disease_score, col="chocolate", main = "Estimate1 disease scores (mean)")
abline(v = theta_disease_score, col="red", lwd=2)
hist(tD2_disease_score, col="cyan", main = "Estimate2 disease scores (mid-range)")
abline(v = theta_disease_score, col="red", lwd=2)



par(mfrow=c(4,2))
hist(tD1_height, col="blue", main = "Estimate1 height (mean): ")
abline(v = theta_height, col="red", lwd=2)
hist(tD2_height, col="grey", main = "Estimate2 height (mid-range)")
abline(v = theta_height, col="red", lwd=2)

hist(tD1_mass, col="green", main = "Estimate1 mass (mean)")
abline(v = theta_mass, col="red", lwd=2)
hist(tD2_mass, col="pink", main = "Estimate2 mass (mid-range)")
abline(v = theta_mass, col="red", lwd=2)

hist(tD1_metabolite_A, col="brown", main = "Estimate1 metabolite A (mean)")
abline(v = theta_metabolite_A, col="red", lwd=2)
hist(tD2_metabolite_A, col="violet", main = "Estimate2 metabolite A (mid-range)")
abline(v = theta_metabolite_A, col="red", lwd=2)

hist(tD1_metabolite_B, col="beige", main = "Estimate1 metabolite B (mean)")
abline(v = theta_metabolite_B, col="red", lwd=2)
hist(tD2_metabolite_B, col="lightblue", main = "Estimate2 metabolite B (mid-range)")
abline(v = theta_metabolite_B, col="red", lwd=2)

hist(tD1_metabolite_C, col="orange", main = "Estimate1 metabolite C (mean)")
abline(v = theta_metabolite_C, col="red", lwd=2)
hist(tD2_metabolite_C, col="yellow", main = "Estimate2 metabolite C (mid-range)")
abline(v = theta_metabolite_C, col="red", lwd=2)

hist(tD1_metabolite_D, col="darkmagenta", main = "Estimate1 metabolite D (mean)")
abline(v = theta_metabolite_D, col="red", lwd=2)
hist(tD2_metabolite_D, col="lightgreen", main = "Estimate2 metabolite D (mid-range)")
abline(v = theta_metabolite_D, col="red", lwd=2)

hist(tD1_metabolite_E, col="coral", main = "Estimate1  metabolite E (mean)")
abline(v = theta_metabolite_E, col="red", lwd=2)
hist(tD2_metabolite_E, col="deeppink", main = "Estimate2 metabolite E(mid-range)")
abline(v = theta_metabolite_E, col="red", lwd=2)

hist(tD1_disease_score, col="chocolate", main = "Estimate1 disease scores (mean)")
abline(v = theta_disease_score, col="red", lwd=2)
hist(tD2_disease_score, col="cyan", main = "Estimate2 disease scores (mid-range)")
abline(v = theta_disease_score, col="red", lwd=2)
##########
#simplified
sample_size<-c(50,50,100,100,100,100,100,100)
colors1<-c("blue","green","brown","beige","orange","darkmagenta","coral","chocolate")
colors2<-c("grey","pink","violet","lightblue","yellow","lightgreen","deeppink","cyan")

for (i in seq_along(continuous_variables)){
  variable_names<-continuous_variables[i]
  variable_vector2<-mydata1[[variable_names]]
  n<-sample_size[i]
  
  if (variable_names=="metabolite_A"){
    variable_vector2<-na.omit(variable_vector2) #OpenAI for removal of NA variables in metabolite_A
  }
  tD1<-tD2<-c() #bootstapping
  theta<-mean(variable_vector2)
  
  for (j in 1:1000) {
    sample_j <- sample(variable_vector2, size=n, replace=TRUE)  # repeated sampling from data
    tD1<-c(tD1, mean(sample_j))
    tD2<-c(tD2, (min(sample_j) + max(sample_j)) / 2)
  }
  par(mfrow=c(4,2))
  hist(tD1, col=colors1[i], main = paste("Estimate1",variable_names, "(mean)"))
  abline(v = theta, col="red", lwd=2)
  hist(tD2, col=colors2[i], main = paste("Estimate2",variable_names, "(mid-range)"))
  abline(v = theta, col="red", lwd=2)
}



##################################
#BIAS AND MEAN SQUARED ERRORS (MSE)
##################################
###BIAS
bias1_height<-mean(tD1_height)-theta_height
bias2_height<-mean(tD2_height)-theta_height

bias1_mass<-mean(tD1_mass)-theta_mass
bias2_mass<-mean(tD2_mass)-theta_mass

bias1_metabolite_A<-mean(tD1_metabolite_A)-theta_metabolite_A
bias2_metabolite_A<-mean(tD2_metabolite_A)-theta_metabolite_A

bias1_metabolite_B<-mean(tD1_metabolite_B)-theta_metabolite_B
bias2_metabolite_B<-mean(tD2_metabolite_B)-theta_metabolite_B

bias1_metabolite_C<-mean(tD1_metabolite_C)-theta_metabolite_C
bias2_metabolite_C<-mean(tD2_metabolite_C)-theta_metabolite_C

bias1_metabolite_D<-mean(tD1_metabolite_D)-theta_metabolite_D
bias2_metabolite_D<-mean(tD2_metabolite_D)-theta_metabolite_D

bias1_metabolite_E<-mean(tD1_metabolite_E)-theta_metabolite_E
bias2_metabolite_E<-mean(tD2_metabolite_E)-theta_metabolite_E

bias1_disease_score<-mean(tD1_disease_score)-theta_disease_score
bias2_disease_score<-mean(tD2_disease_score)-theta_disease_score


###MSE
mse_tD1_height<-mean((tD1_height-theta_height)^2)
mse_tD2_height<-mean((tD2_height-theta_height)^2)

mse_tD1_mass<-mean((tD1_mass-theta_mass)^2)
mse_tD2_mass<-mean((tD2_mass-theta_mass)^2)

mse_tD1_metabolite_A<-mean((tD1_metabolite_A-theta_metabolite_A)^2)
mse_tD2_metabolite_A<-mean((tD2_metabolite_A-theta_metabolite_A)^2)

mse_tD1_metabolite_B<-mean((tD1_metabolite_B-theta_metabolite_B)^2)
mse_tD2_metabolite_B<-mean((tD2_metabolite_B-theta_metabolite_B)^2)

mse_tD1_metabolite_C<-mean((tD1_metabolite_C-theta_metabolite_C)^2)
mse_tD2_metabolite_C<-mean((tD2_metabolite_C-theta_metabolite_C)^2)

mse_tD1_metabolite_D<-mean((tD1_metabolite_D-theta_metabolite_D)^2)
mse_tD2_metabolite_D<-mean((tD2_metabolite_D-theta_metabolite_D)^2)

mse_tD1_metabolite_E<-mean((tD1_metabolite_E-theta_metabolite_E)^2)
mse_tD2_metabolite_E<-mean((tD2_metabolite_E-theta_metabolite_E)^2)

mse_tD1_disease_score<-mean((tD1_disease_score-theta_disease_score)^2)
mse_tD2_disease_score<-mean((tD2_disease_score-theta_disease_score)^2)
#####################################
#simplified
for (x in continuous_variables){
  bias1<-mean(tD1)-theta
  bias2<-mean(tD2)-theta
  mse1<-mean((tD1-theta)^2)
  mse2<-mean((tD2-theta)^2)
  
  print(bias1)
  print(bias2)
  print(mse1)
  print(mse2)
}




#####################################


#BIINARY VARIABLES= sex, smoker, drinker
#considering the fact that all of these variables have 0/1 categorical values:
par(mfrow=c(1,3))
barplot(table(mydata1$sex), col = c("red","blue"), main="Sex data", names.arg = c("Female", "Male")) #female(0), male(1)
barplot(table(mydata1$smoker), col=c("darkgreen","yellow"), main="Smoker data" , names.arg = c("Non-Smoker", "Smoker")) #non-smoker(0), smoker(1)
barplot(table(mydata1$drinker), col = c("aquamarine","brown"), main="Drinker data", names.arg = c("Non-Drinker", "Drinker")) #non-drinker(0), drinker(1)


###################################
#FOLLOWING DAY 4- STATISTICAL TESTING
###################################
###the sample mean is indicated by mu
#One sample t-test
###################################
#Ho height:The sample mean is equal to 179.
#H1 height:The sample mean is not equal to 179.
t.test(mydata1$height, mu=179)

#Ho mass:The sample mean is equal to 74.
#H1 mass:The sample mean is not equal to 74.
t.test(mydata1$mass, mu=74)

#Ho matabolite A: The sample mean is equal to 0.8.
#H1 metabolite A: The sample mean is not equal to 0.8.
t.test(mydata1$metabolite_A, mu=0.8)

#Ho matabolite B: The sample mean is equal to 1.7.
#H1 metabolite B:The sample mean is not equal to 1.7.
t.test(mydata1$metabolite_B,mu=1.7)

#Ho matabolite C:The sample mean is equal to 0.3.
#H1 metabolite C:The sample mean is not equal to 0.3.
t.test(mydata1$metabolite_C, mu=0.3)

#Ho matabolite D: The sample mean is equal to 4.95
#H1 metabolite D: The sample mean is not equal to 4.95
t.test(mydata1$metabolite_D, mu=4.95)

#Ho matabolite E: The sample mean is equal to 10.1
#H1 metabolite E:The sample mean is not equal to 10.1
t.test(mydata1$metabolite_E, mu=10.1)

#Ho disease score: The sample mean is equal to -0.3
#H1 disease score:The sample mean is not equal to -0.3
t.test(mydata1$disease_score, mu=-0.3)


######################################
#Two samples t-test
#####################################
#H0: The binary values do not imply a difference in the mean in the continuous variables. 
#H1: The binary values create differences between the mean of the two groups in the continuous variables. 
#####################################
#HEIGHT
t.test(height~sex, data=mydata1)
t.test(height~smoker, data=mydata1)
t.test(height~drinker, data=mydata1)

#MASS
t.test(mass~sex, data=mydata1)
t.test(mass~smoker, data=mydata1)
t.test(mass~drinker, data=mydata1)

#METABOLITE A
t.test(metabolite_A~sex, data=mydata1)
t.test(metabolite_A~smoker, data=mydata1)
t.test(metabolite_A~drinker, data=mydata1)

#METABOLITE B
t.test(metabolite_B~sex, data=mydata1)
t.test(metabolite_B~smoker, data=mydata1)
t.test(metabolite_B~drinker, data=mydata1)

#METABOLITE C
t.test(metabolite_C~sex, data=mydata1)
t.test(metabolite_C~smoker, data=mydata1)
t.test(metabolite_C~drinker, data=mydata1)

#METABOLITE D
t.test(metabolite_D~sex, data=mydata1)
t.test(metabolite_D~smoker, data=mydata1)
t.test(metabolite_D~drinker, data=mydata1)

#METABOLITE E
t.test(metabolite_E~sex, data=mydata1)
t.test(metabolite_E~smoker, data=mydata1)
t.test(metabolite_E~drinker, data=mydata1)

#DISEASE SCORE
t.test(disease_score~sex, data=mydata1)
t.test(disease_score~smoker, data=mydata1)
t.test(disease_score~drinker, data=mydata1)

######################################
#Chi-squared test- binary variables
######################################
#H0: The binary values are independent.
#H1: The binary values are associated with one another. 
chisq.test(table(mydata1$sex, mydata1$smoker))
chisq.test(table(mydata1$sex, mydata1$drinker))
chisq.test(table(mydata1$smoker, mydata1$drinker))



######################################
#FOLLOWING DAY 5- LINEAR MODELS
######################################
#linear model
model1<-lm(disease_score~height+mass+smoker+sex+drinker+metabolite_A+metabolite_B+metabolite_C+metabolite_D+metabolite_E, data = mydata1)
summary(model1)
par(mfrow=c(2,2))
plot(model1)
anova(model1)

#correlation and visualisation on scatter plots
#checking if the predictors affect the disease score
par(mfrow=c(2,3))
plot1<-plot(mydata1$height, mydata1$disease_score, main = "Height vs Disease Score", xlab = "Height", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])
plot2<-plot(mydata1$mass, mydata1$disease_score, main = "Mass vs Disease Score", xlab = "Mass", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])
plot3<-plot(mydata1$metabolite_A, mydata1$disease_score, main = "Metabolite A vs Disease Score", xlab = "Metabolite A", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])
plot4<-plot(mydata1$metabolite_B, mydata1$disease_score, main = "Metabolite B vs Disease Score", xlab = "Metabolite B", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])
plot5<-plot(mydata1$metabolite_C, mydata1$disease_score, main = "Metabolite C vs Disease Score", xlab = "Metabolite C", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])
plot6<-plot(mydata1$metabolite_D, mydata1$disease_score, main = "Metabolite D vs Disease Score", xlab = "Metabolite D", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])
plot7<-plot(mydata1$metabolite_E, mydata1$disease_score, main = "Metabolite E vs Disease Score", xlab = "Metabolite E", ylab = "Disease Score", col = topo.colors(2)[mydata1$sex + 1])

cor(mydata1$height, mydata1$disease_score)
cor(mydata1$mass, mydata1$disease_score)
cor(mydata1$metabolite_A, mydata1$disease_score, use = "complete.obs") #(OpenAI, 2025) to remove the NA variables
cor(mydata1$metabolite_B, mydata1$disease_score)
cor(mydata1$metabolite_C, mydata1$disease_score)
cor(mydata1$metabolite_D, mydata1$disease_score)
cor(mydata1$metabolite_E, mydata1$disease_score)

  




########################################
#CONSIDERING THE FIRST DATASET (mydata1) IS ASSOCIATED WITH THE DISEASE SCORE OF COLOUR BLINDENESS,
#THE SECOND DATABASE WAS CHOSEN ON THE PRESENCE OF THE VARIABLES OF SEX,HEIGHT, WEIGHT(MASS IN mydataset1), 
#SMOKER BINARY VARIABLE AND EYESIGHT SCORES OF BOTH LEFT AND RIGHT EYE

#THE DATASETS ARE INDEPENDENT FROM ONE ANOTHER, BUT THE EFFECT OF SMOKING ON EYESIGHT AND DISEASE
#ARE CONSIDERED
########################################
mydata2<- read.csv("data/raw/smoking.csv")
View(smoking)
head(mydata) #visualization of data
str(mydata) #structure of the dataset
#######################################
#ALL THE TRUE MEANS ARE CALCULATED FOR THE VARIABLES MENTIONED PREVIOUSLY, AS WELL AS THE ESTIMATED MEANS
calculationforestimate2<-function(a,b) abs(a-b)/2+min(c(a,b))
#HEIGHT#
estimate1_height_smoking<-mean(mydata2$height.cm.) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_height_smoking) #164.65
estimate2_height_smoking<-calculationforestimate2(max(mydata2$height.cm.),min(mydata2$height.cm.))
print(estimate2_height_smoking) #160

#WEIGHT#
estimate1_weight_smoking<-mean(mydata2$weight.kg.) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_weight_smoking) #65.864
estimate2_weight_smoking<-calculationforestimate2(max(mydata2$weight.kg.),min(mydata2$weight.kg.))
print(estimate2_weight_smoking) #82.5

#EYESIGHT LEFT#
estimate1_eyesight_left_smoking<-mean(mydata2$eyesight.left.) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_eyesight_left_smoking) #1.012623
estimate2_eyesight_left_smoking<-calculationforestimate2(max(mydata2$eyesight.left.),min(mydata2$eyesight.left.))
print(estimate2_eyesight_left_smoking) #5

#EYESIGHT RIGHT#
estimate1_eyesight_right_smoking<-mean(mydata2$eyesight.right.) #the first estimate is done based on the central tendency, hence using the mean
print(estimate1_eyesight_right_smoking) #1.007443
estimate2_eyesight_right_smoking<-calculationforestimate2(max(mydata2$eyesight.right.),min(mydata2$eyesight.right.))
print(estimate2_eyesight_right_smoking) #5
#######
#simplified
continuous_variables2<-c("height.cm.","weight.kg.","eysight.left.","eyesight.right.")
for (i in continuous_variables2){
  x<-mydata2[[i]]
  estimate1_2<-mean(x, na.rm=TRUE)
  estimate2_2<-calculationforestimate2(max(x, na.rm=TRUE), min(x, na.rm=TRUE))
  print(estimate1_2) 
  print(estimate2_2)
}



#################################
#Plot visusalisation
#HEIGHT#
tD1_height_smoking<-tD2_height_smoking<-c()
x1<-mydata2$height.cm.
n1<-50
theta_height_smoking <- mean(x1)
for (i in 1:1000) {
  samples1 <- sample(x1, size=n1, replace=TRUE)  # repeated sampling from data
  tD1_height_smoking<- c(tD1_height_smoking, mean(samples1))
  tD2_height_smoking<- c(tD2_height_smoking, (min(samples1) + max(samples1)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_height_smoking, col="blue", main = "Estimate1 height smoking (mean): ")
abline(v = theta_height_smoking, col="red", lwd=2)
hist(tD2_height_smoking, col="grey", main = "Estimate2 height smoking (mid-range)")
abline(v = theta_height_smoking, col="red", lwd=2)

#WEIGHT#
tD1_weight_smoking<-tD2_weight_smoking<-c()
x2<-mydata2$weight.kg.
n2<-50 
theta_weight_smoking <- mean(x2)
for (i in 1:1000) {
  samples2 <- sample(x2, size=n2, replace=TRUE)  # repeated sampling from data
  tD1_weight_smoking<- c(tD1_weight_smoking, mean(samples2))
  tD2_weight_smoking <- c(tD2_weight_smoking, (min(samples2) + max(samples2)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_weight_smoking, col="green", main = "Estimate1 weight smoking (mean)")
abline(v = theta_weight_smoking, col="red", lwd=2)
hist(tD2_weight_smoking, col="pink", main = "Estimate2 weight smoking (mid-range)")
abline(v = theta_weight_smoking, col="red", lwd=2)

#EYESIGHT LEFT#
tD1_eyesight_left_smoking<-tD2_eyesight_left_smoking<-c()
x3<-mydata2$eyesight.left.
n3<-1000
theta_eyesight_left_smoking <- mean(x2)
for (i in 1:1000) {
  samples2 <- sample(x3, size=n3, replace=TRUE)  # repeated sampling from data
  tD1_eyesight_left_smoking<- c(tD1_eyesight_left_smoking, mean(samples2))
  tD2_eyesight_left_smoking <- c(tD2_eyesight_left_smoking, (min(samples2) + max(samples2)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_eyesight_left_smoking, col="brown", main = "Estimate1 eyesight left smoking (mean)")
abline(v = theta_eyesight_left_smoking, col="red", lwd=2)
hist(tD2_eyesight_left_smoking, col="violet", main = "Estimate2 eyesight left smoking (mid-range)")
abline(v = theta_eyesight_left_smoking, col="red", lwd=2)

#EYESIGHT RIGHT#
tD1_eyesight_right_smoking<-tD2_eyesight_right_smoking<-c()
x4<-mydata2$eyesight.right.
n4<-1000
theta_eyesight_right_smoking <- mean(x2)
for (i in 1:1000) {
  samples2 <- sample(x4, size=n4, replace=TRUE)  # repeated sampling from data
  tD1_eyesight_right_smoking<- c(tD1_eyesight_right_smoking, mean(samples2))
  tD2_eyesight_right_smoking <- c(tD2_eyesight_right_smoking, (min(samples2) + max(samples2)) / 2)
}
par(mfrow=c(1,2))
hist(tD1_eyesight_right_smoking, col="beige", main = "Estimate1 eyesight right smoking (mean)")
abline(v = theta_eyesight_right_smoking, col="red", lwd=2)
hist(tD2_eyesight_right_smoking, col="lightblue", main = "Estimate2 eyesight right smoking (mid-range)")
abline(v = theta_eyesight_right_smoking, col="red", lwd=2)
#############################

par(mfrow=c(4,2))
hist(tD1_height_smoking, col="blue", main = "Estimate1 height smoking (mean): ")
abline(v = theta_height_smoking, col="red", lwd=2)
hist(tD2_height_smoking, col="grey", main = "Estimate2 height smoking (mid-range)")
abline(v = theta_height_smoking, col="red", lwd=2)

hist(tD1_weight_smoking, col="green", main = "Estimate1 weight smoking (mean)")
abline(v = theta_weight_smoking, col="red", lwd=2)
hist(tD2_weight_smoking, col="pink", main = "Estimate2 weight smoking (mid-range)")
abline(v = theta_weight_smoking, col="red", lwd=2)

hist(tD1_eyesight_left_smoking, col="brown", main = "Estimate1 eyesight left smoking (mean)")
abline(v = theta_eyesight_left_smoking, col="red", lwd=2)
hist(tD2_eyesight_left_smoking, col="violet", main = "Estimate2 eyesight left smoking (mid-range)")
abline(v = theta_eyesight_left_smoking, col="red", lwd=2)

hist(tD1_eyesight_right_smoking, col="beige", main = "Estimate1 eyesight right smoking (mean)")
abline(v = theta_eyesight_right_smoking, col="red", lwd=2)
hist(tD2_eyesight_right_smoking, col="lightblue", main = "Estimate2 eyesight right smoking (mid-range)")
abline(v = theta_eyesight_right_smoking, col="red", lwd=2)
###########
#simplified 
sample_size<-c(50,50,100,100,100,100,100,100)
colors1<-c("blue","green","brown","beige","orange","darkmagenta","coral","chocolate")
colors2<-c("grey","pink","violet","lightblue","yellow","lightgreen","deeppink","cyan")

for (i in seq_along(continuous_variables2)){
  variable_names<-continuous_variables2[i]
  variable_vector2_2<-mydata2[[variable_names]]
  n<-sample_size[i]
  
  tD1_2<-tD2_2<-c() #bootstapping
  theta_2<-mean(variable_vector2_2)
  
  for (j in 1:1000) {
    sample_j <- sample(variable_vector2_2, size=n, replace=TRUE)  # repeated sampling from data
    tD1_2<-c(tD1_2, mean(sample_j))
    tD2_2<-c(tD2_2, (min(sample_j) + max(sample_j)) / 2)
  }
  par(mfrow=c(4,2))
  hist(tD1_2, col=colors1[i], main = paste("Estimate1",variable_names, "(mean)"))
  abline(v = theta_2, col="red", lwd=2)
  hist(tD2_2, col=colors2[i], main = paste("Estimate2",variable_names, "(mid-range)"))
  abline(v = theta_2, col="red", lwd=2)
}







############################
#T-test
############################
#HYPOTHESES
#HO:  The true mean is equal to the true mean calculated above.
#H1:  The true mean is not equal to the true mean calculated above.
###########
t.test(mydata2$height.cm., mu =164) # p-value<2.2e-16
t.test(mydata2$weight.kg., mu =65) # p-value<2.2e-16
t.test(mydata2$eyesight.left., mu =1.012) #p-value=0.7627
t.test(mydata2$eyesight.right., mu =1.007) #p-value = 0.8298

############################
#Two samples t-test
###########################
#HYPOTHESES
#H0:The true mean is equal to the true mean calculated above.
#H1: The true mean is not equal to the true mean calculated above.
#################
#By gender
t.test(mydata2$height.cm. ~gender, data=mydata2) #p-value<2.2e-16
t.test(mydata2$weight.kg. ~gender, data=mydata2) #p-value<2.2e-16
t.test(mydata2$eyesight.left. ~gender, data=mydata2) #p-value<2.2e-16
t.test(mydata2$eyesight.right. ~gender, data=mydata2) #p-value<2.2e-16

#By smoker status
t.test(mydata2$height.cm. ~smoking, data=mydata2) #p-value<2.2e-16
t.test(mydata2$weight.kg. ~smoking, data=mydata2) #p-value<2.2e-16
t.test(mydata2$eyesight.left. ~smoking, data=mydata2) #p-value<2.2e-16
t.test(mydata2$eyesight.right. ~smoking, data=mydata2) #p-value<2.2e-16

#############################
#Chi-squared test (the two categorical variables)
############################
chisq.test(table(mydata2$gender, mydata2$smoking)) #p-value < 2.2e-16



##########################
#Linear models for interactions
#########################
model_left<- lm(eyesight.left. ~ smoking*gender+height.cm.+weight.kg., data = mydata2)
summary(model_left) #p-value: < 2.2e-16
par(mfrow=c(2,2))
plot(model_left)
anova(model_left)

model_right<- lm(eyesight.right. ~ smoking*gender+height.cm.+weight.kg., data = mydata2)
summary(model_right) #p-value: < 2.2e-16
par(mfrow=c(2,2))
plot(model_right)
anova(model_right)

#######
#data visualisation based on interactions of gender and smoking
######
par(mfrow=c(2,2))
boxplot(eyesight.left.~ smoking+gender, data = mydata2,col = c("lightblue", "pink"), main = "Eyesight left by Smoking Status and Gender")
boxplot(eyesight.right.~ smoking+gender, data = mydata2,col = c("blue", "red"), main = "Eyesight right by Smoking Status and Gender")

ggplot(mydata2, aes(x = interaction(smoking, gender), y = eyesight.left., color = smoking)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~gender) +
  labs(title = "Eyesight vs. Height by Smoking and Gender")

par(mfrow=c(2,2))
plot11<-plot(mydata2$height.cm., mydata2$eyesight.left.  , main = "Height vs Eyesight left", xlab = "Height", ylab = "Eyesight left", col = topo.colors(2)[as.factor(mydata2$gender)])
plot12<-plot(mydata2$height.cm., mydata2$eyesight.right.  , main = "Height vs Eyesight right", xlab = "Height", ylab = "Eyesight right", col = topo.colors(2)[as.factor(mydata2$gender)])
plot13<-plot(mydata2$weight.kg., mydata2$eyesight.left., main = "Weight vs Eyesight left", xlab = "Mass", ylab = "Eyesight left", col = topo.colors(2)[as.factor(mydata2$gender)])
plot14<-plot(mydata2$weight.kg., mydata2$eyesight.right., main = "Weight vs Eyesight right", xlab = "Mass", ylab = "Eyesight right", col = topo.colors(2)[as.factor(mydata2$gender)])
