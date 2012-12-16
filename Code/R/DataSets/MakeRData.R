## This is from the StatLog collection and was downloaded the the UCI MLR.

## IMPORTANT COMMAND: model.matrix, e.g. model.matrix(Class ~. , data=aus);

################################################################################
                            ## Australian Credit ##

## Austrailian - check out .doc file for info ##
aus = read.table("australian.dat", header=FALSE)

## Tell R which columns are categorical.
aus.ctg = c(1,4,5,6,8,9,11,12);
for (i in aus.ctg)
  aus[,i] = factor(aus[,i]);

## The last column is the binary response.
names(aus)[15] = "Class"

## Design matrix and response.
X.aus = as.matrix(model.matrix(Class ~. , data=aus));
y.aus = aus$Class;
  
save(aus, y.aus, X.aus, file="australia.RData");

################################################################################
                              ## German Credit ##

ger = read.table("german-nominal.dat", header=FALSE);

## Tell R which columns are categorical.
ger.ctg = c(1,3,4,6,7,9,10,12,14,15,17,19,20);
for (i in ger.ctg)
  ger[,i] = factor(ger[,i]);

## The last column is the binary response labeled 1=Good, 2=Bad.
ger[,21] = ger[,21]-1;
names(ger)[21] = "Bad"

## Design matrix and response.
X.ger = as.matrix(model.matrix(Bad ~. , data=ger));
y.ger = ger$Bad;

save(ger, y.ger, X.ger, file="germany.RData")

################################################################################
                                ## Heart Data ##

heart = read.table("heart.dat", header=FALSE);

## Tell R which columns are categorical.
heart.ctg = c(11,2,6,9,7,3,13);
for (i in heart.ctg)
  heart[,i] = factor(heart[,i]);

## The last column is the binary response labeled 1=Absence, 2=Presence of heart disease.
heart[,14] = heart[,14]-1;
names(heart)[14] = "heart.disease"

## Design matrix and response.
X.heart = as.matrix(model.matrix(heart.disease ~. , data=heart));
y.heart = heart$heart.disease;

save(heart, y.heart, X.heart, file="heart.RData");

################################################################################
                        ## Pima Indians Diabetes Data ##

diabetes = read.csv("pima-indians-diabetes.dat", header=FALSE);

## Remove missing?
rowSums(diabetes[,2:8]==0)==0 -> not.missing
diabetes = diabetes[not.missing,]

## Predictors are all continuous

## The last column is the binary response labeled 1=Absence, 2=Presence of heart disease.
names(diabetes)[9] = "diabetes"

## Design matrix and response.
X.diabetes = as.matrix(model.matrix(diabetes ~. , data=diabetes));
y.diabetes = diabetes$diabetes;

save(diabetes, y.diabetes, X.diabetes, file="diabetes.RData");

################################################################################
                              ## German Numeric ##

ger.num = read.table("german.dat", header=FALSE);

## The last column is the binary response labeled 1=Good, 2=Bad.
ger.num[,25] = ger.num[,25]-1;
names(ger.num)[25] = "Bad"

## Design matrix and response.
X.ger.num = as.matrix(model.matrix(Bad ~. , data=ger.num));
y.ger.num = ger.num$Bad;

save(ger.num, y.ger.num, X.ger.num, file="german-numeric.RData")
