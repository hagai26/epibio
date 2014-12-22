
#MSet <- mapToGenome(MSet)
#MSet <- preprocessRaw(RGSet)
#beta = getBeta(ratioSet)

N=3
u <- rnorm(N)
x1 <- rnorm(N)
x2 <- rnorm(N)
y <- 1 + x1 + x2 + u
mydat <- data.frame(y,x1,x2)

#list.files(baseDir)
#print(names(targets))

#detP <- detectionP(RGSet)
#failed <- detP>0.01
#colMeans(failed) # Fraction of failed positions per sample
#sum(rowMeans(failed)>0.5) # How many positions failed in >50% of samples?


#head(granges(MSet))
#sex <- getSex(mset)
#plotSex(sex)
#plot(as.matrix(getQC(mset)))

#mean <- relevant_beta
#std <- numeric(length(mean))

#print(sampleNames(ratioSet))
#print("dim(getRed(RGSet)):")
#print(dim(getRed(RGSet)))

#only_level_targets = targets[relevant_targets,]
#cols_num <- nrow(only_level_targets)

#relevant_beta <- beta[,relevant_targets, drop = FALSE]


#study_levels = levels(factor(targets$Study))
#type_levels = levels(factor(targets$Type))
#print("study_levels:")
#print(study_levels)
#print(type_levels)
#all_kinds <- data.frame(kind=character(), number=numeric(), stringsAsFactors=FALSE)
#for(type in type_levels) {
#  for(study in study_levels) {
#    study_no_spaces <- gsub(" ","_", study , fixed=TRUE)
#    type_no_spaces <- gsub(" ","_", type , fixed=TRUE)
#    kind_name <- paste0(study_no_spaces, ".", type_no_spaces)
#    relevant_targets <- subset(targets, Study==study & Type==type)
#    
#    newrow = data.frame(kind=kind_name, number=nrow(relevant_targets))
#    all_kinds = rbind(all_kinds, newrow)
#  }  
#}

get_basename_no_ext <- function(filepath) {
  return(strsplit(basename(filename), "\\.")[[1]][[1]])
}

