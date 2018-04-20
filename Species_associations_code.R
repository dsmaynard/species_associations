##############################################################
# Code for boral co-occurrence analysis and model averaging
#
# Publication: "Species associations overwhelm abiotic conditions to 
# dictate the structure and function of wood-decay fungal communities"
# Ecology, 2018
#
# Authors: Daniel S. Maynard, Kristofer R. Covey, Thomas W. Crowther,
# Noah W. Sokol, Eric W. Morrison, Serita D. Frey, Linda T.A. van Diepen, 
# Mark A. Bradford
#
# Contact: Daniel Maynard, dmaynard@uchicago.edu
#
####################################################################


# read in environmental data
data<-read_csv("environmental_data.csv")

# select environmental variables for boral
environmental.dat<-as.matrix(model.matrix(~-1+w+n+w*n+soil.n+soil.moist+ph+wood.perc.moist+wood.perc.n,data=data))

# cumulative sum scaling to normalize abundances
# Paulson et al 2013 10(12) Nature Methods
# NOTE: this is not needed with the normalized data!
# otutab.raw<-read_csv("otu_table_raw.csv")
# otu.mr<-newMRexperiment(t(otutab))
# p = cumNormStat(otu.mr)
# otu.table<-ceiling(t(MRcounts(otu.mr,norm=TRUE,log=F)))
# 

# or read in cleaned data
out.table<-read_csv("OTU_normalized.csv")

################################ boral analysis ##############################################################

fit.n1<- boral(y = otu.table, X=environmental.dat,num.lv =4,  family = "negative.binomial",row.eff="random",
			mcmc.control = list(n.burnin = 10000, n.iteration = 100000, n.thin = 30, save.model = TRUE,calc.ics = F,
								prior.control = list(hypparams = c(100, 20, 100, 20))))


### get correlations
rcor<-get.residual.cor(use.fit)
ecor<-get.enviro.cor(use.fit)

### get all correlations
resc.full<-resc1.full<-rcor$cor
envc.full<-envc1.full<-ecor$cor

### get significant correlations
envc<-envc1<-envc1.sig<-ecor$sig.cor
resc<-resc1<-resc1.sig<-rcor$sig.cor




#######################################################################
#### MODEL SELECTION #################################################
####################################################################

# dat.o contains the environemntal data for all wood blocks
# dat2 contains the environemntal and biotic data for the sequenced blocks

# subset the data into biotic and abiotic
dat.o<-data[is.na(data$sequence.id),]
dat2<-data[!is.na(data$sequence.id),]

# set the variables to be tested for each model
aform<-"w*n+w+n+soil.n+moist+ph+(1|plot)"
bform<-paste(aform,"+frequency.rank+richness+evenness+e.prop.pos+prop.pos+frequency.rank*evenness+richness*evenness+richness+frequency.rank+prop.pos*w")
bform.r<-paste(aform,"+frequency.rank+e.prop.pos+prop.pos*w+e.prop.pos*w")
bform.e<-paste(aform,"+frequency.rank+e.prop.pos+prop.pos*w+e.prop.pos*w")

detach("package:lmerTest", unload=TRUE)
detach("package:LMERConvenienceFunctions", unload=TRUE)



### these functions all call "model_avg_lmer", appended below, for model selection and averaging

## DENLOSS abiotic
summary(f1<-lmer(paste("perc.denloss.o~",aform),data=dat.o,REML=F,na.action=na.fail))
denloss.a<-model_avg_lmer(f1)


### DENLOSS biotic
summary(f1<-lmer(paste("perc.denloss.o~",bform),data=dat2,REML=F,na.action=na.fail))
denloss.b<-model_avg_lmer(f1)

##### WOOD, abiotic
summary(f1<-lmer(paste("wood.o~",aform),data=dat.o,REML=F,na.action=na.fail))
wood.a<-model_avg_lmer(f1)

#### WOOD, biotic
summary(f1<-lmer(paste("wood.o~",bform),data=dat2,REML=F,na.action=na.fail))
wood.b<-model_avg_lmer(f1)

#### richness, abiotic
summary(f1<-lmer(paste("richness.o~",aform),data=dat2,REML=F,na.action=na.fail))
rich.a<-model_avg_lmer(f1)

#### richness, biotic
summary(f1<-lmer(paste("richness.o~",bform.r),data=dat2,REML=F,na.action=na.fail))
rich.b<-model_avg_lmer(f1)

######## EVEN, abiotic
summary(f1<-lmer(paste("even.o~",aform),data=dat2,REML=F,na.action=na.fail))
even.a<-model_avg_lmer(f1)

#### EVEN biotic
summary(f1<-lmer(paste("even.o~",bform.e),data=dat2,REML=F,na.action=na.fail))
even.b<-model_avg_lmer(f1)

library(lmerTest)

## fit lmer models to compare
summary(denloss.a.lmer<-lmer(add_form(formula(denloss.a$mod.avg)),data=dat.o))
summary(denloss.b.lmer<-lmer(add_form(formula(denloss.b$mod.avg)),data=dat2))
summary(wood.a.lmer<-lmer(add_form(formula(wood.a$mod.avg)),data=dat.o))
summary(wood.b.lmer<-lmer(add_form(formula(wood.b$mod.avg)),data=dat2))
summary(rich.a.lmer<-lmer(add_form(formula(rich.a$mod.avg)),data=dat2))
summary(rich.b.lmer<-lmer(add_form(formula(rich.b$mod.avg)),data=dat2))
summary(even.a.lmer<-lmer(add_form(formula(even.a$mod.avg)),data=dat2))
summary(even.b.lmer<-lmer(add_form(formula(even.b$mod.avg)),data=dat2))


# get model averaged r2
denloss.a$r2.avg
wood.a$r2.avg
rich.a$r2.avg
even.a$r2.avg
denloss.b$r2.avg
wood.b$r2.avg
rich.b$r2.avg
even.b$r2.avg


#############################################################################################
#### this is the main function for model avgeraging. currently set up for parallel processing.
model_avg_lmer<-function(f1,nclust=6,delt.cut=2,fixed.var,mmax=100,cumw=F){

	clu<-makeForkCluster(nclust)

	# dredge
	if(missing(fixed.var)){
		dr<-pdredge(f1,cluster=clu,extra="R^2",m.lim=c(0,mmax))
	}
	else{
		dr<-pdredge(f1,cluster=clu,fixed=fixed.var,extra="R^2",m.lim=c(0,mmax))
	}
	# stop and remove cluster
	stopCluster(clu)
	rm(clu)

	top.mods<-subset(dr, delta < delt.cut)
	
	if(nrow(top.mods)>1){
		# get top models, formulas and summary stats
		if(!cumw){
			top.mods<-subset(dr, delta < delt.cut)
			top.forms<-get.models(dr, delta < delt.cut)
			#model.avg(top.mods,subset=delta<delt.cut)
			my.all<-model.avg(dr)
			my.avg<-model.avg(dr,subset=delta<delt.cut,revised.var=TRUE,fit=TRUE)
		}
		else{
			top.mods<-subset(dr, cumsum(weight) <= .95)
			top.forms<-get.models(dr, cumsum(weight) <= .95)
			#model.avg(top.mods,subset=cumsum(weight) <= .95)
			my.all<-model.avg(dr)
			my.avg<-model.avg(dr,subset=cumsum(weight) <= .95,revised.var=TRUE,fit=TRUE)
		}

		r2.avg<-sum(top.mods$'R^2'*top.mods$weight)/sum(top.mods$weight)
		
		ma1<-summary(my.avg)
		vcv<-vcov(my.avg)
		conf<-confint(my.avg)
		rel<-as.matrix(importance(my.all))
		cma<-as.matrix(ma1$coefmat.subset)
		nea<-as.matrix(ma1$coef.nmod)

		# get fixed and random R2, model average them
		r2.vals<-data.frame(matrix(0,nrow=length(top.forms),ncol=3))
		names(r2.vals)<-c("f","r","weight")
		r2.vals$weight<-as.numeric(top.mods$weight)
		for(i in 1:length(top.forms)){
			r2.vals$f[i]<-r2.mixed(top.forms[[i]])$fR2
			r2.vals$r[i]<-r2.mixed(top.forms[[i]])$rfR2-r2.vals$f[i]
		}
		r2mixed<-data.frame(fixed=sum(r2.vals$f*r2.vals$weight)/sum(r2.vals$weight),random=sum(r2.vals$r*r2.vals$weight)/sum(r2.vals$weight))

		print(showConnections(all=TRUE))

		dr<-NULL
		
		list(top.mods=top.mods,mod.avg=my.avg,summary=ma1,rel.imp=rel,coefmat=cma,neach=nea,r2.mixed=r2mixed,vcv=vcv,conf=conf,r2.avg=r2.avg)
	}
	else{
		
		print("only one model!!")
		
		top.mods<-subset(dr, delta < delt.cut)
		top.forms<-get.models(dr, delta < delt.cut)	
		my.avg<-top.forms[[1]]
		my.all<-my.avg	
		r2.avg<-sum(top.mods$'R^2'*top.mods$weight)/sum(top.mods$weight)

		r2.vals<-r2.mixed(top.forms[[1]])
		r2.vals$weight<-1
		r2mixed<-data.frame(fixed=sum(r2.vals$f*r2.vals$weight)/sum(r2.vals$weight),random=sum(r2.vals$r*r2.vals$weight)/sum(r2.vals$weight))
		
		ma1<-summary(my.avg)
		vcv<-vcov(my.avg)
		conf<-confint(my.avg)
		coefmat<-coef(ma1)
	
		list(top.mods=top.mods,mod.avg=my.avg,summary=ma1,rel.imp=NA,coefmat=coef(ma1),neach=NA,r2.mixed=r2mixed,vcv=vcv,conf=conf,r2.avg=r2.avg)

	}
}






