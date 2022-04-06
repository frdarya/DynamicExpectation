 # Analysis of behavioural data 
  
library(lme4)
library(reshape)
library(effects)
library(scales)
library(phia)
library(emmeans)
library(data.table)
library(ggplot2)
library(ggpubr)
library(car)

##### EXP 1 #####
dataset <- read.csv('Exp1.csv')
dataset$subj <- factor(dataset$subj)
dataset$targ = factor(dataset$targ,
                      levels = c('0','1'),
                      labels = c('miss','hit'))

dataset$targ_cond = factor(dataset$targ_cond,
                           levels = c('1', '2'),
                           labels = c('Expected', 'Unexpected'))
dataset$f1 = factor(dataset$f1,
                    levels = c('0', '1'),
                    labels = c('fa', 'cr'))
dataset$f1_cond = factor(dataset$f1_cond,
                         levels = c('1', '2'),
                         labels = c('Expected', 'Unexpected'))
dataset$f2 = factor(dataset$f2,
                    levels = c('0', '1'),
                    labels = c('fa', 'cr'))
dataset$f2_cond = factor(dataset$f2_cond,
                         levels = c('1', '2'),
                         labels = c('Expected', 'Unexpected'))
dataset$f3 = factor(dataset$f3,
                    levels = c('0', '1'),
                    labels = c('fa', 'cr'))
dataset$f3_cond = factor(dataset$f3_cond,
                         levels = c('1', '2'),
                         labels = c('Expected', 'Unexpected'))
dataset$targ_ord <- rescale(dataset$targ_ord)
dataset$f1_ord <- rescale(dataset$f1_ord)
dataset$f2_ord <- rescale(dataset$f2_ord)
dataset$f3_ord <- rescale(dataset$f3_ord)
dataset$f1_targ <-rescale(dataset$f1_targ)
dataset$f2_targ <-rescale(dataset$f2_targ)
dataset$f3_targ <-rescale(dataset$f3_targ)
dataset$f1_f3 <-rescale(dataset$f1_f3)
dataset$f2_f3 <-rescale(dataset$f2_f3)
dataset$f1_f2 <-rescale(dataset$f1_f2)
##### predicting hits (Target)####

tonly <- glmer(targ~targ_cond+(1|subj),data=dataset,family=binomial) # more hits for expected

t_f1 <- subset(dataset, c(targ_pos > f1_pos & f1 == 'cr'))
t_int <- glmer(targ ~ targ_cond*f1_cond+targ_ord+(1|subj), data = t_f1, family = binomial)
summary(t_int)
Anova(t_int, type = 3)
int<-emmeans(t_int,~ targ_cond*f1_cond)
pairs(int,simple='each')
plot(allEffects(t_int))


t_f2 <- subset(dataset, c(targ_pos > f2_pos, f2=='cr'))
t2_int <- glmer(targ ~ targ_cond*f2_cond+targ_ord+(1|subj), data = t_f2, family = binomial)
summary(t2_int)
Anova(t2_int, type = 3)
int2<-emmeans(t2_int,~ targ_cond*f2_cond)
pairs(int2,simple='each')
plot(allEffects(t2_int))

t_f3 <- subset(dataset, c(targ_pos > f3_pos, f3=='cr'))
t3_int <- glmer(targ ~ targ_cond*f3_cond+targ_ord+(1|subj), data = t_f3, family = binomial)
summary(t3_int)
Anova(t3_int, type = 3)

t_last <- subset(dataset, targ_pos == 4)
tf1_last <- glmer(targ ~ targ_cond*f1_cond+(1|subj), data = t_last, family = binomial)
Anova(tf1_last, type =3)
plot(allEffects(tf1_last))

targ_Fig1 <- data.frame(effect("targ_cond*f1_cond", t_int))
ts <-ggplot(targ_Fig1)+aes(f1_cond, fit, color = targ_cond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous F1 expectation") + ylab("p(hits)") +
  scale_color_manual(values = c("Unexpected" = "#F8766D", "Expected" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14))+
  labs(color='Current target')
plot(ts)

targLast_SuppFig <- data.frame(effect("targ_cond*f1_cond", tf1_last))
tf1ls <-ggplot(targLast_SuppFig)+aes(f1_cond, fit, shape = targ_cond, color = targ_cond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous F1 expectation") + ylab("p(hits)") +
  scale_color_manual(values = c("Unexpected" = "#F8766D", "Expected" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18))
plot(tf1ls)

##### predicting CR1 (F1) #####
f1only <- glmer(f1~f1_cond+(1|subj),data=dataset,family=binomial) # ns

f1_t <- subset(dataset, c(targ_pos < f1_pos & targ=='hit'))
f1_int <- glmer(f1 ~ targ_cond*f1_cond+f1_ord+ (1|subj), data = f1_t, family = binomial)
summary(f1_int)
Anova(f1_int, type = 3)
fint<-emmeans(f1_int,~ targ_cond*f1_cond)
pairs(fint,simple='each')
plot(allEffects(f1_int))

f1_f2 <- subset(dataset, c(f1_pos > f2_pos & f2=='cr'))
f1_2 <- glmer(f1 ~ f2_cond*f1_cond+f1_ord+ (1|subj), data = f1_f2, family = binomial)
summary(f1_2)
Anova(f1_2, type = 3)
plot(allEffects(f1_2))

f1_f3 <- subset(dataset, c(f1_pos > f3_pos & f3=='cr'))
f1_3 <- glmer(f1 ~ f3_cond*f1_cond+f1_ord+ (1|subj), data = f1_f3, family = binomial)
summary(f1_3)
Anova(f1_3, type = 3)
plot(allEffects(f1_3))


f1_last <- subset(dataset, f1_pos == 4)
f1t_last <- glmer(f1 ~ f1_cond*targ_cond+(1|subj), data = f1_last, family = binomial)
Anova(f1t_last, type = 3)
plot(allEffects(f1t_last))

f1_Fig1 <- data.frame(effect("targ_cond*f1_cond", f1_int))
f1s <-ggplot(f1_Fig1)+aes(targ_cond, fit,  color = f1_cond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous target expectation") + ylab("p(CR1)") +
  scale_color_manual(values = c("Unexpected" = "#F8766D", "Expected" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14))+
  labs(color='Current F1')
plot(f1s)

f1tlast_SuppFig <- data.frame(effect("f1_cond*targ_cond", f1t_last))
f1tls <-ggplot(f1tlast_SuppFig)+aes(targ_cond, fit, shape = f1_cond, color = f1_cond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous target expectation") + ylab("p(CR1)") +
  scale_color_manual(values = c("Unexpected" = "#F8766D", "Expected" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18))
plot(f1tls)

##### predicting CR2 & CR3 (F2 and F3) ####
# predicting f2
f2only <- glmer(f2~f2_cond+(1|subj),data=dataset,family=binomial) # more UCR2

f2_t <- subset(dataset, c(targ_pos < f2_pos & targ=='hit'))
f2_int <- glmer(f2 ~ targ_cond*f2_cond+f2_ord+ (1|subj), data = f2_t, family = binomial)
summary(f2_int)
Anova(f2_int, type = 3)
plot(allEffects(f2_int))

f2_f1 <- subset(dataset, c(f2_pos > f1_pos & f1=='cr'))
f2_1int <- glmer(f2 ~ f1_cond*f2_cond+f2_ord+ (1|subj), data = f2_f1, family = binomial)
summary(f2_1int)
Anova(f2_1int, type = 3)
plot(allEffects(f2_1int))

f2_f3 <- subset(dataset, c(f2_pos > f3_pos & f3=='cr'))
f2_3int <- glmer(f2 ~ f3_cond*f2_cond+f2_ord+ (1|subj), data = f2_f3, family = binomial)
summary(f2_3int)
Anova(f2_3int, type = 3)
plot(allEffects(f2_3int))


# predicting F3#
f3only <- glmer(f3~f3_cond+(1|subj),data=dataset,family=binomial) # ns
f3_t <- subset(dataset, c(targ_pos < f3_pos & targ=='hit'))
f3_int <- glmer(f3 ~ targ_cond*f3_cond+f3_ord+ (1|subj), data = f3_t, family = binomial)
summary(f3_int)
Anova(f3_int, type = 3)
plot(allEffects(f3_int))

f3_f1 <- subset(dataset, c(f3_pos > f1_pos & f1=='cr'))
f3_1int <- glmer(f3 ~ f1_cond*f3_cond+f3_ord+ (1|subj), data = f3_f1, family = binomial)
summary(f3_1int)
Anova(f3_1int, type = 3)
plot(allEffects(f3_1int))

f3_f2 <- subset(dataset, c(f3_pos > f2_pos & f2=='cr'))
f3_2int <- glmer(f3 ~ f2_cond*f3_cond+f3_ord+ (1|subj), data = f3_f2, family = binomial)
summary(f3_2int)
Anova(f3_2int, type = 3)


##### combined T & F1 `#####
Tar = subset(dataset, targ_pos > f1_pos)
Tar$targ = factor(Tar$targ,
                  levels = c('miss', 'hit'),
                  labels = c('Inc', 'Cor'))
Tar$f1 = factor(Tar$f1,
                levels = c('fa', 'cr'),
                labels = c('Inc', 'Cor'))

F1 = subset(dataset, f1_pos > targ_pos)
F1$f1 = factor(F1$f1,
               levels = c('fa', 'cr'),
               labels = c('Inc', 'Cor'))
F1$targ = factor(F1$targ,
                 levels = c('miss', 'hit'),
                 labels = c('Inc', 'Cor'))


keeps <- c("subj", "targ","targ_cond","targ_ord","f1","f1_cond","f1_ord","f1_targ")
Tar=Tar[keeps]
F1=F1[keeps]
setnames(Tar, old=c("targ","targ_cond","targ_ord"), new=c("currResp", "currCond","currOrd"))
setnames(Tar, old=c("f1","f1_cond","f1_ord"), new=c("prevResp", "prevCond","prevOrd"))
setnames(F1, old=c("f1","f1_cond","f1_ord"), new=c("currResp", "currCond","currOrd"))
setnames(F1, old=c("targ","targ_cond","targ_ord"), new=c("prevResp", "prevCond","prevOrd"))

All = rbind(Tar,F1)
All.red <- subset(All, prevResp == "Cor")
All.1 <- glmer(currResp ~ currCond*prevCond+(1|subj), data = All.red, family = binomial)
summary(All.1)
Anova(All.1, type=3)
intAll<-emmeans(All.1,~ currCond*prevCond)
pairs(intAll,simple='each')
plot(allEffects(All.1))

concat_Fig1 <- data.frame(effect("currCond*prevCond", All.1))
concs <-ggplot(concat_Fig1)+aes(prevCond, fit, color = currCond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous  expectation") + ylab("p(correct)") +
  scale_color_manual(values = c("Unexpected" = "#F8766D", "Expected" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=14))+
  labs(color='Current expectation')
plot(concs)

ggarrange(ts, f1s, concs, ncol = 3, legend = 'bottom')

All.1 <- glmer(currResp ~ currCond*prevCond+currOrd+(1|subj), data = All.red, family = binomial)
summary(All.1)
Anova(All.1, type=3)

All.sup <-glmer(currResp ~ currCond*prevCond+(1|subj), data = All, family = binomial)
summary(All.sup)
Anova(All.sup, type=3)
intAllsup<-emmeans(All.sup,~ currCond*prevCond)
pairs(intAllsup,simple='each')
plot(allEffects(All.sup))

Allsup_Fig1 <- data.frame(effect("currCond*prevCond", All.sup))
concsup <-ggplot(Allsup_Fig1)+aes(prevCond, fit, shape = currCond, color = currCond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous  expectation") + ylab("p(correct)") +
  scale_y_continuous(limits=(c(0.5,0.8)))+
  scale_color_manual(values = c("Unexpected" = "#F8766D", "Expected" = "#00BFC4"))+
  theme_classic() + theme(text = element_text(size=18))
plot(concsup)

##### standalone effect of expectatoin (first items) #####
t_first <- subset(dataset, targ_pos == 1)
f1_first <- subset(dataset, f1_pos == 1)
f2_first <- subset(dataset, f2_pos == 1)
f3_first <- subset(dataset, f3_pos == 1)

t_first$resp = t_first$targ
t_first$resp_cond = t_first$targ_cond
f1_first$resp = f1_first$f1
f1_first$resp_cond = f1_first$f1_cond
f2_first$resp = f2_first$f2
f2_first$resp_cond = f2_first$f2_cond
f3_first$resp = f3_first$f3
f3_first$resp_cond = f3_first$f3_cond


af1 <- nrow(f1_first[1])
af2 <- nrow(f2_first[1])
af3 <- nrow(f3_first[1])
at <- nrow(t_first[1])


#### add subj
firsts<- data.frame("subj" = cbind(c(t_first$subj,f1_first$subj,f2_first$subj,f3_first$subj)),
                          "item" = cbind(c(replicate(at,"t"),replicate(af1,"f1"),replicate(af2,"f2"),replicate(af3,"f3"))),
                          "resp"= cbind(c(t_first$resp,f1_first$resp,f2_first$resp,f3_first$resp)),
                          "resp_cond" = cbind(c(t_first$resp_cond,f1_first$resp_cond,f2_first$resp_cond,f3_first$resp_cond)))

firsts$resp = factor(firsts$resp,
                          levels = c('1', '2'),
                          labels = c('0', '1'))


firsts$resp_cond = factor(firsts$resp_cond,
                               levels = c('1', '2'),
                               labels = c('expected', 'unexpected'))

first_mod <- glmer(resp~item*resp_cond+(1|subj), data = firsts, family = binomial)
first_target <- subset(firsts,item == 't')
first_foils <- subset(firsts, item != 't')

first_target_mod <- glmer(resp~resp_cond+(1|subj), data = first_target, family = binomial)
summary(first_target_mod)
Anova(first_target_mod, type = 2)
plot(allEffects(first_target_mod))

first_foil_mod <- glmer(resp~item*resp_cond+(1|subj), data = first_foils, family = binomial)
summary(first_foil_mod)
Anova(first_foil_mod, type = 3)
plot(allEffects(first_foil_mod))

intfirsts<-emmeans(first_mod,~ item*resp_cond)
pairs(intfirsts,simple='each',adjust="FDR")

f1sts <- data.frame(effect("item*resp_cond", first_mod))
firstFigs <-ggplot(f1sts,aes(item, fit, fill = resp_cond, shape = resp_cond, color = resp_cond, ymin = lower, ymax = upper))+
  geom_pointrange(position = position_dodge(width = .1)) + xlab("Foil") + ylab("p(correct)") +
  theme_classic() +theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 18))
plot(firstFigs)

firsts_eff <-ggplot(f1sts, aes(x=resp_cond,y=fit,color=item, group=item)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.1)+ geom_point(size=5) + 
  labs(x = 'Condition',y='% Correct') +
  scale_y_continuous(name="p(correct)", limits=c(0.3, 1))+
  theme_classic() +theme(text = element_text(size = 22, face = "bold"))+ 
  scale_color_manual(values = c("t" = "maroon1", "f1" = "palegreen1","f2" = "seashell4","f3" = "sienna2"),
                     name="Set Event",
                     breaks=c("f1", "f2", "f3","t"),
                     labels=c("F1","F2","F3", "Target"))


plot(firsts_eff)



##### EXP 2  #####
dataset2 <- read.csv('Exp2.csv')
dataset2$subj = factor(dataset2$subj)
dataset2$targ = factor(dataset2$targ,
                       levels = c('0', '1'),
                       labels = c('miss', 'hit'))
dataset2$f1 = factor(dataset2$f1,
                     levels = c('0', '1'),
                     labels = c('fa', 'cr'))
dataset2$f2 = factor(dataset2$f2,
                     levels = c('0', '1'),
                     labels = c('fa', 'cr'))
dataset2$targ_ord <- rescale(abs(dataset2$targ_ord))
dataset2$f1_ord <- rescale(abs(dataset2$f1_ord))
dataset2$f2_ord <- rescale(abs(dataset2$f2_ord))
dataset2$targ_f1 <-rescale(abs(dataset2$targ_f1))
dataset2$targ_f2 <-rescale(abs(dataset2$targ_f2))
dataset2$f1_f2 <-rescale(abs(dataset2$f1_f2))

tonly2 <- glmer(targ~targ_cond+(1|subj),data=dataset2,family=binomial) 
f1only2 <- glmer(f1~f1_cond+(1|subj),data=dataset2,family=binomial) 
f2only2 <-glmer(f2~f2_cond+(1|subj),data=dataset2,family=binomial) 

tbf1_2 <- subset(dataset2, targ_pos < f1_pos)
tbf1_2_c <- subset(tbf1_2, targ == "hit")
f1bt_2 <- subset(dataset2, f1_pos < targ_pos)
f1bt_2_c <- subset(f1bt_2, f1 == "cr")

Tar2 = subset(dataset2, targ_pos > f1_pos)
Tar2$targ = factor(Tar2$targ,
                   levels = c('miss', 'hit'),
                   labels = c('Inc', 'Cor'))
Tar2$f1 = factor(Tar2$f1,
                 levels = c('fa', 'cr'),
                 labels = c('Inc', 'Cor'))

F12 = subset(dataset2, f1_pos > targ_pos)
F12$f1 = factor(F12$f1,
                levels = c('fa', 'cr'),
                labels = c('Inc', 'Cor'))
F12$targ = factor(F12$targ,
                  levels = c('miss', 'hit'),
                  labels = c('Inc', 'Cor'))

##### T & F1 ####
keeps <- c("subj", "targ","targ_cond","targ_ord","f1","f1_cond","f1_ord","targ_f1","targ_pos","f1_pos")
Tar2=Tar2[keeps]
F12=F12[keeps]
setnames(Tar2, old=c("targ","targ_cond","targ_ord","targ_pos"), new=c("currResp", "currCond","currOrd",'currPos'))
setnames(Tar2, old=c("f1","f1_cond","f1_ord","f1_pos"), new=c("prevResp", "prevCond","prevOrd","prevPos"))
setnames(F12, old=c("f1","f1_cond","f1_ord","f1_pos"), new=c("currResp", "currCond","currOrd","currPos"))
setnames(F12, old=c("targ","targ_cond","targ_ord","targ_pos"), new=c("prevResp", "prevCond","prevOrd","prevPos"))

All2 = rbind(Tar2,F12)
All2.red <- subset(All2,prevResp == "Cor")
All2.1 <- glmer(currResp ~ currCond*prevCond+(1|subj), data = All2.red, family = binomial)
summary(All2.1)
Anova(All2.1, type = 3)
plot(allEffects(All2.1))
intAll2<-emmeans(All2.1,~ currCond*prevCond)
pairs(intAll2,simple='each')

All2.sup <-glmer(currResp ~ currCond*prevCond+(1|subj), data = All2, family = binomial)
summary(All2.sup)
Anova(All2.sup, type=3)
intAll2sup<-emmeans(All2.sup,~ currCond*prevCond)
pairs(intAll2sup,simple='each')
plot(allEffects(All2.sup))

concat_Fig2 <- data.frame(effect("currCond*prevCond", All2.1))
concs2 <-ggplot(concat_Fig2)+aes(prevCond, fit, shape = currCond, color = currCond,  ymin = lower, ymax = upper)+ 
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous  expectation") + ylab("p(correct)") +
  scale_color_manual(values = c("'Unexpected'" = "#F8766D", "'Expected'" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18))
plot(concs2)

All2sup_Fig <- data.frame(effect("currCond*prevCond", All2.sup))
concsup2 <-ggplot(All2sup_Fig)+aes(prevCond, fit, shape = currCond, color = currCond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous  expectation") + ylab("p(correct)") +
  scale_y_continuous(limits=(c(0.5,0.8)))+
  scale_color_manual(values = c("'Unexpected'" = "#F8766D", "'Expected'" = "#00BFC4"))+
  theme_classic() + theme(text = element_text(size=18))
plot(concsup2)

##### F2 #######
f2_t1 <- subset(dataset2, c(targ_pos < f2_pos & targ=='hit'))
f2_t_f <- glmer(f2~f2_cond*targ_cond+f2_ord+(1|subj), data = f2_t1, family = binomial)
summary(f2_t_f)
Anova(f2_t_f, type = 3)
plot(allEffects(f2_t_f))

f2_f1_1 <- subset(dataset2, c(f1_pos < f2_pos & f1=='cr'))
f2_f1_f <- glmer(f2~f2_cond*f1_cond+f2_ord+(1|subj), data = f2_f1_1, family = binomial)
summary(f2_f1_f)
Anova(f2_f1_f, type = 3)

##### Target #####
tf12 <- glmer(targ~f1_cond*targ_cond+targ_ord+(1|subj), data = f1bt_2_c, family = binomial)
summary(tf12)
Anova(tf12, type =3)
inttf12sup<-emmeans(tf12,~ f1_cond*targ_cond)
pairs(inttf12sup,simple='each')

tf1_sup2Fig <- data.frame(effect("f1_cond*targ_cond", tf12))
tf1s <-ggplot(tf1_sup2Fig)+aes(f1_cond, fit, shape = targ_cond, color = targ_cond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous F1 expectation") + ylab("p(hit)") +
  scale_color_manual(values = c("'Unexpected'" = "#F8766D", "'Expected'" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18))
plot(tf1s)

f2bt_2 <- subset(dataset2, targ_pos > f2_pos & f2=='cr')
tf22 <- glmer(targ~f2_cond*targ_cond+targ_ord+(1|subj), data = f2bt_2, family = binomial)
Anova(tf22,type=3)
inttf22sup<-emmeans(tf22,~ f2_cond*targ_cond)
pairs(inttf22sup,simple='each')
plot(allEffects(tf22))

##### F1 ######
f1t2 <- glmer(f1~targ_cond*f1_cond+f1_ord+(1|subj), data = tbf1_2_c, family = binomial)
summary(f1t2)
Anova(f1t2, type = 3)

f1t_sup2Fig <- data.frame(effect("targ_cond*f1_cond", f1t2))
f1ts <-ggplot(f1t_sup2Fig)+aes(targ_cond, fit, shape = f1_cond, color = f1_cond, ymin = lower, ymax = upper)+
  geom_pointrange(size=2,position = position_dodge(width = .25)) + xlab("Previous target expectation") + ylab("p(CR1)") +
  scale_color_manual(values = c("'Unexpected'" = "#F8766D", "'Expected'" = "#00BFC4"))+
  theme_classic() + theme(plot.title = element_text(hjust = 0.5), text = element_text(size=18))
plot(f1ts)

f2bf2_2 <- subset(dataset2, f1_pos > f2_pos & f2=='cr')
f1f22 <- glmer(f1~f1_cond*f2_cond+f1_ord+(1|subj), data = f2bf2_2, family = binomial)
Anova(f1f22,type=3)
plot(allEffects(tf22))

##### standalone effect of expectatoin (first items) #####
t_first2 <- subset(dataset2, targ_pos == 1)
f1_first2 <- subset(dataset2, f1_pos == 1)
f2_first2 <- subset(dataset2, f2_pos == 1)

t_first2$resp = t_first2$targ
t_first2$resp_cond = t_first2$targ_cond
f1_first2$resp = f1_first2$f1
f1_first2$resp_cond = f1_first2$f1_cond
f2_first2$resp = f2_first2$f2
f2_first2$resp_cond = f2_first2$f2_cond

a2f1 <- nrow(f1_first2[1])
a2f2 <- nrow(f2_first2[1])
a2t <- nrow(t_first2[1])


firsts2<- data.frame("subj" = cbind(c(t_first2$subj,f1_first2$subj,f2_first2$subj)),
                    "item" = cbind(c(replicate(a2t,"t"),replicate(a2f1,"f1"),replicate(a2f2,"f2"))),
                    "resp"= cbind(c(t_first2$resp,f1_first2$resp,f2_first2$resp)),
                    "resp_cond" = cbind(c(t_first2$resp_cond,f1_first2$resp_cond,f2_first2$resp_cond)))

firsts2$resp = factor(firsts2$resp,
                     levels = c('1', '2'),
                     labels = c('0', '1'))


firsts2$resp_cond = factor(firsts2$resp_cond,
                          levels = c('1', '2'),
                          labels = c('expected', 'unexpected'))

first2_target <- subset(firsts2,item == 't')
first2_foils <- subset(firsts2, item != 't')

first2_target_mod <- glmer(resp~resp_cond+(1|subj), data = first2_target, family = binomial)
summary(first2_target_mod)
Anova(first2_target_mod, type = 2)
plot(allEffects(first2_target_mod))

first2_foil_mod <- glmer(resp~item*resp_cond+(1|subj), data = first2_foils, family = binomial)
summary(first2_foil_mod)
Anova(first2_foil_mod, type = 3)
plot(allEffects(first2_foil_mod))

first2_mod <- glmer(resp~item*resp_cond+(1|subj), data = firsts2, family = binomial)
intfirsts2<-emmeans(first2_mod,~ item*resp_cond)
pairs(intfirsts2,simple='each',adjust="FDR")

f1sts2 <- data.frame(effect("item*resp_cond", first2_mod))
first2Figs <-ggplot(f1sts2,aes(item, fit, fill = resp_cond, shape = resp_cond, color = resp_cond, ymin = lower, ymax = upper))+
  geom_pointrange(position = position_dodge(width = .1)) + xlab("Foil") + ylab("p(correct)") +
  theme_classic() +theme(plot.title = element_text(hjust = 0.5),text = element_text(size = 18))
plot(first2Figs)

firsts2_eff <- ggplot(f1sts2, aes(x=resp_cond,y=fit,color=item, group=item)) +
  geom_line() +
  geom_errorbar(aes(ymin=lower,ymax=upper),width=.1)+ geom_point(size=5) + 
  labs(x = 'Condition',y='% Correct') +
  scale_y_continuous(name="p(correct)", limits=c(0.3, 1))+
  theme_classic() +theme(text = element_text(size = 22, face = "bold"))+  
  scale_color_manual(values = c("t" = "maroon1", "f1" = "palegreen1","f2" = "seashell4"))+
  
plot(firsts2_eff)



##### rule  learning - Exp1 and Exp2 ####
rulexp1 <- read.csv('rul1_plot.csv')

rulacc1 <- ggplot(rulexp1, aes(x=Half, y=Accuracy, fill = Half)) + 
  geom_boxplot()+
  geom_dotplot(binaxis = "y",stackdir = "center",binwidth=0.008,position="dodge")+
  ylab("p(correct)")+
  theme_classic()+
  theme(text = element_text(size=18), legend.position="none")+
  scale_fill_brewer(palette="Dark2")
rulacc1  

rulrt1 <- ggplot(rulexp1, aes(x=Half, y=RT, fill = Half)) + 
  geom_boxplot()+
  geom_dotplot(binaxis = "y",stackdir = "center",binwidth=6,position="dodge")+
  ylab("RT (ms)")+
  theme_classic()+ 
  theme(text = element_text(size=18), legend.position = "none")+
  scale_fill_brewer(palette="Dark2")
rulrt1   


rulexp2 <- read.csv('rul2_plot.csv')

rulacc2 <- ggplot(rulexp2, aes(x=Half, y=Accuracy, fill = Half)) + 
  geom_boxplot()+
  geom_dotplot(binaxis = "y",stackdir = "center",binwidth=0.008,position="dodge")+
  ylab("p(correct)")+
  theme_classic()+
  theme(text = element_text(size=18), legend.position="none")+
  scale_fill_brewer(palette="Dark2")
rulacc2  

rulrt2 <- ggplot(rulexp2, aes(x=Half, y=RT, fill = Half)) + 
  geom_boxplot()+
  geom_dotplot(binaxis = "y",stackdir = "center",binwidth=15,position="dodge")+
  ylab("RT (ms)")+
  theme_classic()+ 
  theme(text = element_text(size=18), legend.position = "none")+
  scale_fill_brewer(palette="Dark2")
rulrt2   
