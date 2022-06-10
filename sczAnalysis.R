library(dplyr) #调用select
library(psych) #相关系数检验
library(pacman) #缺失值处理
library(mice)
library("VIM")
library(pheatmap)
library(car)

#####读入数据
# setwd("d:/fx")
patient = read.table(file = "patient_data.csv", header = T , 
                     sep = "," , fill = TRUE , encoding = "UTF-8")
illdata = read.table(file = "ill.csv", header = T , 
                     sep = "," , fill = TRUE , encoding = "UTF-8")
n = nrow(patient)
mydata = select(patient,1,3,27,52:55,56:59,63:68,69:79)

names = c("sex","age","MMSE","PANSS.P","PANSS.N","PANSS.G","PANSS.S",
          "WHOQOL.生","WHOQOL.心","WHOQOL.社","WHOQOL.环",
          "WHODASD1","WHODASD2","WHODASD3","WHODASD4","WHODASD5","WHODASD6",
          "BACS语言记忆尝试1","BACS语言记忆尝试2",
          "BACS语言记忆尝试3","BACS语言记忆尝试4","BACS语言记忆尝试5",
          "BACS数字序列","BACS代币动作测试","BACS语义流畅度动物",
          "BACS语义流畅度国","BACS符号替代测试","BACS伦敦塔")
colnames(mydata) = names
mydata$BMI = patient$身高.cm.*10000/(patient$体重.kg.*patient$体重.kg.)


#####按照crp分组
mydata = mydata[1:45,]
mydata$meancrp = 0
for (i in 1:14){
  illdata[,i*9] = as.numeric(illdata[,i*9])
}
mydata$s = 0
for (i in 1:45){
  for (j in 1:14){
    if ((substr(illdata[i,(j-1)*9+1],1,4)=="2021")&!is.na(illdata[i,(j-1)*9+8])){
      if(length(grep("C反应蛋白",colnames(illdata)[(j-1)*9+8]))==0) {
        print(paste('Error: ', colnames(illdata)[(j-1)*9+8]))
      }
      if (illdata[i,(j-1)*9+8]==0){
        mydata$meancrp[i] = mydata$meancrp[i]+illdata[i,(j-1)*9+9]
      }
        else{
          mydata$meancrp[i] = mydata$meancrp[i]+illdata[i,(j-1)*9+8]
        }
      mydata$s[i] = mydata$s[i]+1
    }
  }
  mydata$meancrp[i] = mydata$meancrp[i]/mydata$s[i]
}

mydata$meancrp

log(mydata$meancrp+1)
library(ggplot2)
ggplot(mydata, aes(x = meancrp)) +
  geom_histogram(binwidth = 5, fill = "lightblue", colour = "black")

p = ggplot(mydata, aes(x=1,y=log(meancrp+1)))
   geom_violin()
   
ggplot(mydata,aes(x=1,y=log(meancrp+1)))+
     geom_violin()+
     theme(panel.grid = element_blank(),
           panel.background = element_blank(),
           axis.line = element_line(),
           axis.text.x = element_text(angle=90,hjust = 1,vjust=0.5),
           plot.title = element_text(hjust=0.5))+
     labs(x=NULL,
          y=NULL,
          title = "Complexity")+
     geom_jitter(size=0.5)

mydata$logcrp <- log(mydata$meancrp+1)
#####回归分析
results1 = matrix(0,nrow=25,ncol=4,
                   dimnames=list(colnames(mydata)[4:28],
                                 c("tstat","pvalue", "mmsep", "bmip")))
#summary(lm(PANSS.S~logcrp+age+sex+BMI,mydata))
for (i in colnames(mydata)[4:28]){
  formula = paste0(i,"~logcrp+age+sex+BMI+MMSE")
  fit1 = summary(lm(formula,mydata))
  results1[i,"tstat"] = fit1$coefficients["logcrp","t value"]
  results1[i,"pvalue"]  = fit1$coefficients["logcrp","Pr(>|t|)"]
  results1[i,"mmsep"]  = fit1$coefficients["MMSE","Pr(>|t|)"]
  results1[i,"bmip"]  = fit1$coefficients["BMI","Pr(>|t|)"]
}
print(results1,digits=4)


p1 <- ggplot(mydata, aes(x=logcrp,y=PANSS.P)) +
     geom_point() +
     geom_smooth(method=lm)
ggsave('PANSS.pdf',p1)

