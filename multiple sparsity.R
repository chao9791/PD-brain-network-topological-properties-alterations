

data.path <- "F:\\meta\\PD-GTA-Meta\\Data extract\\PD\\multiple sparsity handle by fixed effect model.xlsx"
d <- read_excel(path = data.path)


d.cp <- d %>% dplyr::filter(Cp.PD != "NA") %>% mutate(GTA = rep("Cp",nrow(table(Cp.PD != "NA"))))%>%select(Study:N,Cp.PD:Cp.HC.SD,GTA) 
d.lp <- d %>% dplyr::filter(Lp.PD != "NA") %>% mutate(GTA = rep("Lp",nrow(table(Lp.PD != "NA"))))%>%select(Study:N,Lp.PD:Lp.HC.SD,GTA)
d.gamma <- d %>% dplyr::filter(gamma.PD != "NA")%>% mutate(GTA = rep("Gamma",nrow(table(gamma.PD != "NA"))))%>%select(Study:N,gamma.PD:gamma.HC.SD,GTA)
d.lambda <- d %>% dplyr::filter(lambda.PD != "NA")%>% mutate(GTA = rep("Lambda",nrow(table(lambda.PD != "NA"))))%>%select(Study:N,lambda.PD:lambda.HC.SD,GTA)
d.sigma <- d %>% dplyr::filter(sigma.PD != "NA")%>% mutate(GTA = rep("sigma",nrow(table(sigma.PD != "NA"))))%>%select(Study:N,sigma.PD:sigma.HC.SD,GTA)
d.eloc <- d %>% dplyr::filter(Eloc.PD != "NA")%>% mutate(GTA = rep("Eloc",nrow(table(Eloc.PD != "NA"))))%>%select(Study:N,Eloc.PD:Eloc.HC.SD,GTA)
d.eglob <- d %>% dplyr::filter(Eglob.PD != "NA")%>% mutate(GTA = rep("Eglob",nrow(table(Eglob.PD != "NA"))))%>%select(Study:N,Eglob.PD:Eglob.HC.SD,GTA)
d.q <- d %>% dplyr::filter(Modularity.PD != "NA")%>% mutate(GTA = rep("Modularity",nrow(table(Modularity.PD != "NA"))))%>%select(Study:N,Modularity.PD:Modularity.HC.SD,GTA)

colnames(d.cp)[ncol(d.cp)-4:1] = c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.lp)[ncol(d.lp)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.gamma)[ncol(d.gamma)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.lambda)[ncol(d.lambda)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.sigma)[ncol(d.sigma)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.eloc)[ncol(d.eloc)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.eglob)[ncol(d.eglob)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")
colnames(d.q)[ncol(d.q)-4:1]=c("pd.m","pd.sd","hc.m","hc.sd")


d.cp=     escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.cp)
d.lp=     escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.lp)
d.gamma=  escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.gamma)
d.lambda= escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.lambda)
d.sigma=  escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.sigma)
d.eloc=   escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.eloc)
d.eglob=  escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.eglob)
d.q=      escalc(measure="SMD",n1i = PD.N, n2i = HC.N,m1i = pd.m, m2i = hc.m, sd1i = pd.sd, sd2i = hc.sd,data = d.q)

meta.d <- list(d.cp=d.cp,
               d.gamma=d.gamma,
               d.eloc=d.eloc,
               d.lp=d.lp,
               d.lambda=d.lambda,
               d.eglob=d.eglob,
               d.sigma=d.sigma,
               d.q=d.q)

# meta.d <- list(d.cp=d.cp,
#                d.lp=d.lp,
#                d.eglob=d.eglob)

table(d$Study)


for (i in 1:length(meta.d)) {
  #meta.d[[i]]=meta.d[[i]] %>% filter(Study=="Baggio et al. 2014") %>% filter(Subgroup=="PD-NC")
  #meta.d[[i]]=meta.d[[i]] %>% filter(Study=="Baggio et al. 2014") %>% filter(Subgroup=="PD-MCI")
  #meta.d[[i]]=meta.d[[i]] %>% filter(Study=="Gottlich et al. 2013")
  #meta.d[[i]]=meta.d[[i]] %>% filter(Study=="Hou et al. 2018")
  #meta.d[[i]]=meta.d[[i]] %>% filter(Study=="Prajapati et al. 2020")
  meta.d[[i]]=meta.d[[i]] %>% filter(Study=="Sang et al. 2015")
}

output=data.frame(matrix(ncol = 2, nrow = length(meta.d)))
colnames(output)=c("full.d$b","full.d$vb")


for (i in 1:length(meta.d)) {
  if (nrow(meta.d[[i]])==0) {
    print(paste0(names(meta.d[i])," have not availavble data"))
  } else {
    full.d <- rma.uni(yi, vi, slab = Study, data = meta.d[[i]], method = "FE")
    full.d
    forest(full.d,shade = "zebra") 
    output[i,1]=full.d$b
    output[i,2]=full.d$vb
    rownames(output)[i]=names(meta.d[i])
  }
}

