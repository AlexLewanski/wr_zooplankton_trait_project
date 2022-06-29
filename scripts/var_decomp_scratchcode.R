library(tidyverse)
trait_info <- data.frame(species = c(rep('sp1', 6), rep('sp2', 8), rep('sp3', 10), rep('sp4', 16)),
                         trait = c(c(1:6), c(2:9), c(1:10), c(c(tind2*2), c(tind2*2)-2)) 
                         )


trait_info_alt <- data.frame(species = c(rep('sp1', 1), rep('sp2', 8), rep('sp3', 10), rep('sp4', 16)),
                         trait = c(c(1), c(2:9), c(1:10), c(c(tind2*2), c(tind2*2)-2)) 
)


# abund_info <- stack(table(trait_info$species)) %>% 
#   select(ind, values) %>% 
#   mutate(ind = as.character(ind)) %>%
#   rename(species = ind,
#          abund = values)

trait_info1 <- trait_info %>% 
  group_by(species) %>% 
  mutate(rel_freq = 1/n()) %>% 
  ungroup()


abund_info <- data.frame(species = paste0('sp', 1:4),
                         rel_abund = c(0.1, 0.2, 0.3, 0.4))

abund_info <- data.frame(species = paste0('sp', 1:4),
                         rel_abund = c(0.1, 0.2, 0.3, 0.4))

combined_info <- left_join(trait_info1, abund_info, by = 'species') %>% 
  mutate(Pindp = rel_freq*rel_abund) %>% 
  group_by(species) %>% 
  mutate(mean_trait = mean(trait)) %>% 
  ungroup()

xcomp <- combined_info %>% 
  group_by(species) %>% 
  summarize(mean_trait = mean(trait), .groups = "drop") %>% 
  left_join(., abund_info, by = 'species') %>% 
  mutate(xcomp = mean_trait*rel_abund) %>% 
  pull(xcomp) %>% 
  sum()


within_var <- combined_info %>% 
  group_by(species) %>%
  mutate(sst = (trait - mean_trait)^2) %>% 
  summarize(rel_abund = first(rel_abund),
            varx = mean(sst),
            .groups = 'drop') %>%
  mutate(witp = rel_abund*varx) %>% 
  pull(witp) %>% 
  sum()

between_var <- combined_info %>% 
  group_by(species) %>% 
  slice(1) %>% 
  ungroup() %>% 
  mutate(sstp = ((mean_trait - xcomp)^2)*rel_abund) %>% 
  pull(sstp) %>% 
  sum()



for(i in 1:4){
  sstp[i]<-((t[i]-xcomp)^2)
}
betp<-sum(c(sstp[1]*p1, sstp[2]*p2, sstp[3]*p3, sstp[4]*p4))  
  
  xcomp

sscomp<-tind
for(i in 1:length(tind)){
  sscomp[i]<-((tind[i]-xcomp)^2)*Pindp[i]
}
totp<-sum(sscomp)




sst1<-tind1
for(i in 1:length(sst1)){
  sst1[i]<-((tind1[i]-t1)^2)
}





  




var_decomp <- function(abundance_df, trait_df) {
  
  #STEP 1: calculate the fraction of total observations per species represented by each
  #        trait value; aka 1/sample size split by species
  trait_df1 <- trait_df %>% 
    group_by(species) %>% 
    mutate(rel_freq = 1/n()) %>% 
    ungroup()
  
  #STEP 2: - join the processed trait info with the abundance information
  #        - calculate Pindp and mean trait value for each species and add
  #          the info the merged dataframe
  combined_df <- left_join(trait_df1, abundance_df, by = 'species') %>% 
    mutate(Pindp = rel_freq*rel_abund) %>% 
    group_by(species) %>% 
    mutate(mean_trait = mean(trait)) %>% 
    ungroup()
  
  #STEP 3: calculate within species variance
  within_var <- combined_df %>% 
    group_by(species) %>%
    mutate(sst = (trait - mean_trait)^2) %>% 
    summarize(rel_abund = first(rel_abund),
              varx = mean(sst),
              .groups = 'drop') %>%
    mutate(witp = rel_abund*varx) %>% 
    pull(witp) %>% 
    sum()
  
  #STEP 4: calculate betweeen species variance
  #4a: calculate xcomp
  xcomp <- combined_df %>% 
    group_by(species) %>% 
    summarize(mean_trait = mean(trait), .groups = "drop") %>% 
    left_join(., abundance_df, by = 'species') %>% 
    mutate(xcomp = mean_trait*rel_abund) %>% 
    pull(xcomp) %>% 
    sum()
  
  #4b: calculate between species variance
  between_var <- combined_df %>% 
    group_by(species) %>% 
    slice(1) %>% 
    ungroup() %>% 
    mutate(sstp = ((mean_trait - xcomp)^2)*rel_abund) %>% 
    pull(sstp) %>% 
    sum()
  
  ### RETURN VARIANCE INFO ###
  return(
    data.frame(between_var = between_var,
               within_var = within_var,
               between_var_prop = between_var/(between_var + within_var),
               within_var_prop = within_var/(between_var + within_var)
               )
  )
}




var_decomp(abundance_df = abund_info, trait_df = trait_info_alt)

var_decomp(abundance_df = trait_info %>% 
             group_by(species) %>% 
             summarize(abund = n()) %>% 
             ungroup() %>% 
             mutate(rel_abund = abund/nrow(trait_info)) %>% 
             select(species, rel_abund), 
           trait_df = trait_info)


abund_info

trait_info %>% 
  group_by(species) %>% 
  summarize(abund = n()) %>% 
  ungroup() %>% 
  mutate(rel_abund = abund/nrow(trait_info)) %>% 
  select(species, rel_abund)

tind1<-c(1:6)
tind2<-c(2:9)
tind3<-c(1:10)
tind4<-c(c(tind2*2), c(tind2*2)-2)

tind<-c(tind1, tind2, tind3, tind4)

#####relative abundance of individuals, i.e. 1/Nindi#####

Nind1<-1/6
Nind2<-1/8
Nind3<-1/10
Nind4<-1/16

Pind<-c(rep(Nind1, length(tind1)), rep(Nind2, length(tind2)), rep(Nind3, length(tind3)), rep(Nind4, length(tind4)))/4#or c(rep(Nind1/4, length(tind1)), rep(Nind2/4, length(tind2)), rep(Nind3/4, length(tind3)), rep(Nind4/4, length(tind4)))
sum(Pind)

####mean trait per species####
t1<-sum(tind1*Nind1)#or mean(tind1)
t2<-sum(tind2*Nind2)#or mean(tind2)
t3<-sum(tind3*Nind3)#or mean(tind3)
t4<-sum(tind4*Nind4)#or mean(tind4)

t<-c(t1, t2, t3, t4)

####hypothetic relative abundance of species in the community for the cases that pi does not equate 1/Nsp###
p1<-0.1
p2<-0.2
p3<-0.3
p4<-0.4

Pindp<-c(rep(Nind1*p1, length(tind1)), rep(Nind2*p2, length(tind2)), rep(Nind3*p3, length(tind3)), rep(Nind4*p4, length(tind4)))
sum(Pindp)

#mean trait per community when pi=1/Nsp#
xcom<-mean(c(t1, t2, t3, t4))#or sum(c(t1/4, t2/4, t3/4, t4/4))#or   sum(tind*Pind)# or mean(t)

#mean trait per community when pi is not 1/Nsp#
xcomp<-sum(c(t1*p1, t2*p2, t3*p3, t4*p4))#


-------------------
  #variance approach#
  -------------------
  
  #1. #####variance within each species####

#sp1#
sst1<-tind1
for(i in 1:length(sst1)){
  sst1[i]<-((tind1[i]-t1)^2)
}
v1<-mean(sst1)#or sum(sst1)/length(sst1)

#sp2#
sst2<-tind2
for(i in 1:length(sst2)){
  sst2[i]<-((tind2[i]-t2)^2)
}
v2<-mean(sst2)#or sum(sst2)/length(sst2)

#sp3#
sst3<-tind3
for(i in 1:length(sst3)){
  sst3[i]<-((tind3[i]-t3)^2)
}
v3<-mean(sst3)#or sum(sst3)/length(sst3)

#sp4#
sst4<-tind4
for(i in 1:length(sst4)){
  sst4[i]<-((tind4[i]-t4)^2)
}
v4<-mean(sst4)#or sum(sst4)/length(sst4)



#1.2########overall within species variance#######

#overall within species variance for all 4 species, for the cases that pi=1/Nsp# 
wit<-mean(c(v1, v2, v3, v4))#or sum(c(v1, v2, v3, v4))/4

#overall within species variance for all 4 species, for the cases that pi does not equate 1/Nsp# 
witp<-sum(c(v1*p1, v2*p2, v3*p3, v4*p4))

2. ########between species variance#######


#when p1=1/Nsp#
sst<-rep(0, 4)

for(i in 1:4){
  sst[i]<-((t[i]-xcom)^2)
}

#between species variance when p1=1/Nsp#
bet<-mean(sst)#or sum(c(sst[1]*1/4, sst[2]*1/4, sst[3]*1/4, sst[4]*1/4))

#when p1=1/Nsp#
sstp<-rep(0, 4)

for(i in 1:4){
  sstp[i]<-((t[i]-xcomp)^2)
}

#between species variance when p1 is not 1/Nsp#
betp<-sum(c(sstp[1]*p1, sstp[2]*p2, sstp[3]*p3, sstp[4]*p4))


#3. ########total community variance#######

#p1=1/Nsp#
sscom<-tind
for(i in 1:length(tind)){
  sscom[i]<-((tind[i]-xcom)^2)*Pind[i]
}
tot<-sum(sscom)

#p1 is not 1/Nsp#
sscomp<-tind
for(i in 1:length(tind)){
  sscomp[i]<-((tind[i]-xcomp)^2)*Pindp[i]
}
totp<-sum(sscomp)

4. ########total community variance (see 3 )is equal to within species variance (1.2) + between species variance (2)#######

#when pi=1/Nsp#
tot==bet+wit

#when pi is not 1/Nsp#
totp==betp+witp#it gives false only because periodic numbers#
totp
betp+witp
