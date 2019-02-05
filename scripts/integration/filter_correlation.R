
# Load  interaction datasets

# PPIs associated with brain ageing (PBA)
load(file="results/integration/integrated_int.RData")

coexp_integrate_d_int<-integrated_int[integrated_int$interaction_type%in%"coexpression",]
coexp_integrate_d_int <- coexp_integrate_d_int[!coexp_integrate_d_int$data_source%in%"ADN",]
dim(coexp_integrate_d_int)#[1] 289641771         5
brain_coexp_negative<-coexp_integrate_d_int[coexp_integrate_d_int$score<0,]
#dim(brain_coexp_negative) [1] 129910633         5

# summary(brain_coexp_negative)
#ensg.A             ensg.B              score         interaction_type
#Length:129910633   Length:129910633   Min.   :-0.9693   Length:129910633
#Class :character   Class :character   1st Qu.:-0.5820   Class :character
#Mode  :character   Mode  :character   Median :-0.5297   Mode  :character
                                    #  Mean   :-0.5420
                                    #  3rd Qu.:-0.4888
                                    #  Max.   :-0.4060
#data_source
#Length:129910633
#Class :character
#Mode  :character



brain_coexp_positive<-coexp_integrate_d_int[coexp_integrate_d_int$score>0,]
dim(brain_coexp_positive)#[1] 159731138         5

#summary(brain_coexp_positive)
#ensg.A             ensg.B              score        interaction_type
#Length:159731138   Length:159731138   Min.   :0.4060   Length:159731138
#Class :character   Class :character   1st Qu.:0.5117   Class :character
#Mode  :character   Mode  :character   Median :0.5676   Mode  :character
                                      #Mean   :0.5839
                                      #3rd Qu.:0.6415
                                      #Max.   :0.9989
#data_source
#Length:159731138
#Class :character
#Mode  :character

