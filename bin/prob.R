#This will compute the probabilities of seeing all 4 red, and all 4 blue

#this code performs the statistics to explain our ability to identify suppressors

############################################
##### THIS NEEDS TO BE GONE OVER AGAIN #####
############################################

#These are the number of reds in the figure
red_c1 <- 8/34
red_c23 <- 11/34 #2 is 10 and 3 is 12
red_c4 <- 7/34
none_c5 <- 12/34

#probably that we have 4 red in the checklist figure
prob_all_red <- red_c1*red_c23*red_c4
prob_all_red*34

#you would expect 0.532 all reds by chance (ie given these numbers, if reds were distributed randomly we would not expect them
#to occur all on the same row)

pbinom(4, 34, prob_all_red, lower.tail = F)
# 0.000180091
#This tells us that we have more complete reds than would be expected by chance

#probably that we have no reds
prob_no_red <- (1-red_c1)*(1-red_c23)*(1-red_c4)
prob_no_red*34

#this is the probability of seeing this many bacteria w/o a red
pbinom(19, 34, prob_no_red, lower.tail = F)
#value is 0.0278
#This tells us that we are robust at identifying bacteria that aren't capable of
#suppressing the flg22 mediated immune response, and that the reds are not being randomly placed


#prob that we have all blue
b_c1 <- 17/34
b_c23 <- ((12+7)/2)/34
b_c4 <- 27/34
#b_c5 <- 22/34

prob_all_blue <- b_c1*b_c23*b_c4
prob_all_blue*34

pbinom(3,34, prob_no_blue, lower.tail = F) 
#we get around the same number of all blues as we would expect to see by chance

#prob that we dont see a blue
prob_no_blue <- (1-b_c1)*(1-b_c23)*(1-b_c4)
prob_no_blue*34
pbinom(5,34, prob_no_blue, lower.tail = F) 

#We have more bacteria that have no blues than we would expect by chance.