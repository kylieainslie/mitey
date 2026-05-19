library(dplyr)


validation_data <- readRDS("vignettes/articles/validation_data.rds")
validation_data

vink_estimates <-readRDS("vignettes/articles/vink_estimates.rds")
vink_estimates

subset_measles<-validation_data %>% filter (Pathogen == "Measles", Country == "Kenya")
icc_measles <- subset_measles[5] %>% unlist(,use.names = FALSE)
si_estim(icc_measles,n_routes=6)

subset_Influenza_France<-validation_data %>% filter (Pathogen == "Influenza A(H3N2)", Country == "France")
icc_Influenza_France <- subset_Influenza_France[5] %>% unlist(,use.names = FALSE)
si_estim(icc_Influenza_France,n_routes=7)


subset_Influenza_Canada<-validation_data %>% filter (Pathogen == "Influenza A(H1N1)pdm09", Country == "Canada")
icc_Influenza_Canada <- subset_Influenza_Canada[5] %>% unlist(,use.names = FALSE)
si_estim(icc_Influenza_Canada,n_routes=2)

subset_Influenza_USA<-validation_data %>% filter (Pathogen == "Influenza A(H1N1)pdm09", Country == "USA", Author == "France")
icc_Influenza_USA <- subset_Influenza_USA[5] %>% unlist(,use.names = FALSE)
si_estim(icc_Influenza_USA,n_routes=9)

result2<-si_estim(icc_Influenza_Canada,n_routes=2)
result3<-si_estim(icc_Influenza_Canada,n_routes=3)
result4<-si_estim(icc_Influenza_Canada,n_routes=4)
result5<-si_estim(icc_Influenza_Canada,n_routes=5)
result6<-si_estim(icc_Influenza_Canada,n_routes=6)


plot_si_fit_result(result2,icc_Influenza_Canada)
plot_si_fit_result(result3,icc_Influenza_Canada)
plot_si_fit_result(result4,icc_Influenza_Canada)
plot_si_fit_result(result5,icc_Influenza_Canada)
plot_si_fit_result(result6,icc_Influenza_Canada)


subset_Measles_England<-validation_data %>% filter (Pathogen == "Measles", Country == "England", Author == "Fine")
subset_Measles_England
icc_Measles_England <- subset_Measles_England[5] %>% unlist(,use.names = FALSE)
icc_Measles_England
result_new <-si_estim(icc_Measles_England,n_routes=4, wind=7)
plot_si_fit_result(result_new,icc_Measles_England)




icc_Influenza_Canada

result2_new <-si_estim(icc_Influenza_Canada,n_routes=4, wind=2)
result3_new <-si_estim(icc_Influenza_Canada,n_routes=4, wind=3)
plot_si_fit_result(result2_new,icc_Influenza_Canada)
plot_si_fit_result(result3_new,icc_Influenza_Canada)


weights2<- result2_new$wts

weights3<- result3_new$wts




subset<-validation_data %>% filter (Pathogen == "Varicella", Country == "Australia", Author == "Vally")

icc <- subset[5] %>% unlist(,use.names = FALSE)

result_wind1 <-si_estim(ICC_dataset,n_routes=4, wind=1)
plot_si_fit_result(result_wind1,ICC_dataset)
result_wind1$mean
result_wind1$sd

result_wind7 <-si_estim(ICC_dataset,n_routes=4, wind=7)
plot_si_fit_result(result_wind7,ICC_dataset)
result_wind7$mean
result_wind7$sd

result_wind14 <-si_estim(ICC_dataset,n_routes=4, wind=14)
plot_si_fit_result(result_wind14,ICC_dataset)
result_wind14$mean
result_wind14$sd


result_true <-si_estim(ICC_real,n_routes=4, wind=7)
plot_si_fit_result(result_true,ICC_real)
result_true$mean
result_true$sd
