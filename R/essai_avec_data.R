library(dplyr)

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

result2<-si_estim(icc_Influenza_USA,n_routes=2)
result3<-si_estim(icc_Influenza_USA,n_routes=3)
result4<-si_estim(icc_Influenza_USA,n_routes=4)
result5<-si_estim(icc_Influenza_USA,n_routes=5)
result6<-si_estim(icc_Influenza_USA,n_routes=6)


plot_si_fit_result(result2,icc_Influenza_USA)
plot_si_fit_result(result3,icc_Influenza_USA)
plot_si_fit_result(result4,icc_Influenza_USA)
plot_si_fit_result(result5,icc_Influenza_USA)
plot_si_fit_result(result6,icc_Influenza_USA)

