# Create validation data set by combining the individual data sources used in Vink et al 2014 into a single data frame. The data as it was provided in the supplemental material can be found in scripts/vink/data.R.
# We will include the following information as columns in our data frame:
# author: the first author of the paper in which the data was found
# year: the year the paper was published
# pathogen: the virus or bacteria that caused the infection
# country: the country from which the data was collected
# icc_interval: the index case to case intervals as reported in the corresponding paper

# H3N2 Viboud et al 2004
Viboud_2004_InfluenzaAH3N2_France <- c(rep(1,38),rep(2,39),rep(3,30),rep(4,17),rep(5,7))

#pH1N1 Cauchemez et al 2009
Cauchemez_2009_InfluenzaApH1N1_USA <- c(rep(0,15),rep(1,13),rep(2,15),rep(3,15),rep(4,9),rep(5,7),rep(6,2),rep(7,2))

#pH1N1 France et al 2010
France_2010_InfluenzaApH1N1_USA <- c(rep(0,8),rep(1,5),rep(2,16),rep(3,14),rep(4,11),rep(5,7),rep(6,2),rep(7,3), 8,10,11, rep(15,2), 16,19,20,21,22,23)

# pH1N1 Papenburg et al 2010
Papenburg_2010_InfluenzaApH1N1_Canada <- c(rep(1,11),rep(2,7),rep(3,9),rep(4,6),rep(5,5),rep(6,4),7,8,9,10,13,16)

# pH1N1 Morgan, 2010
Morgan_2010_InfluenzaApH1N1_USA <- c(7,2,3,2,1,1,2,9,4,4,1,4,5,6,3,2,3,3,5,5,8,6,2,3,4,2,4,8,3,5,4,5)

# pH1N1 data Hahne et al 2009
Hahne_2010_InfluenzaApH1N1_Netherlands <- c(rep(0,11),rep(1,9),rep(2,13),rep(3,5),rep(4,9),5,6,7)

# Influenza pH1N1 Savage BMC public health 2011
Savage_2011_InfluenzaApH1N1_Canada <- c(0,rep(1,9),rep(2,10),rep(3,15),rep(4,6),rep(5,4),rep(6,3),7,rep(8,2),9,14,16,17,20)

# RSV Crowcroft 2008, confirmed cases
Crowcroft_2008_RSV_England <- c(rep(1,3),rep(2,6),rep(3,3),rep(4,2),rep(5,3),rep(6,5),rep(7,5),rep(8,4),rep(9,3),rep(10,2), rep(11,3), rep(12,3),13,14,rep(15,4),rep(16,2),17,rep(20,2),22,rep(24,2),25,rep(26,2))

# Pertussis Greeff et al. 2010
Greeff_2010_Pertussis_Netherlands <- c(rep(0,4),1,rep(2,4),rep(3,5),rep(4,3),rep(5,7),rep(6,6),rep(7,6),rep(8,5),rep(9,5),rep(10,10), rep(11,4),rep(12,9),rep(13,9),rep(14,6),rep(15,9),rep(16,8),rep(17,6),rep(18,8),rep(19,9),rep(20,11), rep(21,6),rep(22,8),rep(23,8),rep(24,8),rep(25,7),rep(26,3),rep(27,5),rep(28,11),rep(29,6),rep(30,10), rep(31,6),rep(32,4),rep(33,7),rep(34,3),rep(35,4),rep(36,4),rep(37,5),rep(38,2),rep(40,4),rep(41,4), rep(42,7),43,44,45,46,rep(47,5),rep(48,2),rep(49,2),50,rep(51,3),rep(52,3),rep(53,3),rep(54,2),55, rep(56,2),57,58,59,60,61,62,rep(64,2),65,66,rep(67,2),68,70,72,76,77,80,81,84,90,92,93,101,105,108,rep(116,2), 130,155)

# Measles Simpson 1952
Simpson_1952_Measles_England <- c(rep(6,4),rep(7,8),rep(8,14),rep(9,31),rep(10,29),rep(11,42),rep(12,25),rep(13,16),rep(14,16), rep(15,10),rep(16,4),rep(17,2),rep(18,2))

# Measles Bailey 1954
Bailey_1954_Measles_England <- c(rep(0,5),rep(1,13),rep(2,5),rep(3,4),rep(4,3),rep(5,2),rep(6,4),rep(7,11),rep(8,5),rep(9,25), rep(10,37),rep(11,38),rep(12,26),rep(13,12),rep(14,15),rep(15,6),rep(16,3),17,rep(18,3),21)

# Measles Aaby et al. 1990
Aaby_1990_Measles_Kenya <- c(rep(0,42),rep(1,28),rep(2,22),rep(3,14),rep(4,14),rep(5,13),rep(6,25),rep(7,30),rep(8,38), rep(9,24),rep(10,30),rep(11,42),rep(12,26),rep(13,20),rep(14,20),rep(15,10),rep(16,11),rep(17,6), rep(18,3), rep(19,5),rep(20,3),21,22,23,23,24)

# Measles Chapin, 1925
Chapin_1925_Measles_USA <- c(rep(1,531),rep(2,257),rep(3,226),rep(4,189),rep(5,132),rep(6,144),rep(7,254),rep(8,339), rep(9,283),rep(10,401),rep(11,500),rep(12,416),rep(13,463),rep(14,467),rep(15,381),rep(16,188), rep(17,99),rep(18,69),rep(19,48),rep(20,47),rep(21,55),rep(22,47),rep(23,23),rep(24,19),rep(25,14), rep(26,13), rep(27,9),rep(28,15),rep(29,12),rep(30,18))

# Measles Fine 2003, USA
Fine_2003_Measles_USA <- c(9,11,rep(13,3),rep(15,2),rep(17,2),23,rep(24,2),rep(27,4),28,rep(29,6),30,rep(31,6), 32,rep(33,2),34,35,36,rep(37,3),rep(38,15),rep(39,7),rep(40,5),rep(41,13),rep(42,6),rep(43,6),44,48, rep(49,3), rep(50,4),rep(51,2),rep(53,2),54)

# Measles Fine 2003, England
Fine_2003_Measles_England <-c(11,rep(13,2),14,rep(15,3),24,25,rep(26,2),rep(27,4),rep(28,6),38,rep(39,2),40,42,rep(43,2))

# Varicella Simpson, 1952
Simpson_1952_Varicella_England <- c(rep(7,7),rep(8,8),rep(9,5),rep(10,9),rep(11,7),rep(12,11),rep(13,22),rep(14,55), rep(15,19),rep(16,10),rep(17,14),rep(18,6),rep(19,7),rep(20,4))

# Varicella Vally, 2007
Vally_2007_Varicella_Australia <- c(8,10,rep(11,2),rep(12,6),rep(13,4),rep(14,2),rep(15,2),16,17,rep(19,2),25,rep(26,2), rep(27,3),28,29, rep(30,2),rep(32,3),33,rep(34,2),rep(36,2),rep(40,2),47)

# Varicella Lai, 2011
Lai_2011_Varicella_Taiwan <- c(13,13,14,15,16,27,29,29,31,31,rep(32,3),35,45)

# Mumps Simpson, 1952
Simpson_1952_Mumps_England <- c(rep(10,2),rep(12,4),rep(13,3),rep(14,20),rep(15,7),rep(16,15),rep(17,14),rep(18,13), rep(19,12), rep(20,10),rep(21,10),rep(22,8),23,rep(24,6),rep(25,6),rep(26,3),rep(27,3),rep(28,3), 29,32)

# Rubella Aycock, 1946
Aycock_1946_Rubella_Unknown <- c(rep(15,9),rep(16,10),rep(17,25),rep(18,11),rep(19,26),rep(20,8),rep(21,10), rep(22,6), 23,rep(35,3), rep(36,9),rep(37,2),rep(38,2),39,40,42,48)

# Smallpox Tilburg
# Tilburg <- c(22,22,19,21,16,15,21,19,16,33,43,18,19,17,22,25,27,27,38,22,26,16,17,15,15,13,15,15, 30,14,28,15,15,15,17,16,20,15,23,20,23,27,25,33,33)

# Smallpox Fine 2003 Kosovo (originally from: Fenner et al. Smallpox and its eradication chp 23)
Fine_2003_Smallpox_Kosovo <- c(14,15,rep(16,2),rep(17,2),18,rep(19,3),20,rep(28,3),rep(29,3),rep(30,2),rep(31,6),rep(32,9), rep(33,15),rep(34,25),rep(35,12),rep(36,27),rep(37,19),rep(38,9),rep(39,4),rep(40,4),41,43,rep(44,4), rep(45,3),rep(46,4),rep(47,6),50,51,53,54,rep(55,2))

# Smallpox Fine 2003 Germany (originally from: wehrle et al 1970)
Fine_2003_Smallpox_Germany <- c(rep(12,3),13,14,rep(15,2),16,17,rep(18,3),19,rep(21,4),34,37)

# Hepatitis A Brodribb
# Brodribb <- c(20, 21, 21, 23, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 27, 28,28, 28, 28, 29, 29, 30, 31, 31, 32, 32, 32, 47, 48, 48, 50, 50, 51,51, 51, 51, 52, 54, 54, 55, 57, 57, 59, 59, 62, 63, 72, 79, 66, 90, 96)

# Create validation data set
library(dplyr)
library(tidyr)
library(purrr)

# create function to fill vectors so they're all equal lengths
bind_vectors <- function(..., fill_value = NA) {
  # Get the list of vectors
  vec_list <- list(...)
  vec_names <- sapply(substitute(list(...))[-1], deparse)
  # Find the maximum length of the vectors
  max_len <- max(map_int(vec_list, length))
  # Adjust the length of each vector by filling with the specified fill value
  vec_list_padded <- map(vec_list, ~c(.x, rep(fill_value, max_len - length(.x))))
  # Combine the vectors into a data frame with explicit names
  names(vec_list_padded) <- vec_names
  # Combine the vectors into a data frame
  return(as_tibble(vec_list_padded))
}

df <- bind_vectors(Viboud_2004_InfluenzaAH3N2_France,
             Cauchemez_2009_InfluenzaApH1N1_USA,
             France_2010_InfluenzaApH1N1_USA,
             Papenburg_2010_InfluenzaApH1N1_Canada,
             Morgan_2010_InfluenzaApH1N1_USA,
             Hahne_2010_InfluenzaApH1N1_Netherlands,
             Savage_2011_InfluenzaApH1N1_Canada,
             Crowcroft_2008_RSV_England,
             Greeff_2010_Pertussis_Netherlands,
             Simpson_1952_Measles_England,
             Bailey_1954_Measles_England,
             Aaby_1990_Measles_Kenya,
             Chapin_1925_Measles_USA,
             Fine_2003_Measles_USA,
             Fine_2003_Measles_England,
             Simpson_1952_Varicella_England,
             Vally_2007_Varicella_Australia,
             Lai_2011_Varicella_Taiwan,
             Simpson_1952_Mumps_England,
             Aycock_1946_Rubella_Unknown,
             Fine_2003_Smallpox_Kosovo,
             Fine_2003_Smallpox_Germany
             )

val_data <- df %>%
  pivot_longer(cols = everything(),
               names_to = "Group",   # Collect column names into one column
               values_to = "ICC_interval") %>%
  # Split the "Group" column into two columns based on the "_" character
  separate(Group, into = c("Author", "Year", "Pathogen", "Country"), sep = "_") %>%
  arrange(Author) %>%
  filter(!is.na(ICC_interval)) %>%
  mutate(Pathogen = ifelse(Pathogen == "InfluenzaApH1N1", "Influenza A(H1N1)pdm09", Pathogen),
         Pathogen = ifelse(Pathogen == "InfluenzaAH3N2", "Influenza A(H3N2)", Pathogen))

saveRDS(val_data, "vignettes/articles/validation_data.rds")

### Create data frame containing estimates from Vinnk et al (Table 3)
library(tibble)

df <- tribble(
  ~`Author`, ~ `Year`, ~Pathogen, ~Country, ~`Mean`, ~`SD`, ~`95% CI of Mean`,
  "Hahne", 2009, "Influenza A(H1N1)pdm09", "Netherlands", 1.7, 1.22, "1.3, 2.0",
  "Cauchemez", 2009, "Influenza A(H1N1)pdm09", "United States", 2.1, 1.22, "1.8, 2.4",
  "Savage", 2011, "Influenza A(H1N1)pdm09", "Canada", 2.8, 0.82, "2.6, 3.0",
  "Papenburg", 2010, "Influenza A(H1N1)pdm09", "Canada", 2.9, 1.22, "2.5, 3.2",
  "France", 2010, "Influenza A(H1N1)pdm09", "United States", 3.0, 0.92, "2.8, 3.2",
  "Morgan", 2010, "Influenza A(H1N1)pdm09", "United States", 3.7, 1.12, "3.3, 4.1",
  "Viboud", 2004, "Influenza A(H3N2)", "France", 2.2, 0.82, "2.1, 2.4",
  "Aaby", 1990, "Measles", "Kenya", 9.9, 2.42, "9.7, 10.2",
  "Bailey", 1954, "Measles", "England", 10.9, 1.92, "10.6, 11.1",
  "Simpson", 1952, "Measles", "England", 10.9, 2.02, "10.6, 11.2",
  "Chapin", 1925, "Measles", "United States", 11.9, 2.62, "11.8, 12.0",
  "Fine", 2003, "Measles", "England", 13.7, 1.52, "13.1, 14.3",
  "Fine", 2003, "Measles", "United States", 13.8, 2.52, "13.3, 14.3",
  "Simpson", 1952, "Mumps", "England", 18.0, 3.52, "17.4, 18.6",
  "Greeff", 2010, "Pertussis", "Netherlands", 22.8, 6.52, "22.1, 23.5",
  "Crowcroft", 2008, "RSV", "England", 7.5, 2.12, "7.0, 8.1",
  "Aycock", 1946, "Rubella", "Unknown", 18.3, 2.02, "18.0, 18.6",
  "Fine", 2003, "Smallpox", "Germany", 16.7, 3.32, "15.1, 18.3",
  "Fine", 2003, "Smallpox", "Kosovo", 17.3, 1.92, "17.0, 17.6",
  "Vally", 2007, "Varicella", "Australia", 13.1, 2.22, "12.4, 13.8",
  "Simpson", 1952, "Varicella", "England", 14.1, 2.42, "13.7, 14.4",
  "Lai", 2011, "Varicella", "Taiwan", 14.2, 1.32, "13.5, 14.9"
)

saveRDS(df, "vignettes/articles/vink_estimates.rds")
