
library(extraDistr)
# Paramètres
N <- 10000
mu <- 22      # moyenne de la loi normale
sigma <- 3   # écart-type de la loi normale



# Créer la liste P1 : valeurs uniformément réparties entre 0 et 14
P1 <- rdunif(n = N, min = 1, max = 14)

# Créer la liste norm : valeurs suivant une loi normale, en
norm <- trunc(rnorm(n = N, mean = mu, sd = sigma))

# Créer la liste P2 : somme de P1 et norm
P2 <- P1 + norm

I<- P2-P1


breaks1<-seq(0.5,15.5,1)

# Visualisation
par(mfrow = c(1, 3))
hist(P1, main = "Distribution de P1\n(Uniforme)", xlab = "Valeurs", col = "skyblue",breaks = breaks1)
hist(P2, main = "Distribution de P2\n(P1 + norm)", xlab = "Valeurs", col = "lightgreen")
hist(I, main = "Distribution de I = P2 - P1", xlab = "Valeurs", col = "cyan")
par(mfrow = c(1, 1))


fit_mean <-mean(I)
fit_sd<-sd(I)


# Fonction pour transformer les valeurs (6 et 7 modulo 7 deviennent 8 modulo 7)
transformer_valeurs <- function(x) {
  reste <- x %% 7
  # Si le reste est 6, on ajoute 2 pour que 6 -> 8
  x_transforme <- ifelse(reste == 0 & x>3 , x + 1, ifelse(reste == 6 & x>3, x + 2, x))
  return(x_transforme)
}

# Créer les listes P1' et P2' avec la transformation
P1_prime <- transformer_valeurs(P1)
P2_prime <- transformer_valeurs(P2)

I_prime <- P2_prime - P1_prime
# Visualisation
par(mfrow = c(1, 3))
hist(P1_prime, main = "Distribution de P1_prime", xlab = "Valeurs", col = "skyblue",breaks = breaks1)
hist(P2_prime, main = "Distribution de P2_prime", xlab = "Valeurs", col = "lightgreen")
hist(I_prime, main = "Distribution de I_prime = P2_prime - P1-prime", xlab = "Valeurs", col = "cyan")
par(mfrow = c(1, 1))


fit_mean2 <-mean(I_prime)
fit_sd2 <-sd(I_prime)
