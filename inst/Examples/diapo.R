fais un code overleaf pour un support de presentation.
diapo 1 :
  - etude d oubreak in a household ==> interet pour le serial interval
  - definition du serial interval
  - comment il est computed = delta t = t2 - t1

diapo 2 :
   -en faisant ca on fait une hypothese forte : les dates auxquelles ont ete reportees les symptomes sont les vraies date
   - ce n est pas tjr vrai, example
   - situation ds laquelle les cas sont reportes par telephone : les cas doivent appeler un numero public, qui reportera la date du jour
   - ds le cvas classique, premier cas un mardi, il appelle le mardi,  c est reporte le mardi. deuxieme cas le lundi suivant, il appelle le jour meme, c es reporte le jour meme.
   - le vrai icc interval = 6 jours, l icc interval reporte vaut aussi 6 jours, RAS
    - imaginons le cas suivant : premier ca le mardi, tout pareil reporte le meme jour, ;ais le seond cas est le dimanche. il ne peut pas appeler le dimanche, il appelle le lundi, reporte le lundi
   - le vrai icc interval est 5 jours, mais c est rpeorte 6 jours.

diapo 3
 -cet exemple ;ontre un biais ds les donnes : la date de report n est pas tjr la vraie dte d apprition des symptomes
 - il y a umne incertitude autour des dates.

diapo 4
 - une incertitude de n days autour des dates resulte en une incertitude de 2n days autour de l ICC interval. avec une probability density function triangulaire autour de la date reportee.

 diapo 5 ;
  - pr prendre en compte cette incertitude il faut comprendre comment fonctionne l algo.
  - l algo est basee sur des iterations d un algo d expectation maxinization.
  - ds l expectstion step, l algo fait les choses suivantes . pr chaque intervalle reporte il cherche a l assigner a une des routes.
  - fais un graphique avec deux lois normales de param mu, sd et 2mu, sqrt2 sd et un histogramme de largeur de bins = 1, derriere avec des valeurs coherentes avec les deux lois.

 diapo 6
  - pour trouver a quelle route il apartientm l algo calcule la proba que s = d days en faisant l integrale de s-1 a s+1 de la loi normale multipliee avec une trianguklaire centree en s. il fait ceci pour les deux lois et en deduit une appartenance relative a chaque route, de l intervalle.
  - pour prendre en compte l incertitude sur les dates, il faut donc etendre l integrale a d-wind d+wind avec wind le window d incertitude.
  - l algo prendra ainsi en compte l incertitude autour de chaque valeur




