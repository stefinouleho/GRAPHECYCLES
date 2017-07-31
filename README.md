# GRAPHECYCLES
Similarité des molécules avec le graphe de cycles ( version sans j et k)


Première version de calcul de similarité avec le graphe de cycles. Dès que possible je rajouterais la version qui tient compte des conditions  sur la taille des cycles et la distance entre les cycles.


La compilation est pareille qu'avec le code precedent à l'exception :

il y'a un parametre en plus qui prend une valeur 1 ou 2 en fonction de la formule utilisé pour la similarité .

1-) similarité = ( nb sommet + nb aretes du sous graphe commun) *( nb sommet + nb aretes du sous graphe commun) /( nb sommet + nb aretes du graphe 1 ) * ( nb sommet + nb aretes du graphe 2)

2-) similarité = ( nb sommet ) *( nb sommet  du sous graphe commun) /( nb sommet du graphe 1 ) * ( nb sommet du graphe 2)
