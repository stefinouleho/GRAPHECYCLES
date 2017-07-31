#include "lecture_molecule_sdf.h"

struct couple{
	int a1;
	int a2;
}couple;

struct liste_voisins{
	int id_atome;
	int nb_voisins;
	int *id_voisins;
		
}liste_voisins;

struct distance{ 
	
	int *sommets;
	int *predecesseur;
}distance;

struct cycle{
	int poids;
	int *c;
}cycle;

struct vecteur{
	int taille;
	struct couple *sommets;
}vecteur;

struct lesommet{
	int id;
	int taille;
}lesommet;
void creation_fichier_squelette(int mol_courante,int position,struct molecule *M);
void la_clique_max_cycles( struct graphe g,double date);
struct graphe graphe_g12_cycles(struct graphe g12, struct molecule *M, int g1_chebi, int g2_chebi);
void calcul_cl_cycles(struct graphe g,int *dans_clique,int taille_clique,int *candidat,int taille_candidat,double date);
struct graphe graphe_produit_cycles(int g1_chebi, int g2_chebi, struct molecule *M);
float mesure_similarite_cycles (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille);
float mesure_similarite_cycles_2 (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille);
void affichage_graphe_cycle_molecule( int position, struct molecule *M);
void construction_graphe_de_cycles(int position,struct molecule *M);
int position_M( int g1_chebi,struct molecule *M);
struct molecule construction_matrice_mol(struct molecule m);
void affiche_matrice(struct molecule m);
struct couple *construction_couples(struct molecule *M,int pos1, int pos2,int taille);
struct molecule graphe_produit(int g1_chebi,int g2_chebi,struct molecule *M);
void calcul_cl(struct molecule m,int *dans_clique,int taille_clique,int *candidat,int taille_candidat, double date);	
void la_clique_max( struct molecule m,double date);
struct molecule graphe_g12(struct molecule g12, struct molecule *M, int g1_chebi, int g2_chebi);
void liberer_molecule(struct molecule g);
struct molecule * lecture_fichier_chebi();
struct reaction * lecture_fichier_reaction();
float mesure_similarite (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille_limite);
void similarite_all(int g1_chebi,struct molecule *M,double date,int taille_limite);
struct liste_voisins* construction_voisins_mol( int position, struct molecule *M);
void affichage_liste_voisinage(struct liste_voisins* voisins,int position,int nb_sommets,struct molecule *M);
void calcul_distance_sommets(int i,int position,struct molecule *M,struct liste_voisins *v,int ***vecteur_d);
int ***calcul_distance_sommets_all(int position,struct molecule *M,struct liste_voisins *v);
int verification_egalite_tableaux(int *tab1,int *tab2 , int taille);
int intersection_pcc(int *d1, int *d2,int sommet);
void affichage_pcc(int ***vecteur_d, int i , int j);
struct vecteur construction_vecteur(struct molecule *M,int position);
int calcul_poids(int *v, int taille);
int **fonction_xor( int **matrice ,int taille_y, int i , int j);
int premier_un(int **matrice,int taille_y, int i);
int op_xor(int a , int b);
int position_liaison_vecteur( struct vecteur v,int a1,int a2);
int* vecteur_distance(int ***vecteur_d,int a, int b,struct vecteur v);
int* vecteur_cycle(int ***vecteur_d, int a, int b, int sommet, struct vecteur v, struct couple cple);
struct couple *construction_couples_cycles(struct molecule *M,int pos1, int pos2,int taille);
struct molecule construction_matrice_cycles(struct molecule m);
double chrono();
void initialise();
int lecture_type(FILE *F);
int * lecture_molecules(FILE *F, int nbre);
int position_R(int id_reaction,struct reaction *R);
void affiche_reaction_id(int id_reaction,struct reaction *R);
void affiche_reaction_all(struct reaction *R);
void liberer_reaction(struct reaction r);
struct reaction * lecture_fichier_reaction();
int calcul_position_distri(float j);
