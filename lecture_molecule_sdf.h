#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NB_ATOM_NAMES 119
#define NB_MOLECULES 90130


struct nom_at { char c1, c2; };
struct liaison { int A1, A2; int l_type; };


struct couplet{
	int a1;
	int a2;
	int poids;
	int type; //1 s'il s'agit d'une liaison distance et 0 si les deux cycles partagent des liaiason
}couplet;

struct graphe{
	int nb_sommets;
	int nb_arete;
	struct lesommet *som;
	struct couplet *aretes;
	int **matrice_cycles_type;
	int **matrice_cycles_poids;
}graphe;

struct molecule{
	char chebi_name[1024];
	int chebi_id;
	int nb_atomes;
	int nb_hydrogene;
	int nb_liaisons;
	int *liste_atomes;
	int **matrice_liaisons;
	struct liaison *liste_liaisons;
	struct graphe g;
	int g_def;
} molecule;




int atom_num (char *name);
void init_atom_num ();
int lire_num_atome(FILE *F);
int valeur_char (FILE *F);
void ligne_suivante(FILE *F);
int lire_entier_3 (FILE * F);
struct liaison lire_liaison(FILE *F);
int lire_chebi_id(FILE *F);
void lire_chebi_name(FILE *F, struct molecule *M);
void lire_fin_molecule(FILE *F);
void trouver_la_fin_de_M(FILE *F);
struct molecule lire_molecule(FILE *F);
void affiche_mol(struct molecule M);
void tailles_molecules(struct molecule *M);
