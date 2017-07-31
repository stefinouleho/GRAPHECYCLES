#include "graphecycle.h"
#include "string.h"
int last_chrono;
float temps;

int main(int argc, char *argv[])
{
	if( argc != 6){
		fprintf(stdout,"Missing arguments (num chebi1 , num chebi 2 ( 0  for all) , date , lenght limit ,1 similarite classique et 2 (max nb sommets communs / max ( mol1,mol2))\n");
		exit(20);
	}
//on recupère le chebi de la molécule passée en parametre
	int mol_courante,num_chebi2,taille_limite,calcul;
	mol_courante = atoi(argv[1]);
	num_chebi2 = atoi(argv[2]);
	temps 	   = atof(argv[3]); //temps max de calcul en secondes
	taille_limite = atoi(argv[4]); // la taille max du graphe produit
	calcul = atoi(argv[5]);
	printf("Molecule 1 : %d, molecule 2 : %d, temps max : %f ,taille max : %d et type :%d\n",mol_courante,num_chebi2,temps,taille_limite,calcul );	

//Lectures des molécules
	struct molecule *M = lecture_fichier_chebi();

//Trouver la molecule 1 dans la base de données
	int position = position_M(mol_courante,M);

	if(position == -1)
	{	
		printf("La molecule de chebi id %d ne se trouve pas dans la base de donnée actuelle\n",mol_courante );
		exit(21);
	}

	
// Lorsque la molecule est dans la base de donnée
	fprintf(stdout," Molécule  de numéro CHEBI: %d \n" ,mol_courante);

//parcours toutes les molécules de la base et calcule le degré de similarité avec la molécule courante
	float sim;

//On stocke le résultat dans un fichier de type similarite_(id_mol_courante)_temps_taille_type_all.data
	if( num_chebi2 != 0)
	{
		int position_2 = position_M(num_chebi2,M);

		if(position_2 == -1)
		{	
			printf("La molecule de chebi id %d ne se trouve pas dans la base de donnée actuelle\n",num_chebi2);
			exit(22);
		}

		if(calcul == 1)
			{
				sim = mesure_similarite_cycles( mol_courante,num_chebi2 ,M , temps,taille_limite);
				printf("la valeur de similarite entre %d et %d est de %f\n",mol_courante,num_chebi2,sim );
			}
			else if (calcul == 2)
			{
				sim = mesure_similarite_cycles_2( mol_courante,num_chebi2 ,M , temps,taille_limite);
			}
			last_chrono = chrono();
	}
	else
	{
		char src[64];
		int i;
		sprintf(src,"resultat/similarite_%d_%d_%d_%d_all.data",mol_courante,(int)temps,taille_limite,calcul);
		FILE *G = fopen(src,"w");
		
		if( G == NULL)
		{
			fprintf(stdout,"Cannot open the file resultat/similarite_%d_%d_%d_%d_all.data",mol_courante,(int)temps,taille_limite,calcul);
			exit(80);
		}
		
		for( i = 0; i < NB_MOLECULES; i++)
		{
			if( i %1 == 0)
			{ 
				fprintf(stdout,"\r%5d / %d couples %.3lf ",i + 1 ,NB_MOLECULES,chrono());
				fflush(stdout); 
			}

			if(M[i].chebi_id != mol_courante)
			{
				if(calcul == 1)
				{
					sim = mesure_similarite_cycles( mol_courante,M[i].chebi_id ,M , temps,taille_limite);
				}
				else if (calcul == 2)
				{
					sim = mesure_similarite_cycles_2( mol_courante,M[i].chebi_id ,M , temps,taille_limite);
				}
				last_chrono = chrono();
				fprintf(G,"%6d \t %6d \t %.6f\n", mol_courante,M[i].chebi_id ,sim);
				fflush(G);
			}
			 
		}

	}


	printf("\n Libération de la mémoire : %.3lf s\n",chrono());
	int i;
	for(i=0 ; i < NB_MOLECULES; i++) 
	{
		liberer_molecule(M[i]);
	}
		
	free(M);
	


	exit(0);
}
