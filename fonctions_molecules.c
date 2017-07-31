#include "graphecycle.h"
#include <sys/time.h>
#define AUCUNE_LIAISON (-1024)


int taille_clique_max =0;
int *dans_clique_max= NULL;
int liaison_max = 0;
double last_chrono;

double chrono() 
{
	struct timeval tv;
	static double date_deb = 0.0;
	double date_courante;
	gettimeofday(&tv,NULL);
	date_courante = tv.tv_sec + ((double)tv.tv_usec)*1e-6;
	if (date_deb == 0.0) date_deb = date_courante;
	return date_courante-date_deb;
}

int position_M( int g1_chebi,struct molecule *M)
{ // trouve la poition d'une molecule de chebi g1_chebi dans M
	int i;
	for(i=0;i< NB_MOLECULES +1 ;i++)
	{
		if(M[i].chebi_id == g1_chebi)
			return i;
	}
	
	return -1;
}


struct molecule construction_matrice_mol(struct molecule m)
{//construction de la matrice de liaison d'une molecule
	
	int i,j;
	if(m.matrice_liaisons == NULL)
	{
		m.matrice_liaisons =  malloc(m.nb_atomes * sizeof(int *));
		
		for(i=0;i< m.nb_atomes;i++) m.matrice_liaisons[i] =  malloc(m.nb_atomes * sizeof(int));
		
		for(i=0;i< m.nb_atomes;i++)
		{
			for(j=0;j< m.nb_atomes;j++)
			 m.matrice_liaisons[i][j] = AUCUNE_LIAISON;
		}

		for(i =0; i< m.nb_liaisons;i++)
		{
			m.matrice_liaisons[m.liste_liaisons[i].A1-1][m.liste_liaisons[i].A2-1]=m.liste_liaisons[i].l_type;
			m.matrice_liaisons[m.liste_liaisons[i].A2-1][m.liste_liaisons[i].A1-1]=m.liste_liaisons[i].l_type;
		}
		

	}

	return m;
	
}


void affiche_matrice(struct molecule m)
{//affcihe la matrice d'une molecule m
	int i,j;
	printf("Affichage de la matrice m \n");
	
	if( m.matrice_liaisons == 	NULL)
	{
		printf("La matrice de cette molecule n'a pas encore eté défini \n");
		return;
	}
	for(i=0;i< m.nb_atomes;i++)
	{
		for(j=0;j< m.nb_atomes;j++)
			printf("%d ",m.matrice_liaisons[i][j] );
		
		printf("\n");
	}
	
}

struct couple *construction_couples(struct molecule *M,int pos1, int pos2,int taille)
{//Construction des couples d'atomes compatibles
	
	int n = 0,i,j;
	
	struct couple *couple_at;
	couple_at = malloc(taille * sizeof(couple));
	for(i= 0; i < M[pos1].nb_atomes;i++)
	{ 
		for(j= 0; j < M[pos2].nb_atomes;j++)
		{
			if( M[pos1].liste_atomes[i] == M[pos2].liste_atomes[j])
			{
				couple_at[n].a1 = i;
				couple_at[n].a2 = j;
				n++;
			}
		}
	}
	return couple_at;

}

struct molecule graphe_produit(int g1_chebi,int g2_chebi,struct molecule *M)
{ //prend en entrée les chebi id de deux molecules  et contruit le graphe produit
		
		
	struct molecule g12;
	
	//trouve la position des molecules g1 et g2
	int taille= 0,pos1,pos2;
	pos1 =position_M(g1_chebi,M);
	pos2 =position_M(g2_chebi,M);
	
	//calcul de la taille du graphe produit
	int i,j;
	for(i= 0; i < M[pos1].nb_atomes; i++)
	{ 
		for(j= 0; j < M[pos2].nb_atomes;j++)
		{
			if( M[pos1].liste_atomes[i] == M[pos2].liste_atomes[j]) taille ++;
		}
	}
	
	//couple de liaisons entre les nouveaux sommets
	struct couple * couple_atome = construction_couples(M,pos1,pos2,taille);
	
	//construction de la matrice de liaison d'une molecule
	M[pos1]= construction_matrice_mol(M[pos1]);
	M[pos2]=construction_matrice_mol(M[pos2]);
	
	//initialisation de g12
	g12.nb_atomes= taille;
	g12.nb_liaisons = 0;
	g12.liste_atomes = malloc(g12.nb_atomes * sizeof(int));
	g12.liste_liaisons = NULL;
	g12.matrice_liaisons =  malloc(g12.nb_atomes * sizeof(int *));
	
	for(i=0;i< g12.nb_atomes;i++) g12.matrice_liaisons[i] =  malloc(g12.nb_atomes * sizeof(int));

	//remplissage des atomes	
	for(i=0;i<g12.nb_atomes;i++) g12.liste_atomes[i] = M[pos1].liste_atomes[couple_atome[i].a1];
	
	//remplissage des liaisons
	int i1,i2,j1,j2;
	for(i= 0; i < taille ;i++)
	{
		for(j= 0; j < taille ;j++)
			g12.matrice_liaisons[i][j] =0;
		
	}

	
	for(i= 0; i < taille ;i++)
		g12.matrice_liaisons[i][i] = 1;

	for(i= 0; i < taille ;i++)
	{ 
		i1 = couple_atome[i].a1;
		i2 = couple_atome[i].a2;
		for(j= i + 1 ; j < taille ;j++)
		{
			
			j1=couple_atome[j].a1;
			j2=couple_atome[j].a2;
			
			if( M[pos1].matrice_liaisons[i1][j1] == M[pos2].matrice_liaisons[i2][j2] )
			{

				if(((i1 == j1) && (M[pos2].matrice_liaisons[i2][j2] != AUCUNE_LIAISON)) || ((i2 == j2) && (M[pos1].matrice_liaisons[i1][j1] != AUCUNE_LIAISON)) ||( (i1 != j1) && (i2!=j2))||( (i1 ==j1) && (i2==j2)) )
				{
					g12.matrice_liaisons[i][j] = 1;
					g12.matrice_liaisons[j][i] = 1;
					
									
				}
			}
			
		}
	}
	free(couple_atome);
	//affiche_matrice(g12);
	return g12;
	
}


void calcul_cl(struct molecule m,int *dans_clique,int taille_clique,int *candidat,int taille_candidat, double date)
{	//calcul de la clique max recursif
	
	if (date != 0)
	{	//{printf("%f %f\n",chrono(),last_chrono);return;}
		if(chrono() - last_chrono > date) 
			return;
	}
	int i,j;
	if( taille_candidat == 0)
	{
		if(taille_clique > taille_clique_max){
			taille_clique_max = taille_clique;
			for (i = 0 ;  i < m.nb_atomes ; i ++)
				dans_clique_max[i] = dans_clique[i];
		}
		return;
	}
	
	if (taille_candidat + taille_clique <= taille_clique_max)
	{

		return;
	}
	
	//else 
	int taille_candidat_temp;
	int *candidat_temp;
	candidat_temp = malloc( m.nb_atomes * sizeof(int));



	for (i = 0 ;  i < m.nb_atomes ; i ++)
	{
		if ( candidat[i] == 1)
		{
			candidat[i] = 0;
			dans_clique[i] = 1 ;
			taille_candidat_temp = taille_candidat;
			
			for (j = 0 ;  j < m.nb_atomes ; j ++)
			{
				candidat_temp[j] = candidat[j]; 
				if ((candidat[j] == 1) && (m.matrice_liaisons[i][j] == 0))
				{
					candidat_temp[j] = 0;
					taille_candidat_temp--;	
				}	
			}
			
			taille_candidat_temp --;

			calcul_cl(m,dans_clique,taille_clique + 1,candidat_temp,taille_candidat_temp,date);
			dans_clique[i] = 0;
			
		}	
		
	}
	free(candidat_temp);
		 
}

	
void la_clique_max( struct molecule m,double date)
{	//Debut calcul de la clique -- Initialisation
	int i;
	int *candidat;
	int *dans_clique;
	
	dans_clique_max = malloc( m.nb_atomes *sizeof(int));
	if (!dans_clique_max) { fprintf(stderr,"cannot malloc dans_clique_max %d\n",m.nb_atomes); exit(41); }
	dans_clique = malloc( m.nb_atomes *sizeof(int));
	if (!dans_clique) { fprintf(stderr,"cannot malloc dans_clique %d\n",m.nb_atomes); exit(42); }
	candidat = malloc( m.nb_atomes *sizeof(int));
	if (!candidat) { fprintf(stderr,"cannot malloc candidat %d\n",m.nb_atomes); exit(43); }
	
	//initialisation 
	for(i = 0; i < m.nb_atomes ; i++ )
	{
		candidat[i] 	= 1;
		dans_clique[i]	= 0;
		dans_clique_max[i] = 0;
	}
	
	taille_clique_max = 0;
	
	calcul_cl(m,dans_clique,0,candidat,m.nb_atomes,date); // 0 taille de la clique initial  et m.nb_atome = nb sommets candidats
	 
	free(dans_clique);
	free(candidat);
	
}


struct molecule graphe_g12(struct molecule g12, struct molecule *M, int g1_chebi, int g2_chebi)
{ //contruction du graphe commun

	struct molecule clique;
	int nb_at =0,nb_liaisons=0,i,j,i1,j1;
	int pos1,pos2;
	pos1 =position_M(g1_chebi,M);
	pos2 =position_M(g2_chebi,M);
	
	struct couple *couple_atome = construction_couples(M,pos1,pos2,g12.nb_atomes);
	
	int tab[g12.nb_atomes];
	for(i=0;i < g12.nb_atomes ; i++)
	{
		tab[i] = dans_clique_max[i];
	}
	free(dans_clique_max);
	
	for(i=0;i < g12.nb_atomes - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j=i+1;j < g12.nb_atomes ; j++)
			{
				if(tab[j]==1 && (couple_atome[i].a1 == couple_atome[j].a1))
						tab[j] = 0;
			}
		}
	}

	for(i=0;i < g12.nb_atomes ; i++)
	{
		if(tab[i] ==1)
			nb_at++;
	}

	for(i=0;i < g12.nb_atomes - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j = i+1;j < g12.nb_atomes; j++)
			{
				if(tab[j] == 1)
				{
					
					i1 = couple_atome[i].a1;
					j1 = couple_atome[j].a1;
					if(M[pos1].matrice_liaisons[i1][j1] != AUCUNE_LIAISON)
						nb_liaisons ++;
				}
			}
		}
	}
	free(couple_atome);
	
	clique.nb_liaisons = nb_liaisons;
	clique.nb_atomes= nb_at;
	return clique;
	
}


struct graphe graphe_g12_cycles(struct graphe g12, struct molecule *M, int g1_chebi, int g2_chebi)
{ //contruction du graphe commun

	struct graphe clique;
	int nb_at =0,nb_liaisons=0,i,j,i1,j1;
	int pos1,pos2;
	pos1 =position_M(g1_chebi,M);
	pos2 =position_M(g2_chebi,M);
	
	struct couple *couple_atome = construction_couples_cycles(M,pos1,pos2,g12.nb_sommets);
	
	int tab[g12.nb_sommets];
	for(i=0;i < g12.nb_sommets ; i++)
	{
		tab[i] = dans_clique_max[i];
	}
	free(dans_clique_max);
	
	for(i=0;i < g12.nb_sommets - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j=i+1;j < g12.nb_sommets ; j++)
			{
				if(tab[j]==1 && (couple_atome[i].a1 == couple_atome[j].a1))
						tab[j] = 0;
			}
		}
	}

	for(i=0;i < g12.nb_sommets ; i++)
	{
		if(tab[i] ==1)
			nb_at++;
	}

	for(i=0;i < g12.nb_sommets - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j = i+1;j < g12.nb_sommets; j++)
			{
				if(tab[j] == 1)
				{
					
					i1 = couple_atome[i].a1;
					j1 = couple_atome[j].a1;
					if(M[pos1].g.matrice_cycles_poids[i1][j1] != AUCUNE_LIAISON)
						nb_liaisons ++;
				}
			}
		}
	}
	free(couple_atome);
	
	clique.nb_arete = nb_liaisons;
	clique.nb_sommets = nb_at;
	return clique;
	
}



void liberer_molecule(struct molecule g) 
{ //liberation de l'espace memoire d'une molecule
	if (g.liste_atomes) free(g.liste_atomes);
	if (g.liste_liaisons) free(g.liste_liaisons);
	if (g.matrice_liaisons)
	{	
		int i;
		for (i=0 ; i<g.nb_atomes ; i++) 
		{
			free(g.matrice_liaisons[i]);
			
		}
		free(g.matrice_liaisons);
	}

	int i;
	//printf("%d g_def \n", g.g_def);
	if(g.g_def == 1)
	{
		/*if(g.g.nb_sommets !=0)
		{*/
		free(g.g.som);
		free(g.g.aretes);
		if(g.g.matrice_cycles_type != NULL)
		{
			for( i = 0; i < g.g.nb_sommets; i++)
				free(g.g.matrice_cycles_type[i]);
			free(g.g.matrice_cycles_type);
		}
		if(g.g.matrice_cycles_poids != NULL)
		{
			for( i = 0; i < g.g.nb_sommets; i++)
				free(g.g.matrice_cycles_poids[i]);
			free(g.g.matrice_cycles_poids);
		}
		//}
	}
}

// cette fonction est valide pour le classement des molecules par dossier 
/*struct molecule * lecture_fichier_chebi()
{//lecture du fichier chebi.sdf
	
	FILE *F;
	F = fopen("../CHEBI/liste.num.chebi","r");
	
	if ( F == NULL ) 
	{
		 fprintf(stderr,"Cannot open ../CHEBI/liste.num.chebi file\n"); 
		 exit(18); 
	}
	init_atom_num();
	int nb_mol, DEB = 0, FIN = NB_MOLECULES;
	struct molecule *M = malloc((NB_MOLECULES +1)*sizeof(struct molecule));
	
	if (M == NULL)
	{
		fprintf(stderr,"Not enough memory for M\n"); 
		exit(3); 
	}
	struct molecule m;
	printf("1. Lecture des molecules : %.3lf s\n",chrono());
	FILE *G;
	char nom[128];
	int mol_courante;
	for(nb_mol = DEB ; nb_mol < FIN ; nb_mol++)
	{
		fscanf(F,"%d",&mol_courante);
		sprintf(nom,"../CHEBI/%d/%d.sdf",mol_courante,mol_courante);
		G = fopen(nom,"r");
		
		if( G == NULL)
		{
			fprintf(stderr,"Cannot open %s file\n",nom); 
			exit(19); 
			
		}
		if (nb_mol % 1000 == 0) 
		{ 
			fprintf(stdout,"\r %5d / %d",nb_mol,FIN);
			fflush(stdout); 
		}
		m = lire_molecule(G);
		fclose(G);
		M[nb_mol] = m;
	}
	
	fclose(F);
	fprintf(stdout,"\r%5d / %d\n",nb_mol,FIN); 
	
	return M;
	
}*/
struct molecule * lecture_fichier_chebi()
{//lecture du fichier chebi.sdf
	
	FILE *F;
	F = fopen("ChEBI_lite.sdf","r");
	
	if ( F == NULL ) 
	{
		 fprintf(stderr,"Cannot open ChEBI_lite.pdf file\n"); 
		 exit(1); 
	}
	init_atom_num();
	int nb_mol, DEB = 0, FIN = NB_MOLECULES;
	struct molecule *M = malloc(NB_MOLECULES * sizeof(struct molecule));
	
	if (M == NULL)
	{
		fprintf(stderr,"Not enough memory for M\n"); 
		exit(3); 
	}
	struct molecule m;
	printf("1. Lecture des molecules : %.3lf s\n",chrono());

	for(nb_mol = DEB ; nb_mol < FIN ; nb_mol++)
	{
		if (nb_mol % 1000 == 0) 
		{ 
			fprintf(stdout,"\r%5d / %d",nb_mol,FIN);
			fflush(stdout); 
		}
		m = lire_molecule(F);
		M[nb_mol] = m;
	}
	
	fclose(F);
	fprintf(stdout,"\r%5d / %d\n",nb_mol,FIN); 
	printf("Fin de la Lecture des molecules : %.3lf s\n",chrono());
	return M;
	
}

float mesure_similarite (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille_limite)
{//calcul du degré de similarité
	
	float similarite;
	int pos1,pos2;
	pos1 = position_M(g1_chebi,M);
	pos2 = position_M(g2_chebi,M);
	
	struct molecule g12 = graphe_produit(g1_chebi,g2_chebi,M);
	
	if( taille_limite != 0 && ( g12.nb_atomes > taille_limite))
	{
		
		return -2;
	}
	else
	{
		
		la_clique_max(g12,date);
		
		struct molecule clique= graphe_g12(g12,M,g1_chebi,g2_chebi);
		
		
		if(date == 0 || (chrono() - last_chrono <= date))
		{
			float num = (float)((clique.nb_atomes + clique.nb_liaisons)*(clique.nb_atomes + clique.nb_liaisons));
			float denum = (float)((M[pos1].nb_atomes + M[pos1].nb_liaisons)*(M[pos2].nb_atomes + M[pos2].nb_liaisons));
			similarite = num/denum;
		}
		else
		{
			similarite = -1;
		}
	}
	
	//liberer_molecule(g12);

	return similarite;
}

void similarite_all(int g1_chebi,struct molecule *M,double date,int taille_limite)
{
	int i;
	char nom[64];
	char nom2[64];
	sprintf(nom2,"%d",g1_chebi);
	strcpy(nom, "resultats/similatite_");
	strcat(nom,nom2);
	strcat(nom,"_all.data");
	FILE *F;
	F = fopen( nom, "w");
	if(F == NULL){
		fprintf(stdout,"Impossible de creer le fichier de resultat\n");
		exit(45);
	}
	
	float r;
	
	for ( i = 0;  i < NB_MOLECULES;  i++)
	{
		last_chrono = chrono();
		if (i % 1 == 0) 
		{ 
			fprintf(stdout,"\r%5d / %d (%3d atomes) %.3lf ",i,NB_MOLECULES,M[i].nb_atomes,last_chrono);
			fflush(stdout); 
		}
		r = mesure_similarite( g1_chebi, M[i].chebi_id,M, date, taille_limite);
		fprintf(F, "%d %f\n",M[i].chebi_id,r);
		fflush(F);
	}
	
	fclose(F);
}


// fonctions pour horton


struct liste_voisins*  construction_voisins_mol( int position, struct molecule *M){
	
	struct liste_voisins* voisins;
	int i,j,k,a1, a2;
	int nb_sommets = M[position].nb_atomes;
	int tableau_nb_voisins[nb_sommets];
	
//initialisation du tableau à 0
	for(i = 0; i < nb_sommets; i++)
		tableau_nb_voisins[i] = 0;
//remplissage	
	for(i = 0; i < M[position].nb_liaisons; i++)
	{
		a1 = M[position].liste_liaisons[i].A1;
		a2 = M[position].liste_liaisons[i].A2;
		tableau_nb_voisins[ a1 -1] ++;
		tableau_nb_voisins[ a2 -1] ++;
		
	
	}
//allocation memoire
	voisins = malloc( nb_sommets * sizeof(struct liste_voisins));
	if( voisins == NULL)
	{
		fprintf(stdout,"Cannot allocate memory in function construire_voisins_mol\n");
		exit(60);
	}
	for(i = 0; i < nb_sommets; i++)
	{
		voisins[i].id_atome = i + 1 ;
		voisins[i].nb_voisins = tableau_nb_voisins[i];
		
		if( voisins[i].nb_voisins  > 0 )
		{
			voisins[i].id_voisins = malloc(voisins[i].nb_voisins * sizeof(int));
			
			if(voisins[i].id_voisins == NULL)
			{
				fprintf(stdout,"Cannot allocate memory (second call) infunction construire_voisins_mol\n");
				exit(61);
			}
			k = 0;
			
			for(j = 0; j < M[position].nb_liaisons; j++)
			{
				if( M[position].liste_liaisons[j].A1 == voisins[i].id_atome )
				{
					voisins[i].id_voisins[k] = M[position].liste_liaisons[j].A2;
					k++;
				}
				
				else if ( M[position].liste_liaisons[j].A2 == voisins[i].id_atome )
				{
					voisins[i].id_voisins[k] = M[position].liste_liaisons[j].A1;
					k++;
				}
				
				
			}
		} 
	}
	//affichage_liste_voisinage(voisins,position,nb_sommets,M);
	return voisins;
}

void affichage_liste_voisinage(struct liste_voisins* voisins,int position,int nb_sommets,struct molecule *M)
{
	printf("Nombre d'atomes est de %d \n",nb_sommets);
	
	int i,j; 
	for( i = 0; i < M[position].nb_atomes; i++)
	{
			printf("%3d \t ",i+1);
			
			if(voisins[i].nb_voisins > 0)
			{
				for(j = 0; j < voisins[i].nb_voisins; j++)
				{
					printf("%3d ",voisins[i].id_voisins[j]);
				}
			}
			else
			{
				printf(" - ");
			}
			printf(" \n");
	}
	
}

//verifie si l'intersection de deux plus court chemin est egal au sommet i
//retoune 1 si le pcchemin est i et 0 sinon

int intersection_pcc(int *d1, int *d2,int sommet)
{
	int taille1,taille2,i,j;

	if(d1[0] == -1 || d2[0] == -1)
		return 0;
	taille1 = d1[0] + 1;
	taille2 = d2[0] + 1;
	
	int bool = 1;
	//printf("tailles : %d  %d\n",taille1,taille2);
	for( i = 1; i <= taille1  ; i++)
	{
		for(j = 1; j <= taille2; j++)
		{
			//printf("%d  %d %d \n",d1[i],d2[j],sommet);
			if((d1[i] == d2[j]) && ( d1[i] != sommet))
			{
				bool = 0;
				//printf("echec\n");
				break;
			}
		}
		if( bool == 0)
			break;
	}
	
	return bool;
	
}

void calcul_distance_sommets(int i,int position,struct molecule *M,struct liste_voisins *v,int ***vecteur_d){
	
	
	int *d,*predecesseur;
	int taille = M[position].nb_atomes;
	d = calloc(taille ,sizeof(int));
	predecesseur = calloc(taille ,sizeof(int));
	int j;
	for( j = 0; j < taille; j++)
	{
		d[j] = -1;
		predecesseur[j] = -1;
	}
	d[i-1] = 0;
	predecesseur[i-1] = i;

	int sommets[taille];
	int sommets2[taille];
	
	for(j = 0; j < taille ; j++)
	{
		sommets[j] = 0;
		sommets2[j] = 0;
	}
		
	int k,z,t;
	//premiere etape
	for( j = 0; j < v[i-1].nb_voisins; j++)
	{
		k = v[i-1].id_voisins[j] -1 ;
		d[k] = 1;
		predecesseur[k] = i;
		sommets[j] = k +1;
	
	}
	
	
	int *d_last = calloc(taille,sizeof(int));
		//initialisation 
	for ( j = 0; j < taille; j++)
		d_last[j] = 0;

	
	//parcours de tous les autres sommets
	int en_cours = v[i-1].nb_voisins;	
	int profondeur = 1;
	while(verification_egalite_tableaux(d_last,d,taille) !=1)
	{
		//sauvegarde de l'ancienne distance 
		for ( j = 0; j < taille; j++)
			d_last[j] = d[j];
		
		t = 0;
		for(k = 0; k < en_cours; k++)
		{
				for( j = 0; j < v[sommets[k]-1].nb_voisins; j++)
				{
					
					z = v[sommets[k]-1].id_voisins[j] -1 ;
					if(d[z] < 0 )
					{
						d[z] = profondeur + 1;
						predecesseur[z] = sommets[k];
						sommets2[t] = z+1;
						t++;
					}
					else if(d[z] != 0 && (d[z] > profondeur + 1))
					{
						d[z] = profondeur + 1;
						predecesseur[z] = sommets[k];
						sommets2[t] = z+1;
						t++;
					}
				}
		}
		
		for( j = 0; j < taille; j++)
			sommets[j]= sommets2[j];
		en_cours = t;
		profondeur++;
	}

	
	int val;
	int pos;
	int tailles;
	for ( j = 0; j < taille; j++)
	{
		
		
		if(d[j] <0)
			tailles = 2;
		else
			tailles = d[j] +2;

		vecteur_d[i-1][j] = calloc(tailles ,sizeof(int));
		//if(M[position].chebi_id == 3243)
			//printf("taille = %d d[j]= %d\n",taille,d[j] );
		pos = 0;
		k = j;
		val = d[j];
		vecteur_d[i-1][j][pos] = val;
		pos = pos +1;
		vecteur_d[i-1][j][pos] = j+1;
		pos++;
		
		while(predecesseur[k] != i && val > 0)
		{
			vecteur_d[i-1][j][pos] = predecesseur[k];
			k = predecesseur[k] -1;
			val--;	
			pos++;
		}
		
		if(predecesseur[k] == i && val == 1)
		{	vecteur_d[i-1][j][pos] = predecesseur[k];
		}
	}
	free(d);
	free(predecesseur);
	free(d_last);
}

int ***calcul_distance_sommets_all(int position,struct molecule *M,struct liste_voisins *v)
{
	int taille = M[position].nb_atomes;
	
	int i; 
	int ***vecteur_d = malloc(taille * sizeof(int **));
	
	for (i = 0; i < taille ; i++)
		vecteur_d[i] = malloc(taille * sizeof(int *));
	
	for (i = 0; i < taille ; i++)
		calcul_distance_sommets(i+1,position,M,v,vecteur_d);
//affichage des plus chemins
/*
	int j,nb,k;
	for (i = 0; i < taille ; i++)
	{
		
		for ( j = 0; j < taille; j++)
		{
			nb = vecteur_d[i][j][0];
			printf(" le plus court chemin entre %d  et %d de distance %d :\n",i+1,j+1,nb);
			
			for ( k = 1; k <= nb +1; k++)
			{
				printf("%d ",vecteur_d[i][j][k]);
			}
			printf("\n");
		}
	
	}
*/
//fin d'affichage


	return vecteur_d;
}
int* vecteur_distance(int ***vecteur_d,int a, int b,struct vecteur v)
{
	int taille1 = vecteur_d[b][a][0];
	int taille_vecteur = v.taille;
	int *c = malloc(taille_vecteur *sizeof(int));
	int k,i,j;
	for ( k =0; k < taille_vecteur; k++)
		c[k] = 0;
		
	int nb= 0;
	while( nb <= taille1-1)
	{	
		i = vecteur_d[b][a][nb + 1];
		j = vecteur_d[b][a][nb + 2];
		
		for ( k =0; k < taille_vecteur; k++)
		{
			if((v.sommets[k].a1 == i && v.sommets[k].a2 == j) || (v.sommets[k].a1 == j && v.sommets[k].a2 == i))
				break;
		}
		
		if ( k == taille_vecteur)
		{
			fprintf(stdout,"La liaison n'a pas été trouvé fonction vecteur_cycle\n");
			exit(100);
		}
		//printf("%d %d liaison a pos %d\n",i,j,k);
		c[k] = 1;
		nb++;
	}
	
	return c;
}

int* vecteur_cycle(int ***vecteur_d, int a, int b, int sommet,struct vecteur v,struct couple cple)
{
	//printf("sommet %d et %d et de taille %d\n",sommet+1, a+1, vecteur_d[sommet][a][0]);
	int taille1 = vecteur_d[sommet][a][0];
	int taille2 = vecteur_d[sommet][b][0];
	
	int taille_vecteur = v.taille;
	int i,j,k,nb = 0;
	int *c = malloc(taille_vecteur *sizeof(int));
	
	for ( k =0; k < taille_vecteur; k++)
		c[k] = 0;
	
	while( nb <= taille1-1)
	{	
		i = vecteur_d[sommet][a][nb + 1];
		j = vecteur_d[sommet][a][nb + 2];
		
		for ( k =0; k < taille_vecteur; k++)
		{
			if((v.sommets[k].a1 == i && v.sommets[k].a2 == j) || (v.sommets[k].a1 == j && v.sommets[k].a2 == i))
				break;
		}
		
		if ( k == taille_vecteur)
		{
			fprintf(stdout,"La liaison n'a pas été trouvé fonction vecteur_cycle\n");
			exit(100);
		}
		//printf("%d %d liaison a pos %d\n",i,j,k);
		c[k] = 1;
		nb++;
	}
	nb = 0;
	while( nb <= taille2-1)
	{	
		i = vecteur_d[sommet][b][nb + 1];
		j = vecteur_d[sommet][b][nb + 2];
		
		for ( k =0; k < taille_vecteur; k++)
		{
			if((v.sommets[k].a1 == i && v.sommets[k].a2 == j) || (v.sommets[k].a1 == j && v.sommets[k].a2 == i))
				break;
		}
		
		if ( k == taille_vecteur)
		{
			fprintf(stdout,"La liaison n'a pas été trouvé fonction vecteur_cycle\n");
			exit(100);
		}
		//printf("%d %d liaison a pos %d\n",i,j,k);
		c[k] = 1;
		nb++;
	}
	i =cple.a1 ;
	j = cple.a2;
	for ( k =0; k < taille_vecteur; k++)
	{
			if((v.sommets[k].a1 == i && v.sommets[k].a2 == j) || (v.sommets[k].a1 == j && v.sommets[k].a2 == i))
				break;
	}
	c[k] = 1;

	return c;
}


void affichage_pcc(int ***vecteur_d, int i , int j)
{
	int k,nb;
	nb = vecteur_d[i][j][0];
	
	for ( k = 1; k <= nb +1; k++)
	{
		printf("%d ",vecteur_d[i][j][k]);
	}
	//printf("\n");
}

//verifie si deux tableaux sont les memes
int verification_egalite_tableaux(int *tab1,int *tab2 , int taille)
{
	
	int i, verif = 1 ;
	for( i = 0; i < taille ; i++)
	{
		if( tab1[i] != tab2[i])
		{
			verif = 0;
			break;
		}
	}
	
	return verif;
	
}


struct vecteur construction_vecteur(struct molecule *M,int position)
{
	struct vecteur v;
	v.taille = M[position].nb_liaisons;
	v.sommets = malloc(v.taille *sizeof(struct couple));
	int i;
	
	for( i = 0; i < v.taille ; i++) 
	{
		v.sommets[i].a1 = M[position].liste_liaisons[i].A1;
		v.sommets[i].a2 = M[position].liste_liaisons[i].A2;
//printf("%d %d \n",v.sommets[i].a1,v.sommets[i].a2);
	}
	
	return v;
	
}
int position_liaison_vecteur( struct vecteur v,int a1,int a2)
{
	int c = -1;
	int i;
	for ( i = 0; i< v.taille; i++)
	{
		if(((v.sommets[i].a1 == a1) && (v.sommets[i].a2 ==a2)) ||((v.sommets[i].a1 == a2) && (v.sommets[i].a2 ==a1)))
		{
			c = i;
			break;
		}
	}	
	return c;
}

int calcul_poids(int *v, int taille)
{
	int i,poids= 0; 
	
	for (i = 0; i < taille ; i++)
		poids += v[i];
	
	return poids;
}

int op_xor(int a , int b)
{
	if(( a ==0 && b ==0) ||(b ==1 && a==1))
		return 0;
	else
		return 1;
	
}

int **fonction_xor( int **matrice ,int taille_y, int i , int j)
{
	int k;
	
	for (k = 0; k < taille_y; k++)
		matrice[j][k] = op_xor(matrice[i][k] ,matrice[j][k]);
	
	return matrice;
}

int premier_un(int **matrice,int taille_y, int i)
{
	int pos = -1;
	int j;
	
	for ( j  = 0; j < taille_y; j++)
	{
		if(matrice[i][j] == 1)
		{
			pos =j;
			break;
		}
	}
	
	return pos;
	
}

//construction du graphe de cycles

void construction_graphe_de_cycles(int position,struct molecule *M)
{

	
	struct liste_voisins* v;
	v =  construction_voisins_mol( position,M);
	int ***dist;
	dist = calcul_distance_sommets_all(position,M,v);
	struct vecteur vect = construction_vecteur(M,position);
	int nb_sommets = M[position].nb_atomes;
	int *c;
	int sommet;
	int i,x,y,t;
	int pcc;
	struct couple cple;
	struct cycle cycles[16384];
	int k,p;
	int nb_cycles =0;
	int doublon ;
	for(sommet = 2; sommet <= M[position].nb_atomes; sommet++)
	{
		for (i =0 ;  i < M[position].nb_liaisons;i++)
		{
			if( (M[position].liste_liaisons[i].A1 !=sommet) && ( M[position].liste_liaisons[i].A2 !=sommet))
			{
				x = M[position].liste_liaisons[i].A1;
				y = M[position].liste_liaisons[i].A2;
				cple.a1 = M[position].liste_liaisons[i].A1;
				cple.a2 = M[position].liste_liaisons[i].A2;
				doublon = 0;	
				pcc = intersection_pcc(dist[sommet -1][x-1],dist[sommet -1][y-1],sommet);

				if(pcc == 1)
				{
					c =vecteur_cycle(dist, x-1, y-1,sommet - 1,vect,cple) ;
					
					for (k = 0; k < nb_cycles ; k++)
					{
						
						if((verification_egalite_tableaux(cycles[k].c,c,M[position].nb_liaisons) == 1)||(calcul_poids(c,M[position].nb_liaisons) <2))
						{
							doublon =1;
							break;
						}
					}
					
					if(doublon == 0)
					{
						cycles[nb_cycles].poids = calcul_poids(c,M[position].nb_liaisons);
						cycles[nb_cycles].c = malloc(M[position].nb_liaisons * sizeof(int));
						for(p = 0; p < M[position].nb_liaisons; p++)
							cycles[nb_cycles].c[p] = c[p];
						nb_cycles++;
						
					}
					free(c);
				}
			}
			
		}
	}
	int r = 0;

		
	M[position].g.nb_sommets = 0;
	M[position].g.nb_arete = 0;
	if ( nb_cycles >= 1)
	{
		int pos = 0;
		int min = cycles[0].poids;
		int max = cycles[0].poids;
		struct cycle lescycles[nb_cycles];
		
		for (k = 1; k < nb_cycles ; k++)
		{	
			if(cycles[k].poids > max)
			{
				max = cycles[k].poids;
			}
			if(cycles[k].poids < min)
			{
				min = cycles[k].poids;
				pos = k;
			}
		}
		
		while ( r < nb_cycles)
		{	
			lescycles[r].poids = cycles[pos].poids;
			lescycles[r].c = cycles[pos].c;
			cycles[pos].poids = 2 * max;
			r++;
			min = 2*max;
			
			for (k = 0; k < nb_cycles ; k++)
			{
				if(cycles[k].poids < min)
				{
					min = cycles[k].poids;
					pos = k;
				}
			}
		}
		
		
		int j,l;
		int **matrice = malloc(nb_cycles * sizeof(int * ));
		
		for ( i = 0; i < nb_cycles; i++)
			matrice[i] = malloc(M[position].nb_liaisons * sizeof(int));
		
		
		for ( i = 0; i < nb_cycles; i++)
		{
			for( j = 0; j < M[position].nb_liaisons; j++)
			{
				matrice[i][j] = lescycles[i].c[j];
			}
		}
		for ( i = 0; i < nb_cycles; i++)
		{
			pos = premier_un(matrice,M[position].nb_liaisons, i);

			if (pos >= 0)
			{
				for( j = i+1; j <nb_cycles; j++)
				{
					if(matrice[j][pos] == 1)
						matrice = fonction_xor(matrice,M[position].nb_liaisons,i,j);
						
				}
			}

		}
		int cycles_ind[nb_cycles],somme;
		for (i = 0; i < nb_cycles; i++)
		{
			somme  = 0;
			for(j = 0; j< M[position].nb_liaisons;j++)
				somme = somme + matrice[i][j];
			
			if ( somme > 0)
				cycles_ind[i] = 1;
			else 	
				cycles_ind[i] = 0;
		}
		
		int taille_max = 0;
		for (i = 0; i < nb_cycles; i++)
		{
			if(cycles_ind[i] == 1 && lescycles[i].poids > taille_max)
	
				taille_max = lescycles[i].poids ;

		}
		
		for (i = 0; i < nb_cycles; i++)
		{
			if(cycles_ind[i] == 0 && lescycles[i].poids <= taille_max)
			
				cycles_ind[i] = 1;
		}
		
		for ( i = 0; i < nb_cycles; i++)
			free(matrice[i]);
			
		free(matrice);
	//ok here
		int nb_l = M[position].nb_liaisons;
		int arete_cycle[nb_l];
		
		
		for(j = 0; j<nb_l;j++)
			arete_cycle[j] = 0;
		
		
		for( i = 0; i < nb_cycles; i++)
		{
			if(cycles_ind[i] == 1)
			{
				for(j = 0; j< nb_l; j++)
				{
					if( lescycles[i].c[j] == 1)
						arete_cycle[j] = 1;
				}
			}
		}


		
		for (i = 0; i < nb_cycles; i++)
			M[position].g.nb_sommets += cycles_ind[i];
	
		int compt =0;
		M[position].g.som = malloc(M[position].g.nb_sommets * sizeof(lesommet));
		if (M[position].g.som == NULL)
		{
			printf("cannot allocate memory for g.sommets\n");
			exit(101);
		}
		
		for (i = 0; i < nb_cycles; i++)
		{
			
			if(cycles_ind[i] == 1)
			{
				M[position].g.som[compt].id = i + 1;
				M[position].g.som[compt].taille = lescycles[i].poids;
				//printf("%d %d %d\n",compt,M[position].g.som[compt].id,M[position].g.som[compt].taille);
				compt++;
			}
		}
		//printf("end \n");

	int a1, a2,b1,b2,w;
		int *c;
		int bool;
		struct couplet aretes_graphe[16384];
		int nb_actu =0;
		for( i = 0; i < nb_cycles; i++)
		{
			for( j = i+1; j < nb_cycles; j++)
			{
				if( (cycles_ind[i] == cycles_ind[j]) &&(cycles_ind[j] == 1))
				{
					pcc = -1;
					
					for( k = 0; k < M[position].nb_liaisons; k++) //k ieme liaison dans i
					{
						for( l = 0; l < M[position].nb_liaisons; l++) //l ieme laison dans j
						{ 
							
							
							if( (lescycles[i].c[k] == 1) && (lescycles[i].c[k] == lescycles[j].c[l]))
							{
								a1 = vect.sommets[k].a1;
								a2 = vect.sommets[k].a2;
								b1 = vect.sommets[l].a1;
								b2 = vect.sommets[l].a2;
								
								if( ((pcc ==-1)&&(dist[a1-1][b1-1][0] >=0)) ||((pcc>=0 )||(dist[a1-1][b1-1][0] < pcc ) ))
								{
									
									c = vecteur_distance(dist,a1-1,b1-1,vect);
									bool = 1;
									for ( t = 0; t < vect.taille; t++)
									{
										if((c[t] ==1) &&(c[t] == arete_cycle[t]))
										{
											bool =0;
											break;
										}
									}
									if( bool == 1)
									{
										pcc = dist[a1-1][b1-1][0];
									}
									free(c);
								}
								//a1 et b2
								if( ((pcc ==-1)&&(dist[a1-1][b2-1][0] >=0)) ||((pcc>=0 )||(dist[a1-1][b2-1][0] < pcc ) ))
								{
									
									c = vecteur_distance(dist,a1-1,b2-1,vect);
									bool = 1;
									for ( t = 0; t < vect.taille; t++)
									{
										if((c[t] ==1) &&(c[t] == arete_cycle[t]))
										{
											bool =0;
											break;
										}
									}
									if( bool == 1)
									{
										//printf("%d %d %d\t",a1,b2,dist[a1-1][b2-1][0]);
										pcc = dist[a1-1][b2-1][0];
									}
								
									free(c);
								}
								//a2 et b1
								if( ((pcc ==-1)&&(dist[a2-1][b1-1][0] >=0)) ||((pcc>=0 )||(dist[a2-1][b1-1][0] < pcc ) ))
								{
									
									c = vecteur_distance(dist,a2-1,b1-1,vect);
									bool = 1;
									for ( t = 0; t < vect.taille; t++)
									{
										if((c[t] ==1) &&(c[t] == arete_cycle[t]))
										{
											bool =0;
											break;
										}
									}
									if( bool == 1)
									{
										//printf("%d %d %d\t",a2,b1,dist[a2-1][b1-1][0]);
										pcc = dist[a2-1][b1-1][0];
									}
									free(c);
								}
								//a2 et b2
								if( ((pcc ==-1)&&(dist[a2-1][b2-1][0] >=0)) ||((pcc>=0 )||(dist[a2-1][b2-1][0] < pcc ) ))
								{
									
									c = vecteur_distance(dist,a2-1,b2-1,vect);
									bool = 1;
									for ( t = 0; t < vect.taille; t++)
									{
								
										if((c[t] ==1) &&(c[t] == arete_cycle[t]))
										{
											bool =0;
											break;
										}
									}
									if( bool == 1)
									{
										//printf("%d %d %d\t",a2,b2,dist[a2-1][b2-1][0]);
										pcc = dist[a2-1][b2-1][0];
									}
									//printf("\n");
									free(c);
								}
								
								
							}
						}
					}
					if(pcc >= 0)
					{
						M[position].g.nb_arete ++;
						aretes_graphe[nb_actu].a1 = i+1;
						aretes_graphe[nb_actu].a2 = j+1;
						
						somme = 0;
						for (w = 0; w < nb_l; w++)
						{
							if((lescycles[i].c[w] == 1)&&(lescycles[i].c[w] == lescycles[j].c[w]))
								somme++;
						}
						//printf("la somme est %d\n",somme);
						if(somme > 0)
						{
							aretes_graphe[nb_actu].poids = somme;
							aretes_graphe[nb_actu].type = 0;

						}
						else
						{
							if(somme ==0 && pcc == 0)
							{
								//printf("mol %d\n", M[position].chebi_id);
								aretes_graphe[nb_actu].poids = pcc;
								aretes_graphe[nb_actu].type = 0;							
							}
							else
							{
								aretes_graphe[nb_actu].poids = pcc;
								aretes_graphe[nb_actu].type = 1;
							}
						}
						nb_actu++;
					
					
					
					}
					//fin de la fonction
				}
				
			}
		}
		
		M[position].g.aretes = malloc( M[position].g.nb_arete * sizeof(couplet));
	
		for (i = 0; i< M[position].g.nb_arete; i++)
		{
			M[position].g.aretes[i].a1 = aretes_graphe[i].a1;
			M[position].g.aretes[i].a2 = aretes_graphe[i].a2;
			M[position].g.aretes[i].type = aretes_graphe[i].type;
			M[position].g.aretes[i].poids = aretes_graphe[i].poids;
		}
		M[position].g_def = 1;	

	}	
	//liberation de la memoire
	//if(M[position].chebi_id == 3385)
		//printf("nb sommets %d\n",nb_sommets );
	for (i = 0; i < nb_sommets; i++)
	{
		if(v[i].nb_voisins > 0)
			free(v[i].id_voisins);
	}
	free(v);
	//int max,j;
	int j;
	for (i = 0; i < nb_sommets ; i++)
	{

		for (j = 0; j < nb_sommets; j++)
		{
			free(dist[i][j]);
		}
		free(dist[i]);
	}
		
	free(dist);
	free(vect.sommets);
	for (k = 0; k < nb_cycles ; k++)
		free(cycles[k].c);
	
}

void affichage_graphe_cycle_molecule( int position, struct molecule *M)
{
	printf("affichage du graphe de molecule\n");
	printf("nb sommet %d et nb aretes %d\n",M[position].g.nb_sommets,M[position].g.nb_arete);
	printf("\n Liste des sommets et leurs poids\n {");
	int i;
	for( i = 0; i < M[position].g.nb_sommets; i++)
		printf("(%d,%d),",M[position].g.som[i].id, M[position].g.som[i].taille);
	printf("}\n");
	printf("\n Liste des aretes et leur etiquettes et types \n {");	
	for( i = 0; i < M[position].g.nb_arete; i++)
		printf("(%d,%d,%d,%d),",M[position].g.aretes[i].a1,M[position].g.aretes[i].a2,M[position].g.aretes[i].type,M[position].g.aretes[i].poids);
	printf("}\n");
}


float mesure_similarite_cycles (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille)
{
	last_chrono = chrono();
	float similarite = 0;
	int pos1,pos2;
	int i;
	pos1 = position_M(g1_chebi,M);
	pos2 = position_M(g2_chebi,M);
	
	if (pos1 ==-1 || pos2 == -1)
		return -3;
	if( M[pos1].g_def == 0)
	{
		construction_graphe_de_cycles(pos1,M);
		M[pos1].g_def = 1 ;
	}
	
	if(M[pos2].g_def == 0)
	{
		construction_graphe_de_cycles(pos2,M);
		M[pos2].g_def = 1;
	}
	//construction du graphe produit

	//printf("%d %d \n",M[pos1].g.nb_sommets,M[pos2].g.nb_sommets );
	if (M[pos1].g.nb_sommets == 0 || M[pos2].g.nb_sommets == 0 )
	{
		similarite = -1;
		M[pos2].g_def = 0;
	}
	else
	{
		struct graphe g12 = graphe_produit_cycles(g1_chebi,g2_chebi,M);
		if(g12.nb_sommets > taille)
		{
			similarite = -3;
		}
		else
		{
			la_clique_max_cycles(g12,date);
			if(date == 0 || (chrono() - last_chrono <= date))
			{
			
				struct graphe clique = graphe_g12_cycles(g12,M,g1_chebi,g2_chebi);
				//printf("nb sommets %d et nb_liaison %d\n",clique.nb_sommets,clique.nb_arete);
		
				float num = (float)((clique.nb_sommets + clique.nb_arete)*(clique.nb_sommets+ clique.nb_arete));
				float denum = (float)((M[pos1].g.nb_sommets + M[pos1].g.nb_arete)*(M[pos2].g.nb_arete + M[pos2].g.nb_sommets));
				similarite = num/denum;
			}
			else
			{
				similarite = -2;
			}
		}
		free(g12.som);
		free(g12.aretes);
		
		for (i = 0; i < g12.nb_sommets; i++) free(g12.matrice_cycles_type[i]);
		free(g12.matrice_cycles_type);
		free(M[pos2].g.som);
		free(M[pos2].g.aretes);


		for(i=0;i< M[pos2].g.nb_sommets ;i++)
			free(M[pos2].g.matrice_cycles_type[i]);
		free(M[pos2].g.matrice_cycles_type);

		for(i=0;i< M[pos2].g.nb_sommets ;i++)
			free(M[pos2].g.matrice_cycles_poids[i]);
		free(M[pos2].g.matrice_cycles_poids);
		
		M[pos2].g_def = 0;
	}

	return similarite;
}
float mesure_similarite_cycles_2 (int g1_chebi,int g2_chebi,struct molecule *M,double date,int taille)
{
	last_chrono = chrono();
	//printf("g2 chebi %d\n", g2_chebi);
	float similarite = 0;
	int pos1,pos2;
	int i;
	pos1 = position_M(g1_chebi,M);
	pos2 = position_M(g2_chebi,M);
	
	if (pos1 ==-1 || pos2 == -1)
		return -3;
	if( M[pos1].g_def == 0)
	{
		construction_graphe_de_cycles(pos1,M);
		M[pos1].g_def = 1 ;
	}
	
	if(M[pos2].g_def == 0)
	{
		construction_graphe_de_cycles(pos2,M);
		M[pos2].g_def = 1;
	}
	//construction du graphe produit

	//printf("%d %d \n",M[pos1].g.nb_sommets,M[pos2].g.nb_sommets );
	if (M[pos1].g.nb_sommets == 0 || M[pos2].g.nb_sommets == 0 )
	{
		similarite = -1;
		M[pos2].g_def = 0;
	}
	else
	{
		struct graphe g12 = graphe_produit_cycles(g1_chebi,g2_chebi,M);
		if(g12.nb_sommets > taille)
		{
			similarite = -3;
		}
		else
		{
			la_clique_max_cycles(g12,date);
			if(date == 0 || (chrono() - last_chrono <= date))
			{
			
				struct graphe clique = graphe_g12_cycles(g12,M,g1_chebi,g2_chebi);
				//printf("nb sommets %d et nb_liaison %d\n",clique.nb_sommets,clique.nb_arete);
		
				float num = (float)(clique.nb_sommets );
				int max;
				if (M[pos1].g.nb_sommets > M[pos2].g.nb_sommets)
					max = M[pos1].g.nb_sommets;
				else 
					max = M[pos2].g.nb_sommets;
				float denum = (float)(max);
				similarite = num/denum;
			}
			else
			{
				similarite = -2;
			}
		}
		free(g12.som);
		free(g12.aretes);
		
		for (i = 0; i < g12.nb_sommets; i++) free(g12.matrice_cycles_type[i]);
		free(g12.matrice_cycles_type);
		free(M[pos2].g.som);
		free(M[pos2].g.aretes);


		for(i=0;i< M[pos2].g.nb_sommets ;i++)
			free(M[pos2].g.matrice_cycles_type[i]);
		free(M[pos2].g.matrice_cycles_type);

		for(i=0;i< M[pos2].g.nb_sommets ;i++)
			free(M[pos2].g.matrice_cycles_poids[i]);
		free(M[pos2].g.matrice_cycles_poids);
		
		M[pos2].g_def = 0;
	}

	return similarite;
}
struct graphe graphe_produit_cycles(int g1_chebi, int g2_chebi, struct molecule *M)
{
	
	int i,j;
	int taille = 0;
	int pos1,pos2;
	pos1 =position_M(g1_chebi,M);
	pos2 =position_M(g2_chebi,M);
	
	for(i= 0; i < M[pos1].g.nb_sommets; i++)
	{ 
		for(j= 0; j < M[pos2].g.nb_sommets; j++)
		{
			if( M[pos1].g.som[i].taille == M[pos2].g.som[j].taille) 
				taille ++;
		}
	}
	//couple de cycles entre les nouveaux sommets
	struct couple * couple_cycles = construction_couples_cycles(M,pos1,pos2,taille);
	

	M[pos2]= construction_matrice_cycles(M[pos2]);
	M[pos1]= construction_matrice_cycles(M[pos1]);

	
	struct graphe g12;
	g12.nb_sommets = taille;
	g12.nb_arete = 0;

	//les sommets
	g12.som = malloc(g12.nb_sommets * sizeof(struct lesommet));

	if( g12.som == NULL)
	{
		fprintf(stderr, "Cannot allocate memory\n");	
		exit(84);
	}
	for(i= 0; i < taille ;i++)
	{
		g12.som[i].id = couple_cycles[i].a1;
		g12.som[i].taille = couple_cycles[i].a2;
	}
	
	//les aretes 
	//printf("affichage des aretes\n");
	/*for (i = 0; i < M[pos1].g.nb_sommets; i++)
	{
		for (j = 0; j < M[pos1].g.nb_sommets; j++)
	{
			printf("%d %d\n",M[pos1].g.matrice_cycles_poids[i][j],M[pos1].g.matrice_cycles_poids[i][j]);
		}
	}*/
	//printf("end \n");
	//compter le nombre d'aretes
	int i1,i2,j1,j2;
	int nb_ar = 0;
	
	//printf("not me taille = %d\n",taille);
	for(i= 0; i < taille -1 ;i++)
	{
		
		i1 = couple_cycles[i].a1;
		i2 = couple_cycles[i].a2;
		//printf("not me %d %d \n",i1,i2);
		for(j= i + 1 ; j < taille ;j++)
		{
			
			j1=couple_cycles[j].a1;
			j2=couple_cycles[j].a2;
			//printf("not me %d %d \n",j1,j2);
			//printf("%d %d %d %d \n",M[pos1].g.som[i1].id, M[pos1].g.som[j1].id,M[pos2].g.som[i2].id , M[pos2].g.som[j2].id);//,M[pos1].g.matrice_cycles_poids[i1][j1],M[pos2].g.matrice_cycles_poids[i2][j2],M[pos1].g.matrice_cycles_type[i1][j1],M[pos2].g.matrice_cycles_type[i2][j2]);
			if( (M[pos1].g.som[i1].id != M[pos1].g.som[j1].id) && (M[pos2].g.som[i2].id != M[pos2].g.som[j2].id) &&(M[pos1].g.matrice_cycles_poids[i1][j1] == M[pos2].g.matrice_cycles_poids[i2][j2]) && (M[pos1].g.matrice_cycles_type[i1][j1] == M[pos2].g.matrice_cycles_type[i2][j2]) )
			{

				nb_ar++;
			}
			
		}
				//printf("toto is ok here %d\n",i);
	}

	//printf("nb ar =%d\n",nb_ar);
	g12.aretes = malloc(nb_ar * sizeof(struct couplet));
	if( g12.aretes == NULL)
	{
		fprintf(stderr, "Canot allocate memory\n");	
		exit(85);
	}
	nb_ar = 0;
	
	//constructiion matrice de liaison de g
	g12.matrice_cycles_type = malloc(taille * sizeof(int *));
	for (i = 0;i < taille ; i++) 
		g12.matrice_cycles_type[i] = malloc(taille * sizeof(int ));
	

	for (i = 0;i < taille ; i++) 
	{
		for( j = 0; j< taille ; j++)
		g12.matrice_cycles_type[i][j] = 0; 
	}
	for( i= 0; i < taille -1 ;i++)
	{
		i1 = couple_cycles[i].a1;
		i2 = couple_cycles[i].a2;
		for(j= i + 1 ; j < taille ;j++)
		{
			
			j1 = couple_cycles[j].a1;
			j2 = couple_cycles[j].a2;
			
			if( (M[pos1].g.som[i1].id != M[pos1].g.som[j1].id) && (M[pos2].g.som[i2].id != M[pos2].g.som[j2].id) &&(M[pos1].g.matrice_cycles_poids[i1][j1] == M[pos2].g.matrice_cycles_poids[i2][j2]) && (M[pos1].g.matrice_cycles_type[i1][j1] == M[pos2].g.matrice_cycles_type[i2][j2]) )
			{

				g12.aretes[nb_ar].a1 = i;
				g12.aretes[nb_ar].a2 = j;
				g12.aretes[nb_ar].type = M[pos1].g.matrice_cycles_type[i1][j1];
				g12.aretes[nb_ar].poids = M[pos1].g.matrice_cycles_poids[i1][j1];
				nb_ar++;
				g12.matrice_cycles_type[i][j] = 1; 
				//printf("(%d, %d),",i,j);
			}
			
		}
	}
	//printf("\n");
	

	//stocker le nombre d'aretes
	free(couple_cycles);
	return g12;
}

void la_clique_max_cycles( struct graphe g,double date)
{	
	//Debut calcul de la clique -- Initialisation
	int i;
	int *candidat;
	int *dans_clique;
	
	dans_clique_max = malloc( g.nb_sommets *sizeof(int));
	if (!dans_clique_max) { fprintf(stderr,"cannot malloc dans_clique_max %d\n",g.nb_sommets); exit(41); }
	dans_clique = malloc( g.nb_sommets *sizeof(int));
	if (!dans_clique) { fprintf(stderr,"cannot malloc dans_clique %d\n",g.nb_sommets); exit(42); }
	candidat = malloc( g.nb_sommets *sizeof(int));
	if (!candidat) { fprintf(stderr,"cannot malloc candidat %d\n",g.nb_sommets); exit(43); }
	
	//initialisation 
	for(i = 0; i < g.nb_sommets;i++)
	{
		candidat[i] 	= 1;
		dans_clique[i]	= 0;
		dans_clique_max[i] = 0;
	}
	
	taille_clique_max = 0;
	
	calcul_cl_cycles(g,dans_clique,0,candidat,g.nb_sommets,date); // 0 taille de la clique initial  et m.nb_atome = nb sommets candidats
	 
	free(dans_clique);
	free(candidat);
	
}


void calcul_cl_cycles(struct graphe g,int *dans_clique,int taille_clique,int *candidat,int taille_candidat,double date)
{	//calcul de la clique max recursif
	if (date != 0)
	{	//{printf("%f %f\n",chrono(),last_chrono);return;}
		if(chrono() - last_chrono > date) {
		
		return;
		}
			
	}

	int i,j;
	if( taille_candidat == 0)
	{
		if(taille_clique > taille_clique_max){
			taille_clique_max = taille_clique;
			for (i = 0 ;  i < g.nb_sommets ; i ++)
				dans_clique_max[i] = dans_clique[i];
		}
		return;
	}
	
	if (taille_candidat + taille_clique <= taille_clique_max)
	{

		return;
	}
	
	//else 
	int taille_candidat_temp;
	int *candidat_temp;
	candidat_temp = malloc( g.nb_sommets * sizeof(int));



	for (i = 0 ;  i < g.nb_sommets ; i ++)
	{
		if ( candidat[i] == 1)
		{
			candidat[i] = 0;
			dans_clique[i] = 1 ;
			taille_candidat_temp = taille_candidat;
			
			for (j = 0 ;  j < g.nb_sommets ; j ++)
			{
				candidat_temp[j] = candidat[j]; 
				if ((candidat[j] == 1) && (g.matrice_cycles_type[i][j] == 0))
				{
					candidat_temp[j] = 0;
					taille_candidat_temp--;	
				}	
			}
			
			taille_candidat_temp --;

			calcul_cl_cycles(g,dans_clique,taille_clique + 1,candidat_temp,taille_candidat_temp,date);
			dans_clique[i] = 0;
			
		}	
		
	}
	free(candidat_temp);
		 
}


struct molecule construction_matrice_cycles(struct molecule m)
{//construction de la matrice de liaison d'une molecule
	
	int i,j;
	if(m.g.matrice_cycles_type == NULL)
	{
		m.g.matrice_cycles_type =  malloc(m.g.nb_sommets * sizeof(int *));
		
		for(i=0;i< m.g.nb_sommets ;i++) m.g.matrice_cycles_type[i] =  malloc(m.g.nb_sommets * sizeof(int));
		
		for(i=0;i< m.g.nb_sommets ;i++)
		{
			for(j=0;j< m.g.nb_sommets ;j++)
			{
				m.g.matrice_cycles_type[i][j] = AUCUNE_LIAISON;
			 }
		}

		for(i =0; i< m.g.nb_arete;i++)
		{
			m.g.matrice_cycles_type[m.g.aretes[i].a1-1][m.g.aretes[i].a2-1] = m.g.aretes[i].type;
			m.g.matrice_cycles_type[m.g.aretes[i].a2-1][m.g.aretes[i].a1-1] = m.g.aretes[i].type;
		}
		

	}
	
	if(m.g.matrice_cycles_poids == NULL)
	{
		m.g.matrice_cycles_poids =  malloc(m.g.nb_sommets * sizeof(int *));
		
		for(i=0;i< m.g.nb_sommets ;i++) m.g.matrice_cycles_poids[i] =  malloc(m.g.nb_sommets * sizeof(int));
		
		for(i=0;i< m.g.nb_sommets ;i++)
		{
			for(j=0;j< m.g.nb_sommets ;j++)
			{
				m.g.matrice_cycles_poids[i][j] = AUCUNE_LIAISON;
			 }
		}

		for(i =0; i< m.g.nb_arete;i++)
		{
			m.g.matrice_cycles_poids[m.g.aretes[i].a1-1][m.g.aretes[i].a2-1] = m.g.aretes[i].poids;
			m.g.matrice_cycles_poids[m.g.aretes[i].a2-1][m.g.aretes[i].a1-1] = m.g.aretes[i].poids;
		}
		

	}
	/*printf("matrice cycles %d \n",m.g.nb_sommets );
	for(i=0;i< m.g.nb_sommets ;i++)
		{
			for(j=0;j< m.g.nb_sommets ;j++)
			{
				printf("%d ",m.g.matrice_cycles_poids[i][j]);
			 }
			 printf("\n");
		}*/

	return m;
	
}

struct couple *construction_couples_cycles(struct molecule *M,int pos1, int pos2,int taille)
{//Construction des couples de cycles compatibles
	
	int n = 0,i,j;
	
	struct couple *couple_at;
	couple_at = malloc(taille * sizeof(couple));
	for(i= 0; i <  M[pos1].g.nb_sommets;i++)
	{ 
		for(j= 0; j <  M[pos2].g.nb_sommets;j++)
		{
			if( M[pos1].g.som[i].taille == M[pos2].g.som[j].taille)
			{
				couple_at[n].a1 = i;
				couple_at[n].a2 = j;
				n++;
			}
		}
	}
	return couple_at;

}
//fin des fonctions horton


int lecture_type(FILE *F)
{
	char c,c1,c3;
	c = fgetc(F); //lecture ";"
	c1 = fgetc(F);
	c = fgetc(F);
	c3 = fgetc(F);
	
	if( c1 == ' ')
		return 0;
	else
	{
		if( c3 ==' ')
			return 1;
		else 
			return 2;
		
	}
}


//procédure permattant de créer le fichier dot qui servira pour le dessin du graphe de cycles de la molécule
void creation_fichier_squelette(int mol_courante, int position, struct molecule *M)
{
	char src[64];
	sprintf(src,"%d_squelette.dot",mol_courante);
	//affichage_graphe_cycle_molecule(position,M);

	if ( M[position].g_def != 0) // le squelette a été construit
	{
		FILE *D = fopen(src, "w");

		if ( D == NULL)
		{
			fprintf(stderr, "Cannot create the file %s \n", src);
			exit(45);
		}

		fprintf(D, "graph G {\n");
		fprintf(D, "node [fixedsize = true,width = 0.2,height = 0.2,fontsize = 5, shape=circle];\n edge [fontsize = 8];\n");
		
		int i;
		for ( i = 0 ; i < M[position].g.nb_sommets ; i++)
		{

			fprintf(D, "%d [ label=\"%d\" ];\n",M[position].g.som[i].id,M[position].g.som[i].taille);
		}

		for ( i = 0 ; i < M[position].g.nb_arete ; i++)
		{
			fprintf(D,"%d -- %d",M[position].g.aretes[i].a1,M[position].g.aretes[i].a2);
			
			if ( M[position].g.aretes[i].type >= 0)
			{
				fprintf(D,"[color = blue,label = \"%d\"];\n",M[position].g.aretes[i].poids);
			}
			else
			{
				fprintf(D,"[color = red,label = \"%d\"];\n",M[position].g.aretes[i].poids);
			}
		}
			
		
		fprintf(D,"}");
		fclose(D);
	}
	else
	{
		fprintf(stderr, "Le graphe du squelette de la molecule %d n'a pas été calculé \n",mol_courante );
		exit(44);
	}
	
}

int calcul_position_distri( float j)
{
	int pos;
	if( j < 0.5)
	{
		if(j >=0 && j < 0.1)
			pos = 0;
		else
		{
			if(j >=0.1 && j < 0.2)
				pos =1;
				else
				{
					if(j >=0.2 && j < 0.3)
						pos =2;
						else
						{
							if(j >=0.3 && j < 0.4)
								pos =3;
								else

										pos =4;
				
						}
				}
		}
	}
	else 
	{
		if(j >=0.5 && j < 0.6)
			pos = 5;
		else
		{
			if(j >=0.6 && j < 0.7)
				pos =6;
				else
				{
					if(j >=0.7 && j < 0.8)
						pos =7;
						else
						{
							if(j >=0.8 && j < 0.9)
								pos =8;
								else

										pos =9;
				
						}
				}
		}		
		
	}
	
	return pos;
}

void construction_fichier_distribution(int pos)
{
	int i;
	float valeur;
	char src[64];
	sprintf(src,"distri_%d.data",pos);
	FILE *F = fopen(src,"r");
	if ( F == NULL)
	{
		fprintf(stdout,"Cannot open the file distri.data\n");
		exit(83);
	}
	int tableau[10];
	
	for( i = 0; i < 10; i++)
		tableau[i] = 0;
		

	int position;
	

	for( i = 0; i < NB_MOLECULES ; i++)
	{
		fscanf(F,"%f",&valeur);
		position = calcul_position_distri(valeur);
		tableau[position]++;
	}
	fclose(F);

	sprintf(src,"fichier_distri_%d.data",pos);
	F= fopen(src,"w");
	if ( F == NULL)
	{
		fprintf(stdout,"Cannot open the file distri.data\n");
		exit(82);
	}
	for( i = 0; i < 10 ; i++)
	{
		fprintf(F,"%.1f %d\n",0.1 * i, tableau[i]);
	}
	fclose(F);

}
