#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define twopi 6.283185307
#define pi 3.141592653
#define halfpi 1.570796325



//Para os Neuronios
#define Rede 10
#define Excitadores 10
#define Inibidores 0
#define Interacoes 1000

#define T1 50
#define T2 200
#define T3 300


#define POP 30           //popn size
#define LEN 100           //genotypes
#define MUT 0.05          // mutn prob at each locus
#define REC 0.6          //recomb prob at each locus 0.6
#define writeBestGene 1 //Save the best genotype of the simulation in a file
#define writePopulation 1 //Write the last population of a simulation in a file
#define StepSize 0.1
#define TotIterations 1000

bool loadPopulationFromFile; //Load all individuals genotype from file
bool logDynamic; //Save the dynamic of a individual (usuallly the best)
bool loadBestGeneFromFile; //Load only the best  from the last simulation
int  trials;        //Trials for each individual
int  end;
int  start;		//start value for logGA (column "iteration")
float fitness[POP];      //Fitness of each individual
float gene[POP][LEN];    //global array containing the binary genotypes of popn

float target1= 50;
float target2= 70;

//Protótipos
void init_pop();
void write_pop(int best);
float evaluate (int n );
void read_conf_file();

int main(){

	srand(time(NULL));
	int i,t,x,y,W,L,best=0;
	time_t t1,t2; (void) time(&t1);
  //Calcula o tempo de execu
	srand48(time(0)); srand (time(0));

  //inicia o fitness da população
	for (int k=0;k<POP;k++)
		fitness[k]=0;

	read_conf_file(); //Leitura do arquivo de configuração da simulação
	init_pop();//inicializa a população

	FILE * file;
	file = fopen("logGA.txt","w");//Armazena evolução do GA
	t=0;
	//Escolhe os dois primeiros robôs
	x=POP*drand48();
	y=POP*drand48();
	for (t=start;t<end+start;t++){
		//Avalia os 2 robs
		//fitness[x]=evaluate(x);
		//fitness[y]=evaluate(y);
		fitness[x]=0;
		fitness[y]=0;
		int t2;
		for (t2=0;t2<30;t2++){
			fitness[x]+=evaluate(x);
		}
		for (t2=0;t2<30;t2++){
			fitness[y]+=evaluate(y);
		}

		//Verifica qual foi o vencedor
		if (fitness[x] > fitness[y])  {W=x; L=y;}
		else                          {W=y; L=x;}
		//Altera os parâmetros (gene) do perdedor
		for ( i=0; i<LEN; i++){
			if (drand48()<REC)
			gene[L][i]=gene[W][i];
			if (drand48()<MUT){
				gene[L][i]+=(0.3*drand48())-0.15;
				if (gene[L][i]>1) gene[L][i]-=1;
				else if (gene[L][i]<0) gene[L][i]+=1;
			}
		}
		//Escolhe outros 2 indivíduos
		x=POP*drand48();
		do{
			y=POP*drand48();
		}while(y==x);

		if (t%100 ==0 || t+1==end+start){
			float fitpopulation;
			fitpopulation=fitness[0];
			best=0;
			//Atualiza o melhor indivíduo e o fitness da população.
			for (int k=1;k<POP;k++){
				if (fitness[best]<fitness[k])
					best=k;
				fitpopulation+=fitness[k];
			}
			fprintf(file,"%6d %2d %.2f %.2f\n",t,best,fitness[best],fitpopulation/POP);
			printf("%6d %2d %5.2f %5.2f %5.2f\n",t,best,fitness[best],fitpopulation/POP,fitness[W]);
			write_pop(best);
		}
	}
	logDynamic=true;
	printf("\nBest individual's (%d) internal dynamic saved - fitness: %9.5f",best, evaluate(best));
	write_pop(best);
	fclose(file);
	(void) time(&t2);
	printf("\nRunning time: %ld\n",t2-t1);
	return 0;

}

//Avalia o fitness de um determinado robô da população.
float evaluate (int n ){
	int m=0;
	float cont=0;
	float fitness =0.0; //Fitness do robô
	float worst = 1.0; //Fitness do pior robô
  float freq;
  //Declara��o da rede
  float a[Rede], b[Rede], c[Rede], d[Rede], S[Rede][Rede];
  float v[Rede], u[Rede], I[Rede];
  float re[Excitadores], ri[Inibidores];

  //Vetores de Disparos
  float disparado[Rede], disparando[Rede];

  //Auxiliares de manipula��o de la�os
  int i, j, t, corrente;


  /**********************/
  /*   INICIALIZA��ES   */
  /**********************/


  //for (i = 0; i < Excitadores; i++) {re[i] = (rand())/RAND_MAX;}
  //for (i = 0; i < Inibidores; i++)  {ri[i] = (rand())/RAND_MAX;}

  //Peso dos Exicitadores

  a[0]=0.4;
  b[0]=-0.1;
  c[0]=-55;
  d[0]=4;

  for (i = 1; i < Excitadores; i++){
    a[i] = 0.4;
    b[i] = -0.1;
    c[i] = -55;
    d[i] = 4;

    //a[i] = 0.02;
    //b[i] = 0.2;
    //c[i] = -65 + 15 * re[i];
    //d[i] = 8 - 6 * re[i];
  }

  //Peso dos Inibidores
  for (i = Excitadores; i < Rede; i++){
    a[i] = 0.02 + 0.08 * ri[(i - Excitadores - 1)];
    b[i] = 0.25 - 0.05 * ri[(i - Excitadores - 1)];
    c[i] = -65;
    d[i] = 2;
  }
//for media
  //inicializa a Rede
  for (i = 0; i < Rede; i++){
    //inicializa aqueles que dispararam
    disparando[i] = 0;
    //potencial inicial dos neuronios
    v[i] = -60;
    //taxa de restitui��o da membrana
    u[i] = (v[i] * b[i]);
    //Peso das liga��es dos Exicitadores
    for (j = 0; j < Excitadores; j++)// {S[i][j] = gene[n][m];m++;}
		{
				if(i!=j){
					S[i][j] = gene[n][m] * 20;
					m++;
				}
				else{
					S[i][j] = 0;
					m++;
				}

		}
    //Peso das liga��es dos Inibidores em S
    for (; j < Rede; j++) {S[i][j] = gene[n][m];m++;}
    //Como um neuronio nao liga com si mesmo,
    // o peso dessas liga��es ser� 0

  }


  /**********************/
  /*      SIMULA��O     */
  /**********************/

  //SIMULA��O DA COLUNA DA MATRIZ S

	float inputa=rand()%10+1;
	float inputb=rand()%10+1;

	while (inputa==inputb){
		inputb=rand()%10+1;
	}
	inputa=inputa*10;
	inputb=inputb*10;
	//printf("%f,%f",inputa,inputb);
  for (t = 0; t < Interacoes; t++){
    //Resetar o vetor de disparandos e passar seu valor para disparados
    for (i = 0; i < Rede; i++){
      disparado[i] = disparando[i];
      if (disparado[i]){
        //fprintf(saida, "%d\t%d\n", t, i);
      }
      disparando[i] = 0;
    }

    //CRIA UM INPUT ALEATORIO
    for (i = 0; i < Excitadores; i++) {I[i] = (0);}
		if(t>=T1 && t<=T1+100){
				I[0]=inputa;
		}
		if(t>=T2 && t<=T2+100){
			I[0]=inputb;
		}

    for ( ; i < Rede; i++) {I[i] = (0);}


    //SIMULA��O DA LINHA DA MATRIZ S
    for (corrente = 0; corrente < Rede; corrente++){
      //ATUALIZA OS NOVOS VALORES DE ENTRADA
      for (i = 0; i < Rede; i++){
        if (disparado[i]) {I[corrente] += S[corrente][i];}
      }

      //ATUALIZA OS VALORES DE V e U
      v[corrente] += 0.25 * ((0.04 * (v[corrente] * v[corrente])) + 4.1 * v[corrente] + 108 - u[corrente] + I[corrente]);
      v[corrente] += 0.25 * ((0.04 * (v[corrente] * v[corrente])) + 4.1 * v[corrente] + 108 - u[corrente] + I[corrente]);
      u[corrente] += 0.5 * a[corrente] * (b[corrente] * v[corrente] - u[corrente]);

      //VERIFICA SE A POSI��O CORRENTE DISPAROU

      if (v[corrente] >= 30){
        if(corrente==5 && (t>=T3 && t<=T3+500)){
            cont++;
        }
        disparando[corrente] = 1;
        v[corrente]  = c[corrente];
        u[corrente] += d[corrente];
      }else{
      }
    }
  }

  freq=cont/0.5;
	//printf("C:%f   F:%f ",cont, freq);
	if (freq==0){
		return 0;
	}

	if(inputa>inputb){

		return 100-fabs(target1-freq);
		/*if (freq!=target1){
			return 1-fabs(target1-freq);
			fitness -= 0.01/(fabs(target1-freq));
		}
		else{
				return 0.01;
				//fitness += 0.01;
		}
		if(logDynamic)
			printf("\nFitness do robo registrado no arquivo log.txt: %.2f \n",fitness);
		 //the worst trial will be returned
		if (worst > fitness)
			worst = fitness;*/
	}else{
		return 100-fabs(target2-freq);
		/*if (freq!=target2){
			return 100-fabs(target1-freq);
			fitness -= 0.01/(fabs(target2-freq));
		}
		else{
				return 0.01;
				fitness +=0.01;
		}
		if(logDynamic)
			printf("\nFitness do robo registrado no arquivo log.txt: %.2f \n",fitness);
		 //the worst trial will be returned
		if (worst > fitness)
			worst = fitness;*/
	}


		//if (worst==0) {break;} //do not complete the trials
	return (worst);
}


/*Initialise the population at random or loading from file.*/
void init_pop() {
	int i,j;
	if (loadPopulationFromFile) //Load the population from file
	{ FILE * popfile = fopen("logPopulation.txt","r");
		for (int i=0;i<POP;i++)
		{    for (int j=0;j<LEN;j++)
					 fscanf(popfile,"%f",&(gene[i][j]));
		}
		fclose(popfile);
	}
	else
	{
		for (i=0;i<POP;i++)
		{    for (j=0;j<LEN;j++)
				   gene[i][j]=drand48(); //Number between 0 and 1
				fitness[i]=0;

		}
		if (loadBestGeneFromFile) //Load only the best gene from another simulation
		{	 FILE * bestGenefile = fopen("logBestGene.txt","r");
			 for(int n=0;n<LEN;n++)
				 fscanf(bestGenefile,"%f",&(gene[0][n])); //Store in position 0
				//printf(" %7.5f ",gene[0][n]);
				fclose(bestGenefile);
		}
	}
}

void write_pop(int best){
	if (writeBestGene)
	{ FILE * bestGenefile = fopen("logBestGene.txt","w");
		for ( int i=0; i<LEN; i++)
			fprintf(bestGenefile,"%8.6f ",gene[best][i]);
		fclose(bestGenefile);
	}
	if (writePopulation)
	{ FILE * popfile = fopen("logPopulation.txt","w");
		for (int i=0;i<POP;i++){
		    for (int j=0;j<LEN;j++){
			     fprintf(popfile,"%8.6f ",gene[i][j]);
			  }
	     fprintf(popfile,"\n");
	  }
		fclose(popfile);
	}
}

void read_conf_file(){
  char *cmd;
  char buf[200];
  FILE *arq;
  arq = fopen("config.txt","r");
 if (arq == NULL){
    printf("\nErro ao abrir arquivo de configuracao\n");
    return;
  }
  fgets(buf,200,arq);
  while (!feof(arq))
  {
    cmd = strtok(buf,"=");
    int val = atoi(strtok(NULL,";"));
		if (strcmp(cmd,"evolution") == 0)
			if (val==0)
			{
				loadPopulationFromFile=false;
				 trials=1;
				 loadBestGeneFromFile=true;
				 logDynamic=true;
				 break;
			}
		if (strcmp(cmd,"loadPopulationFromFile") == 0)
			loadPopulationFromFile=val;
		if (strcmp(cmd,"trials") == 0)
			trials=val;
		if (strcmp(cmd,"logDynamic")==0)
			logDynamic=val;
		if (strcmp(cmd,"loadBestGeneFromFile")==0)
			loadBestGeneFromFile=val;
		if (strcmp(cmd,"end")==0)
			end=val;
		if (strcmp(cmd,"start")==0)
			start=val;
    fgets(buf,200,arq);

  }
  //printf("%d %d %d %d\n",loadPopulationFromFile,trials,logDynamic,loadBestGeneFromFile);
  fclose(arq);
}
