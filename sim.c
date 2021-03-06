#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/*******************************/
/*     CONSTANTES DEFINIDAS    */
/*******************************/

//Para os Neuronios
#define Rede 10
#define Excitadores 10
#define Inibidores (Rede - Excitadores)
#define Interacoes 1000
#define POP 30           //popn size
#define LEN 100           //genotypes
#define T1 50
#define T2 200
#define T3 300

float cont=0;
float gene[POP][LEN];

void imprime (float ms[Rede][Rede]){
  int i;
  int j;
  for (i=0;i<Rede;i++){
    for(j=0;j<Rede;j++){
      printf("%.2f ",ms[i][j]);
    }
    printf("\n");
  }
}

void imprimevetor (float m[]){
  int i;
  for (i=0;i<sizeof(m);i++){
      printf("%.2f ",m[i]);
  }
}

int main(int argc, char const *argv[]) {

  srand(time(NULL));



  // Arquivo saida
  FILE *saida = fopen ("out.txt", "w");
  FILE *saida2 = fopen ("ativacao.txt","w");
  //Declaracao da rede
  float a[Rede], b[Rede], c[Rede], d[Rede], S[Rede][Rede];
  float v[Rede], u[Rede], I[Rede];
  float re[Excitadores], ri[Inibidores];

  //Vetores de Disparos
  float disparado[Rede], disparando[Rede];

  //Auxiliares
  int i, j, t, corrente;


  /**********************/
  /*   INICIALIZACOES   */
  /**********************/




  //Peso da entrada

  a[0]=0.4;
  b[0]=-0.1;
  c[0]=-55;
  d[0]=4;

  //Peso dos Excitadores
  for (i = 1; i < Excitadores; i++){
    a[i] = 0.4;
    b[i] = -0.1;
    c[i] = -55;
    d[i] = 4;
  }

  //Peso dos Inibidores
  for (i = Excitadores; i < Rede; i++){
    a[i] = 0.4;
    b[i] = -0.1;
    c[i] = -55;
    d[i] = 4;
  }


  int m=0;
  int n=0;

  FILE * bestGenefile = fopen("logBestGene.txt","r");
	if(bestGenefile!=NULL){
    for(n=0;n<LEN;n++){
      fscanf(bestGenefile,"%f",&(gene[0][n])); //Store in position 0
    }
  }


  n=0;
  if(bestGenefile!=NULL){
    //inicializa a Rede
    for (i = 0; i < Rede; i++){
      //inicializa aqueles que dispararam
      disparando[i] = 0;
      //potencial inicial dos neuronios
      v[i] = -60;
      //taxa de restituicaoo da membrana
      u[i] = (v[i] * b[i]);
      //Peso das ligacoes dos Exicitadores
      for (j = 0; j < Excitadores; j++) {S[i][j] = gene[n][m]*20;m++;}
      //Peso das ligacoes dos Inibidores em S
      for (; j < Rede; j++) {S[i][j] = gene[n][m];m++;}
      //Como um neuronio nao liga com si mesmo,
      // o peso dessas ligacoes deve ser 0
      S[i][i] = 0;
    }
  }else{
    for (i = 0; i < Rede; i++){
      disparando[i] = 0;
      v[i] = -60;
      u[i] = (v[i] * b[i]);
      for (j = 0; j < Excitadores; j++) {S[i][j] = 20 * drand48();}
      for (; j < Rede; j++) {S[i][j] = -20*drand48();}
      S[i][i] = 0;
    }
  }


  //fclose(bestGenefile);

  float inputa;
  float inputb;
  printf("Inputa\n");
  scanf("%f",&inputa);
  printf("Inputb\n");
  scanf("%f",&inputb);

  /**********************/
  /*      SIMULACAO     */
  /**********************/

  //SIMULACAO DA COLUNA DA MATRIZ S

  for (t = 0; t < Interacoes; t++){
    //Resetar o vetor de disparandos e passar seu valor para disparados
    for (i = 0; i < Rede; i++){
      disparado[i] = disparando[i];
      if (disparado[i]){
        fprintf(saida, "%d\t%d\n", t, i);
      }
      disparando[i] = 0;
    }



    //Reinicia os valores de INPUT
    for (i = 0; i < Excitadores; i++) {I[i] = (0 * (rand())/RAND_MAX);}
    for ( ; i < Rede; i++) {I[i] = (0 * (rand())/RAND_MAX);}

    //inputs de entrada
    if(t>=T1 && t<=T1+100){
      I[0]=inputa;
    }
    if(t>=T2 && t<=T2+100){
      I[0]=inputb;
    }


    //printf("\n");
    //imprimevetor(I);
    //printf("\n");
    //imprimevetor(v);
    //printf("\n");

    //SIMULACAO DA LINHA DA MATRIZ S
    for (corrente = 0; corrente < Rede; corrente++){
      //ATUALIZA OS NOVOS VALORES DE ENTRADA
      for (i = 0; i < Rede; i++){
        if (disparado[i]) {I[corrente] += S[corrente][i];}
      }

      //ATUALIZA OS VALORES DE V e U



      v[corrente] += 0.25 * ((0.04 * (v[corrente] * v[corrente])) + 4.1 * v[corrente] + 108 - u[corrente] + I[corrente]);
      v[corrente] += 0.25 * ((0.04 * (v[corrente] * v[corrente])) + 4.1 * v[corrente] + 108 - u[corrente] + I[corrente]);
      u[corrente] += 0.5 * a[corrente] * (b[corrente] * v[corrente] - u[corrente]);

      //VERIFICA SE A POSICAO CORRENTE DISPAROU
      if (v[corrente] >= 30){

        if (t>=T1 && t<=T3+500){
          //printf("1 ");
          fprintf(saida2,"1 ");
        }
        if(corrente==5 && (t>=T3 && t<=T3+500)){
            cont++;
        }
        disparando[corrente] = 1;
        v[corrente]  = c[corrente];
        u[corrente] += d[corrente];
      }else{
        if (t>=T1 && t<=T3+500){
          //printf("0 ");
          fprintf(saida2,"0 ");
        }

      }
    }

    if (t>=T1 && t<=T3+500){
      //printf("%d \n",t);
      fprintf(saida2,"%d \n",t);
    }

  }

    imprime(S);
    printf("\n\n%f ",cont);
    printf("\n\nFrequencia: %f\n\n",cont/0.5);


  return 0;
}
