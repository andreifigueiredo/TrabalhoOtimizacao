#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define Q 80				// Influencia na quantidade de feromonios depositados
#define RHO 0.95			// Influencia na taxa de evaporação do feronomio
#define ALPHA 1.0			    // Influencia do feromonio de um caminho
#define BETA 4.0			    // Influencia da distancia de um caminho
#define numCrit 3            //numero de criterios
#define MIN_T 1.0
#define MAX_T 6.0

#define N_ITENS 300		    // número de itens


typedef struct item {
	int weight;
	int profit[numCrit];
	double pheromone;
} ITEM;

typedef struct itemList{
	ITEM *item;
	struct itemList *next, *prev;
} ITEM_NODE;

typedef struct solution{
	ITEM **itens;
	int nItens;
	int profit;
} SOLUTION;


void addItem(ITEM *item, ITEM_NODE *head){
    ITEM_NODE *newNode = (ITEM_NODE*)malloc(sizeof(ITEM_NODE));

    if(newNode == NULL){
        fprintf(stderr, "Unable to allocate memory for new node\n");
        exit(-1);
    }

    newNode->item = item;
    newNode->next = NULL;
    newNode->prev = NULL;

    //printf(" - add test\n");
	if(head->next == NULL){
        head->next = newNode;
        //printf(" - first added\n");
	} else {
	    ITEM_NODE *currentNode = head;
	    while(1) {
            //printf(" - run while\n");
            if(currentNode->next == NULL){
                currentNode->next = newNode;
                currentNode->next->prev = currentNode;
                break;
            }
            currentNode = currentNode->next;
        }
        //printf(" - other added\n");
	}
}
void printItemList(ITEM_NODE *head){
    printf("\n ------------ ITENS NA LISTA: ------------\n");
    if(head->next != NULL){
	    ITEM_NODE *currentNode = head->next;
	    int countItens = 0;
	    while(1) {
            //printf(" - run while\n");
            countItens++;
            printf(" %d - p0: %d, p1: %d, p2: %d, w: %d, r: %f | cand: %d\n", countItens, currentNode->item->profit[0], currentNode->item->profit[1], currentNode->item->profit[2], currentNode->item->weight, currentNode->item->pheromone, currentNode);
            if(currentNode->next == NULL){
                break;
            }
            currentNode = currentNode->next;
        }
        printf(" ------------ Total de itens: %d ------------\n\n", countItens);
    }
}
void printItens(ITEM *itens, int nItens){

    int i;

    printf("\n ------------| TODOS OS ITENS: |------------\n");
    for(i = 0; i < nItens; i++){
        printf(" | %d - p0: %d, p1: %d, p2: %d,\tw: %d,\tr: %lf\t| cand: %d\n", i, itens[i].profit[0], itens[i].profit[1], itens[i].profit[2], itens[i].weight, itens[i].pheromone, &itens[i]);
    }
    printf(" ------------| Total de itens: %d |------------\n\n", i);
}
int sizeListItens(ITEM_NODE *head){
    if(head->next != NULL){
	    ITEM_NODE *currentNode = head->next;
	    int countItens = 0;
	    while(1) {
            //printf(" - run while\n");
            countItens++;
            if(currentNode->next == NULL){
                break;
            }
            currentNode = currentNode->next;
        }
        return countItens;
    } else{
        return 0;
    }
}

void clearItemList(ITEM_NODE **head) {
    ITEM_NODE *deleteMe;
    deleteMe = *head;
    while ((deleteMe = *head)) {
            printf(" -- teste 1 -- next: %d\n", deleteMe->next->item);
    //printf(" -- deleteme: %d", deleteMe->next);
    //printf(" -- deleteme: %s", deleteMe->next?"yes":"no");
        *head = deleteMe->next;
            printf(" -- teste\n");
        free(deleteMe);
            printf(" -- teste 2 -- next: %d\n", *head);
    }
}
/*
void addItemSolution(ITEM *item, SOLUTION *solution){
    solution->nItens = solution->nItens + 1;
    if(solution->itens == NULL)
        solution->itens = (ITEM**)malloc(sizeof(ITEM*));

    if(solution->itens == NULL){
        printf("\n erro de alocacao\n");
        free(solution->itens);
        exit(-1);
    }
    else{
        solution->itens = (ITEM**)realloc(solution->itens, solution->nItens*sizeof(ITEM*));

        if(solution->itens == NULL){
            printf("\n erro de realocacao: %d\n",solution->itens);
            free(solution->itens);
            exit(1);
        }
    }
    solution->itens[solution->nItens-1] = item;
}
*/
void addItemSolution(ITEM *item, SOLUTION *solution){
    solution->itens[solution->nItens] = item;
    solution->nItens = solution->nItens + 1;
}

void printSolution(SOLUTION *solution){
    int i;
    printf(" ------------------ Itens na solucao ---------------------\n");
    for(i = 0; i < solution->nItens; i++){
        printf("-: - p: %d,\tw: %d,\tr: %lf\t| cand: %d\n", solution->itens[i]->profit, solution->itens[i]->weight, solution->itens[i]->pheromone, solution->itens[i]);
    }
    printf(" -------------- n de itens na solucoes: %d ---------------\n", i);
}


int createNeighbourhood(ITEM_NODE *neighbourhood, ITEM *itens, int nItens, int currentC){
	int i, sizeNeighbourhood = 0;
	ITEM_NODE *currentNode;
	currentNode = neighbourhood;
	for(i = 0; i < nItens; i++){
        if(itens[i].weight <= currentC){
            addItem(&itens[i], neighbourhood);
            sizeNeighbourhood++;
        } else{
            nItens++;
        }
	}
	return sizeNeighbourhood;
}

int filterNeighbourhood(ITEM_NODE *neighbourhood, ITEM *candidate, int currentC){
    int countRemoved[3] = {0,0,0};
    int countItens = 0;
	if(currentC >= 0){

		ITEM_NODE *currentItem;
		struct itemList *auxCurrentItem = NULL; //auxCurrentItem tem que ser ponteiro?
		currentItem = neighbourhood->next;

		while(currentItem != NULL){
            //printf("entrou no while");
            //printf("p: %d, w: %d, r: %f\n", currentItem->item.profit, currentItem->item.weight, currentItem->item.pheromone);

			//if(&currentItem->item == candidate)
                //printf(" current Item: %d = candidate: %d\n", &currentItem->item, candidate);
            countItens++;
            auxCurrentItem = currentItem->next;
			if(currentItem->item->weight > currentC || currentItem->item == candidate){
                //printf(" e no if");
                //countItens--;
				if(currentItem->next != NULL && currentItem->prev != NULL){
                    countRemoved[1]++;
					currentItem->next->prev = currentItem->prev;
					currentItem->prev->next = currentItem->next;
                    //printf(" e removeu do meio");

				} else if(currentItem->prev != NULL){
                    countRemoved[2]++;
					currentItem->prev->next = NULL;
                    //printf(" e removeu do fim");

				} else if(currentItem->next != NULL){
                    countRemoved[0]++;
					currentItem->next->prev = NULL;
					neighbourhood->next = currentItem->next;
                    //printf(" e removeu do inicio");
				} else{
                    countRemoved[0]++;
					neighbourhood->next = NULL;
                    //printf(" e removeu o ultimo");
				}
	        	free(currentItem);
			}
            //printf("\n");
			currentItem = auxCurrentItem;
		}
                                    //printf(" - REMOVIDOS DA LISTA - Inicio: %d , Meio: %d , Final: %d | Total: %d de %d (sobra %d)\n", countRemoved[0], countRemoved[1], countRemoved[2], (countRemoved[0]+countRemoved[1]+countRemoved[2]), countItens, (countItens - countRemoved[0] - countRemoved[1] - countRemoved[2]));
	}
	return (countItens - countRemoved[0] - countRemoved[1] - countRemoved[2]); //retorna numero de itens na vizinhança
}

double calculateDelta(int bestProfit, int globalProfit){
    double quotient = (double) (globalProfit - bestProfit) / globalProfit;
    //double quotient = (double) (globalProfit - bestProfit);
	return 1 / (1 + quotient);
}


int updatePheromones(SOLUTION *cicleSolution, SOLUTION *BEST, ITEM *itens, int nItens){
    int i;
    double delta;
    //printf(" - cicleSolution nItens: %d\n", cicleSolution->nItens);
    /**/
    for(i = 0; i < nItens; i++){
        itens[i].pheromone = itens[i].pheromone * RHO;
        /*if(i == 11){*/
            //printf(" - Item %d -> p: %d | w: %d: | r: %f      -     ends.: %d - %d\n", i, cicleSolution->itens[i]->profit, cicleSolution->itens[i]->weight,itens[i].pheromone, cicleSolution->itens[i], &itens[i]);
        /*}*/
        if(itens[i].pheromone < MIN_T)
            itens[i].pheromone = MIN_T;
    }


    delta = calculateDelta(cicleSolution->profit, BEST->profit);
    //printf("    - Delta: %lf | best: %d | global: %d \n", delta, cicleSolution->profit, BEST->profit);

    for (i = 0; i < cicleSolution->nItens; i++) {
        cicleSolution->itens[i]->pheromone = cicleSolution->itens[i]->pheromone + (50  * delta);
        //printf("%d - %lf \t", &cicleSolution->itens[i]->pheromone, cicleSolution->itens[i]->pheromone);
        if(cicleSolution->itens[i]->pheromone >= MAX_T)
            cicleSolution->itens[i]->pheromone = MAX_T;
    }
    //printf("\n");
}
/*
int updateTrails(ITEM* itens, int nItens, SOLUTION cicleSolution){
	int i;
	double auxSum;
	for (i = 0; i < nItens; ++i) {
        auxSum = itens[i].pheromone + calculateDelta(itens[i].profit, cicleSolution.profit);
        if(auxSum <= MAX_T)
            itens[i].pheromone = auxSum;
	}
}

int evaporate(ITEM *itens, int nItens){
	if(itens != NULL && nItens > 0){
		int i;
        double auxProduct;
		for(i = 0; i < nItens; i++){
            auxProduct = itens[i].pheromone * RHO;
            if(auxProduct >= MIN_T)
                itens[i].pheromone = auxProduct;
			if(i == 4){
                printf(" - pheromore %d -> p: %d | w: %d: | r: %f\n", i, itens[i].profit, itens[i].weight, itens[i].pheromone);
			}
		}
		return 1;
	}
	return 0;
}
*/

double probProduct(ITEM item){
    double product;
	double quotient = ((double)((item.profit[0] + item.profit[1] + item.profit[2]) / numCrit) / (double)item.weight);
	//printf(" - p: %f, w: %f, q: %f\n", (double)item.profit, (double)item.weight, quotient);
	//printf(" - quotient: %f\n", (pow(item.pheromone, ALPHA) * pow(quotient, BETA)));
	product = (pow(item.pheromone, ALPHA) * pow(quotient, BETA));
	return product;
}

ITEM* selectCandidate(ITEM_NODE *neighbourhood) {
	double denominator = 0.0, c = 0.0, r;
	int i;
	//time_t t;
	ITEM_NODE *currentNode;
    currentNode = neighbourhood->next;
    //printf(" - p: %d , w: %d , a: %d\n", (*neighbourhood).item.profit, (*neighbourhood).item.weight, &(*neighbourhood).item.weight);
    //printf(" - p: %d , w: %d , a: %d\n", (*currentNode).item.profit, (*currentNode).item.weight, &(*currentNode).item.weight);


	r = (double)rand()/(double)RAND_MAX;
	//printf(" - random: %f\n", r);
    //printf("dfdfsdfsdf\n\n");
	while(currentNode != NULL) {
        //printf(" - p: %f, w: %f\n", (double)currentNode->item.profit, (double)currentNode->item.weight);
		denominator += probProduct(*currentNode->item);
        //printf(" - probproduct: %d | %d | %f\n", currentNode->item.profit, currentNode->item.weight, probProduct(currentNode->item));
        //printf(" - denominator: %f\n", denominator);
		currentNode = currentNode->next;
	}
    currentNode = neighbourhood->next;

	if(denominator > 0.0) {
		while(currentNode != NULL) {
			c += probProduct(*currentNode->item)/denominator;
			//printf("confere select: r: %lf <= c: %lf    oq e: %s     | currentNode: %d , next: %d\n", r, c, (r==c)?"igual":"diferente", currentNode, currentNode->next);
			if((r <= c) || (currentNode->next == NULL)){
                //printf(" - entrou if\n");
                //printf(" - selected adress: %d\n", currentNode->item);
                //printf(" - selected: p: %d, w: %d\n", currentNode->item.profit, currentNode->item.weight);
                break;
			}

            currentNode = currentNode->next;
		}
		return currentNode->item;
	} else {
		return NULL;
	}
}

void generateItens(ITEM * itens, int nItens, int minP, int maxP, int minW, int maxW){
    int i, j;
    int arrP[nItens], arrW[nItens];
    ITEM newItem;
    for(i = 0; i < nItens; i++){
        arrP[i] = (rand() % (maxP - minP)) + minP;
        arrW[i] = (rand() % (maxW - minW)) + minW;
        j = 0;
        while(j < i){
            //printf("%d - entrou %d\t|  %d - %d\t|  %d - %d\n", j, i, arrP[j], arrP[i], arrW[j], arrW[i]);
            if((arrP[i] == arrP[j]) && (arrW[i] == arrW[j])){
                arrP[i] = (rand() % (maxP - minP)) + minP;
                arrW[i] = (rand() % (maxW - minW)) + minW;
                j = 0;
            } else{
                j++;
            }
        }
        newItem.profit[0] = arrP[i];
        newItem.profit[1] = arrP[i];
        newItem.profit[2] = arrP[i];
        newItem.weight = arrW[i];
        newItem.pheromone = MAX_T;
        //printf(" - max %d | %d\n", maxW, minW);
        //printf(" - inseriu um\n");
        //printf(" - p: %d, w: %d\n", newItem.profit, newItem.weight);
        itens[i] = newItem;
    }
}


int main(void){
	int nCycles = 50;						//numero de ciclos
	int nAnts = 10;							//numero de formigas
	int capacity = 3000;					//capacidade total da mochila

	int currentC = capacity;				//capacidade total da mochila

	int countCycles = 0, countAnts = 0, countItens = 0, nNeighbours, auxCount;
	ITEM *candidate;							//item selecionado na iteração atual
	SOLUTION solution;		 				//lista da solução parcial
	SOLUTION cicleSolution; 					//lista da melhor solução (por formiga)
	SOLUTION BEST; 				//lista da melhor solução global (por ciclo)

    solution.itens = (ITEM**)malloc(N_ITENS*sizeof(ITEM*));
    cicleSolution.itens = (ITEM**)malloc(N_ITENS*sizeof(ITEM*));
    BEST.itens = (ITEM**)malloc(N_ITENS*sizeof(ITEM*));

    solution.nItens = 0;
    solution.profit = 0;

    cicleSolution.profit = 0;
    cicleSolution.nItens = 0;

    BEST.profit = 0;
    BEST.nItens = 0;


	ITEM itens[N_ITENS];					//objetos
	generateItens(itens, N_ITENS, 1, 10, 1, 100);
	printf(" - Itens generated\n");


	//srand(time(NULL));
	srand(15);

    int test = 0;
	//printf("p - %d | w - %d | r - %d\n", itens[test].profit, itens[test].weight, itens[test].pheromone);

	ITEM_NODE neighbourhood;
	neighbourhood.next = NULL;
	neighbourhood.prev = NULL;
	nNeighbours = createNeighbourhood(&neighbourhood, itens, N_ITENS, currentC);
	printf(" - neighbourhood created | size: %d -\n", sizeListItens(&neighbourhood));


	while(countCycles < nCycles){
		while(countAnts < nAnts){
			while(currentC > 0 && nNeighbours != 0){

                //printf(" - p: %d , w: %d\n", neighbourhood.item.profit, neighbourhood.item.weight);

                candidate = selectCandidate(&neighbourhood); //seleciona candidato e o remove da vizinhança

                //printf(" -- n teste: %d --\n", test);
                //test++;

				addItemSolution(candidate, &solution);
                //printf(" --- solution - w: %d\n", solution.head.next->item.weight);
                                //printf("\n - Capaciade antes: %d |", currentC);
				currentC -= candidate->weight;
                                //printf(" Capacidade depois : %d\n\n", currentC);
				solution.profit += (candidate->profit[0] + candidate->profit[1] + candidate->profit[2]) / numCrit;
                                //printf("        - Melhor solucao global ate agora: %d -\n\n", BEST.profit);
                                //printf("        - CICLO %d | FORMIGA %d | ITEM %d -\n", countCycles, countAnts, countItens);
                                //printf("       - Candidato selecionado: %d -\n", candidate);
                                //printf(" - profit: %d , weight: %d , pheromone: %f\n", candidate->profit, candidate->weight, candidate->pheromone);
                /*
				printf(" - itens adicionados: %d | itens na solucao: %d\n", countItens+1, solution.nItens);
				printf(" - itens adicionados: %d | itens na Best solucao: %d\n", countItens+1, cicleSolution.nItens);
				printf(" - itens adicionados: %d | itens na Global solucao: %d\n", countItens+1, BEST.nItens);
                */

                //printItemList(&neighbourhood); // IMPRIME OS ITENS QUE PODEM SER SELECIONADOS
                				//printf("            -- capacity: %d --\n", currentC);
                //printItemList(&neighbourhood);
                //printf(" - neighbourhood size: %d\n", sizeListItens(&neighbourhood));
                nNeighbours = filterNeighbourhood(&neighbourhood, candidate, currentC); //retorna numero de vizinhos (= itens que podem ser selecionados)
                                //printf("\n - numero de itens na vizinhanca: %d | capaciade posterior: %d\n", nNeighbours, currentC);

                countItens++;
                                //printf("_________________________________________________________________\n\n");

			}
			countItens = 0;
            currentC = capacity;
            nNeighbours = createNeighbourhood(&neighbourhood, itens, N_ITENS, currentC);

			if(solution.profit > cicleSolution.profit){
				//printf(" - teste do segmentation %d > %d\n",solution.profit, cicleSolution.profit );
                for(auxCount = 0; auxCount < solution.nItens; auxCount++){
            		//printf(" -- teste %d -- \n", auxCount);
                    cicleSolution.itens[auxCount] = solution.itens[auxCount];
                }
				cicleSolution.nItens = solution.nItens;
				cicleSolution.profit = solution.profit;
				//printf("     end. solution itens: %d - %d \n end. cicleSolution itens: %d - %d\n", solution.itens[0], solution.itens[0]->profit, cicleSolution.itens[0], cicleSolution.itens[0]->profit);
			}
			//printSolution(&solution);

            //printf("teste %d\n", *(solution.head.next));
            //printf(" - solution item profit: %d\n", solution.profit);
            //free(solution.itens);
            //solution.itens = NULL;
			solution.nItens = 0;
			solution.profit = 0;
			countAnts++;
		}
		countAnts = 0;

        if(cicleSolution.profit > BEST.profit){
            for(auxCount = 0; auxCount < cicleSolution.nItens; auxCount++){
                BEST.itens[auxCount] = cicleSolution.itens[auxCount];
            }
            BEST.nItens = cicleSolution.nItens;
            BEST.profit = cicleSolution.profit;
        }
		//evaporate(itens, N_ITENS);
		//updateTrails(itens, N_ITENS, cicleSolution);
		updatePheromones( &cicleSolution, &BEST, itens, N_ITENS );

        printItens(itens, N_ITENS);
        //printSolution(&cicleSolution);
        //printSolution(&BEST);

		cicleSolution.profit = 0;
		cicleSolution.nItens = 0;
		//cicleSolution.itens = NULL;

		countCycles++;
		printf(" Best solution till cycle %d:\t PROFIT:%d\t Quantidade de Criterios:%d\n",countCycles, BEST.profit, numCrit);
	}
	free(solution.itens);
}
/*
C – is the total knapsack load capacity,
Vc – is the current knapsack load capacity, V C C g w g Si = −∑ ∈ ( ) ,
Si – is a partial solution,
g S ( ) wg ∑ ∈ i – is the weight of all objects which were included in the partial
solution Si,
wj – is the weight of selected object j,
zj – is the profit of selected object j,
µj – is the attractiveness of selecting an object j.

begin
	while (a cycle exists) do
		while (an ant k, which has not yet worked, exists) do
			while (Vc ≥ 0) do
				select a next object oj from neighbourhood with probability pj
				add a selected object to a partial solution S = S + {oj}
				update the current knapsack load capacity VC = VC – wj
				update the profit Z = Z + zj
				update the neighbourhood of the current state Ni = {oi : wi ≤ VC}
			end
			remember the best solution if a better solution has been found
		end
		remember a global best solution if a better solution has been found
		use an evaporation mechanism τ = ρτ
		update pheromone trails τ = τ + ∆τ
	end
end.

*/