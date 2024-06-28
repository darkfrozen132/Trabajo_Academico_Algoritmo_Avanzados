#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;


struct Individual {
    vector<int> chromosome;
    int fitness;

    Individual(vector<int> chromo) : chromosome(chromo), fitness(0) {}
};


double calculate_standard_deviation(vector<int> bin_usage, int max_value_per_bin) {
    double med = 0.0;
    int num_bins = bin_usage.size();
    vector<double> unused_space(num_bins);

    // Calcular el espacio no utilizado en cada mochila
    for (int i = 0; i < num_bins; ++i) {
        unused_space[i] = max_value_per_bin - bin_usage[i];
        med += unused_space[i];
    }
    med =med/ num_bins;

    // Calcular la suma de las diferencias al cuadrado respecto a la media
    double variance = 0.0;
    for (int i = 0; i < num_bins; ++i) {
        variance += (unused_space[i] - med) * (unused_space[i] - med);
    }
    variance =variance/ num_bins;

    // Retornar la desviación estándar
    return sqrt(variance);
}

void calculate_fitness(Individual& ind,  vector<int> itemPool,  vector<int>bin_capacities, int max_value_per_bin) {
    vector<int> bin_usage(bin_capacities.size(), 0);

    for (int i = 0; i < ind.chromosome.size(); ++i) {
        int bin_index = ind.chromosome[i];
        bin_usage[bin_index] += itemPool[i];
    }

    int total_overflow = 0;
    int full_bins = 0;

    for (int i = 0; i < bin_capacities.size(); ++i) {
        if (bin_usage[i] > max_value_per_bin) {
            ind.fitness = 0; // Si mochila sobrepasa el maximo valor, fitness automaticamente 0
            return;
        } else if (bin_usage[i] == max_value_per_bin) {
            full_bins += 1;//Calcula las mochilas completamente llenas
        }
    }

    // Calcular la deviacion estnadar para elegir las mochilas con espacio sobrantes mas lejos uno de otros
    double std_dev = calculate_standard_deviation(bin_usage, max_value_per_bin);
    ind.fitness = full_bins * 10000 + std_dev; //Asignar un gran valor a cada mochila llena
}

void evaluate_population(vector<Individual>& population,  vector<int>& itemPool,  vector<int>& bin_capacities, int max_value_per_bin) {
    for (Individual& ind : population) {
        calculate_fitness(ind, itemPool, bin_capacities, max_value_per_bin);
    }
}

pair<Individual, Individual> crossover_onepoint(Individual& parent1, Individual& parent2) {
    int size = parent1.chromosome.size();
    int crossover_point = rand() % size;
    vector<int> child1_chromo(size);
    vector<int> child2_chromo(size);

    for (int i = 0; i < size; ++i) {
        child1_chromo[i] = (i < crossover_point) ? parent1.chromosome[i] : parent2.chromosome[i];
        child2_chromo[i] = (i < crossover_point) ? parent2.chromosome[i] : parent1.chromosome[i];
    }

    return make_pair(Individual(child1_chromo), Individual(child2_chromo));
}

pair<Individual, Individual> crossover_uniform(Individual& parent1, Individual& parent2) {
    int size = parent1.chromosome.size();
    vector<int> child1_chromo(size);
    vector<int> child2_chromo(size);

    for (int i = 0; i < size; ++i) {
        if (rand() % 2) {
            child1_chromo[i] = parent1.chromosome[i];
            child2_chromo[i] = parent2.chromosome[i];
        } else {
            child1_chromo[i] = parent2.chromosome[i];
            child2_chromo[i] = parent1.chromosome[i];
        }
    }

    return make_pair(Individual(child1_chromo), Individual(child2_chromo));
}

void mutation_swap(Individual& ind, int num_bins) {
    int index1 = rand() % ind.chromosome.size();
    int index2 = rand() % ind.chromosome.size();
    ind.chromosome[index1] = rand() % num_bins;  // Asegurarse de que los índices sean válidos
    ind.chromosome[index2] = rand() % num_bins;  // Asegurarse de que los índices sean válidos
}

void mutation_insercion(Individual& ind){
    int index1 = rand() % ind.chromosome.size();
    int index2 = rand() % ind.chromosome.size();
    if(index1==index2) return;
    if(index2>index1){
        int a;
        a=index1;
        index1=index2;
        index2=a;
    }
    if(index2==index1+1) return;
    int aux;
    aux=ind.chromosome[index1+1];
    ind.chromosome[index1+1]=ind.chromosome[index2];
    for(int i=index1+2;i<=index2;i++){
        int aux2;
        aux2=ind.chromosome[i];
        ind.chromosome[i]=aux;
        aux=aux2;
    }
}

Individual tournament_selection(vector<Individual>& population, int tournament_size) {
    vector<Individual> tournament;
    for (int i = 0; i < tournament_size; ++i) {
        int random_index = rand() % population.size();
        tournament.push_back(population[random_index]);
    }

    return *max_element(tournament.begin(), tournament.end(), [](Individual& a, Individual& b) {
        return a.fitness < b.fitness;
    });
}



vector<Individual> select_survivors_ranking(vector<Individual>& population, vector<Individual>& offspring_population, int numsurvivors) {
    vector<Individual> next_population;
    population.insert(population.end(), offspring_population.begin(), offspring_population.end());
    sort(population.begin(), population.end(), [](Individual& a, Individual& b) {
        return a.fitness > b.fitness;
    });
    for (int i = 0; i < numsurvivors; ++i) {
        next_population.push_back(population[i]);
    }
    return next_population;
}

void genetic_algorithm(vector<Individual>& population, vector<int>& itemPool,vector<int>& bin_capacities, int max_value_per_bin, int generations, double mutation_rate, int tournament_size){
    int popsize = population.size();
    evaluate_population(population, itemPool, bin_capacities, max_value_per_bin); 
    vector<int> bestfitness;
    auto best_individual = max_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness < b.fitness;
    });
    bestfitness.push_back(best_individual->fitness);

    cout << "Poblacion inicial, best_fitness = " << best_individual->fitness << endl;
    cout << "Mejor individuo en la primera generacion: [";
    for (int gene : best_individual->chromosome) {
        cout << gene << " ";
    }
    cout << "] con fitness: " << best_individual->fitness << endl;
   

    // loop de generaciones
    for (int g = 0; g < generations; ++g) {
        // crea las parejas a cruzarse (mating pool)
        vector<pair<Individual, Individual>> mating_pool;
        for (int i = 0; i < popsize / 2; ++i) {
            mating_pool.push_back(make_pair(tournament_selection(population, tournament_size), tournament_selection(population, tournament_size)));
        }

        // cruza las parejas del mating pool. Cada cruzamiento genera 2 hijos
        vector<Individual> offspring_population;
        for (auto& parents : mating_pool) {  // por cada pareja del mating pool
            pair<Individual, Individual> children = crossover_uniform(parents.first, parents.second);  // cruzamiento uniforme
            //pair<Individual, Individual> children = crossover_onepoint(parents.first, parents.second);  // cruzamiento uniforme

            if ((double)rand() / RAND_MAX < mutation_rate) { // intenta mutar el hijo 1 de acuerdo a la tasa de mutacion
                mutation_swap(children.first, bin_capacities.size());
                //mutation_insercion(children.first);
                
            }
            if ((double)rand() / RAND_MAX < mutation_rate) { // intenta mutar el hijo 2 de acuerdo a la tasa de mutacion
                mutation_swap(children.second, bin_capacities.size());
               // mutation_insercion(children.second);
                
            }
            offspring_population.push_back(children.first);  // agrega el hijo 1 a la poblacion descendencia
            offspring_population.push_back(children.second); // agrega el hijo 2 a la poblacion descendencia
        }

        evaluate_population(offspring_population, itemPool, bin_capacities, max_value_per_bin); // evalua poblacion descendencia
        population = select_survivors_ranking(population, offspring_population, popsize); // selecciona sobrevivientes por ranking

        // obtiene el mejor individuo de la poblacion sobreviviente
        best_individual = max_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness < b.fitness;
        });
        bestfitness.push_back(best_individual->fitness);

//        if (g % 2 == 0) { // reporta cada 2 generaciones
//            cout << "Generation " << g << ", Best fitness: " << best_individual->fitness << endl;
//        }
    }
    
    cout << "Mejor individuo en la ultima generacion: [";
    for (int gene : best_individual->chromosome) {
        cout << gene << " ";
    }
    cout << "] con fitness: " << best_individual->fitness << endl;

}

vector<Individual> init_population(int population_size, int num_items, int num_bins) {
    vector<Individual> population;
    for (int i = 0; i < population_size; ++i) {
        vector<int> chromosome(num_items);
        for (int j = 0; j < num_items; ++j) {
            chromosome[j] = rand() % num_bins;
        }
        population.push_back(Individual(chromosome));
    }
    return population;
}

int main() {
   srand(time(0));

    const int NUM_ITEMS = 10;
    const int MAX_ITEM_WEIGHT = 150;
    const int NUM_BINS = 7;
    
    const int MAX_VALUE_PER_BIN = 150;


    // crea el pool de items
    vector<int> itemPool;
    for (int i = 0; i < NUM_ITEMS; ++i) {
        int weight = (rand() % (MAX_ITEM_WEIGHT / 10) + 1) * 10;
        itemPool.emplace_back(weight);
        cout << "Item " << i << ": Weight = " << itemPool[i] << endl; // Imprime ítem
    }

    vector<int> bin_capacities(NUM_BINS);

    const int POPSIZE = 100;
    const int GENERATIONS = 200;
    const double MUTATION_RATE = 0.8;
    const int TOURNAMENT_SIZE = 3;

    // inicializa poblacion
    

    for(int i=0;i<10;i++){
        vector<Individual> population = init_population(POPSIZE, NUM_ITEMS, NUM_BINS);
     genetic_algorithm(population, itemPool, bin_capacities, MAX_VALUE_PER_BIN, GENERATIONS, MUTATION_RATE, TOURNAMENT_SIZE);
    }
    return 0;
}