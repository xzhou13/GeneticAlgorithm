#include <string>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <fstream>
using namespace std;

class GA {
	public:
		// initializes the population P
		GA(int nb, int n, int g, double mp, double cp, int sd);
		// prints out the population
		void print(vector <string>& P);
		// prints out the new generation
		void newPrint();
		// calculates and prints the fitness, normalized fitness, and running total normalized fitness
		void fitness(vector <string>& P);
		// randomly selects two parents from population
		vector <int> selectParents();
		// mate 2 parents to produce 2 offsprings
		vector <string> mate(int p1, int p2, vector <string>& po);
		// perform mutation on an offspring
		vector <string> mutate(vector <string>& offspring);
		// generate another generation
		void secondGeneration(vector <string>& po);
		// calulate and returns average fitness value in one generation
		double avgFitness();
		// find the index of the best fitness in one generation
		int bestFitness();
		// returns the best fitness value in one generation
		double bestFitnessValue(int bf);
		// returns the number of ones in the most fit individual
		int numOnes(int bf, vector <string>& P);
		// output avg fitness, best fitness, and number of correct bits (standard output)
		void results(vector <string>& po);
		// run for number of generations and store results in vectors
		void run();
	protected:
		int numBits;				    // number of gene bits
		int N, G;						// N=population size; G=number of generations
		double pm, pc;					// pm=mutation probability; pc=crossover probability
		int seed;						// seed for random number generator
		vector <string> P;				// initial population
		vector <string> newP;			// new generation	
		vector <double> Pfit;			// fitness values for initial population
		vector <double> newPfit;		
		vector <double> nPfit;			// normalized fitness for initial population
		vector <double> ntotFit;		// running normalized total fitness
		vector <double> totalFitness;	// totalFitness values for all generations
		vector <double> af;				// avg fitness values for all generations
		vector <double> bf;				// best fitness for all generations
		vector <int> ones;				// number of ones in best fitness for all generations
};

// initializes the population P
GA::GA(int nb, int n, int g, double mp, double cp, int sd){
	numBits = nb;
	N = n;
	G = g;
	pm = mp;
	pc = cp;
	seed = sd;
	
	// initialize population P size
	P.resize(N);
	newP.resize(N);
	for (int i=0; i < N; i++){
		P[i].resize(numBits);
	}

	// initialize population P to '0' and '1' binary bits
	char bits[2] = {'0', '1'};
	srand(seed);
	for (int i=0; i < N; i++){
		for (int j=0; j < numBits; j++){
			P[i][j] = bits[rand() % 2]; 		
		}
	}

	// initialize sizes for vectors Pfit and nPfit
	Pfit.resize(N);
	nPfit.resize(N);
	ntotFit.resize(N);

	// initialize Pfit values to 0
	for (int i=0; i < N; i++){
		Pfit[i] = 0.00;
	}
	
}

// prints the population P
void GA::print(vector <string>& P){
	for (int i=0; i < N; i++){
		for (int j=0; j < numBits; j++){
			cout << P[i][j]; 		
		}
		cout << endl;
	}
	cout << endl;
}

// calculate fitness, normalized fitness, and running total normalized fitness
void GA::fitness(vector <string>& P){
	double tf = 0.00;
	double x;
	double ntf = 0.00;
	
	for (int i=0; i < N; i++){
		x = 0.00;
		for (int j=0; j < numBits; j++){
			if (P[i][j] == '1'){
				x += pow(2, (numBits - 1 - j));
			}
		}
		Pfit[i] = pow(x/pow(2, numBits), 10);
		
		// calculate running total fitness
		tf += Pfit[i];
	}

	cout << "Total fitness: " << tf << endl;

	for (int i=0; i < N; i++){
		// calculate normalized population fitness
		nPfit[i] = Pfit[i]/tf;
		ntf += nPfit[i];

		// calculate running normalized total fitnesses
		ntotFit[i] = ntf;
	}

	printf("%s  %s  %s  %s  \n", "Individual", "Fitness value", "Normalized fitness value", "Running total");
	for (int i=0; i < N; i++){
		printf("%*s%03d     %*s%f   %*s%f           %*s%f    \n", 4, "", i, 4, "", Pfit[i], 7, "", nPfit[i], 3, "", ntotFit[i]);
	}
	printf("\n");

}


// selects parents
vector <int> GA::selectParents(){
	double randNum1, randNum2;
	vector <int> parents;
	int parent1, parent2;
	
	parents.resize(2);

	randNum1 = (double) rand()/RAND_MAX;
	
	// pick 1st parent
	for (int i=0; i < N; i++){
		if (ntotFit[i] == randNum1){
			parent1 = i;
			parents[0] = parent1;
			break;
		}
		if (ntotFit[i] > randNum1){
			parent1 = i;
			parents[0] = parent1;
			break;
		}
	}
	
	// pick 2nd parent
	while(1){
		randNum2 = (double) rand()/RAND_MAX;
		//cout << "second random number is " << randNum2 << endl;
		for (int i=0; i < N; i++){
			if (ntotFit[i] == randNum2){
				parent2 = i;
			}
			if ((ntotFit[i] > randNum2) && (ntotFit[i-1] < randNum2)){
				parent2 = i;
			}
		}
		if (parent1 != parent2){
			parents[1] = parent2;
			break;
		}
	}
	
	return parents;
}

// mate 2 parents to produce 2 offsprings
vector <string> GA::mate(int p1, int p2, vector <string>& po){
	double cross;
	int crossSpot;
	vector <string> offsprings;
	string offspring1, offspring2;

	offsprings.resize(2);
	offspring1.resize(numBits);
	offspring2.resize(numBits);

	//determine if to cross or not
	cross = drand48();
	if ((cross <= pc) || (pc == 1.00)){
		// determine site of crossover and perform crossover
		crossSpot = rand() % (numBits-2)+1;
		offspring1 = po[p1].substr(0, crossSpot+1) + po[p2].substr(crossSpot+1, string::npos);
		offspring2 = po[p2].substr(0, crossSpot+1) + po[p1].substr(crossSpot+1, string::npos);
	}else{
		// no crossover
		offspring1 = po[p1];
		offspring2 = po[p2];
	}
	offsprings[0] = offspring1;
	offsprings[1] = offspring2;
	//cout << "two parents have mated to produce " << offspring1 << " and " << offspring2 <<  endl << endl;

	return offsprings;
}

// perform mutation on offsprings
vector <string> GA::mutate(vector <string>& offsprings){
	double mutate1, mutate2;
	for (int i=0; i < numBits; i++){
		mutate1 = drand48();
		mutate2 = drand48();
		if ((mutate1 <= pm) || (mutate1 == 1.00)){
			if (offsprings[0][i] == '0'){
				offsprings[0][i] = '1';
			}else{
				offsprings[0][i] = '0';
			}
		}
		if ((mutate2 <= pm) || (mutate2 == 1.00)){
			if (offsprings[1][i] == '0'){
				offsprings[1][i] = '1';
			}else{
				offsprings[1][i] = '0';
			}
		}
	}
	return offsprings;
}

// calculates another generation
void GA::secondGeneration(vector <string>& po){
	vector <int> parents;
	vector <string> ofs, ofs2;

	for (int i=0; i < N/2; i++){
		parents = selectParents();
		ofs = mate(parents[0], parents[1], po);
		ofs2 = mutate(ofs);
		newP[2*i] = ofs2[0];
		newP[2*i+1] = ofs2[1];
	}
}

// calulate and returns average fitness value in one generation
double GA::avgFitness(){
	double sum = 0.00;
	for (int i=0; i < N; i++){
		sum += Pfit[i];
	}
	return (double) sum/N;
}

// find the index of the best fitness in one generation
int GA::bestFitness(){
	double max = Pfit[0];
	int maxIndex = 0;
	for (int i=1; i < N; i++){
		if (Pfit[i] > max){
			max = Pfit[i];
			maxIndex = i;
		}
	}
	return maxIndex;
}

// returns the best fitness value in one generation
double GA::bestFitnessValue(int bf){
	return Pfit[bf];
}

// returns the number of ones in the most fit individual
int GA::numOnes(int bf, vector <string>& P){
	int num = 0;
	for (int j=0; j < P[bf].size(); j++){
		if (P[bf][j] == '1'){
			num++;
		}
	}
	return num;
}

// output avg fitness, best fitness, and number of 1s in best fit individual (standard output)
void GA::results(vector <string>& po){
	fitness(po);
	cout << "Avg fitness: " << avgFitness() << endl;
	af.push_back(avgFitness());
	int index = bestFitness();
	cout << "Best fitness: " << bestFitnessValue(index) << endl;
	bf.push_back(bestFitnessValue(index));
	cout << "Number of ones in the best fitness: " << numOnes(index, po) << endl << endl;
	ones.push_back(numOnes(index, po));
}

// run for number of generations and store results in vectors
void GA::run(){
	cout << "1st Generation: " << endl;
	print(P);
	
	FILE * fout;
	stringstream ss;
	ss << seed;
	string str = ss.str()+"results.csv";
	fout = fopen(str.c_str(), "w+");
	
	cout << "1st generation results: " << endl;
	results(P);
	fprintf(fout, "%s%s%s%s\n", "Generation,", "Average fitness,", "Best fitness value,", "# of ones,");
	fprintf(fout, "%d%s%f%s%f%s%d%s\n", 1, ",", af[0], ",", bf[0], ",", ones[0], ",");

	secondGeneration(P);
	cout << "2nd Generation: " << endl;
	print(newP);
	cout << "2nd generation results: " << endl;
	results(newP);
	fprintf(fout, "%d%s%f%s%f%s%d%s\n", 2, ",", af[1], ",", bf[1], ",", ones[1], ",");

	for (int i=0; i < G-2; i++){
		secondGeneration(newP);
		cout << i+3 << " generation: " << endl;
		print(newP);
		cout << i+3 << " generation results: " << endl;
		results(newP);
	}

	for (int i=0; i < G-2; i++){
		fprintf(fout, "%d%s%f%s%f%s%d%s\n", i+3, ",", af[i+2], ",", bf[i+2], ",", ones[i+2], ",");
	}

	fclose(fout);
}

int main(int argc, char** argv){
	int nb, n, g, sd;
	double mp, cp;
	istringstream ss;

	if (argc != 7){
		fprintf(stderr, "usage: number of genes in genetic string, population size, number of generations, mutation probability, crossover probability, seed for random number generator\n");
		return -1;
	}
	
	// read in and error check command line arguments entered
	ss.clear();
	ss.str(argv[1]);
	if ( !(ss >> nb) || (nb <=0) ){
		fprintf(stderr, "usage: number of genes in genetic string, must an integer > 0\n");
		return -2;
	}
	ss.clear();
	ss.str(argv[2]);
	if ( !(ss >> n) || (n <=0) ){
		fprintf(stderr, "usage: population size must an integer > 0\n");
		return -3;
	}
	ss.clear();
	ss.str(argv[3]);
	if ( !(ss >> g) || (g <=0) ){
		fprintf(stderr, "usage: number of generations must an integer > 0\n");
		return -4;
	}
	ss.clear();
	ss.str(argv[4]);
	if ( !(ss >> mp) || (mp < 0) || (mp > 1) ){
		fprintf(stderr, "usage: mutation probability must a decimal >= 0 and <= 1\n");
		return -5;
	}
	ss.clear();
	ss.str(argv[5]);
	if ( !(ss >> cp) || (cp < 0) || (cp > 1) ){
		fprintf(stderr, "usage: crossover probability must an decimal >= 0 and <=1\n");
		return -6;
	}
	ss.clear();
	ss.str(argv[6]);
	if ( !(ss >> sd) || (sd < 0) ){
		fprintf(stderr, "usage: seed for random number generator must an integer >= 0\n");
		return -7;
	}

	cout << "Number of genes(bits) in genetic string: " << nb << endl;
	cout << "Population size N: " << n << endl;
	cout << "Number of generations G: " << g << endl;
	cout << "Mutation probability pm: " << mp << endl;
	cout << "Crossover probability pc: " << cp << endl;
	cout << "Seed for random number generator: " << sd << endl << endl;
	
	// initialize genetic algorithm and run
	GA ga(nb, n, g, mp, cp, sd);
	ga.run();
	
	return 0;
}
