import random
import numpy as np
from generation import Generation
from chromosome import Chromosome
from cluster import Clustering

random.seed(1)


class Genetic:
    def __init__(self, numberOfIndividual, Ps, Pm, Pc, budget, data, generationCount, kmax):
        self.numberOfIndividual = numberOfIndividual
        self.Ps = Ps
        self.Pm = Pm
        self.Pc = Pc
        self.budget = budget
        self.data = data
        self.generationCount = generationCount
        self.kmax = kmax

    def geneticProcess(self, generation):
        budget = self.budget
        Ps = self.Ps
        Pm = self.Pm
        Pc = self.Pc
        numOfInd = self.numberOfIndividual

        print("------------Generation:",
              self.generationCount, "-----------------")
        generation.sortChromosomes()

        # ------------------------simple ranking selection------------------------

        generation = self.selection(generation)

        #  ------------------------------Crossover---------------------------------

        generation = self.crossover(generation)

        #  ------------------------------Mutation---------------------------------

        generation = self.mutation(generation)

        self.generationCount += 1
        return generation, self.generationCount

    def selection(self, generation):
        numOfInd = self.numberOfIndividual
        Ps = self.Ps

        # replace the worst Ps*numOfInd individual with the best Ps*numOfInd individual
        for i in range(0, int(Ps * numOfInd)):
            generation.chromosomes[numOfInd -
                                   1 - i] = generation.chromosomes[i]

        # sort chromosomes after ranking selection
        generation.sortChromosomes()
        return generation

    def crossover(self, generation):
        numOfInd = self.numberOfIndividual
        Pc = self.Pc

        index = random.sample(
            range(0, numOfInd - 1), int(Pc * numOfInd))

        for i in range(int(len(index) / 2),+2):  # do how many time
            generation = self.doCrossover(
                generation, i, index)

        generation.sortChromosomes()

        return generation

    def doCrossover(self, generation, i, index):

        chromo = generation.chromosomes
        length = chromo[0].length
        cut = random.randint(1, length - 1)
        parent1 = chromo[index[i]]
        parent2 = chromo[index[i + 1]]
        genesChild1 = parent1.genes[0:cut] + parent2.genes[cut:length]
        genesChild2 = parent1.genes[cut:length] + parent2.genes[0:cut]
        child1 = Chromosome(genesChild1, len(genesChild1))
        child2 = Chromosome(genesChild2, len(genesChild2))

        # ----clustering----
        clustering = Clustering(generation, self.data, self.kmax)
        child1 = clustering.calcChildFit(child1)
        child2 = clustering.calcChildFit(child2)
        # -------------------

        listA = []
        listA.append(parent1)
        listA.append(parent2)
        listA.append(child1)
        listA.append(child2)
        # sort parent and child by fitness / dec
        listA = sorted(listA, reverse=True,
                       key=lambda elem: elem.fitness)

        generation.chromosomes[index[i]] = listA[0]
        generation.chromosomes[index[i + 1]] = listA[1]

        return generation

    def mutation(self, generation):
        numOfInd = self.numberOfIndividual
        fitnessList = []
        generationAfterM = Generation(numOfInd, generation.generationCount)
        flagMutation = (np.zeros(numOfInd)).tolist()

        for i in range(numOfInd):
            temp = generation.chromosomes[i]
            fitnessList.append(temp.fitness)

        for i in range(numOfInd):
            if i == 0:  # Ibest doesn't need mutation
                generationAfterM.chromosomes.append(generation.chromosomes[0])
                flagMutation[0] = 0
            else:
                generationAfterM = self.doMutation(
                    generation.chromosomes[i],	generationAfterM, flagMutation, fitnessList, i)

        generationAfterM.sortChromosomes()
        return generationAfterM

    def doMutation(self, chromosomeBeforeM, generationAfterM, flagMutation, fitnessList, i):
        Pm = self.Pm
        dice = []
        length = len(chromosomeBeforeM.genes)
        chromosome = Chromosome([], length)
        geneFlag = []

        for j in range(length):
            dice.append(float('%.2f' % random.uniform(0.0, 1.0)))
            if dice[j] > Pm:
                chromosome.genes.append(chromosomeBeforeM.genes[j])
                geneFlag.append(0)

            if dice[j] <= Pm:
                chromosome.genes.append(
                    float('%.2f' % random.uniform(0.0, 1.0)))
                geneFlag.append(1)

        check = sum(geneFlag)

        if check == 0:
            flagMutation[i] = 0
            chromosome.fitness = fitnessList[i]
        else:
            flagMutation[i] = 1

            #---clustering----
            clustering = Clustering(chromosomeBeforeM, self.data, self.kmax)
            chromosome = clustering.calcChildFit(
                chromosome)
            #------------------

        generationAfterM.chromosomes.append(chromosome)
        return generationAfterM
