from chromosome import Chromosome


class Generation:
    def __init__(self, numberOfIndividual, generationCount):
        self.numberOfIndividual = numberOfIndividual
        self.chromosomes = []
        self.generationCount = generationCount

    def sortChromosomes(self):
        self.chromosomes = sorted(
            self.chromosomes, reverse=True, key=lambda elem: elem.fitness)
        return self.chromosomes

    def randomGenerateChromosomes(self, lengthOfChromosome):
        for i in range(0, self.numberOfIndividual):
            chromosome = Chromosome([], lengthOfChromosome)
            chromosome.randomGenerateChromosome()
            self.chromosomes.append(chromosome)
