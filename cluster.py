import pandas as pd
import random
import math
import json
import numpy as np
#from generation import Generation
#from chromosome import Chromosome
pd.options.mode.chained_assignment = None
random.seed(1)


class Cluster:
    def __init__(self, generation, data, kmax):
        self.generation = generation
        self.data = data
        self.dim = data.shape[1]
        self.penalty = 1000000
        self.kmax = kmax

    def calcDistance(self, chromosome):
        kmax = self.kmax
        dim = self.dim
        data = self.data
        dis = 0
        disSet = []
        tmp = []
        allIndex = []
        allDis = []

        for z in range(0, data.shape[0]):
            for i in range(0, kmax):
                if chromosome.genes[i] >= 0.5:		# tag=1 -> used
                    # euclidian distance
                    for j in range(0, dim):  # j is # of dim
                        # print("i",i,"j",j)
                        # print("gene",kmax + dim * i + j)
                        # print("size of chromosome",chromosome.length)
                        # print("chromosome:",chromosome.genes[kmax + dim * i + j])
                        square = pow(
                            chromosome.genes[kmax + dim * i + j] - data.loc[z][j], 2)
                        tmp.append(square)		# save (center-point)^2 in tmp

                    for t in range(0, dim):
                        dis += tmp[t]

                    dis = math.sqrt(dis)
                    disSet.append(dis)
                    dis = 0
                    tmp = []

                elif chromosome.genes[i] < 0.5:  # tag=0 -> not used
                    disSet.append(self.penalty)

            allDis, allIndex = self.findMin(
                disSet, chromosome, allIndex, allDis)
            disSet = []  # clear disSet	# calculate distance

        return allDis, allIndex

    # --SHOULD_BE_UNDERSTAND-- find minimum
    def findMin(self, disSet, chromosome, allIndex, allDis):
        tmp = []
        dis = 0

        if not disSet == []:
            n = disSet.index(min(disSet))  # n is index
            minDis = disSet[n]
            allIndex.append(n + 1)

        allDis.append(minDis)

        # handle no center is used -> assign to center1
        for k, item in enumerate(allDis):
            if item == 1000000:
                for z in range(0, self.data.shape[0]):
                    for j in range(0, self.dim):  # j is # of dim
                        i = 0
                        square = pow(
                            chromosome.genes[self.kmax + self.dim * i + j] - self.data.loc[z][j], 2)
                        tmp.append(square)
                    for t in range(self.dim):
                        dis += tmp[t]
                    dis = math.sqrt(dis)
                allDis[k] = dis
                tmp = []

        return allDis, allIndex

    def calcInter(self, allDis, allIndex, countInter, interSum):
        kmax = self.kmax
        data = self.data
        allSum = 0

        for i in range(0, data.shape[0]):  # norm_data.shape[0]
            for k in range(0, kmax):		# k: # of cluster(Kmax)
                centroid = k + 1

                if allIndex[i] == centroid:
                    countInter[centroid] += 1
                    interSum[centroid] += allDis[i]

        for i in range(0, kmax):
            if countInter[i + 1] == 0:  # 0 can't be Denominator
                countInter[i + 1] = 1  # set Denominator = 1

        for i in range(0, kmax):  # average distance / number for each cluster
            allSum += interSum[i + 1] / countInter[i + 1]

        # average
        interScore = allSum / (kmax)

        return interScore

    def calcIntra(self, chromosome):
        kmax = self.kmax
        dim = self.dim
        centerTemp = []
        centerList = []
        count = 0

        for i in range(0, kmax):  # catch center point in chromo to centerList
            if chromosome.genes[i] > 0.5:
                count += 1
                for j in range(0, dim):  # j is # of dim
                    centerTemp.append(chromosome.genes[kmax + dim * i + j])
                centerList.append(centerTemp)		#
                centerTemp = []
        i = 0

        if count == 0:
            intraScore = 0
        elif count == 1:
            intraScore = 0
        elif count >= 2:
            for j in range(i + 1, count):
                for i in range(0, count - 1):
                    v1 = centerList[i]
                    v2 = centerList[j]
                    # Combination : count get 2
                    vec = list(map(lambda x: x[0] - x[1], zip(v1, v2)))
                    # square of vec
                    vecSquare = list(map(lambda x: pow(x, 2), vec))
                    sumIntra = sum(vecSquare)
            intraScore = math.sqrt(sumIntra) / count

        return intraScore

    def calcFitness(self, intraScore, interScore, chromosome, countFitTime):
        if interScore != 0:
            fitness = float(intraScore) / float(interScore)
        else:
            fitness = 0

        chromosome.fitness = fitness
        countFitTime += 1
        return fitness, chromosome, countFitTime

    # childChromosome, kmax, countFitTime
    def calcChildFit(self, childChromosome, countFitTime):
        kmax = self.kmax

        countInter = (np.zeros(kmax + 1)).tolist()
        interSum = (np.zeros(kmax + 1)).tolist()
        allDis, allIndex = self.calcDistance(childChromosome)
        interScore = self.calcInter(
            allDis, allIndex, countInter, interSum)
        intraScore = self.calcIntra(childChromosome)

        # cal fitness
        if interScore != 0:
            fitness = float(intraScore) / float(interScore)
        else:
            fitness = 0

        childChromosome.fitness = fitness
        countFitTime += 1

        return childChromosome, countFitTime

    def calcChromosomesFit(self):
        kmax = self.kmax
        generation = self.generation
        numOfInd = generation.numberOfIndividual
        data = self.data
        chromo = generation.chromosomes
        countFitTime = 0

        for i in range(0, numOfInd):
            countInter = (np.zeros(kmax + 1)).tolist()
            interSum = (np.zeros(kmax + 1)).tolist()

            allDis, allIndex = self.calcDistance(chromo[i])

            interScore = self.calcInter(
                allDis, allIndex, countInter, interSum)
            # exit(0)
            intraScore = self.calcIntra(chromo[i])
            fitness, chromosome, countFitTime = self.calcFitness(
                intraScore, interScore, chromo[i], countFitTime)
            generation.chromosomes[i] = chromosome

        return generation, countFitTime

    def printIBest(self, iBest, countFitTime):
        kmax = self.kmax

        countInter = (np.zeros(kmax + 1)).tolist()
        interSum = (np.zeros(kmax + 1)).tolist()
        allDis, allIndex = self.calcDistance(iBest)
        interScore = self.calcInter(
            allDis, allIndex, countInter, interSum)
        intraScore = self.calcIntra(iBest)

        if interScore != 0:
            fitness = float(intraScore) / float(interScore)
        else:
            fitness = 0

        print("Cal fit time:", countFitTime)
        print("iBest Fitness:", fitness)
        print("Cluster Index:", allIndex)
        print()
        print("")

    def output_result(self, iBest, _data):
        print("Saving the result...")
        kmax = self.kmax
        countInter = (np.zeros(kmax + 1)).tolist()
        interSum = (np.zeros(kmax + 1)).tolist()
        allDis, allIndex = self.calcDistance(iBest)

        # decode cluster center
        rule_index = []

        for i in range(self.kmax):
            if iBest.genes[i] >= 0.5:
                rule_index.append(1)
            else:
                rule_index.append(0)

        rDict = dict()
        tmp = iBest.genes[self.kmax:]
        for i in range(self.kmax):
            if rule_index[i] == 1:
                begin_index = i * self.dim
                rDict[i + 1] = tmp[begin_index: begin_index + self.dim]

        with open('result/cluster_center.json', 'w') as outfile:
            json.dump(rDict, outfile, sort_keys=True,
                      indent=4, separators=(',', ': '))

        # rename df header
        col_name = list()
        for i in range(_data.shape[1]):
            col_name.append("f{0}".format(i))
        _data.columns = col_name

        # insert cluster result
        _data['Cluster Index'] = pd.Series(allIndex, index=_data.index)
        _data.to_csv('result/result.csv', index=None)
        print("Done.")