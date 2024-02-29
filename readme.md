# Genetic Algorithm on K-Means Clustering

This Project is mainly based on the [Genetic-Kmeans-Algorithm-GKA-](https://github.com/charlie60507/Genetic-Kmeans-Algorithm-GKA-)

## The approaches which I used
- Min-max normalization for standardization
- Davies–Bouldin index for evaluation of each cluster
- IN GENETIC :
  * Rank-based selection
  * One-point crossover

## Requirements
- Panda
- NumPy

## Getting Started
```
python __main__.py
```

## Input
- The data that I analyzed is from Iris
  * ```data/iris.csv``` have 3 column and ```data/iris2.csv``` have 4 column and ```data/isis_with_header.csv``` with header
- ```config.txt``` contain control parameters
  * kmax: maximum number of clusters
  * budget: budget of how many times run GA
  * numOInd: number of Individual
  * Ps: the probability of ranking Selection
  * Pc: the probability of crossover
  * Pm: the probability of mutation

## Output
- ```norm_data.csv``` is normalization data
- ```cluster_json``` is centroid of each cluster
- ```result.csv``` is data with labeled to each cluster

## Analysis
- the accuracy of GA on K-means: 88%
- the accuracy of k-means++: 83%

## Reference
- [Genetic-Kmeans-Algorithm-GKA-](https://github.com/charlie60507/Genetic-Kmeans-Algorithm-GKA-)
