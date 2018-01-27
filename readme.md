# Genetic Algorithm on K-Means Clustering
## The approaches which I used
- Minmax normalization for standardization
- Davies–Bouldin index for evaluation of each cluster
- IN GENETIC :
  * Rank based selection
  * One point crossover

## Requirements
- panda
- numpy

## Getting Started
```
python __main__.py
```

## Input
- data which I analysis them is Iris
  * ```data/iris.csv``` have 3 column and ```data/iris2.csv``` have 4 column and ```data/isis_with_header.csv``` with header
- ```config.txt``` contain control parameters
  * kmax : maximum number of clusters
  * budget : budget of how many times run GA
  * numOInd : number of Individual
  * Ps : probability of ranking Selection
  * Pc : probability of crossover
  * Pm : probability of mutation

## Output
- ```norm_data.csv``` is normalization data
- ```cluster_json``` is centroid of each cluster
- ```result.csv``` is data with labeled to each cluster

## Analysis
- the accuracy of GA on K-means : 88%
- the accuracy of k-means++ : 83%
