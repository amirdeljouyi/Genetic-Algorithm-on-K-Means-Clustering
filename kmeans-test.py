import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
import numpy as np
from sklearn.cluster import KMeans
from sklearn import preprocessing
from scipy.spatial import distance
import pandas as pd

# correct = 0
# for i in range(len(x)):
#     predict_me = np.array(x[i].astype(float))
#     predict_me = predict_me.reshape(-1, len(predict_me))
#     prediction = clf.predict(predict_me)
#     if prediction[0] == y[i]:
#         correct+=1


def main():
    data = pd.read_csv('result/norm_data.csv', header=None)
    clf = KMeans(n_clusters=3)
    clf.fit_predict(data)
    print(clf.cluster_centers_)
    print(clf.labels_)

    centroids = clf.cluster_centers_
    # 10 clusters
    labels = clf.labels_
    correct_answer = 0
    for i in range(0,50):
        if labels[i] == 0:
            correct_answer+=1
    for i in range(50,100):
        if labels[i] == 1:
            correct_answer+=1
    for i in range(100,150):
        if labels[i] == 2:
            correct_answer+=1
        
    accuracy = (correct_answer/150)*100
    print(accuracy)


if __name__ == "__main__":
    main()
