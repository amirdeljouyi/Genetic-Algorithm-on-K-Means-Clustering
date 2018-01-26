import matplotlib.pyplot as plt
from matplotlib import style
style.use('ggplot')
import numpy as np
from sklearn.cluster import KMeans
from sklearn import preprocessing
import pandas as pd

data = pd.read_csv('result/norm_data.csv', header=None)

clf = KMeans(n_clusters=3)
clf.init='k-means++'
clf.fit_predict(data)
print(clf.cluster_centers_)
print(clf.labels_)


# correct = 0
# for i in range(len(x)):
#     predict_me = np.array(x[i].astype(float))
#     predict_me = predict_me.reshape(-1, len(predict_me))
#     prediction = clf.predict(predict_me)
#     if prediction[0] == y[i]:
#         correct+=1

# print(correct/len(x))