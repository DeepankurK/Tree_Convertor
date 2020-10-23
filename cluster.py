from sklearn_extra.cluster import KMedoids
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn import metrics
from sklearn.cluster import Birch
from sklearn.cluster import KMeans
from sklearn.cluster import OPTICS
from sklearn.cluster import SpectralClustering

def affinity(distmat):
    labels=[]
    cluster = AffinityPropagation().fit(distmat)
    print(cluster.labels_)
    #avg_score = metrics.silhouette_score(distmat, cluster.labels_)
    labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels


def agglo(distmat,k):
    max_score=0
    labels=[-1]
    for i in range(2,k):
        #print(i)
        cluster =  AgglomerativeClustering(n_clusters=i).fit(distmat)
        if -1 in cluster.labels_:
            continue
        avg_score = metrics.silhouette_score(distmat, cluster.labels_)
        if max_score<avg_score:
            max_score=avg_score
            labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels

def birch(distmat,k):
    max_score=0
    labels=[-1]
    for i in range(2,k):
        cluster =  Birch(n_clusters=i).fit(distmat)
        if -1 in cluster.labels_:
            continue
        avg_score = metrics.silhouette_score(distmat, cluster.labels_)
        if max_score<avg_score:
            max_score=avg_score
            labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels


def kmeans(distmat,k):
    max_score=0
    labels=[]
    for i in range(2,k):
        cluster = KMeans(n_clusters=i).fit(distmat)
        if -1 in cluster.labels_:
            continue
        avg_score = metrics.silhouette_score(distmat, cluster.labels_)
        if max_score<avg_score:
            max_score=avg_score
            labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels

def kmediods(distmat,k):
    max_score=0
    labels=[]
    for i in range(2,k):
        cluster = KMedoids(n_clusters=i, random_state=0).fit(distmat)
        if -1 in cluster.labels_:
            continue
        avg_score = metrics.silhouette_score(distmat, cluster.labels_)
        if max_score<avg_score:
            max_score=avg_score
            labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels


def optics(distmat):
    labels=[]
    cluster = OPTICS(min_samples=2).fit(distmat)
    #print(cluster.labels_)
    labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels

def spectral(distmat,k):
    max_score=0
    labels=[]
    for i in range(2,k):
        cluster =  SpectralClustering(n_clusters=i).fit(distmat)
        if -1 in cluster.labels_:
            continue
        avg_score = metrics.silhouette_score(distmat, cluster.labels_)
        if max_score<avg_score:
            max_score=avg_score
            labels=cluster.labels_
    if -1 in labels:
        print("Wrong Clustering selected")
        exit()
    return labels
