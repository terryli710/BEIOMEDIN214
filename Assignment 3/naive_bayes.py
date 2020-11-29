''' This file use naive bayes model to predict protein binding site using features.csv '''

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import BernoulliNB
from sklearn.metrics import plot_confusion_matrix
import matplotlib.pyplot as plt
from sklearn.feature_selection import mutual_info_classif
import numpy as np


def load_data(file_name):
    ''' load and split data'''
    data = pd.read_csv(file_name)
    x_binary = data.iloc[:,:-1]=='y'
    y = data.iloc[:,-1]
    xtr, xts, ytr, yts = train_test_split(x_binary, y, test_size=0.2, random_state=1115)
    return xtr, xts, ytr, yts


def nBModel():
    ''' train and evaluate bayes model '''
    bnb = BernoulliNB()
    xtr, xts, ytr, yts = load_data('features.csv')
    bnb.fit(xtr, ytr)
    preds = bnb.predict(xts, yts)
    plot_confusion_matrix(bnb, xts, yts)
    plt.show()
    importace = mutual_info_classif(pd.concat([xtr, xts], axis=0), pd.concat([ytr, yts], axis=0),
                                    discrete_features=True)
    np.argsort(importace)