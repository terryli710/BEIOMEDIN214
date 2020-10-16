"""
This file implements KNN method for classification of each sample as healthy or patient
Usage: python knn.py expfile sampfile
"""

import os
from typing import List
import sys
import numpy as np
import pandas as pd


class KNN:
    def __init__(self, expfile: str = None, sampfile: str = None) -> None:
        # Initializing knn
        self.expfile = None
        self.sampfile = None
        self.load_data(expfile, sampfile)

    def load_data(self, expfile: str = None, sampfile: str = None) -> None:
        """
        take the paths to the expression and sample file, read them in, and store within your KNN class
        :param expfile: expression file path
        :param sampfile: sample file path
        """
        # if either nothing inputs, or can't find file, pretend nothing happened
        if (expfile is None) or (sampfile is None) or \
                (not os.path.isfile(expfile)) or (not os.path.isfile(sampfile)):
            # print("File not found")
            return
        self.sampfile = pd.read_csv(sampfile, delimiter='\t', header=None)
        self.expfile = pd.read_csv(expfile, delimiter='\t', header=0, index_col=0)
        # reindex the file
        self.expfile = self.expfile[self.sampfile[0].values]
        return

    def get_assignments(self, k: int, fn: float) -> List[int]:
        """
        return the class assignments for all samples for given values of $k$ and $fn$ as a list of integer 0s and 1s
        :param k: k (parameter) from knn
        :param fn: positive ratio
        :return: list of assignments of class
        """
        return [self.predict(i, k, fn) for i in range(self.expfile.shape[1])]

    def predict(self, index: int, k: int, fn: float) -> int:
        """
        make prediction from k neighbors for the index-th sample
        :param index: index to make prediction for
        :param k: k (parameter) from knn
        :param fn: positive ratio
        :return: class {0,1} this belongs to
        """
        # calculate distance
        distances = self.calcDistance(index)
        # get nearest neighbors (indexes)
        kn_indexes_np = np.argsort(distances)[:k + 1]
        kn_indexes = list(kn_indexes_np[kn_indexes_np != index])
        # make prediction
        pred_value = np.mean(self.sampfile.iloc[kn_indexes, 1])
        # pred_value = np.mean(self.sampfile[self.sampfile[0].isin(self.expfile.columns[kn_indexes])][1])
        return int(pred_value > fn)

    def calcDistance(self, index: int) -> List[float]:
        """
        calculate the distance of point[index] to every point
        :param index: index to calculate the distance for
        :return: list of distances
        """
        target_exp = self.expfile.iloc[:, index]
        return [np.linalg.norm(target_exp - self.expfile.iloc[:, i]) for i in range(self.expfile.shape[1])]

    def calc_metrics(self, k: int, fn: float) -> List[float]:
        """
        should return a list of float values [sensitivity,specificity] of a KNN classifier using the given values of $k$ and $fn$.
        This method should run get_assignments at some point to initiate assignments used for performance metrics.
        :param k: k (parameter) from knn
        :param fn: positive ratio
        :return: a list of float values [sensitivity,specificity]
        """
        # get predictions
        y_preds = self.get_assignments(k, fn)
        y_true = list(self.sampfile[1])
        assert len(y_true) == len(y_preds), "y and yhat length diff"
        # calculate metrics
        tp = sum(1 for pred, true in zip(y_preds, y_true) if pred == true == 1)
        FN = sum(1 for pred, true in zip(y_preds, y_true) if pred == 0 and true == 1)
        tn = sum(1 for pred, true in zip(y_preds, y_true) if pred == true == 0)
        fp = sum(1 for pred, true in zip(y_preds, y_true) if pred == 1 and true == 0)
        # sensitivity = TP / (TP + FN)
        sens = tp / (tp + FN)
        # specificity = TN / (TN + FP)
        spec = tn / (tn + fp)
        return [round(sens, 2), round(spec, 2)]


def main(expfile: str, sampfile: str) -> None:
    """
    Main function for command line usage (should be modified to use)
    :param expfile: expression file path
    :param sampfile: sample file path
    :return: None
    """
    knn = KNN(expfile, sampfile)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
    # knn = KNN('GSE25628_filtered_expression.txt', 'GSE25628_samples.txt')
    # knn.calc_metrics(k=3, fn=0.5)
