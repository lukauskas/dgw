__author__ = 'saulius'

from sklearn import mixture
import pandas as pd

def gmm_cluster(data, n_components_range, return_gmm=False):

    min_gmm = None
    min_bic = float('inf')

    for n in n_components_range:
        gmm = mixture.GMM(n_components=n, covariance_type='diag')
        gmm.fit(data.values)

        bic = gmm.bic(data.values)

        if bic < min_bic:
            min_bic = bic
            min_gmm = gmm

    predictions = min_gmm.predict(data.values)

    predictions_pd = pd.Series(predictions, index=data.index)

    if return_gmm:
        return predictions_pd, min_gmm
    else:
        return predictions_pd
