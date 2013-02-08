__author__ = 'saulius'
import pandas as pd

class AggregatedAlignmentsPanel(pd.Panel):
    def mean(self, axis='items', skipna=True):
        # Override axis parameter in the pd.Panel mean function
        return super(AggregatedAlignmentsPanel, self).mean(axis=axis, skipna=skipna)

