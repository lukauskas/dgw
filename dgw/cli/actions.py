__author__ = 'saulius'
import argparse
import os

class StoreUniqueFilenameAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        seen = set()
        for value in values:
            if value in seen:
                raise argparse.ArgumentError(self, 'Dataset {0!r} is twice in the list'.format(value))
            if not os.path.isfile(value):
                raise argparse.ArgumentError(self, 'File not found {0!r}'.format(value))

            seen.add(value)

        setattr(namespace, self.dest, values)

class StoreFilenameAction(argparse.Action):

    def __call__(self, parser, namespace, value, option_string=None):
        if not os.path.isfile(value):
            raise argparse.ArgumentError(self, 'File not found {0!r}'.format(value))

        setattr(namespace, self.dest, value)

