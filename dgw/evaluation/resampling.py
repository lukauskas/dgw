from collections import defaultdict
import numpy as np
from ..data.containers import AlignmentsData
from ..dtw.distance import no_nans_len
import random
import pandas as pd

def extend_point(sequence, pos, extended_length):
    sequence = np.asarray(sequence)
    length = no_nans_len(sequence)

    new_seq = []
    for i, v in enumerate(sequence):
        if i >= length:
            break

        if i == pos:
            new_seq.extend([v] * extended_length)
        else:
            new_seq.append(v)

    return np.asarray(new_seq)


def shrink_to_a_single_point(sequence, pos, shrink_subset_len):
    sequence = np.asarray(sequence)
    length = no_nans_len(sequence)

    new_seq = []
    shrinked_part = []
    for i, v in enumerate(sequence):
        if i >= length:
            break
        if pos <= i < pos + shrink_subset_len:
            shrinked_part.append(v)
        elif i == pos + shrink_subset_len:
            new_seq.append(np.mean(shrinked_part, axis=0))
            shrinked_part = []
            new_seq.append(v)
        else:
            new_seq.append(v)

    if shrinked_part:
        new_seq.append(np.mean(shrinked_part, axis=0))

    return np.asarray(new_seq)

def mutate_sequence(sequence, prob):
    sequence = np.asarray(sequence)
    length = no_nans_len(sequence)

    prob_adjusted = prob / length

    for i in range(length):
        mutate = np.random.choice([True, False], p=[prob_adjusted, 1 - prob_adjusted])

        if mutate:
            sequence[i] = sequence[random.randint(0, length)]

    return sequence

def randomly_warp_sequence(sequence, max_number_of_extensions=10, max_number_of_shrinks=10,
                           max_extend_length=4, max_shrink_len=4, may_flip=True, mutation_prob=0.01):
    sequence = np.copy(sequence)

    n_ext = 0
    n_shrinks = 0

    number_of_extensions = random.randint(0, max_number_of_extensions)
    number_of_shrinks = random.randint(0, max_number_of_shrinks)

    while n_ext < number_of_extensions or n_shrinks < number_of_shrinks:
        if mutation_prob > 0:
            mutate_sequence(sequence, mutation_prob)

        length = no_nans_len(sequence)

        available_actions = []
        if n_ext < number_of_extensions:
            available_actions.append('extend')
        if n_shrinks < number_of_shrinks:
            available_actions.append('shrink')

        if len(available_actions) == 1:
            action = available_actions[0]
        else:
            action = random.choice(available_actions)

        if action == 'extend':
            pos = random.randint(0, length)
            k = random.randint(2, max_extend_length)
            sequence = extend_point(sequence, pos, k)
            n_ext += 1
        else:
            pos = random.randint(0, length)
            k = random.randint(2, max_shrink_len)
            sequence = shrink_to_a_single_point(sequence, pos, k)
            n_shrinks += 1

    if may_flip:
        flip = random.choice([True, False])
        if flip:
            sequence = sequence[::-1]

    return sequence

def generate_new_dataset(starting_dataset, desired_size, *args, **kwargs):
    new_dataset = {}

    item_counters = defaultdict(lambda: 0)
    for i in range(desired_size):
        random_ix = random.choice(starting_dataset.items)
        random_item = starting_dataset.ix[random_ix]

        transformed_item = pd.DataFrame(randomly_warp_sequence(random_item, *args, **kwargs))

        item_counters[random_ix] += 1
        item_id = item_counters[random_ix]

        new_dataset['{0}-{1}'.format(random_ix, item_id)] = transformed_item

    new_dataset = pd.Panel(new_dataset)
    return AlignmentsData(new_dataset)




