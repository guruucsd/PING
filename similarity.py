"""
Script for running similarity matrix comparisons PING data.
"""
from matplotlib import pyplot as plt

from ping.similarity import compare_all_similarity_matrices


if __name__ == '__main__':
    compare_all_similarity_matrices()
    plt.show()
