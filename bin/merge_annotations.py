import pandas as pd
import sys


def main(context, phenotypes, mutation_rates):
    return context.merge(phenotypes, how='left').merge(mutation_rates, how='left')


if __name__ == '__main__':
    context = pd.read_table(sys.argv[1])
    mutation_rates = pd.read_table(sys.argv[2])
    phenotypes = pd.read_table(sys.argv[3])
    
    main(context, mutation_rates, phenotypes).to_csv(sys.argv[4], sep='\t', index=False)
