import pandas as pd
import sys


def main(context, mutation_rates):
    return context.merge(mutation_rates, how='left')


if __name__ == '__main__':
    context = pd.read_table(sys.argv[1])
    mutation_rates = pd.read_table(sys.argv[2])
    
    main(context, mutation_rates).to_csv(sys.argv[3], sep='\t', index=False)
