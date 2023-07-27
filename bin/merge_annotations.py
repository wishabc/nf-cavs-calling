import pandas as pd
import sys
from tqdm import tqdm

tqdm.pandas()


mutations = ['C>A', 'C>T', 'C>G', 'T>A']
_comp = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }

def revcomp(s):
    return ''.join([_comp.get(x) for x in s[::-1]])

def palindromic(ref, alt):
    return {ref, alt} == {'A', 'T'} or {ref, alt} == {'G', 'C'}


def get_mutation_stats(row):
    ref = row['ref']
    alt = row['alt']
    sequence = row['sequence']
    preceding, ref_at_seq, following = sequence[:20], sequence[20], sequence[21:]
    assert ref_at_seq == ref

    flank_len = len(preceding)
    if palindromic(ref, alt):
        initial_fwd = f'{ref}>{alt}' in mutations
        ref_orient = True
        fwd = initial_fwd  # if cycle dosen't break
        palindromic_res = [True] + [False] * flank_len
        for i in range(flank_len):
            left = preceding[-i - 1]
            right = following[i]
            if not palindromic(left, right):
                fwd = 'T' not in {left, right} and not (left == right == 'G')
                ref_orient = not (fwd ^ initial_fwd)
                break
            palindromic_res[i + 1] = True
        if initial_fwd:
            sub = f'{ref}>{alt}'
        else:
            sub = f'{alt}>{ref}'
        
    else:        
        palindromic_res = [False] * (flank_len + 1)
        for sub, ref_orient, fwd in [
            (f'{ref}>{alt}', True, True),
            (f'{_comp.get(ref)}>{_comp.get(alt)}', True, False),
            (f'{alt}>{ref}', False, True),
            (f'{_comp.get(alt)}>{_comp.get(ref)}', False, False),
        ]:
            if sub in mutations:
                break
            
    if not fwd:
        preceding, following = revcomp(following), revcomp(preceding)
        
    return pd.Series(f"{preceding[-1]}[{sub}]{following[0]}" + [sub, fwd, ref_orient])


def main(context, mutation_rates):
    result = context.merge(
        mutation_rates, how='left'
    )
    result[['signature1', 'sub', 'fwd', 'ref_orient']] = result.progress_apply(
        get_mutation_stats, axis=1
    )
    result['cpg'] = ((result['sub'] != 'A>T') & (result['1'] == 'G')) | ((result['sub'] == 'C>G') & (result['-1'] == 'C'))

    return result


if __name__ == '__main__':
    context = pd.read_table(sys.argv[1])
    mutation_rates = pd.read_table(sys.argv[2])
    
    main(context, mutation_rates).to_csv(
        sys.argv[3], sep='\t', index=False
    )
