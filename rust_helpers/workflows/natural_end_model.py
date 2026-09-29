"""Explicit, synthetic library-recovery model; not fitted to experimental data."""
import numpy as np

LENGTHS = np.array([27,28,29,30,31])
LENGTH_PROB = np.array([.10,.25,.35,.20,.10])
MODAL_OFFSETS = np.array([14,15,15,16,16])
JITTER = np.array([-1,0,1])
JITTER_PROB = np.array([.08,.84,.08])
SCENARIOS = ['control','five_prime','three_prime','both_ends']


def dinucleotide_weights(kmers, end):
    if end == 'five':
        first = {'A':.95,'C':.85,'G':1.35,'T':1.0}
        second = {'A':1.05,'C':1.10,'G':.95,'T':.90}
        interaction = {'GG':1.15,'GC':1.12,'TA':.80}
    elif end == 'three':
        first = {'A':.95,'C':1.05,'G':1.10,'T':.90}
        second = {'A':.90,'C':1.40,'G':1.0,'T':.85}
        interaction = {'TC':1.20,'GC':1.10,'AA':.90}
    else:
        raise ValueError('end must be five or three')
    if any(len(k) != 2 or set(k) - set('ACGT') for k in kmers):
        raise ValueError('Expected DNA dinucleotides in biological orientation')
    return np.array([first[k[0]] * second[k[1]] * interaction.get(k,1.) for k in kmers])


def noisy_site_weights(base, rng):
    base = np.asarray(base, dtype=float)
    if not np.all(np.isfinite(base) & (base > 0)):
        raise ValueError('Base weights must be finite and positive')
    patches = np.zeros(len(base))
    innovations = rng.normal(0,.18,len(base))
    for i in range(len(base)):
        patches[i] = (.85 * patches[i-1] if i else 0) + innovations[i]
    weights = base * np.exp(patches) * rng.gamma(shape=6,scale=1/6,size=len(base))
    weights[rng.random(len(base)) < .07] *= .05
    weights[rng.random(len(base)) < .003] *= 5
    return weights


def fragment_candidates(sequence, sites, site_weights):
    starts=[]; positions=[]; widths=[]; offsets=[]; weights=[]
    for length, length_prob, modal in zip(LENGTHS,LENGTH_PROB,MODAL_OFFSETS):
        for jitter, jitter_prob in zip(JITTER,JITTER_PROB):
            offset = int(modal+jitter)
            candidate_start = sites - offset
            valid = (candidate_start >= 0) & (candidate_start+length <= len(sequence))
            starts.extend(candidate_start[valid]); positions.extend(sites[valid])
            widths.extend([int(length)] * int(valid.sum()))
            offsets.extend([offset] * int(valid.sum()))
            weights.extend(site_weights[valid] * length_prob * jitter_prob)
    if not starts:
        raise ValueError('No complete fragment fits')
    starts=np.asarray(starts,dtype=int); positions=np.asarray(positions,dtype=int)
    widths=np.asarray(widths,dtype=int); offsets=np.asarray(offsets,dtype=int)
    sequences=[sequence[s:s+l] for s,l in zip(starts,widths)]
    first=[s[:2] for s in sequences]; last=[s[-2:] for s in sequences]
    return {'start':starts,'site':positions,'length':widths,'offset':offsets,
            'weight':np.asarray(weights), 'sequence':sequences,'five_kmer':first,'three_kmer':last,
            'five_weight':dinucleotide_weights(first,'five'),
            'three_weight':dinucleotide_weights(last,'three')}


def sample_library(candidates, count, scenario, rng):
    if scenario not in SCENARIOS:
        raise ValueError('Unknown scenario')
    if count < 0 or int(count) != count:
        raise ValueError('Read count must be a nonnegative integer')
    probability = candidates['weight'].copy()
    if scenario in ['five_prime','both_ends']:
        probability *= candidates['five_weight']
    if scenario in ['three_prime','both_ends']:
        probability *= candidates['three_weight']
    return rng.multinomial(int(count), probability/probability.sum())
