#!/usr/bin/env python3
import argparse
import math
import sys
from typing import Callable, List

import numpy as np


DEFAULT_ALPHA = 1.0
DEFAULT_BETA = 0.25

DEFAULT_STEP = 0.01
DEFAULT_MAX_T = 1.0


def pK(a: str, b: str, t: float, /,
       alpha: float = DEFAULT_ALPHA, beta: float = DEFAULT_BETA):
    def is_transition(n1: str, n2: str) -> bool:
        return (
                n1 != n2
                and
                (
                    (n1 in {'A', 'T'} and n2 in {'A', 'T'})
                    or
                    (n1 in {'C', 'G'} and n2 in {'C', 'G'})
                )
        )

    def is_transversion(n1: str, n2: str) -> bool:
        return (
                n1 != n2
                and
                not is_transition(n1, n2)
        )

    same = 0.25 + 0.25 * math.e ** (-4 * beta * t) + 0.5 * math.e ** (-2 * (alpha + beta) * t)
    transition = 0.25 + 0.25 * math.e ** (-4 * beta * t) - 0.5 * math.e ** (-2 * (alpha + beta) * t)
    transversion = 0.25 - 0.25 * math.e ** (-4 * beta * t)

    same_n = 0
    transi_n = 0
    transv_n = 0
    for a1, a2 in zip(a, b):
        if a1 == '-' or a2 == '-':
            continue
        elif a1 == a2:
            same_n += 1
        elif is_transition(a1, a2):
            transi_n += 1
        elif is_transversion(a1, a2):
            transv_n += 1
        else:
            print(f'[ERROR] Should not be here, a1 = {a1}, a2 = {a2}')
            sys.exit(1)

    probability = (same ** same_n) * (transition ** transi_n) * (transversion ** transv_n)
    return probability


def pJC(a: str, b: str, t: float, /,
        alpha: float = DEFAULT_ALPHA):
    same = 0.25 + 0.75 * math.e ** (-4 * alpha * t)
    trans = 0.25 - 0.25 * math.e ** (-4 * alpha * t)

    same_n = 0
    trans_n = 0
    for a1, a2 in zip(a, b):
        if a1 == '-' or a2 == '-':
            continue
        elif a1 == a2:
            same_n += 1
        elif a1 != a2:
            trans_n += 1
        else:
            print(f'[ERROR] Should not be here, a1 = {a1}, a2 = {a2}')
            sys.exit(1)

    probability = (same ** same_n) * (trans ** trans_n)
    return probability


def plot(data: list, step: float, /,
         title: str, output_file: str = ''):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    def format_fn(tick_val, tick_pos):
        return tick_val * step
    fig, ax = plt.subplots()
    ax.xaxis.set_major_formatter(FuncFormatter(format_fn))
    ax.plot(data)
    plt.title(title)
    if output_file:
        fig.savefig(output_file, dpi=fig.dpi)
    plt.show()


def plot_K(a: str, b: str, /,
           alpha: float = DEFAULT_ALPHA, beta: float = DEFAULT_BETA, *,
           max_t: float = 10.0, step: float = DEFAULT_STEP,
           output_file: str = ''):
    props = [pK(a, b, t, alpha, beta) for t in np.arange(0, max_t, step)]
    tmax = optT_p(props, step)
    plot(props, step, title=f'Kimura model: alpha = {alpha}, beta = {beta}, tmax = {tmax}', output_file=output_file)


def plot_JC(a: str, b: str, /,
            alpha: float = DEFAULT_ALPHA, *,
            max_t: float = DEFAULT_MAX_T, step: float = DEFAULT_STEP,
            output_file: str = ''):
    props = [pJC(a, b, t, alpha) for t in np.arange(0, max_t, step)]
    tmax = optT_p(props, step)
    plot(props, step, title=f'Jukes-Cantor model: alpha = {alpha}, tmax = {tmax}', output_file=output_file)


def optT_p(props: List[float], step: float) -> float:
    tmax_index = props.index(max(props))
    return tmax_index * step


def optT(a: str, b: str, p: Callable, /, *,
         max_t: float = DEFAULT_MAX_T, step: float = DEFAULT_STEP):
    props = [p(a, b, t) for t in np.arange(0, max_t, step)]
    return optT_p(props, step)


def read_sequences(filename: str) -> List[str]:
    from Bio import SeqIO
    seqs = [record.seq for record in SeqIO.parse(filename, 'fasta')]
    return seqs


def get_two_seqs(filename: str) -> List[str]:
    seqs = read_sequences(filename)
    if len(seqs) < 2:
        print(f'[ERROR] Needs at least 2 sequences, found = {len(seqs)}')
        sys.exit(1)
    return seqs[:2]


def simmatrix(filename: str, p: Callable, /, *,
              output_file: str = ''):
    seqs = read_sequences(filename)
    smatrix = np.array([
            np.array([
                optT(seq1, seq2, p)
                for seq2 in seqs
            ])
            for seq1 in seqs
    ])
    if output_file:
        np.save(output_file, smatrix)
    return smatrix


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Probability of sequences transforms with evolution model in time')
    parser.add_argument('file', type=str, help='fasta file with aligned sequences, 2 sequences for \'tmax_plot\'')
    parser.add_argument('mode', type=str, choices=['tmax_plot', 'simmatrix'], help='mode')
    parser.add_argument('-m', '--model', type=str, default='jukes-cantor', choices=['kimura', 'jukes-cantor'], help='evolution model')
    parser.add_argument('-a', '--alpha', type=float, default=DEFAULT_ALPHA, help='alpha param')
    parser.add_argument('-b', '--beta', type=float, default=DEFAULT_BETA, help='beta param (for Kimura model)')
    parser.add_argument('-t', '--time', type=float, default=DEFAULT_MAX_T, help='max t')
    parser.add_argument('-s', '--step', type=float, default=DEFAULT_STEP, help='step')
    parser.add_argument('-o', '--output', type=str, help='output file')
    args = parser.parse_args()

    output_file = args.output if args.output else ''

    if args.mode == 'tmax_plot':
        seqs = get_two_seqs(args.file)
        if args.model == 'jukes-cantor':
            plot_JC(seqs[0], seqs[1], args.alpha, max_t=args.time, step=args.step, output_file=output_file)
        elif args.model == 'kimura':
            plot_K(seqs[0], seqs[1], args.alpha, args.beta, max_t=args.time, step=args.step, output_file=output_file)
        else:
            print(f'[ERROR] Should not be here, model = {args.model}')
            sys.exit(1)

    if args.mode == 'simmatrix':
        if args.model == 'jukes-cantor':
            evo_func = pJC
        elif args.model == 'kimura':
            evo_func = pK
        else:
            print(f'[ERROR] Should not be here, model = {args.model}')
            sys.exit(1)
        smatrix = simmatrix(args.file, evo_func, output_file=output_file)
        print(smatrix)
