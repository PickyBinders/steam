"""A multi-round cascade starts with seed counts, then uses combined scoring."""
import argparse
import os
from pathlib import Path
import random
import re
import shlex
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steam', type=Path, required=True)
    args = parser.parse_args()
    steam = args.steam.resolve()
    environment = {k: v for k, v in os.environ.items()
                   if not k.startswith(('STEAM_', 'MMSEQS_'))}

    def run(*arguments, success=True):
        result = subprocess.run([str(steam), *map(str, arguments)], env=environment,
                                capture_output=True, text=True)
        log = result.stdout + result.stderr
        assert (result.returncode == 0) == success, log
        return log

    with tempfile.TemporaryDirectory(prefix='cluster-seed-count-',
                                     dir=os.environ.get('TMPDIR')) as temporary:
        work = Path(temporary)
        rng = random.Random(42)
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        names = {f'seq{i}' for i in range(64)}
        for kind in ('tea', 'aa'):
            fasta = ''.join(f'>seq{i}\n' + ''.join(rng.choices(alphabet, k=100)) + '\n'
                            for i in range(64))
            (work / f'{kind}.fasta').write_text(fasta)
        db = work / 'db'
        run('createdb', work / 'tea.fasta', work / 'aa.fasta', db, '--threads', 2)
        for rounds, diagonals, caps in ((1, [1], [1000]),
                                       (2, [0, 1], [100, 1000]),
                                       (3, [0, 1, 1], [100, 200, 1000])):
            output = work / f'clu{rounds}'
            log = run('cluster', db, output, work / f'tmp{rounds}',
                      '--cluster-steps', rounds, '--threads', 2, '-e', '0.01', '-c', 0)
            seed_matchers = [shlex.split(line) for line in log.splitlines()
                             if line.startswith('kmermatcher ')]
            assert seed_matchers and all(
                cmd[cmd.index('--kmer-per-seq') + 1] == '300'
                for cmd in seed_matchers), 'Linclust must retain 300 seeds per sequence'
            children = [shlex.split(line) for line in log.splitlines()
                        if line.startswith('prefilter ') and re.search(r'/pref_step\d+ ', line)]
            actual_diagonals = [int(cmd[cmd.index('--diag-score') + 1]) for cmd in children]
            actual_caps = [int(cmd[cmd.index('--max-seqs') + 1]) for cmd in children]
            assert actual_diagonals == diagonals, ('cascade diagonal modes', rounds, actual_diagonals)
            assert actual_caps == caps, ('cascade candidate caps', rounds, actual_caps)
            assert log.count('Combined TEA+AA ungapped diagonal scoring enabled') == sum(diagonals), log
            tsv = work / f'clu{rounds}.tsv'
            run('createtsv', db, db, output, tsv, '--threads', 2)
            members = [line.split('\t')[1] for line in tsv.read_text().splitlines()]
            assert len(members) == len(names) and set(members) == names
            print(f'PASS {rounds}-round cascade: diagonal modes {diagonals}, caps {caps}', flush=True)
        # Standalone combined prefilter still requires diagonal scoring.
        rejected = run('prefilter', db, db, work / 'invalid', '--diag-score', 0,
                       '--threads', 2, success=False)
        assert 'requires --diag-score 1' in rejected, rejected
        print('PASS standalone combined-scoring validation remains enforced', flush=True)


if __name__ == '__main__':
    main()
