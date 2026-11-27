"""Native observed words drive both seed membership and Markov correction."""
import argparse
import os
from pathlib import Path
import random
import subprocess
import tempfile
from test_profile_prefilter import profile, ALPHABET
from test_prefilter_setup import rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steam', type=Path, required=True)
    steam = parser.parse_args().steam.resolve()
    environment = {k: v for k, v in os.environ.items() if not k.startswith('STEAM_')}
    def run(*arguments, success=True):
        result = subprocess.run([str(steam), *map(str, arguments)], env=environment,
            capture_output=True, text=True)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
    with tempfile.TemporaryDirectory(prefix='steam-observed-', dir=os.environ.get('TMPDIR')) as temporary:
        work = Path(temporary)
        rng = random.Random(71)
        aa, first, observed = [''.join(rng.choices(ALPHABET, k=72)) for _ in range(3)]
        (work / 'aa.fa').write_text(f'>first\n{aa}\n>observed\n{aa}\n')
        (work / 'tea.fa').write_text(f'>first\n{first}\n>observed\n{observed}\n')
        target, query = work / 'target', work / 'query'
        run('createdb', work / 'tea.fa', work / 'aa.fa', target, '--threads', 2, '-v', 1)
        profile(query, first)
        profile(work / 'query_aa', aa, positive=80)
        catalogue = ''.join(f'0\t{i}\t{observed[i:i+7]}\n' for i in range(len(observed)-6))
        (work / 'query.words').write_text(catalogue)
        lookup = {int(line.split()[0]): line.split()[1] for line in
            (work / 'target.lookup').read_text().splitlines()}
        candidate_rows = []
        for sensitivity in (1, 8):
            pref = work / f'pref-{sensitivity}'
            run('prefilter', query, target, pref, '-k', 6, '--spaced-kmer-pattern', '1110111',
                '-s', sensitivity, '--max-seqs', 100, '--threads', 2, '-v', 1)
            candidate_rows.append(rows(pref))
            hits = candidate_rows[-1][0].rstrip(b'\0\n').splitlines()
            assert len(hits) == 1 and lookup[int(hits[0].split()[0])] == 'observed', hits
        assert candidate_rows[0] == candidate_rows[1]
        aligned = []
        for correction in (0, 1):
            output = work / f'align-{correction}'
            run('align', query, target, work / 'pref-1', output,
                '--seed-pattern', '1110111', '--seed-correction', correction, '-a', 1,
                '--threads', 2, '-v', 1)
            aligned.append(rows(output)[0].rstrip(b'\0\n').split())
        assert int(aligned[1][1]) > int(aligned[0][1]), aligned
        assert aligned[0][0] == aligned[1][0] and aligned[0][2:] == aligned[1][2:]
        (work / 'query.words').write_text(catalogue + catalogue)
        run('align', query, target, work / 'pref-1', work / 'duplicates',
            '--seed-pattern', '1110111', '--seed-correction', 1, '-a', 1, '--threads', 2, '-v', 1)
        assert rows(work / 'duplicates') == rows(work / 'align-1')
        (work / 'query.words').write_text('')
        run('prefilter', query, target, work / 'empty', '-k', 6,
            '--spaced-kmer-pattern', '1110111', '--threads', 2, '-v', 1)
        assert not rows(work / 'empty')[0].rstrip(b'\0\n')
        run('prefilter', query, target, work / 'wrong-pattern', '-k', 5,
            '--spaced-kmer-pattern', '1101101', '--threads', 2, '-v', 1, success=False)
    print('PASS observed membership, sensitivity invariance, Markov correction, duplicates and empty catalogue')


if __name__ == '__main__':
    main()
