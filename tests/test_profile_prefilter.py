"""Paired PSSMs, not consensus sequences, must determine candidate ranking."""
import argparse
import os
from pathlib import Path
import random
import struct
import subprocess
import tempfile

from test_prefilter_setup import rows


ALPHABET = 'ACDEFGHIKLMNPQRSTVWY'


def profile(path, preferred, positive=32):
    # MMseqs rows: 20 signed scores in eighth-bits, query/consensus letters,
    # Neff and two gap bytes. Deliberately fixed metadata is NOT the consensus.
    data = bytearray()
    for letter in preferred:
        scores = [positive if aa == letter else -32 for aa in ALPHABET]
        data.extend(struct.pack('20b5B', *scores, 0, 0, 0, 0, 0))
    data.append(0)
    path.write_bytes(data)
    Path(str(path) + '.index').write_text(f'0\t0\t{len(data)}\n')
    Path(str(path) + '.dbtype').write_bytes(struct.pack('<I', 2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steam', type=Path, required=True)
    steam = parser.parse_args().steam.resolve()
    environment = {k: v for k, v in os.environ.items() if not k.startswith('STEAM_')}

    def run(*args, success=True):
        result = subprocess.run([str(steam), *map(str, args)], env=environment,
                                capture_output=True, text=True)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result.stdout + result.stderr

    with tempfile.TemporaryDirectory(prefix='steam-profiles-', dir=os.environ.get('TMPDIR')) as tmp:
        work = Path(tmp)
        rng = random.Random(7)
        tea = ''.join(rng.choices(ALPHABET, k=48))
        alternate = tea[:24] + ''.join(rng.choices(ALPHABET, k=24))
        aa = ['AC' * 24, 'CA' * 24]
        for channel in ('aa', 'tea'):
            target = work / (channel + '-target')
            tea_targets = [tea, tea] if channel == 'aa' else [tea, alternate]
            aa_targets = aa if channel == 'aa' else [aa[0], aa[0]]
            for track, sequences in [('tea', tea_targets), ('aa', aa_targets)]:
                (work / track).write_text(''.join(
                    f'>target{i}\nGGG{seq}\n' for i, seq in enumerate(sequences)))
            run('createdb', work / 'tea', work / 'aa', target, '--threads', 1, '-v', 1)
            lookup = {int(line.split()[0]): line.split()[1] for line in
                      Path(str(target) + '.lookup').read_text().splitlines()}
            for winner in (0, 1):
                query = work / f'{channel}-query{winner}'
                profile(query, tea_targets[winner])
                profile(Path(str(query) + '_aa'), aa_targets[winner], positive=24)
                for splits in (1, 2):
                    output = work / f'{channel}-result{winner}-{splits}'
                    run('prefilter', query, target, output, '--max-seqs', 1,
                        '--split', splits, '--split-mode', 0, '--threads', 1, '-v', 1)
                    hits = rows(output)[0].rstrip(b'\0\n').splitlines()
                    assert len(hits) == 1, hits
                    key, score, diagonal = map(int, hits[0].split())
                    assert lookup[key] == f'target{winner}', (channel, hits)
                    assert score > 255 and diagonal == -3, (channel, hits)
                print(f'PASS {channel} profile chooses target{winner}; splits agree', flush=True)
            for options in [('--exact-kmer-matching', 1), ('--ungapped-tea-aa', 0),
                            ('--aa-weight', 3), ('--tea-scale', 2)]:
                log = run('prefilter', query, target, work / 'rejected', *options,
                          '--threads', 1, success=False)
                assert 'channel weights must already be included' in log, log
        print('PASS profile-only option checks', flush=True)
        log = run('prefilter', query, target, work / 'reduced', '--alph-size', 13,
                  '--threads', 1, success=False)
        assert 'full 21-symbol alphabet' in log, log
        run('createindex', target, work / 'index-tmp', '--alph-size', 13,
            '-k', 5, '--threads', 1, '-v', 1)
        log = run('prefilter', query, Path(str(target) + '.idx'), work / 'reduced-index',
                  '--threads', 1, success=False)
        assert 'full 21-symbol target index' in log, log
        print('PASS reduced alphabet and index rejection', flush=True)


if __name__ == '__main__':
    main()
