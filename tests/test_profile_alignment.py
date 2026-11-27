"""Paired profile scores and seed coordinates, including a shared gap path."""
import argparse
import math
import os
from pathlib import Path
import random
import struct
import subprocess
import tempfile

from test_profile_prefilter import ALPHABET, profile
from test_prefilter_setup import rows


def reference_score(aa, tea, target_aa, target_tea):
    """Scalar local affine-gap reference in the native half-bit score units."""
    width = len(target_aa) + 1
    previous = [0] * width
    deletion = [-10**9] * width
    best = 0
    for q in range(len(aa)):
        current = [0] * width
        insertion = -10**9
        for t in range(1, width):
            deletion[t] = max(previous[t] - 24, deletion[t] - 2)
            insertion = max(current[t - 1] - 24, insertion - 2)
            score = (6 if aa[q] == target_aa[t - 1] else -8)
            score += 8 if tea[q] == target_tea[t - 1] else -8
            current[t] = max(0, previous[t - 1] + score, deletion[t], insertion)
            best = max(best, current[t])
        previous = current
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steam', type=Path, required=True)
    steam = parser.parse_args().steam.resolve()
    environment = {k: v for k, v in os.environ.items() if not k.startswith('STEAM_')}

    def run(*args, success=True):
        result = subprocess.run([str(steam), *map(str, args)], env=environment,
                                text=True, capture_output=True)
        assert (result.returncode == 0) == success, result.stdout + result.stderr
        return result.stdout + result.stderr

    with tempfile.TemporaryDirectory(prefix='steam-profile-align-', dir=os.environ.get('TMPDIR')) as tmp:
        work = Path(tmp)
        rng = random.Random(31)
        aa = ''.join(rng.choices(ALPHABET, k=48))
        tea = ''.join(rng.choices(ALPHABET, k=48))
        query = work / 'query'
        profile(query, tea)
        profile(Path(str(query) + '_aa'), aa, positive=24)
        for gap in ('', 'GGG'):
            target = work / f'target{len(gap)}'
            targets = ['GGG' + seq[:24] + gap + seq[24:] + 'FFF' for seq in (aa, tea)]
            for track, sequence in zip(('aa', 'tea'), targets):
                (work / track).write_text(f'>target\n{sequence}\n')
            run('createdb', work / 'tea', work / 'aa', target, '--threads', 1, '-v', 1)
            prefilter = work / 'prefilter'
            prefilter.write_bytes(b'0\t672\t-3\n\0')
            Path(str(prefilter) + '.index').write_text(f'0\t0\t{prefilter.stat().st_size}\n')
            Path(str(prefilter) + '.dbtype').write_bytes(struct.pack('<I', 7))
            output = work / f'alignment{len(gap)}'
            run('align', query, target, prefilter, output, '--threads', 1,
                '--max-seq-len', 32, '--alignment-mode', 3, '-a', 1, '-v', 1)
            fields = rows(output)[0].rstrip(b'\0\n').split()
            assert int(fields[1]) == reference_score(aa, tea, *targets), fields
            assert math.isinf(float(fields[3])), fields
            assert list(map(int, fields[4:10])) == [0, 47, 48, 3, 50 + len(gap), 54 + len(gap)], fields
            assert fields[10] == (b'24M3D24M' if gap else b'48M'), fields
            print(f'PASS paired profile score, coordinates and traceback with gap={len(gap)}', flush=True)
        # Correction-only mode preserves ranking-only E-values. Deliberately
        # bogus consensus metadata must not determine the supported words.
        corrected = work / 'corrected'
        run('align', query, target, prefilter, corrected, '--threads', 1,
            '--seed-correction', 1, '-s', 10, '-a', 1, '-v', 1)
        corrected_fields = rows(corrected)[0].rstrip(b'\0\n').split()
        assert int(corrected_fields[1]) > int(fields[1]), (fields, corrected_fields)
        assert corrected_fields[4:] == fields[4:]
        assert math.isinf(float(corrected_fields[3]))
        # Reuse sequence coefficients explicitly as an experimental filter.
        calibrated = work / 'calibrated'
        run('align', query, target, prefilter, calibrated, '--threads', 1,
            '--seed-correction', 1, '-s', 10, '-e', 'inf', '-a', 1, '-v', 1)
        calibrated_fields = rows(calibrated)[0].rstrip(b'\0\n').split()
        score = float(calibrated_fields[1])
        assert calibrated_fields[1] == corrected_fields[1]
        expected_e = 10**(0.2839183138420103 - 0.005633978159801299*score
                          + (0.005633978159801299-0.0014312498238873948)*max(0, score-586))
        assert math.isclose(float(calibrated_fields[3]), expected_e, rel_tol=0.01)
        for scale in (0.5, 2.):
            filtered = work / f'filtered-{scale}'
            run('align', query, target, prefilter, filtered, '--threads', 1,
                '--seed-correction', 1, '-s', 10, '-e', expected_e*scale, '-a', 1, '-v', 1)
            assert (rows(filtered)[0] != b'\0') == (scale > 1)
        print('PASS profile correction, finite E-values and filtering', flush=True)
        for options in [('--aa-weight', 3),
                        ('--min-seq-id', 0.5), ('--min-tfident', 0.5)]:
            log = run('align', query, target, prefilter, work / 'rejected', *options,
                      '--threads', 1, success=False)
            assert 'sequence weighting and identity' in log, log
        profile(query, tea, positive=-32)
        profile(Path(str(query) + '_aa'), aa, positive=-32)
        output = work / 'impossible'
        run('align', query, target, prefilter, output, '--threads', 1, '-a', 1, '-v', 1)
        assert rows(output)[0] == b'\0'
        print('PASS incompatible filters and impossible alignment', flush=True)


if __name__ == '__main__':
    main()
