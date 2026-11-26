"""Exercise paired Linclust2 through the direct, easy and cascaded commands."""
import os
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile


def main(binary):
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(('STEAM_', 'MMSEQS_'))}

    def run(*args):
        result = subprocess.run([str(binary), *(str(arg) for arg in args)], env=env,
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout + result.stderr

    with tempfile.TemporaryDirectory() as directory:
        work = Path(directory)
        rng = random.Random(42)
        alphabet = 'ACDEFGHIKLMNPQRSTVWY'
        tea_sequence, related, unrelated = [''.join(rng.choices(alphabet, k=180)) for _ in range(3)]
        tea = work / 'tea.fasta'
        aa = work / 'aa.fasta'
        tea.write_text(''.join(f'>{name}\n{tea_sequence}\n' for name in ('one', 'two', 'three')))
        aa.write_text(f'>one\n{related}\n>two\n{unrelated}\n>three\n{related}\n')
        db = work / 'paired'
        run('createdb', tea, aa, db, '--threads', 2)
        settings = ('--linclust-version', 2, '--min-seq-id', '0.95', '-c', '0.9',
                    '-e', 'inf', '--cov-mode', 0, '--threads', 2)
        for command in ('linclust', 'easy-linclust', 'cluster'):
            result = work / command
            inputs = (tea, aa) if command == 'easy-linclust' else (db,)
            log = run(command, *inputs, result, work/f'tmp-{command}', *settings)
            assert log.count('Linclust2 paired TEA+AA alignment') == 2, log
            assert re.search(r'k-mers per sequence[ \t]+21\b', log), log
            assert re.search(r"^kmermatcher .* --kmer-per-seq 21 .* -k 0 ", log, re.MULTILINE), log
            if command == 'easy-linclust':
                table = Path(str(result) + '_cluster.tsv')
            else:
                table = work / f'{command}.tsv'
                run('createtsv', db, db, result, table, '--threads', 2)
            assignment = {member: representative for representative, member in
                          (row.split('\t') for row in table.read_text().splitlines())}
            assert assignment['one'] == assignment['three'] != assignment['two'], assignment

        keys = {name: int(key) for key, name, _ in
                (row.split('\t') for row in Path(str(db)+'.lookup').read_text().splitlines())}
        one, two, three = (keys[name] for name in ('one', 'two', 'three'))
        first_pass = work/'first-pass.tsv'
        first_pass.write_text(''.join(f'{q}\t{t}\n' for q, t in
                                      sorted([(one, one), (three, three), (three, two)])))
        first_clusters = work/'first-clusters'
        run('tsv2db', first_pass, first_clusters, '--output-dbtype', 6)
        reps = work/'reps'
        run('createsubdb', first_clusters, db, reps, '--subdb-mode', 1)
        candidates = work/'candidates.tsv'
        candidates.write_text(''.join(f'{q}\t{t}\t100\t0\n' for q, t in
                                      sorted([(one, one), (one, three), (three, three)])))
        pre = work/'pre'
        run('tsv2db', candidates, pre, '--output-dbtype', 7)
        for use_gate, expected in ((False, 1), (True, 2)):
            result = work/f'refine-{use_gate}'
            filters = ('--filter-cludb-file', first_clusters, '--filter-seqdb-file', db) if use_gate else ()
            run('align2clust', reps, pre, result, '--min-seq-id', '0.95', '-c', '0.9',
                '-e', 'inf', '--cov-mode', 0, '--cluster-mode', 0, '--threads', 2, *filters)
            actual = len(Path(str(result)+'.index').read_text().splitlines())
            assert actual == expected, (use_gate, actual, expected)

        # A central insertion breaks coverage on either ungapped seed diagonal.
        # The seeded block alignment must recover the same paired score and CIGAR.
        tea_insert = tea_sequence[:90] + 'HRLK' + tea_sequence[90:]
        aa_insert = related[:90] + 'MAVP' + related[90:]
        tea.write_text(f'>q\n{tea_sequence}\n>t\n{tea_insert}\n')
        aa.write_text(f'>q\n{related}\n>t\n{aa_insert}\n')
        indel_db = work/'indel'
        run('createdb', tea, aa, indel_db, '--threads', 2)
        lookup = Path(str(indel_db)+'.lookup')
        keys = {name: int(key) for key, name, _ in (row.split('\t') for row in lookup.read_text().splitlines())}
        q, t = keys['q'], keys['t']
        candidates = work/'indel-pref.tsv'
        candidates.write_text(''.join(f'{a}\t{b}\t100\t0\n' for a, b in sorted([(q,q), (q,t), (t,t)])))
        pre = work/'indel-pref'
        run('tsv2db', candidates, pre, '--output-dbtype', 7)
        clusters, exact = work/'indel-clusters', work/'indel-exact'
        options = ('-a', 1, '--alignment-mode', 3, '-c', '0.9', '-e', '0.01', '--threads', 2)
        log = run('align2clust', indel_db, pre, clusters, '--include-align-files', 1, *options)
        assert 'banded accepted=1' in log, log
        run('align', indel_db, indel_db, pre, exact, *options)
        rows = []
        for i, db_path in enumerate((Path(str(clusters)+'_aln'), exact)):
            hits = work/f'indel-{i}.tsv'
            run('convertalis', indel_db, indel_db, db_path, hits, '--threads', 2,
                '--format-output', 'query,target,raw,cigar,evalue')
            rows.append(next(row for row in hits.read_text().splitlines() if row.startswith('q\tt\t')))
        assert rows[0] == rows[1], rows


if __name__ == '__main__':
    main(Path(sys.argv[1]).resolve())
