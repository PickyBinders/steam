"""Sequence-only search, indexed inputs and TEA-identity clustering checks."""

import argparse
from pathlib import Path
import random
import re
import shlex
import subprocess
import tempfile


def run(*args):
    result = subprocess.run([str(arg) for arg in args], capture_output=True, text=True)
    if result.returncode:
        raise AssertionError(f"{' '.join(map(str, args))}\n{result.stdout}\n{result.stderr}")
    return result


def sequence(path):
    return "".join(path.read_text().splitlines()[1:])


def rejected(*args, message):
    result = subprocess.run(list(map(str, args)), capture_output=True, text=True)
    assert result.returncode and message in result.stdout + result.stderr, result


def records(prefix):
    parts = [prefix] if prefix.exists() else sorted(
        (p for p in prefix.parent.glob(prefix.name + '.*') if p.suffix[1:].isdigit()),
        key=lambda p: int(p.suffix[1:]))
    data = b''.join(p.read_bytes() for p in parts)
    return {int(key): data[int(offset):int(offset) + int(length)]
            for key, offset, length in (line.split() for line in
                Path(str(prefix) + '.index').read_text().splitlines())}


def check_cascade_calibration(steam, work):
    # Unrelated synthetic sequences keep representatives in every cascade stage.
    rng = random.Random(42)
    names = {f"sequence{i}" for i in range(16)}
    inputs = []
    for alphabet in ("aa", "tea"):
        fasta = work / f"cascade-{alphabet}.fasta"
        fasta.write_text("".join(
            f">{name}\n{''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=128))}\n"
            for name in sorted(names)))
        inputs.append(fasta)
    db = work / "cascade-db"
    run(steam, "createdb", inputs[1], inputs[0], db, "--threads", 2)
    cases = [(cutoff, 1, 0) for cutoff in ("0.001", "0.01", "0.1", "1")]
    cases += [("0.01", 1, 0.9), ("inf", 0, 0)]
    for cutoff, correction, identity in cases:
        output = work / f"cascade-{cutoff}-{identity}"
        extra = ["--comp-bias-corr", 0] if correction == 0 else []
        result = run(steam, "cluster", db, output, work / f"cascade-tmp-{cutoff}-{identity}",
                     "-c", 0, "-e", cutoff, "--min-seq-id", identity,
                     "--threads", 2, "-v", 3, *extra)
        log = result.stdout + result.stderr
        matcher = [line for line in log.splitlines() if line.startswith("kmermatcher ")]
        assert matcher and all(" -k 0 " in line for line in matcher), log
        alignments = [line for line in log.splitlines()
                      if line.startswith("align ") and "aln_step" in line]
        assert len(alignments) == (1 if identity >= 0.7 else 3), log
        assert all(f"--comp-bias-corr {correction} " in line for line in alignments), log
        prefilter = [line for line in log.splitlines()
                     if line.startswith("prefilter ") and "pref_step0" in line]
        assert len(prefilter) == 1 and "--comp-bias-corr 0 " in prefilter[0], log
        prefilters = [line for line in log.splitlines()
                      if line.startswith("prefilter ") and "pref_step" in line]
        caps = [int(re.search(r"--max-seqs (\d+)", line)[1]) for line in prefilters]
        assert caps == ([1000] if identity >= 0.7 else [100, 200, 1000]), log
        assert all("--diag-score 1 " in line for line in prefilters), log
        assert len(re.findall(r"Combined TEA\+AA ungapped diagonal scoring enabled \(AA weight 3(?:\.0+)?\)", log)) == len(prefilters), log
        tsv = work / f"cascade-{cutoff}-{identity}.tsv"
        run(steam, "createtsv", db, db, output, tsv, "--threads", 2)
        members = [line.split("\t")[1] for line in tsv.read_text().splitlines()]
        assert len(members) == len(names) and set(members) == names, members
    output = work / "easy-cascade"
    run(steam, "easy-cluster", inputs[1], inputs[0], output, work / "easy-cascade-tmp",
        "-c", 0, "--threads", 2)
    members = [line.split("\t")[1] for line in
               Path(str(output) + "_cluster.tsv").read_text().splitlines()]
    assert len(members) == len(names) and set(members) == names, members

    # Explicit limits remain upper bounds; the last stage always uses the
    # requested cap, including one- and two-stage cascades.
    for steps, limit, caps in ((1, 37, [37]), (2, 137, [100, 137]),
                               (3, 50, [50, 50, 50]), (4, 700, [100, 200, 200, 700])):
        result = run(steam, "cluster", db, work / f"caps-{steps}-{limit}",
                     work / f"caps-tmp-{steps}-{limit}", "-c", 0,
                     "--cluster-steps", steps, "--max-seqs", limit, "--threads", 2)
        prefilters = [line for line in (result.stdout + result.stderr).splitlines()
                      if line.startswith("prefilter ") and "pref_step" in line]
        assert [int(re.search(r"--max-seqs (\d+)", line)[1]) for line in prefilters] == caps

    # Custom scoring options must reach all child alignments and the combined
    # diagonal scorer, not fall back to the default matrices/AA weight.
    custom_matrix = work / "custom.out"
    source = Path(__file__).parents[1] / "data" / "matcha.out"
    custom_matrix.write_bytes(source.read_bytes())
    result = run(steam, "cluster", db, work / "custom", work / "custom-tmp",
                 "--matcha", custom_matrix, "--aa-weight", 4, "-e", "inf", "-c", 0,
                 "--threads", 2)
    log = result.stdout + result.stderr
    assert len(re.findall(r"Combined TEA\+AA ungapped diagonal scoring enabled \(AA weight 4(?:\.0+)?\)", log)) == 3, log
    alignments = [line for line in log.splitlines() if line.startswith("align ")]
    for line in alignments:
        arguments = shlex.split(line)
        assert arguments[arguments.index("--matcha") + 1] == str(custom_matrix), line
        assert float(arguments[arguments.index("--aa-weight") + 1]) == 4, line
    for mode in ("--single-step-clustering", "--ungapped-tea-aa"):
        value = 1 if mode == "--single-step-clustering" else 0
        result = run(steam, "cluster", db, work / f"mode-{mode}", work / f"tmp-{mode}",
                     mode, value, "-c", 0, "--threads", 2)
        log = result.stdout + result.stderr
        assert ("Combined TEA+AA ungapped diagonal scoring enabled" in log) == (value == 1), log
        if value == 1:
            assert "clusthash " not in log, log
    for method in ("linclust", "easy-linclust"):
        output = work / method
        arguments = (inputs[1], inputs[0]) if method.startswith("easy") else (db,)
        run(steam, method, *arguments, output, work / f"{method}-tmp", "--threads", 2)

    # Force a stale, invalid cascaded assignment, then exercise the actual
    # reassignment branch (including its merged paired target database).
    reassign_tmp = work / "reassign-tmp"
    run(steam, "cluster", db, work / "before-reassign", reassign_tmp,
        "-c", 0, "--threads", 2)
    cached = (reassign_tmp / "latest").resolve()
    invalid = work / "invalid-clusters.tsv"
    invalid.write_text("".join(f"0\t{i}\n" for i in range(16)))
    run(steam, "tsv2db", invalid, cached / "clu", "--output-dbtype", 6)
    reassigned = work / "reassigned"
    result = run(steam, "cluster", db, reassigned, reassign_tmp, "-c", 0,
                 "--cluster-reassign", 1, "--force-reuse", 1, "--threads", 2)
    assert "seq_seeds.merged.keys" in result.stdout + result.stderr
    assert (cached / "seq_seeds.merged_aa.dbtype").exists()
    tsv = work / "reassigned.tsv"
    run(steam, "createtsv", db, db, reassigned, tsv, "--threads", 2)
    members = [line.split("\t")[1] for line in tsv.read_text().splitlines()]
    assert len(members) == len(names) and set(members) == names, members

    # Identical TEA strings must not bypass an AA-identity requirement through
    # the upstream TEA-only clusthash shortcut.
    aa_input, tea_input = work / "distinct-aa.fasta", work / "identical-tea.fasta"
    rng = random.Random(123)
    aa_input.write_text("".join(f">{name}\n{''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=128))}\n"
                                for name in ("one", "two")))
    letters = ''.join(rng.choices('ACDEFGHIKLMNPQRSTVWY', k=128))
    tea_input.write_text(f">one\n{letters}\n>two\n{letters}\n")
    run(steam, "createdb", tea_input, aa_input, work / "distinct", "--threads", 2)
    run(steam, "cluster", work / "distinct", work / "distinct-clu", work / "distinct-tmp",
        "--single-step-clustering", 1, "--min-seq-id", 0.99, "--threads", 2)
    assert len(Path(str(work / "distinct-clu") + ".index").read_text().splitlines()) == 2

    for k in (5, 6):
        result = run(steam, "cluster", db, work / f"k{k}", work / f"k{k}-tmp",
                     "-k", k, "--threads", 2)
        assert f"STEAM k-mer selection: k={k}, pattern=" in result.stdout + result.stderr


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--steam", type=Path, required=True)
    args = parser.parse_args()
    steam = args.steam.resolve()
    fixture = Path(__file__).parent / "data" / "ungapped_tea_aa"
    aa = sequence(fixture / "query_aa.fasta")
    tea = sequence(fixture / "query_tea.fasta")
    assert len(aa) == len(tea)
    changed = list(tea)
    for pos in range(5, len(changed), 40):
        changed[pos] = "A" if changed[pos] != "A" else "C"
    tea_variant = "".join(changed)

    with tempfile.TemporaryDirectory(prefix="steam-tfident-") as directory:
        work = Path(directory)
        check_cascade_calibration(steam, work)
        aa_fasta = work / "aa.fasta"
        tea_fasta = work / "tea.fasta"
        aa_fasta.write_text(f">first\n{aa}\n>second\n{aa}\n")
        tea_fasta.write_text(f">first\n{tea}\n>second\n{tea_variant}\n")
        db = work / "db"
        run(steam, "createdb", tea_fasta, aa_fasta, db, "--threads", 2)

        for method in ("cluster",):
            for cutoff, expected_clusters in ((0, 1), (1, 2)):
                result = work / f"{method}-{cutoff}"
                extra = []
                run(steam, method, db, result, work / f"{method}-tmp-{cutoff}",
                    *extra, "--min-tfident", cutoff, "-e", "inf", "--threads", 2)
                tsv = work / f"{method}-{cutoff}.tsv"
                run(steam, "createtsv", db, db, result, tsv, "--threads", 2)
                rows = [line.split("\t") for line in tsv.read_text().splitlines()]
                representatives = {row[0] for row in rows}
                assert len(representatives) == expected_clusters, (method, cutoff, rows)

        for method in ("easy-cluster",):
            output = work / method
            run(steam, method, tea_fasta, aa_fasta, output, work / f"{method}-tmp",
                "--min-tfident", 1,
                "-e", "inf", "--threads", 2)
            rows = [line.split("\t") for line in Path(str(output) + "_cluster.tsv").read_text().splitlines()]
            assert len({row[0] for row in rows}) == 2, (method, rows)

        # Retired iteration controls must not silently become one-pass runs.
        rejected(steam, "search", db, db, work / "iterations", work / "iterations-tmp",
                 "--num-iterations", 2, message='Unrecognized parameter "--num-iterations"')
        # Upstream utilities remain available, but STEAM must reject their
        # profile databases rather than interpret binary PSSMs as sequences.
        for track, sequence_text in (("aa", aa), ("tea", tea)):
            stockholm = work / f"{track}.sto"
            stockholm.write_text(
                f"# STOCKHOLM 1.0\n#=GF ID paired\nfirst {sequence_text}\n"
                f"second {sequence_text}\n//\n"
                f"# STOCKHOLM 1.0\n#=GF ID paired2\nfirst {sequence_text}\n"
                f"second {sequence_text}\n//\n"
            )
            msa = work / f"{track}.msa"
            profile = work / ("query_aa" if track == "aa" else "query")
            run(steam, "convertmsa", stockholm, msa)
            run(steam, "msa2profile", msa, profile, "--match-mode", 0, "--threads", 2)

        rejected(steam, "search", work / "query", db, work / "profile-hits", work / "profile-tmp",
                 "-e", "inf", message="STEAM search requires paired sequence databases")
        pref = work / "pref"
        pref.write_bytes(b'\0')
        Path(str(pref) + '.index').write_text('0\t0\t1\n')
        Path(str(pref) + '.dbtype').write_bytes((7).to_bytes(4, 'little'))
        rejected(steam, "align", work / "query", db, pref, work / "profile-align",
                 "-e", "inf", message="STEAM alignment requires paired sequence databases")

        for k, pattern in ((5, '1101101'), (6, '1110111')):
            indexed = work / f'indexed-{k}'
            run(steam, 'createdb', tea_fasta, aa_fasta, indexed, '--threads', 2)
            before = work / f'raw-{k}'
            run(steam, 'search', db, indexed, before, work / f'raw-tmp-{k}',
                '-k', k, '-e', 'inf', '-a', 1, '--threads', 2)
            run(steam, 'createindex', indexed, work / f'index-tmp-{k}', '-k', k,
                '--seed-sub-mat', 'VTML80.out', '--spaced-kmer-pattern', pattern,
                '--mask', 0, '--comp-bias-corr', 1, '--comp-bias-corr-scale', 0,
                '--threads', 2)
            for selected_k in (k, 0):
                result = work / f'index-{k}-selected-{selected_k}'
                run(steam, 'search', db, Path(str(indexed) + '.idx'), result,
                    work / f'index-search-tmp-{k}-{selected_k}',
                    '-k', selected_k, '-e', 'inf', '-a', 1, '--threads', 2)
                assert records(before) == records(result), (k, selected_k)
        print("PASS: finite-E cascades retain alignment settings; "
              "cluster splits by tfident; easy-cluster passes it through; "
              "profile/iteration inputs rejected; k5/k6 explicit and automatic index search unchanged")


if __name__ == "__main__":
    main()
