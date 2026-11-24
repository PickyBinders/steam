"""Check indexed search, AA output, and exact double-hit filtering."""

import argparse
from pathlib import Path
import re
import subprocess
import tempfile

from test_prefilter_setup import rows


def check_indexed_search(run, work):
    fixture = Path(__file__).parent / "data/ungapped_tea_aa"
    query, target = work / "query", work / "target"
    query_fastas = [fixture / f"query_{alphabet}.fasta" for alphabet in ("tea", "aa")]
    run("createdb", *query_fastas, query)
    target_fastas = []
    repeated = "".join(f">repeat-{i}\n{'A' * 60}\n" for i in range(32))
    for alphabet in ("tea", "aa"):
        fasta = work / f"target_{alphabet}.fa"
        fasta.write_text((fixture / f"target_{alphabet}.fasta").read_text().rstrip() + "\n" + repeated)
        target_fastas.append(fasta)
    run("createdb", *target_fastas, target)
    seeds = ["-k", 5, "--spaced-kmer-pattern", "1101101", "--seed-sub-mat", "VTML80.out",
             "--mask", 0, "--comp-bias-corr", 1, "--comp-bias-corr-scale", 0]
    columns = ["--format-output", "query,target,fident,qseq,tseq,qaln,taln,ppos,qteaaln,tfident"]
    aln = work / "alignments"
    run("search", query, target, aln, work / "search", *seeds, "-e", "inf", "-a")
    expected = work / "expected.tsv"
    run("convertalis", query, target, aln, expected, *columns)
    assert expected.stat().st_size > 0
    native_pref = work / "native-prefilter"
    native_log = run("prefilter", query, target, native_pref, *seeds, "-v", 3)

    run("createindex", query, work / "query-index", *seeds)
    query_idx = Path(f"{query}.idx")
    target_idx = Path(f"{target}.idx")
    run("createindex", target, work / "target-index", *seeds, "--split", 2)
    assert Path(f"{target_idx}.0").is_file()
    indexed_pref = work / "indexed-prefilter"
    log = run("prefilter", query, target_idx, indexed_pref, *seeds, "-v", 3)
    assert "across target index splits" in log, log
    threshold = r"\[steam-bucket-filter\] global coverage=[^\n%]+% -> bucket threshold (\d+)"
    native_count, indexed_count = (re.search(threshold, output) for output in (native_log, log))
    assert native_count and indexed_count, (native_log, log)
    assert int(native_count[1]) >= 10 and indexed_count[1] == native_count[1], (native_log, log)
    assert rows(native_pref) == rows(indexed_pref)

    indexed = work / "indexed.tsv"
    run("convertalis", query_idx, target_idx, aln, indexed, *columns)
    assert indexed.read_bytes() == expected.read_bytes()

    easy = work / "easy.tsv"
    run("easy-search", *query_fastas, target_idx, easy, work / "easy", *seeds, "-e", "inf", *columns)
    assert easy.read_bytes() == expected.read_bytes()

    tea_only = work / "easy-tea.tsv"
    run("easy-search", *query_fastas, target_idx, tea_only, work / "easy-tea", *seeds,
        "-e", "inf", "--format-output", "qteaaln,tfident")
    assert tea_only.read_text() == "".join(
        "\t".join(line.split("\t")[-2:]) + "\n" for line in expected.read_text().splitlines())

    Path(f"{target}_aa.dbtype").rename(work / "hidden-aa.dbtype")
    run("convertalis", query, target_idx, aln, work / "missing-aa.tsv", *columns,
        error="requires the amino acid companion database")


def check_double_hits(run, work):
    word = "ACDEFGHI"
    query_fasta, target_fasta = work / "double-query.fa", work / "double-target.fa"
    query_fasta.write_text(f">query\n{'X' * 270}{word}\n")
    targets = {f"{kind}-{diagonal}": f"{'X' * (270 - diagonal)}{motif}"
               for diagonal in (0, 10, 255, 256)
               for kind, motif in (("single", word[:-1]), ("double", word))}
    target_fasta.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in targets.items()))
    query, target = work / "double-query", work / "double-target"
    run("createdb", query_fasta, query_fasta, query)
    run("createdb", target_fasta, target_fasta, target)
    lookup = dict(line.split()[:2] for line in Path(f"{target}.lookup").read_text().splitlines())
    for diag_score in (0, 1):
        pref = work / f"double-prefilter-{diag_score}"
        log = run("prefilter", query, target, pref, "-k", 6, "--spaced-kmer-pattern", "1110111",
            "--mask", 0, "--exact-kmer-matching", 1, "--diag-score", diag_score,
            "--ungapped-tea-aa", diag_score, "--min-ungapped-score", 1)
        result, = rows(pref).values()
        names = {lookup[line.split(b"\t", 1)[0].decode()] for line in result.rstrip(b"\0").splitlines()}
        assert names == {name for name in targets if name.startswith("double-")}, (names, log)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steam", type=Path, required=True)
    steam = parser.parse_args().steam.resolve()

    def run(*args, error=None):
        result = subprocess.run([str(steam), *map(str, args), "--threads", "2"],
                                capture_output=True, text=True)
        log = result.stdout + result.stderr
        if error is None:
            assert result.returncode == 0, log
        else:
            assert result.returncode > 0 and error in log, log
        return log

    with tempfile.TemporaryDirectory(prefix="steam-indexed-search-") as directory:
        work = Path(directory)
        check_indexed_search(run, work)
        print("PASS indexed split search, AA output, companion error and automatic backtrace")
        check_double_hits(run, work)
        print("PASS diagonal zero and modulo-256 double-hit filtering")


if __name__ == "__main__":
    main()
