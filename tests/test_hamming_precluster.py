"""Hamming pre-clustering filters identity/coverage, not alignment E-values."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def records(prefix):
    parts = [prefix] if prefix.exists() else sorted(
        (p for p in prefix.parent.glob(prefix.name + ".*") if p.suffix[1:].isdigit()),
        key=lambda p: int(p.suffix[1:]))
    data = b"".join(p.read_bytes() for p in parts)
    return {int(key): data[int(offset):int(offset) + int(size)]
            for key, offset, size in (line.split() for line in
                Path(str(prefix) + ".index").read_text().splitlines())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steam", type=Path, required=True)
    args = parser.parse_args()
    steam = args.steam.resolve()
    fixture = Path(__file__).parent / "data/ungapped_tea_aa"
    environment = {k: v for k, v in os.environ.items()
                   if not k.startswith(("STEAM_", "MMSEQS_"))}

    def run(*arguments):
        result = subprocess.run([str(steam), *map(str, arguments)], env=environment,
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stdout + result.stderr

    with tempfile.TemporaryDirectory(prefix="hamming-precluster-",
                                     dir=os.environ.get("TMPDIR")) as temporary:
        work = Path(temporary)
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        for kind in ("aa", "tea"):
            sequence = "".join((fixture / f"query_{kind}.fasta").read_text().splitlines()[1:])
            different = ("".join(alphabet[(alphabet.index(c) + 1) % len(alphabet)]
                                 for c in sequence) if kind == "aa" else sequence)
            (work / f"{kind}.fasta").write_text(
                f">original\n{sequence}\n>copy\n{sequence}\n"
                f">different_aa\n{different}\n>short\n{sequence[:30]}\n")
        db = work / "db"
        run("createdb", work / "tea.fasta", work / "aa.fasta", db, "--threads", 2)
        keys = {name: int(key) for key, name, _ in
                (line.split() for line in Path(str(db) + ".lookup").read_text().splitlines())}
        (work / "pref.tsv").write_text("".join(
            f"{q}\t{t}\t100\t0\n" for q in sorted(keys.values())
            for t in sorted(keys.values())))
        run("tsv2db", work / "pref.tsv", work / "pref", "--output-dbtype", 7)
        outputs = []
        for label, cutoff in (("finite", "0.001"), ("infinite", "inf")):
            prefix = work / label
            run("rescorediagonal", db, db, work / "pref", prefix,
                "--rescore-mode", 0, "--filter-hits", 0, "--min-seq-id", 0.9,
                "-c", 0.8, "-e", cutoff, "--threads", 2)
            outputs.append(records(prefix))
        assert outputs[0] == outputs[1], "Hamming membership changed with E-value cutoff"
        accepted = {key: {int(line.split(b"\t")[0])
                          for line in value.rstrip(b"\x00\n").splitlines()}
                    for key, value in outputs[0].items()}
        expected = {keys["original"], keys["copy"]}
        for name in ("original", "copy"):
            assert accepted[keys[name]] == expected, accepted
        assert accepted[keys["different_aa"]] == {keys["different_aa"]}, accepted
        assert accepted[keys["short"]] == {keys["short"]}, accepted
        run("linclust", db, work / "clu", work / "tmp", "--min-seq-id", 0.9,
            "-c", 0.8, "-e", "0.001", "--threads", 2)
        preliminary = records((work / "tmp/latest").resolve() / "pre_clust")
        groups = [set(map(int, value.rstrip(b"\x00\n").split()))
                  for value in preliminary.values()]
        assert expected in groups, "Identical pair did not collapse in Hamming pre-clustering"
        print("PASS: finite-E Hamming pre-clustering merges identical pairs; "
              "AA-identity and coverage filters remain enforced")


if __name__ == "__main__":
    main()
