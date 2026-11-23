"""Standalone prefilter must match search's native candidate rows."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def rows(prefix):
    files = [prefix] if prefix.exists() else sorted(
        (p for p in prefix.parent.glob(prefix.name+'.*') if p.suffix[1:].isdigit()),
        key=lambda p:int(p.suffix[1:]))
    data = b''.join(path.read_bytes() for path in files)
    return {int(key):data[int(offset):int(offset)+int(size)]
            for key,offset,size in (line.split() for line in
                Path(str(prefix)+'.index').read_text().splitlines())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--steam',type=Path,required=True)
    args = parser.parse_args()
    steam = args.steam.resolve()
    root = Path(__file__).resolve().parents[1]
    fixture = root/'tests/data/ungapped_tea_aa'
    environment = {k:v for k,v in os.environ.items() if not k.startswith('STEAM_')}

    def run(*arguments,success=True):
        result = subprocess.run([str(steam),*map(str,arguments)],env=environment,
            capture_output=True,text=True)
        assert (result.returncode == 0) == success, result.stdout+result.stderr
        return result.stdout+result.stderr

    with tempfile.TemporaryDirectory(prefix='steam-prefilter-',dir=os.environ.get('TMPDIR')) as temporary:
        work = Path(temporary)
        query,target = work/'query',work/'target'
        for prefix,role in [(query,'query'),(target,'target')]:
            run('createdb',fixture/(role+'_tea.fasta'),fixture/(role+'_aa.fasta'),prefix,
                '--threads',2,'-v',1)
        cases = [
            ('defaults',[],5,'combined_best'),
            ('tea-only',['--ungapped-tea-aa',0],5,'tea_best'),
            ('aa-weight',['--aa-weight',1.4],5,'combined_best'),
            ('aa-matrix',['--sub-mat','PAM30.out'],5,None),
            ('split',['--split',2,'--split-mode',0],5,'combined_best'),
            ('k6',['-k',6],6,None),
            ('pattern-selects-k',['--spaced-kmer-pattern','1110111'],6,None),
            ('auto-k',['-k',0],5,'combined_best'),
            ('custom-matrix',['--matcha',root/'tests/data/legacy_matcha.out'],5,None)]
        lookup = {int(line.split()[0]):line.split()[1]
                  for line in Path(str(target)+'.lookup').read_text().splitlines()}
        for name,options,k,winner in cases:
            pref = work/(name+'-pref')
            common = [*options,'--max-seqs',1,'--threads',2,'-v',3]
            direct_log = run('prefilter',query,target,pref,*common)
            search_tmp = work/(name+'-search-tmp')
            search_log = run('search',query,target,work/(name+'-align'),search_tmp,
                *common,'-e','inf','--seed-correction',0)
            assert f'STEAM k-mer selection: k={k}, pattern=' in direct_log,direct_log
            assert f'STEAM k-mer selection: k={k}, pattern=' in search_log,search_log
            assert rows(pref) == rows(search_tmp/'pref'),name
            assert rows(pref),name
            combined = name != 'tea-only'
            message = 'Combined TEA+AA ungapped diagonal scoring enabled'
            assert (message in direct_log) == combined,name
            if winner is not None:
                payload = next(iter(rows(pref).values()))
                assert lookup[int(payload.split(b'\t',1)[0])] == winner,(name,payload)
            print('PASS',name,flush=True)
        run('createindex',target,work/'index-tmp','-k',6,
            '--seed-sub-mat','VTML80.out','--spaced-kmer-pattern','1110111',
            '--mask',0,'--comp-bias-corr',1,'--comp-bias-corr-scale',0,'--threads',2)
        indexed = Path(str(target)+'.idx')
        for options in ([],['-k',0]):
            pref = work/('indexed-auto' if options else 'indexed-default')
            log = run('prefilter',query,indexed,pref,*options,
                '--max-seqs',1,'--threads',2,'-v',3)
            assert 'STEAM k-mer selection: k=6, pattern=1110111' in log,log
            assert rows(pref) == rows(work/'k6-pref')
        print('PASS cached-index automatic k/pattern and native rows',flush=True)
        rejected = run('prefilter',query,target,work/'bad-scale','--tea-scale',1,
            '--threads',2,success=False)
        assert 'requires --tea-scale 2' in rejected,rejected
        Path(str(query)+'_aa.dbtype').unlink()
        rejected = run('prefilter',query,target,work/'missing-aa','--threads',2,success=False)
        assert 'requires paired TEA/AA' in rejected,rejected
        run('prefilter',query,target,work/'unpaired-tea-only','--ungapped-tea-aa',0,'--threads',2)
        print('PASS scale and companion checks; TEA-only compatibility',flush=True)


if __name__ == '__main__':
    main()
