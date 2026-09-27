"""Standalone native-helper tests; build tests/seed_evidence_bridge.cpp first."""
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


BRIDGE=Path(os.environ.get('STEAM_SEED_BRIDGE',Path(__file__).resolve().parents[1]/
    'artifacts/experiments/seed-linear-v1/seed_evidence_bridge'))


@unittest.skipUnless(BRIDGE.is_file(),'Compile the standalone seed-evidence bridge first')
class AlignmentSeedEvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.model=Path(self.temporary.name)/'model'
        self.text='STEAM_SEED_MARKOV_V1\nACDEFGHIKLMNPQRSTVWY\n100\n'+(' '.join(['0.05']*20)+'\n')*21
        self.model.write_text(self.text)
        # Invalid-input tests must not pass merely because the bridge or its
        # compiler runtime cannot load in the test environment.
        self.run_bridge('AAA AAA 0 0 MMM 100')

    def run_bridge(self,row,pattern='1101101',success=True):
        result=subprocess.run([str(BRIDGE),str(self.model),pattern],input=row+'\n',text=True,capture_output=True)
        if not success:
            self.assertNotEqual(result.returncode,0)
            self.assertTrue(result.stderr)
            return
        self.assertEqual(result.returncode,0,result.stderr)
        evidence,score,denominator=result.stdout.split()
        return float(evidence),int(score),int(denominator)

    def test_repeat_deduplication_and_case(self):
        e,s,n=self.run_bridge('aaaaaaaaa AAAAAAAAA 0 0 MMMMMMMMM 100')
        self.assertAlmostEqual(e,5*math.log(20)-math.log(100))
        self.assertEqual(n,1)
        self.assertEqual(s,math.floor(100*(1+e)+.5))

    def test_zero_support_keeps_score(self):
        for row in ['AAA AAA 0 0 MMM 100','AAAAAAA CCCCCCC 0 0 MMMMMMM 100',
                    'AAAAAAA AAAAAA 0 0 MMMIMMM 100','- - 0 0 - 100']:
            self.assertEqual(self.run_bridge(row)[:2],(0.,100))

    def test_unused_positions_and_k6(self):
        self.assertGreater(self.run_bridge('AAXAAXA AAAAAAA 0 0 MMMMMMM 100')[0],0)
        self.assertEqual(self.run_bridge('XAAAAAA AAAAAAA 0 0 MMMMMMM 100')[0],0)
        e,_,_=self.run_bridge('AAAXAAA AAAAAAA 0 0 MMMMMMM 100','1110111')
        self.assertAlmostEqual(e,6*math.log(20)-math.log(100))

    def test_invalid_path_and_overflow_fail(self):
        for row in ['AAAAAAA AAAAAAA 1 0 MMMMMMM 100','AAAAAAA AAAAAAA 0 0 MZ 100',
                    'AAAAAAA AAAAAAA 0 0 MMMMMMM 2147483647','AAA AAA 0 0 MMM -1']:
            self.run_bridge(row,success=False)

    def test_bad_models_and_patterns_fail(self):
        row='AAAAAAA AAAAAAA 0 0 MMMMMMM 100'
        self.run_bridge(row,'11111',success=False)
        for text in [self.text.replace('0.05','0',1),self.text.replace('ACDEFGHIKLMNPQRSTVWY','ACDEFGHIKLMNPQRSTVWX'),
                     self.text+'extra',self.text.rsplit('\n',2)[0],self.text.replace('\n100\n','\n-1\n')]:
            self.model.write_text(text)
            self.run_bridge(row,success=False)


if __name__=='__main__':unittest.main()
