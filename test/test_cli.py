"""Exercise the CLI with a small synthetic tree and real Tabix indexes."""

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import pysam

ROOT = Path(__file__).resolve().parents[1]
TREE = ROOT / 'test' / 'fixtures' / 'tree.json'

class CommandLineTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.work = Path(self.temporary.name)
        self.database = self.work / 'database with spaces; literal'
        result = self.run_script('createTreeInTabix.py', TREE, '--output-dir', self.database, '--index')
        self.assertEqual(result.returncode, 0, result.stderr)

    def run_script(self, script, *args, input=None, env=None):
        return subprocess.run([sys.executable, str(ROOT / script), *map(str, args)],
                              input=input, text=True, capture_output=True, cwd=self.work, env=env)

    def predict(self, text, *args, env=None):
        return self.run_script('clade_finder.py', '-', '--database', self.database,
                               *args, input=text, env=env)

    def test_file_prediction_is_json(self):
        result = self.run_script('clade_finder.py', ROOT / 'test' / 'sample_for_clade_finder.txt',
                                 '--database', self.database)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout), {
            'clade': 'J-YP879', 'score': 2.0,
            'nextPrediction': {'clade': 'J-Z6065', 'score': 1.0}})

    def test_stdin_aliases_and_all_predictions(self):
        result = self.predict(' Z6065+\nALIAS879+, L21-\n', '--all')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([(r['clade'], r['score']) for r in json.loads(result.stdout)],
                         [('J-YP879', 2.0), ('J-Z6065', 1.0)])

    def test_ties_are_reported_and_repeatable(self):
        outputs = []
        for seed in ('1', '2', '3'):
            result = self.predict('L21+, M81+', env=dict(os.environ, PYTHONHASHSEED=seed))
            self.assertEqual(result.returncode, 0, result.stderr)
            outputs.append(json.loads(result.stdout))
        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(outputs[0], outputs[2])
        self.assertIn('share the highest score', outputs[0]['warning'])

    def test_bad_or_unusable_calls_fail_without_json(self):
        for calls in ('', 'Z6065', '+', 'Z6065+ L21-', 'UNKNOWN+', 'Z6065-', 'YP879+, ALIAS879-'):
            with self.subTest(calls=calls):
                result = self.predict(calls)
                self.assertEqual(result.returncode, 2, result.stderr)
                self.assertEqual(result.stdout, '')
                self.assertNotIn('Traceback', result.stderr)

    def test_unknown_calls_with_all_fail(self):
        result = self.predict('UNKNOWN+', '--all')
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, '')

    def test_missing_input_and_index(self):
        result = self.run_script('clade_finder.py', self.work / 'missing.txt')
        self.assertEqual(result.returncode, 2)
        (self.database / 'SNPclade.bgz.tbi').unlink()
        result = self.predict('Z6065+')
        self.assertEqual(result.returncode, 2)
        self.assertIn('Missing database file', result.stderr)
        self.assertEqual(result.stdout, '')

    def test_custom_paths_and_old_flag_names(self):
        clades = self.work / 'custom clades.tsv'
        snps = self.work / 'custom snps.tsv'
        result = self.run_script('createTreeInTabix.py', TREE, '--cladeSNPFilePath', clades,
                                 '--SNPcladeFilePath', snps, '--create_tabix_files')
        self.assertEqual(result.returncode, 0, result.stderr)
        result = self.run_script('clade_finder.py', '-', '--clade-snps', clades.with_suffix('.bgz'),
                                 '--snp-clades', snps.with_suffix('.bgz'), input='Z6065+')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)['clade'], 'J-Z6065')

    def test_marker_outputs_products_and_ignored_aliases(self):
        markers = self.work / 'markers.tsv'
        markers.write_text('100\tYP879\tG\n')
        products = self.work / 'products.tsv'
        products.write_text('ALIAS879\t123\n')
        ignore = self.work / 'ignore.txt'
        ignore.write_text('ALIAS879\n')
        result = self.run_script('createTreeInTabix.py', TREE, '--output-dir', self.database,
                                 '--hg19-markers', markers, '--hg38-markers', markers,
                                 '--products', products, '--ignore', ignore, '--index')
        self.assertEqual(result.returncode, 0, result.stderr)
        for build in ('hg19', 'hg38'):
            with pysam.TabixFile(str(self.database / (build + 'positionMarkers.bgz'))) as database:
                self.assertEqual(list(database.fetch('100', 0, 1)), ['100\t1\t1\tYP879\tG'])
        with pysam.TabixFile(str(self.database / 'SNPclade.bgz')) as database:
            self.assertNotIn('ALIAS879', database.contigs)
            self.assertEqual(list(database.fetch('YP879/ALIAS879', 2, 3)),
                             ['YP879/ALIAS879\t3\t3\t123\t.'])

    def test_legacy_positional_builder(self):
        markers = self.work / 'markers.tsv'
        markers.write_text('100\tZ6065\tG\n')
        products = self.work / 'products.tsv'
        products.write_text('Z6065\t123\n')
        ignore = self.work / 'ignore.txt'
        ignore.write_text('')
        outputs = [self.work / name for name in ('clades', 'snps', 'hg19', 'hg38')]
        result = self.run_script('createTreeInTabix.py', TREE, markers, markers,
                                 *outputs, products, ignore)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(all(path.is_file() for path in outputs))
        self.assertIn('Z6065\t3\t3\t123', outputs[1].read_text())

    def test_output_conflicts_do_not_overwrite_inputs(self):
        tree = self.work / 'tree.json'
        original = TREE.read_text()
        tree.write_text(original)
        result = self.run_script('createTreeInTabix.py', tree, '--clade-snps', tree)
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(tree.read_text(), original)
        result = self.run_script('createTreeInTabix.py', tree, '--clade-snps', self.work / 'same',
                                 '--snp-clades', self.work / 'same')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.work / 'same').exists())

    def test_marker_output_requires_input(self):
        result = self.run_script('createTreeInTabix.py', TREE, '--hg19-output', self.work / 'unused')
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('requires --hg19-markers', result.stderr)

    def test_help_without_running_classifier(self):
        result = self.run_script('clade_finder.py', '--help')
        self.assertEqual(result.returncode, 0)
        self.assertIn('--database', result.stdout)
        self.assertNotIn('csv', result.stdout)

if __name__ == '__main__':
    unittest.main()
