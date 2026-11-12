"""Classify a file of positive and negative Y-DNA SNP calls."""

import argparse
import contextlib
import json
from pathlib import Path
import re
import sys

def parse_calls(text):
    calls = [call.strip() for call in re.split(r'[,\r\n]', text) if call.strip()]
    if not calls:
        raise ValueError('No SNP calls found in the input.')
    for call in calls:
        if len(call) < 2 or call[-1] not in '+-' or any(c.isspace() for c in call):
            raise ValueError('Invalid SNP call {!r}; expected a name followed by + or -.'.format(call))
    return calls

def predict(text, clade_snps, snp_clades, all_predictions=False):
    calls = parse_calls(text)
    for path in (clade_snps, snp_clades):
        for required in (Path(path), Path(str(path) + '.tbi')):
            if not required.is_file():
                raise ValueError('Missing database file: {}'.format(required))

    import tabix
    from Common.CommonMethods import getEncodedPositivesNegatives, getJSONObject
    positives, negatives = getEncodedPositivesNegatives(calls)
    # The original classifier can print diagnostics during database lookups.
    try:
        with contextlib.redirect_stdout(sys.stderr):
            result = getJSONObject('all,score', positives, negatives, str(clade_snps), str(snp_clades), None)
    except tabix.TabixError as error:
        raise ValueError('Cannot read the Tabix database: {}'.format(error)) from error
    if isinstance(result, dict) and 'error' in result:
        raise ValueError(result['error'])
    if not result:
        raise ValueError('No matching clade found for the supplied SNP calls.')
    result.sort(key=lambda candidate: (-candidate['score'], candidate['clade']))
    if all_predictions:
        return result
    best = result[0]
    if len(result) > 1:
        best['nextPrediction'] = {'clade': result[1]['clade'], 'score': result[1]['score']}
        if result[1]['score'] == best['score']:
            warning = 'Multiple clades share the highest score; use --all to inspect them.'
            best['warning'] = (best.get('warning', '') + ' ' + warning).strip()
    return best

def main(argv=None):
    parser = argparse.ArgumentParser(description='Predict a Y-DNA clade from signed SNP calls.')
    parser.add_argument('input', help='Text file containing SNP calls, or - to read stdin.')
    parser.add_argument('--database', type=Path, default=Path('tabix_output'),
                        help='Directory containing cladeSNP.bgz and SNPclade.bgz (default: tabix_output).')
    parser.add_argument('--clade-snps', type=Path, help='Override the clade-to-SNP database path.')
    parser.add_argument('--snp-clades', type=Path, help='Override the SNP-to-clade database path.')
    parser.add_argument('--all', action='store_true', help='Return all candidate clades as a JSON array.')
    args = parser.parse_args(argv)
    try:
        if args.input == '-':
            text = sys.stdin.read()
        else:
            text = Path(args.input).read_text(encoding='utf-8-sig')
        result = predict(text, args.clade_snps or args.database / 'cladeSNP.bgz',
                         args.snp_clades or args.database / 'SNPclade.bgz', args.all)
    except (OSError, ValueError, UnicodeError) as error:
        parser.exit(2, '{}: error: {}\n'.format(parser.prog, error))
    print(json.dumps(result))
    return 0

if __name__ == '__main__':
    sys.exit(main())
