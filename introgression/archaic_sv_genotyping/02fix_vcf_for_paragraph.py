#!/usr/bin/env python3
"""
Fix SURVIVOR-merged VCF for Paragraph compatibility:
1. INS: convert indel-style to symbolic <INS> + SEQ in INFO
2. DEL: ensure INFO has END field (END = POS + REF_len - 1)
3. Remove FORMAT and all sample GT columns

Usage:
  python3 fix_vcf_for_paragraph.py 5caller_paragraph_chm13_fixed.ins_del.vcf.gz 5caller_paragraph_ready.vcf.gz
  python3 fix_vcf_for_paragraph.py input.vcf.gz output.vcf.gz
"""
import sys, gzip

IN_VCF = sys.argv[1]
OUT_VCF = sys.argv[2]

MIN_FLANK = 150  # Paragraph needs at least read_length bp upstream of variant

total = 0; ins_count = 0; del_count = 0; skipped_flank = 0

with gzip.open(IN_VCF, 'rt') as fin:
    import subprocess
    bgzip = subprocess.Popen(['bgzip', '-c'], stdin=subprocess.PIPE, stdout=open(OUT_VCF, 'wb'))

    header_done = False
    for line in fin:
        if line.startswith('##'):
            bgzip.stdin.write(line.encode())
            # Add SEQ INFO definition before #CHROM
            if not header_done and line.startswith('##INFO') and 'ID=END' in line:
                bgzip.stdin.write(b'##INFO=<ID=SEQ,Number=1,Type=String,Description="Inserted sequence for <INS>">\n')
                header_done = True
            continue

        if line.startswith('#CHROM'):
            # Ensure SEQ header is written if not already (before any END INFO line)
            if not header_done:
                bgzip.stdin.write(b'##INFO=<ID=SEQ,Number=1,Type=String,Description="Inserted sequence for <INS>">\n')
            parts = line.strip().split('\t')
            bgzip.stdin.write(('\t'.join(parts[:8]) + '\n').encode())
            continue

        total += 1
        parts = line.strip().split('\t')
        if len(parts) < 8:
            continue  # skip empty/malformed lines

        chrom, pos, vid, ref, alt, qual, filt, info = parts[0:8]
        pos_int = int(pos)

        # Skip variants too close to chromosome start
        if pos_int <= MIN_FLANK:
            skipped_flank += 1
            continue

        # Parse INFO
        info_items = [x for x in info.split(';') if not x.startswith('END=')]
        svtype = ''
        for item in info_items:
            if item.startswith('SVTYPE='):
                svtype = item.split('=')[1]
                break

        pos_int = int(pos)

        if svtype == 'INS':
            # Convert to <INS> + SEQ
            ins_seq = alt[len(ref):] if alt.startswith(ref) else alt[1:]
            alt = '<INS>'
            info_items.append(f'SEQ={ins_seq}')
            info_items.append(f'END={pos_int}')
            ins_count += 1

        elif svtype == 'DEL':
            end = pos_int + len(ref) - 1
            info_items.append(f'END={end}')
            del_count += 1

        else:
            # Other types: ensure END if not present
            has_end = any(x.startswith('END=') for x in info_items)
            if not has_end:
                end = pos_int + len(ref) - 1
                info_items.append(f'END={end}')

        info = ';'.join(info_items)
        bgzip.stdin.write(f'{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\n'.encode())

    bgzip.stdin.close()
    bgzip.wait()

# Index
import subprocess
subprocess.run(['tabix', '-p', 'vcf', OUT_VCF])

print(f'Done. Total: {total}, INS: {ins_count}, DEL: {del_count}, Skipped (POS<={MIN_FLANK}): {skipped_flank}', file=sys.stderr)
