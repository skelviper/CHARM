import argparse
import gzip
import sys
import os

def process_line(line):
    parts = line.strip().split('\t')
    read_id = parts[0]
    segments = parts[1:]
    
    # Find all R1 and R2 segments
    r1_segments = []
    r2_segments = []
    for segment in segments:
        if segment.endswith('R1'):
            r1_segments.append(segment)
        elif segment.endswith('R2'):
            r2_segments.append(segment)

    r1_haplotype = "."
    for r1_segment in r1_segments:
        fields = r1_segment.split('!')
        if fields[4] != '.':
            r1_haplotype = fields[4]
            break

    output_lines = []
    for r2_segment in r2_segments:
        r2_fields = r2_segment.split('!')
        if r2_fields[4] == '.':
            r2_fields[4] = r1_haplotype
        bed_format = f"{r2_fields[0]}\t{r2_fields[1]}\t{r2_fields[2]}\t{r2_fields[4]}\t{r2_fields[5]}\t{r2_fields[3]}"
        output_lines.append(bed_format)
    
    return output_lines

def convert_seg_to_bed(input_file, output_file):
    if input_file.endswith('.gz'):
        infile = gzip.open(input_file, 'rt')
    else:
        infile = open(input_file, 'r') if input_file != '-' else sys.stdin

    outfile = open(output_file, 'w') if output_file != '-' else sys.stdout

    with infile, outfile:
        for line in infile:
            if line.strip(): 
                bed_lines = process_line(line)
                for bed_line in bed_lines:
                    outfile.write(bed_line + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert seg file to bed file')
    parser.add_argument('input_file', type=str, help='Input seg file, use "-" for stdin')
    parser.add_argument('output_file', type=str, help='Output bed file, use "-" for stdout')
    args = parser.parse_args()
    convert_seg_to_bed(args.input_file, args.output_file)
