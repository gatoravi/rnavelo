import sys
from simplesam import Reader, Writer

# Modified from https://github.com/velocyto-team/velocyto.py/issues/130
def main():
    if len(sys.argv) != 2:
        print("Usage: python3 add_cb_ub_tags.py  BAM")
        sys.exit(1)
    in_file = open(sys.argv[1], 'r')
    sample = sys.argv[1].split(".bam", 1)[0]
    print("Sample: " + sample)
    out_file = sample + "_withtags.sam"
    in_sam = Reader(in_file)
    x = next(in_sam)
    print(x.tags)
    barcode_tag = 'CB'
    umi_tag = 'UB'
    with Reader(open(sys.argv[1])) as in_bam:
        with Writer(open(out_file, 'w'), in_bam.header) as out_sam:
            for read in in_bam:
                #print(read.qname)
                #read[umi_tag] = read.qname.split(":")[2] # add the umi tag
                #read[barcode_tag] = read.qname.split(":")[1] # add the barcode tag
                read[umi_tag] = "dummy_umi" # add the umi tag
                read[barcode_tag] = sample # add the barcode tag
                out_sam.write(read)

main()
