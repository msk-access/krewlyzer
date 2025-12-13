import pysam
import os

def create_test_data():
    # 1. Create Reference FASTA
    with open("test.fasta", "w") as f:
        f.write(">chr1\n")
        f.write("A" * 1000 + "\n")
        f.write(">chr2\n")
        f.write("G" * 1000 + "\n")
    pysam.faidx("test.fasta")

    # 2. Create BAM
    header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1000, 'SN': 'chr1'},
                   {'LN': 1000, 'SN': 'chr2'}] }

    with pysam.AlignmentFile("test.bam", "wb", header=header) as outf:
        # PAIRED READS (Valid)
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence="A"*100
        a.flag = 99
        a.reference_id = 0
        a.reference_start = 100
        a.mapq = 60
        a.cigar = ((0,100),)
        a.next_reference_id = 0
        a.next_reference_start=200
        a.template_length=167
        a.is_proper_pair = True
        outf.write(a)

        b = pysam.AlignedSegment()
        b.query_name = "read1"
        b.query_sequence="T"*100
        b.flag = 147
        b.reference_id = 0
        b.reference_start = 200
        b.mapq = 60
        b.cigar = ((0,100),)
        b.next_reference_id = 0
        b.next_reference_start=100
        b.template_length=-167
        b.is_proper_pair = True
        outf.write(b)

    pysam.index("test.bam")
    print("Created test.bam, test.bam.bai, test.fasta, test.fasta.fai")

if __name__ == "__main__":
    create_test_data()
