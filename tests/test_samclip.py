import pytest
from snippy_ng.utils.samclip import samclip_filter_lines, SamclipError

CONTIGS = {"ref": 1000}

def make_sam(qname, flag, ref="ref", pos=500, cigar="100M", seq="A"*100):
    return f"{qname}\t{flag}\t{ref}\t{pos}\t60\t{cigar}\t=\t{pos}\t0\t{seq}\t{'I'*len(seq)}\n"

def test_samclip_filter_basic():
    header = ["@HD\tVN:1.6\tSO:queryname\n"]
    lines = [
        make_sam("r1", 0, cigar="100M"),      # Good
        make_sam("r2", 0, cigar="10S90M"),    # Bad (10S > 5)
        make_sam("r3", 0, cigar="5S95M"),     # Good (5S <= 5)
    ]
    
    out = list(samclip_filter_lines(header + lines, CONTIGS, max_clip=5))
    
    # Should keep header + r1 + r3
    assert len(out) == 3
    assert out[0] == header[0]
    assert "r1" in out[1]
    assert "r3" in out[2]

def test_samclip_invert():
    header = ["@HD\tVN:1.6\tSO:queryname\n"]
    lines = [
        make_sam("r1", 0, cigar="100M"),      # Good
        make_sam("r2", 0, cigar="10S90M"),    # Bad
    ]
    
    out = list(samclip_filter_lines(header + lines, CONTIGS, max_clip=5, invert=True))
    
    # Should keep header + r2 (the bad one)
    assert len(out) == 2
    assert "r2" in out[1]

def test_samclip_fix_mate_logic():
    header = ["@HD\tVN:1.6\tSO:queryname\n"]
    # Pair: r1/1 (bad), r1/2 (good)
    # r1/2 should get MUNMAP (0x8) added to its flag (0x1 | 0x2 = 3 -> 3 | 8 = 11)
    # Note: PAIRED=1. Initial flag 3 = PAIRED(1) + PROPER_PAIR(2).
    lines = [
        make_sam("r1", 99, cigar="10S90M"),   # Bad (read1)
        make_sam("r1", 147, cigar="100M"),    # Good (read2)
    ]
    
    out = list(samclip_filter_lines(header + lines, CONTIGS, max_clip=5, fix_mate=True))
    
    # Should keep header + r1 (second line only)
    assert len(out) == 2
    
    line2 = out[1]
    fields = line2.split("\t")
    qname = fields[0]
    flag = int(fields[1])
    
    assert qname == "r1"
    # Check MUNMAP bit (0x8) is set
    assert (flag & 0x8) != 0, f"Flag {flag} should have MUNMAP (8) set"

def test_samclip_no_fix_mate():
    header = ["@HD\tVN:1.6\tSO:queryname\n"]
    lines = [
        make_sam("r1", 99, cigar="10S90M"),   # Bad
        make_sam("r1", 147, cigar="100M"),    # Good
    ]
    
    # fix_mate=False
    out = list(samclip_filter_lines(header + lines, CONTIGS, max_clip=5, fix_mate=False))
    
    # Should keep header + r1 (second line)
    # But NO flag modification
    assert len(out) == 2
    
    line2 = out[1]
    flag = int(line2.split("\t")[1])
    assert (flag & 0x8) == 0, "Should not set MUNMAP when fix_mate=False"

def test_samclip_allows_unsorted_if_no_fix_mate():
    # SO:coordinate allowed if fix_mate=False
    header = ["@HD\tVN:1.6\tSO:coordinate\n"]
    lines = [
        make_sam("r1", 0),
    ]
    
    out = list(samclip_filter_lines(header + lines, CONTIGS, fix_mate=False))
    assert len(out) == 2

def test_samclip_enforces_sorted_header_if_fix_mate():
    header = ["@HD\tVN:1.6\tSO:coordinate\n"]
    lines = []
    
    with pytest.raises(SamclipError, match="Input must be queryname-sorted"):
        list(samclip_filter_lines(header + lines, CONTIGS, fix_mate=True))

def test_samclip_handles_out_of_order_stream_gracefully():
    # Since we removed the check, it should just process them.
    # If fix_mate=True, it groups by ID. If they are split, they are just treated as separate groups.
    header = ["@HD\tVN:1.6\tSO:queryname\n"]
    lines = [
        make_sam("r2", 0),
        make_sam("r1", 0),
    ]
    
    out = list(samclip_filter_lines(header + lines, CONTIGS, fix_mate=True))
    
    # Should get both back
    qnames = [l.split("\t")[0] for l in out[1:]]
    assert qnames == ["r2", "r1"]

def test_samclip_start_end_clipping():
    header = ["@HD\tVN:1.6\tSO:queryname\n"]
    # 5S at start, 6S at end -> total clip = 11? 
    # Logic in _bad: (L > max) or (R > max).
    # L = hl + sl. R = hr + sr.
    # 5S90M6S: L=5, R=6. Max=5. R > Max -> bad.
    
    lines = [
        make_sam("r1", 0, cigar="5S90M5S"), # L=5, R=5. Good.
        make_sam("r2", 0, cigar="6S90M"),   # L=6. Bad.
        make_sam("r3", 0, cigar="90M6S"),   # R=6. Bad.
    ]
    
    out = list(samclip_filter_lines(header + lines, CONTIGS, max_clip=5))
    assert len(out) == 2 # header + r1
    assert "r1" in out[1]
