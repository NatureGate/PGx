import os
import sys


MTRNR1_variants = {
	"rs267606618":"m.1095T>C",
	"rs267606619":"m.1494C>T",
	"rs267606617":"m.1555A>G",
	"rs28358569":"m.827A>G",
    "Reference":"Reference"
}

if __name__ == "__main__":
    MTRNR1_FILE = sys.argv[1]
    OUTCALL_FILE = sys.argv[2]

    # Check if MTRNR1.txt exists and is not empty

    if os.path.getsize(MTRNR1_FILE) == 0:
        with open(MTRNR1_FILE,'w') as f:
            f.write("Reference")
        

    with open(MTRNR1_FILE) as f:
        key = f.readline().strip()


    append_newline = False
    if os.path.exists(OUTCALL_FILE) and os.path.getsize(OUTCALL_FILE) > 0:
        if key in MTRNR1_variants.keys():
            with open(OUTCALL_FILE, 'rb') as f:
                # 另外起一行，将key和对应的MTRNR1_variants[key]写入outcall.tsv
                f.seek(-1, os.SEEK_END)
                last_char = f.read(1).decode()
                if last_char != '\n':
                    append_newline = True
                else:
                    append_newline = False

            with open(OUTCALL_FILE, 'a') as out:
                if append_newline:
                    out.write(f"\nMT-RNR1\t{MTRNR1_variants[key]}")
                else:
                    out.write(f"MT-RNR1\t{MTRNR1_variants[key]}")
    else:
        with open(OUTCALL_FILE, 'w') as out:
            out.write(f"MT-RNR1\t{MTRNR1_variants[key]}")