import re

infile  = "tavg_bflux.dat"
outfile = "tavg_bflux_corrected.dat"

with open(infile, "r") as f:
    lines = f.readlines()

out = []
for line in lines:
    stripped = line.rstrip("\n")

    # Header line: starts with '#Re'
    if stripped.startswith("#Re"):
        cols = stripped.split()
        # cols[0]='#Re', cols[1]='B', cols[2]='Pe', ...
        cols[1], cols[2] = cols[2], cols[1]
        out.append("    ".join(cols) + "    \n")

    # Data lines: not starting with '#'
    elif not stripped.startswith("#"):
        cols = stripped.split()
        if len(cols) >= 3:
            cols[1], cols[2] = cols[2], cols[1]
        out.append("    ".join(cols) + "    \n")

    # Comment lines (# Index ..., # Re = ..., etc.): pass through unchanged
    else:
        out.append(line)

with open(outfile, "w") as f:
    f.writelines(out)

print(f"Written {outfile}")
