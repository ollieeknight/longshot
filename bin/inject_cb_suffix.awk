#!/usr/bin/awk -f
# Appends a deterministic per-library suffix to every CB:Z: tag in a SAM
# stream, so cell barcodes stay unique when libraries are pooled for joint
# experiment-level isoform discovery. Usage: awk -v suf="_01" -f inject_cb_suffix.awk
{
    if ($0 ~ /^@/) {
        print
    } else {
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^CB:Z:/) { $i = $i suf }
        }
        print $0
    }
}
