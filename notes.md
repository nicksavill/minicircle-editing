more realistic probabilities
test all gRNAs in "binding region" and calculate prob based on mfe

# COX3
min anchor 10

# CR3
min anchor 10

# CR4
min anchor 10
- CR4-II-165: max anchor length of 9, but anchors to G-rich region which might cause strong binding
- CR4-Orphan-302: max anchor length 10, but transcripts of gRNAs have untemplated poly-A tail which extends anchor
- CR4-II-426: max anchor length 8, but transcripts of gRNAs have untemplated poly-A tail which extends anchor

what is the proportion of correctly edited mRNA?
try with cr3 min anchor length of 10
need to use all expressed genes, but as we don't know where a non-canonical gene might start
we need to have sequence from start of initiation site and allow the anchor to bind after 
the start if the gene. gene domain = initiation site -> +42