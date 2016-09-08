version: 23c0dee

Description
===========

This issue is with converting coordinates with no clear equivalent in the reference (coordinates in insertions or deletions).

The "intervals" in the script each represent a single alignment between the  reference and the (merged) assembly. The endpoints of these alignments are thus always at coordinates which are defined in both the reference and the assembly. But the script sometimes splits them, which can result in an interval whose endpoints are defined in reference coordinates, but fall in an insertion or deletion in the assembly. The conversion of those reference coordinates to assembly coordinates is ambiguous, and can result in a final, "cleaned" sequence which has small, but distinct differences from the consensus in the input assembly.

N.B., the issue only comes up when there is an indel that falls right on the edge of a repeat region (a region of the reference present more than once in the assembly). Otherwise, indels are preserved perfectly.


Insertions:

If there is an insertion in the assembly right at the point where an interval ends up being split in clean-merge.py, the whole insertion will be lost. Instead, the output will contain the reference allele.

You can see this phenomenon in the "ins2" files (ins2{.lav,-cleaned.fa,-merge.fa,-raw.fa}). There is a 5bp insertion ("GATTA") in the assembly (ins2-raw.fa, ins2-merge.fa), relative to the reference (ref.fa). ins2-merge.fa is the same as ins2-raw.fa, except a ~250bp portion of the sequence has been copied from the middle and pasted to the end. This segment ends in the middle of the 5bp insertion. Then several SNPs were introduced into the original copy of this segment, so that the copy pasted onto the end is a better match to the original, raw assembly. This causes the algorithm to split the main, full-length alignment. There will now be a section before the repeat, where it takes the original sequence, a section for the repeat, where it uses the copy at the end, and a section after the repeat, where it uses the original sequence again.

N.B., no actual conversion error is encountered, so this cannot even be detected using the "throw" failure option. This is because each section aligns properly to the reference, but the insertion simply lies outside the alignment. So no reference coordinate is encountered which does not exist in the assembly.

The "ins" files are a simple example showing how an insertion can be preserved, if it stays in the middle of an interval. 

Here are the commands to replicate the effect:
$ lastz ref.fa ins2-merge.fa > ins2.lav
$ clean-merge.py ins2.lav ins2-merge.fa ins2-raw.fa > ins2-cleaned.fa


Deletions:

If the inverse happens, and there is a deletion right where an interval is split, it is likely that the deletion will be preserved, but with a 1bp insertion in it. You can see this in the "del" files.

This occurs because when converting coordinates that exist in the reference but not in the assembly (because of a deletion), the "tryharder" conversion algorithm ends up using the coordinate on the closest edge of the deletion. If the interval endpoint isn't in the exact center of the deletion, then both the end of the first interval and the start of the following one will be closest to the same edge of the deletion, and convert to the same coordinate. So you'll end up "repeating" one base.

I say "repeating" in quotes because it won't necessarily be the same base. To understand, let's look at an example where it *is* the same base: "del". Basically, here it's the same base because the deletion is identical in both copies of the repeated interval. In detail: clean-merge.py wants to include the second copy (at the end of del-merge.fa), because it's a better match to the raw assembly. The end of this interval exists in both the reference and the assembly, so it just adds the sequence of the copy at the end to the growing output. All is correct so far. Then, we return to the original, long sequence, which was interrupted by the repeat. This interval has been split, and we're currently looking at the second half, which starts at the deletion. But the way the interval was split results in the starting coordinate for this second half being just "where the repeat aligned on the reference + 1". The repeat aligned with its end coordinate being the base before the deletion. 1bp after this is in sequence that's deleted in the assembly; it doesn't exist. So it takes the closest coordinate that does exist: the base before the deletion. So it pastes in sequence from the assembly that starts at that base, and ends at the end of the "normal" sequence (before the final repeat). Problem is, it's already included the base before the deletion when it was adding the sequence from the repeat. Now, it's included that base again.

Now, "del2" is an example where it "repeats" 1bp, but a different base. This is just because the deletion in the second copy of the repeat is a bit different. It's a little smaller, so the final base comes 1bp later than the final base in the first copy. In the former case, it's a "T", and in the latter, it's an "A". So it's still repeating "the last base before the deletion", but it's a different base in each copy.



Solution
========

It'd be relatively easy to catch this issue in the case of assembly deletions, by catching coordinate conversion errors (fail='throw'). But catching it for insertions would be hard. It might require pre-processing the LAV to build a table of indel locations, and checking it every time an interval is split to make sure the new endpoints don't fall in one.
