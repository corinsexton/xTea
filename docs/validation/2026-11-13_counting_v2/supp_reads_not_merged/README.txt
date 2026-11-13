Supplementary junction reads as split-read evidence (AF_SR): implemented on branch counting-v2 in
6325b4d (labeled B16 there) and reverted in 7840dfd. Not merged: on the 8 validation datasets it
cost 44 TP / 13 FP Alu calls through is_polyA_dominant_one_side, gave no GT gain (errors 120 -> 123)
and moved het SR medians above 0.5. No bug ID is reserved for it. See the open issue in
docs/BUG_LOG.md.
  supp_reads_6325b4d.patch  the code change and its test, as committed
  test_gntp_supp_clips.py   the tests (not collected: outside tests/)
The measurement (before the fix) is in ../measure/ and ../summary_hclip*.txt; the validation of
the reverted arm (with its own refit, 5089d75) is in ../supp_reads_arm/.
