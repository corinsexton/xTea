B16 (supplementary junction reads as split-read evidence) was implemented in 6325b4d and reverted:
on the 8 validation datasets it cost 44 TP / 13 FP Alu calls through is_polyA_dominant_one_side,
gave no GT gain (errors 120 -> 123) and moved het SR medians above 0.5. Kept here for later:
  b16_6325b4d.patch         the code change and its test, as committed
  test_gntp_supp_clips.py   the tests (not collected: outside tests/)
The measurement (before the fix) is in ../measure/ and ../summary_hclip*.txt, and the B16 arm's
validation outputs in ../b16_arm/.
