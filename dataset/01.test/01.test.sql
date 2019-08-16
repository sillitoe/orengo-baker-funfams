SELECT
  cs.superfamily_id AS superfamily_id,
  cs.superfamily_id || '-ff-' || cs.funfam_number AS funfam_id,
  cs.num_members_in_seed_aln,
  ff.seed_dops_score,
  cs.funfam_name
FROM
  cluster_summary cs, funfam ff
WHERE
  cs.superfamily_id = ff.superfamily_id
  AND
  cs.funfam_number = ff.funfam_number
  AND
  cs.rep_source_id = 'cath'
  AND
  cs.num_ec4_codes = 1
  AND
  cs.num_members_in_seed_aln > 1000 AND cs.num_members_in_seed_aln < 5000
  AND
  cs.superfamily_id NOT LIKE '1.20.5.%' AND cs.superfamily_id NOT LIKE '4.%'
ORDER BY
  seed_dops_score DESC;
