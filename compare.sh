cmp --silent scaffolds_test.txt scaffolds.txt && echo '### SUCCESS: scaffolds Are Identical! ###' || echo '### WARNING: scaffolds Are Different! ###'
cmp --silent staples_test.txt staples.txt && echo '### SUCCESS: staples Are Identical! ###' || echo '### WARNING: staples Are Different! ###'

