word=$(ps | grep test_regrid.sh)
kill "${word:0:8}"
word=$(ps | grep oskar)
kill "${word:0:8}"