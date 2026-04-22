#!/bin/bash

array=(test-out/test-*)
max=0

for test_file in "${array[@]}"; do
    test_file_split=(${test_file//-/ })
    last_num="${test_file_split[-1]}"

    if [[ $last_num -gt $max ]]; then
        max=$last_num
    fi
done

test_num=$((max+1))

test_dir="test-out/test-${test_num}"

mkdir -p $test_dir

mv ./*.ms $test_dir
mv ./*.out $test_dir
mv ./*.fits $test_dir

echo "Moved files to ${test_dir}."
