for file in 11_13*.txt; do

    echo "$file"
    python collect_and_analyze_zebrafish_behavior_data.py $file
done
