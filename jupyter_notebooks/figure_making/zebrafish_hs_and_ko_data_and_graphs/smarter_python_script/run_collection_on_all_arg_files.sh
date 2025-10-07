for file in *.txt; do

    echo "$file"
    python collect_and_analyze_zebrafish_behavior_data.py $file
done
