diamond="/vol/software/bin/diamond-2.0.0"


# add all proteomes to the database
files="proteomes/*"
for f in $files; do
    echo "Processing $f"
    $diamond makedb --in "$f" -d test
done

