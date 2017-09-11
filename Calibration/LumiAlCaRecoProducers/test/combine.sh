#BASH script that will combine all the database files in a given directory. 
for file in correction*.db
do
    echo "Combining file $file to Global"
    #below, TestCorrections is the current tag and must be updated to match 
    conddb -y --db $file copy TestCorrections TestCorrections --destdb GlobalCorr.db
done
