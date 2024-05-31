prepare_readgroup() {
    forward_read=$1
    samplename=$(basename "$forward_read")
    if [[ $forward_read == *"fastq.gz" ]]; then
        firstLine=$(gzip -dc "$forward_read" | head -n 1)
        if [[ $firstLine == *":"* ]]; then
            id_name=$(echo "$firstLine" | awk -F':' '{print $2}' | awk '{print $1}')
        elif [[ $firstLine == *"/"* ]]; then
            id_name=$(echo "$firstLine" | awk -F'/' '{print $2}' | awk '{print $1}')
        else
            id_name="1"
        fi
        split_field="@RG\\tID:${id_name}\\tSM:${samplename}\\tLB:1\\tPL:Illumina"
        echo "$split_field"
    else
        echo "Sequence read file extension not recognized."
        exit 1
    fi
}

# Call the function to prepare read group
split_field=$(prepare_readgroup "$1")
echo "split_field: $split_field"

# Print arguments for debugging
echo "Arguments:"
echo "1: $1"
echo "2: $2"
echo "3: $3"
echo "4: $4"
echo "5: $5" 

# create file name   
bwa mem -M -R "$split_field" -t "$4" "$3" "$1" "$2" > "$5"